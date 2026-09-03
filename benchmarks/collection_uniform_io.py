#!/usr/bin/env python3
"""Paired collection legacy/uniform timing, RSS, and scientific-parity gate.

The recommended corpus is a routed collection with many sidecar sections so directory loading is
visible rather than lost in process-start noise. The script alternates arm order within every block,
discards warmups, sends measured stdout to /dev/null, and samples Linux VmRSS for each child.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import statistics
import subprocess
import sys
import time
from pathlib import Path
from typing import Any


def file_identity(path: Path) -> dict[str, Any]:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return {
        "path": str(path.resolve()),
        "bytes": path.stat().st_size,
        "sha256": digest.hexdigest(),
    }


def median_absolute_deviation(values: list[float]) -> float:
    median = statistics.median(values)
    return statistics.median(abs(value - median) for value in values)


def command(args: argparse.Namespace, uniform: bool, top: int) -> list[str]:
    cmd = [
        str(args.binary),
        "collection",
        args.kind,
        str(args.collection),
        args.locus,
        "--top",
        str(top),
    ]
    if args.kind == "junction":
        cmd.extend(["--min-support", str(args.min_support)])
    cmd.extend(["--format", "json"] if uniform else ["--json"])
    return cmd


def capture(cmd: list[str]) -> dict[str, Any]:
    completed = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return json.loads(completed.stdout)


def named_table(envelope: dict[str, Any], name: str) -> dict[str, Any]:
    for table in envelope["data"]["tables"]:
        if table["name"] == name:
            return table
    raise AssertionError(f"uniform result lacks table {name!r}")


def table_dicts(table: dict[str, Any]) -> list[dict[str, Any]]:
    names = [field["name"] for field in table["schema"]["fields"]]
    return [dict(zip(names, row, strict=True)) for row in table["rows"]]


def assert_scientific_parity(kind: str, legacy: dict[str, Any], uniform: dict[str, Any]) -> None:
    assert uniform["$schema"] == "gravlax.result-envelope.v1"
    summary = uniform["data"]["summary"]
    samples = table_dicts(named_table(uniform, "samples"))
    cells = table_dicts(named_table(uniform, "cells"))
    legacy_samples = {row["sample"]: row for row in legacy["samples"]}
    uniform_samples = {row["sample"]: row for row in samples}
    assert legacy_samples.keys() == uniform_samples.keys()

    if kind == "junction":
        assert legacy["totals"] == {"umis": summary["umis"], "cells": summary["cells"]}
        sample_fields = ["present", "supporting_children", "umis", "cells"]
        cell_fields = ["barcode", "umis"]
        legacy_cell_key = "top_cells"
    elif kind == "region":
        assert legacy["totals"] == {
            "molecules": summary["molecules"],
            "umis": summary["umis"],
            "cells": summary["cells"],
        }
        sample_fields = ["present", "molecules", "umis", "cells"]
        cell_fields = ["barcode", "umis"]
        legacy_cell_key = "top_cells"
    else:
        raise AssertionError(f"unsupported parity kind {kind}")

    for sample, old in legacy_samples.items():
        new = uniform_samples[sample]
        for field in sample_fields:
            assert old[field] == new[field], (sample, field, old[field], new[field])
        old_cells = [{field: row[field] for field in cell_fields} for row in old[legacy_cell_key]]
        new_cells = [
            {field: row[field] for field in cell_fields}
            for row in cells
            if row["sample"] == sample
        ]
        assert old_cells == new_cells, (sample, old_cells, new_cells)


def read_rss_kib(pid: int) -> int:
    try:
        with open(f"/proc/{pid}/status", encoding="ascii") as handle:
            for line in handle:
                if line.startswith("VmRSS:"):
                    return int(line.split()[1])
    except FileNotFoundError:
        pass
    return 0


def measured_run(cmd: list[str]) -> tuple[float, int]:
    started = time.perf_counter()
    with open(os.devnull, "wb") as sink:
        process = subprocess.Popen(cmd, stdout=sink, stderr=sink)
        peak_rss = 0
        while process.poll() is None:
            peak_rss = max(peak_rss, read_rss_kib(process.pid))
            time.sleep(0.001)
        return_code = process.wait()
    elapsed = time.perf_counter() - started
    if return_code != 0:
        raise subprocess.CalledProcessError(return_code, cmd)
    return elapsed, peak_rss


def paired_case(args: argparse.Namespace, top: int) -> dict[str, Any]:
    legacy_cmd = command(args, False, top)
    uniform_cmd = command(args, True, top)
    for block in range(args.warmups):
        arms = (legacy_cmd, uniform_cmd) if block % 2 == 0 else (uniform_cmd, legacy_cmd)
        for cmd in arms:
            measured_run(cmd)

    times: dict[str, list[float]] = {"legacy": [], "uniform": []}
    rss: dict[str, list[int]] = {"legacy": [], "uniform": []}
    for block in range(args.iterations):
        arms = (
            (("legacy", legacy_cmd), ("uniform", uniform_cmd))
            if block % 2 == 0
            else (("uniform", uniform_cmd), ("legacy", legacy_cmd))
        )
        for name, cmd in arms:
            elapsed, peak = measured_run(cmd)
            times[name].append(elapsed)
            rss[name].append(peak)

    legacy_median = statistics.median(times["legacy"])
    uniform_median = statistics.median(times["uniform"])
    legacy_rss = max(rss["legacy"])
    uniform_rss = max(rss["uniform"])
    ratio = uniform_median / legacy_median
    rss_allowance = max(4096, int(legacy_rss * 0.05))
    return {
        "top": top,
        "commands": {
            "legacy": legacy_cmd,
            "uniform": uniform_cmd,
        },
        "legacy_seconds": times["legacy"],
        "uniform_seconds": times["uniform"],
        "legacy_median_seconds": legacy_median,
        "uniform_median_seconds": uniform_median,
        "legacy_mad_seconds": median_absolute_deviation(times["legacy"]),
        "uniform_mad_seconds": median_absolute_deviation(times["uniform"]),
        "legacy_rss_kib": rss["legacy"],
        "uniform_rss_kib": rss["uniform"],
        "uniform_over_legacy": ratio,
        "legacy_peak_rss_kib": legacy_rss,
        "uniform_peak_rss_kib": uniform_rss,
        "rss_allowance_kib": rss_allowance,
        "rss_gate_pass": uniform_rss <= legacy_rss + rss_allowance,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--collection", type=Path, required=True)
    parser.add_argument("--kind", choices=("junction", "region"), default="junction")
    parser.add_argument("--locus", required=True)
    parser.add_argument("--min-support", type=int, default=0)
    parser.add_argument("--scan-top", type=int, default=5)
    parser.add_argument("--output-top", type=int, default=0)
    parser.add_argument("--warmups", type=int, default=2)
    parser.add_argument("--iterations", type=int, default=8)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()

    if not args.binary.is_file() or not args.collection.is_file():
        parser.error("--binary and --collection must name existing files")
    legacy = capture(command(args, False, args.output_top))
    uniform = capture(command(args, True, args.output_top))
    assert_scientific_parity(args.kind, legacy, uniform)

    scan = paired_case(args, args.scan_top)
    output = paired_case(args, args.output_top)
    scan["time_gate"] = 1.03
    scan["time_gate_pass"] = scan["uniform_over_legacy"] <= 1.03
    output["time_gate"] = 1.10
    output["time_gate_pass"] = output["uniform_over_legacy"] <= 1.10
    result = {
        "schema": "gravlax.benchmark.collection-uniform-io.v2",
        "harness": file_identity(Path(__file__)),
        "binary": file_identity(args.binary),
        "collection": file_identity(args.collection),
        "invocation": sys.argv,
        "environment": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "cwd": str(Path.cwd()),
        },
        "kind": args.kind,
        "locus": args.locus,
        "min_support": args.min_support,
        "scan_top": args.scan_top,
        "output_top": args.output_top,
        "warmups": args.warmups,
        "iterations": args.iterations,
        "scientific_parity": True,
        "scan_case": scan,
        "output_case": output,
        "pass": all(
            case[gate]
            for case in (scan, output)
            for gate in ("time_gate_pass", "rss_gate_pass")
        ),
    }
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.out:
        args.out.write_text(encoded, encoding="utf-8")
    print(encoded, end="")
    if not result["pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
