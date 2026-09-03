#!/usr/bin/env python3
"""Paired scientific-artifact, wall-time, and RSS gate for cohort end reports.

Run this on a representative cohort design. Every timed invocation receives a fresh
scientific output directory. Legacy and uniform runs alternate AB/BA order; uniform
JSON is streamed to the null device so the gate measures encoding without retaining
the report in the harness. A separate report-file run verifies envelope completion.
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
import tempfile
import time
from pathlib import Path
from typing import Any


ARTIFACTS = {
    "transcript-ends": (
        "sites.tsv",
        "genes.tsv",
        "genes.polyasite.tsv",
        "polyasite-mixture-sites.tsv",
        "polyasite-mixture-genes.tsv",
        "fragment-kernel.tsv",
    ),
    "polyasite-mixture": ("sites.tsv", "genes.tsv", "fragment-kernel.tsv"),
}

RESULT_SCHEMAS = {
    "transcript-ends": "gravlax.cohort.transcript-ends.result.v1",
    "polyasite-mixture": "gravlax.cohort.polyasite-mixture.result.v1",
}

TABLES = {
    "transcript-ends": (
        "samples",
        "sites",
        "site_counts",
        "mixture_sites",
        "mixture_site_counts",
        "genes",
        "gene_usages",
        "fragment_kernel",
        "artifacts",
    ),
    "polyasite-mixture": (
        "samples",
        "sites",
        "site_counts",
        "genes",
        "gene_usages",
        "fragment_kernel",
        "heldout_kernel",
        "artifacts",
    ),
}


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as reader:
        while chunk := reader.read(1024 * 1024):
            value.update(chunk)
    return value.hexdigest()


def file_identity(path: Path) -> dict[str, Any]:
    path = path.resolve()
    return {"path": str(path), "bytes": path.stat().st_size, "sha256": digest(path)}


def design_inputs(design: Path) -> list[dict[str, Any]]:
    lines = design.read_text().splitlines()
    if not lines or lines[0].rstrip("\r") != "sample\tcondition\tarchive\tgroups":
        raise ValueError("design has the wrong header")
    base = design.parent
    inputs: list[dict[str, Any]] = []
    for line_number, line in enumerate(lines[1:], 2):
        fields = line.rstrip("\r").split("\t")
        if len(fields) != 4:
            raise ValueError(f"design line {line_number} does not have four fields")
        archive = Path(fields[2])
        groups = Path(fields[3])
        if not archive.is_absolute():
            archive = base / archive
        if not groups.is_absolute():
            groups = base / groups
        inputs.append(
            {
                "sample": fields[0],
                "condition": fields[1],
                "archive": file_identity(archive),
                "groups": file_identity(groups),
            }
        )
    return inputs


def median_summary(values: list[float]) -> dict[str, float]:
    median = statistics.median(values)
    return {
        "median": median,
        "mad": statistics.median(abs(value - median) for value in values),
        "min": min(values),
        "max": max(values),
    }


def rss_kib(pid: int) -> int:
    try:
        with open(f"/proc/{pid}/status", encoding="ascii") as status:
            for line in status:
                if line.startswith("VmRSS:"):
                    return int(line.split()[1])
    except FileNotFoundError:
        return 0
    return 0


def measured_run(command: list[str]) -> tuple[float, int]:
    started = time.perf_counter()
    with open(os.devnull, "wb") as sink:
        process = subprocess.Popen(command, stdout=sink, stderr=sink)
        peak = 0
        while process.poll() is None:
            peak = max(peak, rss_kib(process.pid))
            time.sleep(0.001)
        return_code = process.wait()
    if return_code:
        raise subprocess.CalledProcessError(return_code, command)
    return time.perf_counter() - started, peak


def base_command(args: argparse.Namespace) -> list[str]:
    command = [
        str(args.binary),
        "cohort",
        args.kind,
        "--design",
        str(args.design),
        "--gtf",
        str(args.annotation),
        "--genome",
        str(args.genome),
        "--polyasite",
        str(args.polyasite),
        "--group-contrast",
        args.group_contrast,
    ]
    for option, value in (
        ("--site-gap", args.site_gap),
        ("--tail-extend", args.tail_extend),
        ("--min-site-umis", args.min_site_umis),
        ("--min-site-samples", args.min_site_samples),
        ("--min-group-gene-umis", args.min_group_gene_umis),
        ("--min-samples", args.min_samples),
        ("--min-distal-umis", args.min_distal_umis),
        ("--max-sites", args.max_sites),
    ):
        command.extend([option, str(value)])
    if args.kind == "transcript-ends":
        command.extend(["--motif-min-samples", str(args.motif_min_samples)])
    if args.shuffle_seed is not None:
        command.extend(["--shuffle-seed", str(args.shuffle_seed)])
    return command


def command(
    args: argparse.Namespace,
    output: Path,
    uniform: bool,
    report_output: Path | None = None,
) -> list[str]:
    result = [*base_command(args), "--out-dir", str(output)]
    if uniform:
        result.extend(["--report-format", "json"])
        if report_output is not None:
            result.extend(["--report-output", str(report_output)])
    return result


def artifact_parity(left: Path, right: Path, kind: str) -> dict[str, bool]:
    return {
        name: digest(left / name) == digest(right / name)
        for name in ARTIFACTS[kind]
    }


def envelope_probe(path: Path, kind: str) -> dict[str, Any]:
    size = path.stat().st_size
    with path.open("rb") as reader:
        prefix = reader.read(min(size, 2 * 1024 * 1024))
        if size > 64:
            reader.seek(-64, os.SEEK_END)
        suffix = reader.read()
    schema = RESULT_SCHEMAS[kind].encode()
    missing_tables = [
        name
        for name in TABLES[kind]
        if f'"name":"{name}"'.encode() not in prefix
    ]
    # Schemas for late tables can begin beyond the prefix on a large result. A streaming
    # substring pass validates names without deserializing the result into harness memory.
    if missing_tables:
        remaining = set(missing_tables)
        with path.open("rb") as reader:
            overlap = b""
            while remaining and (chunk := reader.read(1024 * 1024)):
                haystack = overlap + chunk
                for name in tuple(remaining):
                    if f'"name":"{name}"'.encode() in haystack:
                        remaining.remove(name)
                overlap = haystack[-256:]
        missing_tables = sorted(remaining)
    return {
        "bytes": size,
        "result_schema_present": schema in prefix,
        "missing_tables": missing_tables,
        "complete_json_suffix": suffix.endswith(b"]}}\n"),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument(
        "--kind", required=True, choices=("transcript-ends", "polyasite-mixture")
    )
    parser.add_argument("--design", required=True, type=Path)
    parser.add_argument("--annotation", required=True, type=Path)
    parser.add_argument("--genome", required=True, type=Path)
    parser.add_argument("--polyasite", required=True, type=Path)
    parser.add_argument("--group-contrast", required=True)
    parser.add_argument("--site-gap", type=int, default=24)
    parser.add_argument("--tail-extend", type=int, default=2_000)
    parser.add_argument("--min-site-umis", type=int, default=10)
    parser.add_argument("--min-site-samples", type=int, default=3)
    parser.add_argument("--motif-min-samples", type=int, default=4)
    parser.add_argument("--min-group-gene-umis", type=int, default=20)
    parser.add_argument("--min-samples", type=int, default=6)
    parser.add_argument("--min-distal-umis", type=int, default=20)
    parser.add_argument("--max-sites", type=int, default=1_000_000)
    parser.add_argument("--shuffle-seed", type=int)
    parser.add_argument("--repetitions", type=int, default=7)
    parser.add_argument("--warmups", type=int, default=1)
    parser.add_argument("--time-ratio-limit", type=float, default=1.10)
    parser.add_argument("--rss-ratio-limit", type=float, default=1.05)
    parser.add_argument("--min-report-bytes", type=int, default=1_000_000)
    parser.add_argument("--scratch-root", type=Path)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    for path in (
        args.binary,
        args.design,
        args.annotation,
        args.genome,
        args.polyasite,
    ):
        if not path.is_file():
            parser.error(f"not a file: {path}")
    if args.repetitions < 1 or args.warmups < 0:
        parser.error("repetitions must be positive and warmups nonnegative")
    if args.scratch_root is not None:
        args.scratch_root.mkdir(parents=True, exist_ok=True)

    legacy_times: list[float] = []
    uniform_times: list[float] = []
    legacy_rss: list[float] = []
    uniform_rss: list[float] = []
    with tempfile.TemporaryDirectory(
        prefix="gravlax-end-uniform-", dir=args.scratch_root
    ) as directory:
        root = Path(directory)
        for block in range(args.repetitions + args.warmups):
            order = ("legacy", "uniform") if block % 2 == 0 else ("uniform", "legacy")
            for label in order:
                output = root / f"timed-{block}-{label}"
                invocation = command(args, output, label == "uniform")
                elapsed, rss = measured_run(invocation)
                if block < args.warmups:
                    continue
                (legacy_times if label == "legacy" else uniform_times).append(elapsed)
                (legacy_rss if label == "legacy" else uniform_rss).append(float(rss))

        parity_legacy = root / "parity-legacy"
        parity_uniform = root / "parity-uniform"
        report = root / "report.json"
        subprocess.run(
            command(args, parity_legacy, False),
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        subprocess.run(
            command(args, parity_uniform, True, report),
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        parity = artifact_parity(parity_legacy, parity_uniform, args.kind)
        envelope = envelope_probe(report, args.kind)

    legacy_time = statistics.median(legacy_times)
    uniform_time = statistics.median(uniform_times)
    legacy_memory = statistics.median(legacy_rss)
    uniform_memory = statistics.median(uniform_rss)
    time_ratio = uniform_time / legacy_time
    rss_ratio = uniform_memory / legacy_memory
    gates = {
        "time": time_ratio <= args.time_ratio_limit,
        "rss": rss_ratio <= args.rss_ratio_limit,
        "scientific_artifact_bytes": all(parity.values()),
        "result_schema": envelope["result_schema_present"],
        "named_tables": not envelope["missing_tables"],
        "complete_report": envelope["complete_json_suffix"],
        "representative_report_size": envelope["bytes"] >= args.min_report_bytes,
    }
    result = {
        "$schema": "gravlax.benchmark.cohort-end-uniform-io.v1",
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "cpu_affinity": sorted(os.sched_getaffinity(0))
            if hasattr(os, "sched_getaffinity")
            else None,
        },
        "inputs": {
            "binary": file_identity(args.binary),
            "design": file_identity(args.design),
            "annotation": file_identity(args.annotation),
            "genome": file_identity(args.genome),
            "polyasite": file_identity(args.polyasite),
            "design_rows": design_inputs(args.design),
            "command": base_command(args),
        },
        "measurement": {
            "repetitions": args.repetitions,
            "warmups": args.warmups,
            "alternating_order": True,
            "legacy_seconds": legacy_times,
            "uniform_seconds": uniform_times,
            "legacy_time": median_summary(legacy_times),
            "uniform_time": median_summary(uniform_times),
            "time_ratio": time_ratio,
            "legacy_peak_rss_kib": legacy_rss,
            "uniform_peak_rss_kib": uniform_rss,
            "legacy_rss": median_summary(legacy_rss),
            "uniform_rss": median_summary(uniform_rss),
            "rss_ratio": rss_ratio,
        },
        "semantic": {
            "scientific_artifact_byte_parity": parity,
            "report": envelope,
        },
        "limits": {
            "time_ratio": args.time_ratio_limit,
            "rss_ratio": args.rss_ratio_limit,
            "min_report_bytes": args.min_report_bytes,
        },
        "gates": gates,
        "passed": all(gates.values()),
    }
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.out is None:
        print(encoded, end="")
    else:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(encoded)
    if not result["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
