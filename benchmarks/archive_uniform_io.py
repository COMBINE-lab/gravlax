#!/usr/bin/env python3
"""Paired time/RSS and byte-parity gate for archive uniform I/O.

The scan gate compares legacy and uniform ``inspect-archive`` on the same archive.
The output gate compares legacy and uniformly reported ``replay-rows`` runs in fresh
directories.  Runs alternate order to reduce warm-cache and scheduler bias.
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


def file_identity(path: Path) -> dict[str, Any]:
    resolved = path.resolve()
    return {
        "path": str(resolved),
        "bytes": resolved.stat().st_size,
        "sha256": digest(resolved),
    }


def median_absolute_deviation(values: list[float]) -> float:
    median = statistics.median(values)
    return statistics.median(abs(value - median) for value in values)


def read_rss_kib(pid: int) -> int:
    try:
        with open(f"/proc/{pid}/status", encoding="ascii") as status:
            for line in status:
                if line.startswith("VmRSS:"):
                    return int(line.split()[1])
    except FileNotFoundError:
        pass
    return 0


def measured_run(command: list[str]) -> tuple[float, int]:
    started = time.perf_counter()
    with open(os.devnull, "wb") as sink:
        process = subprocess.Popen(command, stdout=sink, stderr=sink)
        peak_rss = 0
        while process.poll() is None:
            peak_rss = max(peak_rss, read_rss_kib(process.pid))
            time.sleep(0.001)
        return_code = process.wait()
    elapsed = time.perf_counter() - started
    if return_code:
        raise subprocess.CalledProcessError(return_code, command)
    return elapsed, peak_rss


def paired_measure(
    legacy_command: list[str],
    uniform_command: list[str],
    repetitions: int,
    warmups: int,
) -> dict[str, Any]:
    legacy_times: list[float] = []
    uniform_times: list[float] = []
    legacy_rss: list[int] = []
    uniform_rss: list[int] = []
    for block in range(repetitions + warmups):
        order = (
            (("legacy", legacy_command), ("uniform", uniform_command))
            if block % 2 == 0
            else (("uniform", uniform_command), ("legacy", legacy_command))
        )
        for label, command in order:
            elapsed, rss = measured_run(command)
            if block < warmups:
                continue
            if label == "legacy":
                legacy_times.append(elapsed)
                legacy_rss.append(rss)
            else:
                uniform_times.append(elapsed)
                uniform_rss.append(rss)
    legacy_median = statistics.median(legacy_times)
    uniform_median = statistics.median(uniform_times)
    legacy_rss_median = statistics.median(legacy_rss)
    uniform_rss_median = statistics.median(uniform_rss)
    return {
        "commands": {
            "legacy": legacy_command,
            "uniform": uniform_command,
        },
        "legacy_seconds": legacy_times,
        "uniform_seconds": uniform_times,
        "legacy_median_seconds": legacy_median,
        "uniform_median_seconds": uniform_median,
        "legacy_mad_seconds": median_absolute_deviation(legacy_times),
        "uniform_mad_seconds": median_absolute_deviation(uniform_times),
        "time_ratio": uniform_median / legacy_median,
        "legacy_peak_rss_kib": legacy_rss,
        "uniform_peak_rss_kib": uniform_rss,
        "legacy_median_peak_rss_kib": legacy_rss_median,
        "uniform_median_peak_rss_kib": uniform_rss_median,
        "rss_ratio": uniform_rss_median / legacy_rss_median,
    }


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as reader:
        while chunk := reader.read(1024 * 1024):
            value.update(chunk)
    return value.hexdigest()


def replay_command(
    binary: Path,
    archive: Path,
    annotation: Path,
    barcodes: Path,
    output: Path,
    report: Path | None,
) -> list[str]:
    command = [
        str(binary),
        "replay-rows",
        str(archive),
        "--gtf",
        str(annotation),
        "--barcodes",
        str(barcodes),
        "--out-dir",
        str(output),
    ]
    if report is not None:
        command.extend(["--report-format", "json", "--report-output", str(report)])
    return command


def replay_gate(args: argparse.Namespace, root: Path) -> dict[str, Any]:
    # Replay cannot reuse output paths: metadata and reports are deliberately no-clobber.
    legacy_times: list[float] = []
    uniform_times: list[float] = []
    legacy_rss: list[int] = []
    uniform_rss: list[int] = []
    parity_pair: tuple[Path, Path] | None = None
    for block in range(args.replay_repetitions + args.warmups):
        order = ("legacy", "uniform") if block % 2 == 0 else ("uniform", "legacy")
        outputs: dict[str, Path] = {}
        for label in order:
            output = root / f"replay-{block}-{label}"
            output.mkdir()
            outputs[label] = output
            report = root / f"replay-{block}-report.json" if label == "uniform" else None
            command = replay_command(
                args.binary,
                args.archive,
                args.annotation,
                args.barcodes,
                output,
                report,
            )
            elapsed, rss = measured_run(command)
            if block < args.warmups:
                continue
            if label == "legacy":
                legacy_times.append(elapsed)
                legacy_rss.append(rss)
            else:
                uniform_times.append(elapsed)
                uniform_rss.append(rss)
        if block == args.warmups:
            parity_pair = outputs["legacy"], outputs["uniform"]

    assert parity_pair is not None
    files = ("matrix.mtx", "features.tsv", "barcodes.tsv")
    parity = {
        name: digest(parity_pair[0] / name) == digest(parity_pair[1] / name)
        for name in files
    }
    metadata = json.loads((parity_pair[1] / "metadata.json").read_text())
    legacy_median = statistics.median(legacy_times)
    uniform_median = statistics.median(uniform_times)
    legacy_rss_median = statistics.median(legacy_rss)
    uniform_rss_median = statistics.median(uniform_rss)
    return {
        "command_template": {
            "legacy": replay_command(
                args.binary,
                args.archive,
                args.annotation,
                args.barcodes,
                Path("<fresh-output-directory>"),
                None,
            ),
            "uniform": replay_command(
                args.binary,
                args.archive,
                args.annotation,
                args.barcodes,
                Path("<fresh-output-directory>"),
                Path("<fresh-report.json>"),
            ),
        },
        "legacy_seconds": legacy_times,
        "uniform_seconds": uniform_times,
        "legacy_median_seconds": legacy_median,
        "uniform_median_seconds": uniform_median,
        "legacy_mad_seconds": median_absolute_deviation(legacy_times),
        "uniform_mad_seconds": median_absolute_deviation(uniform_times),
        "time_ratio": uniform_median / legacy_median,
        "legacy_peak_rss_kib": legacy_rss,
        "uniform_peak_rss_kib": uniform_rss,
        "legacy_median_peak_rss_kib": legacy_rss_median,
        "uniform_median_peak_rss_kib": uniform_rss_median,
        "rss_ratio": uniform_rss_median / legacy_rss_median,
        "scientific_artifact_byte_parity": parity,
        "metadata_schema": metadata.get("result_schema"),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--archive", required=True, type=Path)
    parser.add_argument("--annotation", required=True, type=Path)
    parser.add_argument("--barcodes", required=True, type=Path)
    parser.add_argument("--scan-repetitions", type=int, default=15)
    parser.add_argument("--replay-repetitions", type=int, default=5)
    parser.add_argument("--warmups", type=int, default=1)
    parser.add_argument("--scratch-root", type=Path)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    for path in (args.binary, args.archive, args.annotation, args.barcodes):
        if not path.is_file():
            parser.error(f"not a file: {path}")

    legacy_inspect = [str(args.binary), "inspect-archive", str(args.archive)]
    uniform_inspect = [
        str(args.binary),
        "inspect-archive",
        str(args.archive),
        "--format",
        "json",
    ]
    scan = paired_measure(
        legacy_inspect,
        uniform_inspect,
        args.scan_repetitions,
        args.warmups,
    )
    if args.scratch_root is not None:
        args.scratch_root.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(
        prefix="gravlax-archive-io-",
        dir=args.scratch_root,
    ) as directory:
        replay = replay_gate(args, Path(directory))

    time_scan_limit = 1.03
    time_output_limit = 1.10
    rss_limit = 1.10
    gates = {
        "inspect_time": scan["time_ratio"] <= time_scan_limit,
        "inspect_rss": scan["rss_ratio"] <= rss_limit,
        "replay_time": replay["time_ratio"] <= time_output_limit,
        "replay_rss": replay["rss_ratio"] <= rss_limit,
        "replay_bytes": all(replay["scientific_artifact_byte_parity"].values()),
        "replay_metadata_schema": replay["metadata_schema"]
        == "gravlax.replay.mex-artifact.v1",
    }
    result = {
        "$schema": "gravlax.benchmark.archive-uniform-io.v2",
        "harness": file_identity(Path(__file__)),
        "invocation": sys.argv,
        "environment": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "cwd": str(Path.cwd()),
        },
        "inputs": {
            "binary": file_identity(args.binary),
            "archive": file_identity(args.archive),
            "annotation": file_identity(args.annotation),
            "barcodes": file_identity(args.barcodes),
        },
        "configuration": {
            "scan_repetitions": args.scan_repetitions,
            "replay_repetitions": args.replay_repetitions,
            "warmups": args.warmups,
            "scratch_root": str(args.scratch_root.resolve())
            if args.scratch_root is not None
            else None,
        },
        "limits": {
            "inspect_time_ratio": time_scan_limit,
            "replay_time_ratio": time_output_limit,
            "rss_ratio": rss_limit,
        },
        "inspect": scan,
        "replay": replay,
        "gates": gates,
        "pass": all(gates.values()),
    }
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.out:
        args.out.write_text(encoded)
    print(encoded, end="")
    if not result["pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
