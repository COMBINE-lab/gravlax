#!/usr/bin/env python3
"""Paired end-to-end performance/RSS gate for `aie extend` uniform reports."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import shutil
import statistics
import subprocess
import sys
import tempfile
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def identity(path: Path) -> dict[str, object]:
    path = path.resolve()
    return {"path": str(path), "bytes": path.stat().st_size, "sha256": sha256(path)}


def parse_metrics(path: Path) -> tuple[float, int]:
    values = {}
    for line in path.read_text().splitlines():
        key, value = line.split("=", 1)
        values[key] = value
    return float(values["elapsed_seconds"]), int(values["max_rss_kib"])


def run_once(
    *,
    mode: str,
    binary: Path,
    archive: Path,
    gtf: Path,
    root: Path,
    repetition: int,
) -> dict[str, object]:
    out_gtf = root / f"{repetition}-{mode}.gtf"
    legacy_report = root / f"{repetition}-{mode}.tsv"
    uniform_report = root / f"{repetition}-{mode}.json"
    metrics = root / f"{repetition}-{mode}.time"
    stderr = root / f"{repetition}-{mode}.stderr"
    command = [
        str(binary),
        "extend",
        str(archive),
        "--gtf",
        str(gtf),
        "--out-gtf",
        str(out_gtf),
        "--report",
        str(legacy_report),
    ]
    if mode == "uniform":
        command.extend(["--report-format", "json", "--report-output", str(uniform_report)])
    timed = [
        "/usr/bin/time",
        "-f",
        "elapsed_seconds=%e\nmax_rss_kib=%M",
        "-o",
        str(metrics),
        *command,
    ]
    with open(os.devnull, "wb") as stdout, stderr.open("wb") as errors:
        subprocess.run(timed, stdout=stdout, stderr=errors, check=True)
    elapsed, rss = parse_metrics(metrics)
    result = {
        "mode": mode,
        "repetition": repetition,
        "elapsed_seconds": elapsed,
        "max_rss_kib": rss,
        "gtf_bytes": out_gtf.stat().st_size,
        "gtf_sha256": sha256(out_gtf),
        "legacy_report_bytes": legacy_report.stat().st_size,
        "legacy_report_sha256": sha256(legacy_report),
        "stderr_bytes": stderr.stat().st_size,
    }
    if mode == "uniform":
        document = json.loads(uniform_report.read_bytes())
        if document.get("result_schema") != "gravlax.extend.result.v1":
            raise RuntimeError("uniform report has the wrong result schema")
        result["uniform_report_bytes"] = uniform_report.stat().st_size
        result["genes_extended"] = document["data"]["summary"]["genes_extended"]
    out_gtf.unlink()
    legacy_report.unlink()
    if uniform_report.exists():
        uniform_report.unlink()
    metrics.unlink()
    stderr.unlink()
    return result


def median(records: list[dict[str, object]], field: str) -> float:
    return statistics.median(float(record[field]) for record in records)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate", required=True, type=Path)
    parser.add_argument("--baseline", required=True, type=Path)
    parser.add_argument("--archive", required=True, type=Path)
    parser.add_argument("--gtf", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--work-dir", type=Path)
    parser.add_argument("--repetitions", type=int, default=3)
    parser.add_argument("--legacy-time-ratio-limit", type=float, default=1.03)
    parser.add_argument("--uniform-time-ratio-limit", type=float, default=1.10)
    parser.add_argument("--rss-ratio-limit", type=float, default=1.10)
    args = parser.parse_args()
    if args.repetitions < 2:
        parser.error("--repetitions must be at least 2 for alternating order")
    if args.legacy_time_ratio_limit <= 0:
        parser.error("--legacy-time-ratio-limit must be positive")
    if args.uniform_time_ratio_limit <= 0:
        parser.error("--uniform-time-ratio-limit must be positive")
    if args.rss_ratio_limit <= 0:
        parser.error("--rss-ratio-limit must be positive")
    for path in [args.candidate, args.baseline, args.archive, args.gtf]:
        if not path.is_file():
            parser.error(f"missing input: {path}")

    temp_parent = str(args.work_dir) if args.work_dir else None
    root = Path(tempfile.mkdtemp(prefix="gravlax-extend-uniform-", dir=temp_parent))
    try:
        records: list[dict[str, object]] = []
        for repetition in range(args.repetitions):
            order = ["baseline", "candidate_legacy", "uniform"]
            if repetition % 2:
                order.reverse()
            for mode in order:
                binary = args.baseline if mode == "baseline" else args.candidate
                records.append(
                    run_once(
                        mode=mode,
                        binary=binary,
                        archive=args.archive,
                        gtf=args.gtf,
                        root=root,
                        repetition=repetition,
                    )
                )

        reference_gtf = {(row["gtf_bytes"], row["gtf_sha256"]) for row in records}
        reference_report = {
            (row["legacy_report_bytes"], row["legacy_report_sha256"]) for row in records
        }
        artifact_parity = {
            "extended_gtf_bytes_identical": len(reference_gtf) == 1,
            "legacy_report_bytes_identical": len(reference_report) == 1,
        }
        by_mode = {
            mode: [record for record in records if record["mode"] == mode]
            for mode in ["baseline", "candidate_legacy", "uniform"]
        }
        medians = {
            mode: {
                "elapsed_seconds": median(rows, "elapsed_seconds"),
                "max_rss_kib": median(rows, "max_rss_kib"),
            }
            for mode, rows in by_mode.items()
        }
        ratios = {
            "candidate_legacy_vs_baseline_elapsed": (
                medians["candidate_legacy"]["elapsed_seconds"]
                / medians["baseline"]["elapsed_seconds"]
            ),
            "candidate_legacy_vs_baseline_rss": (
                medians["candidate_legacy"]["max_rss_kib"]
                / medians["baseline"]["max_rss_kib"]
            ),
            "uniform_vs_candidate_legacy_elapsed": (
                medians["uniform"]["elapsed_seconds"]
                / medians["candidate_legacy"]["elapsed_seconds"]
            ),
            "uniform_vs_candidate_legacy_rss": (
                medians["uniform"]["max_rss_kib"]
                / medians["candidate_legacy"]["max_rss_kib"]
            ),
        }
        limits = {
            "candidate_legacy_vs_baseline_elapsed_ratio_max": args.legacy_time_ratio_limit,
            "uniform_vs_candidate_legacy_elapsed_ratio_max": args.uniform_time_ratio_limit,
            "rss_ratio_max": args.rss_ratio_limit,
        }
        gates = {
            "candidate_legacy_elapsed": (
                ratios["candidate_legacy_vs_baseline_elapsed"]
                <= args.legacy_time_ratio_limit
            ),
            "candidate_legacy_rss": (
                ratios["candidate_legacy_vs_baseline_rss"] <= args.rss_ratio_limit
            ),
            "uniform_elapsed": (
                ratios["uniform_vs_candidate_legacy_elapsed"]
                <= args.uniform_time_ratio_limit
            ),
            "uniform_rss": (
                ratios["uniform_vs_candidate_legacy_rss"] <= args.rss_ratio_limit
            ),
            **artifact_parity,
        }
        report = {
            "schema": "gravlax.benchmark.extend-uniform-io.v2",
            "harness": identity(Path(__file__)),
            "invocation": [str(Path(__file__).resolve()), *sys.argv[1:]],
            "fixture": {
                "archive": identity(args.archive),
                "gtf": identity(args.gtf),
                "baseline_binary": identity(args.baseline),
                "candidate_binary": identity(args.candidate),
            },
            "environment": {
                "platform": platform.platform(),
                "python": platform.python_version(),
                "cpu_count": os.cpu_count(),
            },
            "method": {
                "repetitions": args.repetitions,
                "order": "alternating baseline,candidate-legacy,uniform / reverse",
                "timing": "/usr/bin/time elapsed wall seconds",
                "rss": "/usr/bin/time maximum resident set size",
                "command_parameters": "extend defaults, legacy per-gene --report in every mode",
                "artifact_parity": "SHA-256 and byte length for extended GTF and legacy TSV",
            },
            "records": records,
            "medians": medians,
            "ratios": ratios,
            "limits": limits,
            "gates": gates,
            "artifact_parity": artifact_parity,
            "passed": all(gates.values()),
        }
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
        if not report["passed"]:
            raise SystemExit(1)
    finally:
        shutil.rmtree(root)


if __name__ == "__main__":
    main()
