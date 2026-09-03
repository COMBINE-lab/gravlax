#!/usr/bin/env python3
"""Paired semantic, legacy-compatibility, time, and RSS gate for query uniform I/O."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import statistics
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Case:
    name: str
    legacy: tuple[str, ...]
    uniform: tuple[str, ...]
    time_limit: float


def run(argv: list[str], capture: bool) -> subprocess.CompletedProcess[bytes]:
    return subprocess.run(
        argv,
        check=True,
        stdout=subprocess.PIPE if capture else subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def path_identity(path: Path) -> dict[str, str | int]:
    hasher = hashlib.sha256()
    with path.open("rb") as reader:
        while chunk := reader.read(1024 * 1024):
            hasher.update(chunk)
    return {
        "path": str(path.resolve()),
        "size_bytes": path.stat().st_size,
        "sha256": hasher.hexdigest(),
    }


def historical_path_identity(path: Path) -> dict[str, object]:
    """Record a mutable artifact's measured location without claiming it is current."""
    result: dict[str, object] = path_identity(path)
    result["path_at_measurement"] = result.pop("path")
    result["path_semantics"] = (
        "Historical path at measurement time; size_bytes and sha256 are authoritative."
    )
    return result


def sample_summary(values: list[float]) -> dict[str, float]:
    median = statistics.median(values)
    return {
        "median": median,
        "mad": statistics.median(abs(value - median) for value in values),
        "min": min(values),
        "max": max(values),
    }


def paired_times(
    left: list[str], right: list[str], iterations: int, warmups: int
) -> tuple[list[float], list[float]]:
    samples = ([], [])
    for index in range(iterations + warmups):
        order = ((0, left), (1, right))
        if index % 2:
            order = tuple(reversed(order))
        for side, argv in order:
            started = time.perf_counter()
            run(argv, False)
            elapsed = time.perf_counter() - started
            if index >= warmups:
                samples[side].append(elapsed)
    return samples


def max_rss_kib(argv: list[str], record: Path) -> int:
    subprocess.run(
        ["/usr/bin/time", "-f", "%M", "-o", str(record), *argv],
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    return int(record.read_text().strip())


def paired_rss(
    left: list[str], right: list[str], repetitions: int, scratch: Path, label: str
) -> tuple[list[float], list[float]]:
    samples = ([], [])
    record = scratch / f"{label}.rss"
    for index in range(repetitions):
        order = ((0, left), (1, right))
        if index % 2:
            order = tuple(reversed(order))
        for side, argv in order:
            samples[side].append(float(max_rss_kib(argv, record)))
    return samples


def table(result: dict, name: str) -> list[dict]:
    value = next(item for item in result["data"]["tables"] if item["name"] == name)
    fields = [field["name"] for field in value["schema"]["fields"]]
    return [dict(zip(fields, row, strict=True)) for row in value["rows"]]


def semantic_batch(legacy: dict, uniform: dict) -> bool:
    normalized = []
    counts = []
    for query in legacy["queries"]:
        normalized.append(
            {
                "id": query["id"],
                "kind": query["kind"],
                "chrom": query["chrom"],
                "start": query["start"],
                "end": query["end"],
                "present": query.get("present", True),
                "archive_supporting_children": query.get("supporting_children"),
                "archive_posting_chunks": query.get("posting_chunks"),
                "molecules": query.get("molecules"),
                "anchor_semantics": query.get("anchor_semantics"),
                "umis": query["umis"],
                "cells": query["cells"],
            }
        )
        counts.extend(
            (query["id"], "cell", row["barcode"], row["umis"])
            for row in query["cell_rows"]
        )
    actual = [
        {key: row[key] for key in normalized[0]}
        for row in table(uniform, "queries")
    ]
    actual_counts = [
        (row["query_id"], row["aggregation"], row["entity"], row["umis"])
        for row in table(uniform, "counts")
    ]
    return normalized == actual and sorted(counts) == sorted(actual_counts)


def semantic_junctions(legacy: dict, uniform: dict) -> bool:
    expected = [
        (
            row["chrom"], row["donor"], row["acceptor"], row["supporting_children"],
            row["posting_chunks"], row.get("umis"), row.get("cells"),
        )
        for row in legacy["junctions"]
    ]
    actual_rows = table(uniform, "junctions")
    actual = [
        (
            row["chrom"], row["donor"], row["acceptor"],
            row["archive_supporting_children"], row["archive_posting_chunks"],
            row["umis"], row["cells"],
        )
        for row in actual_rows
    ]
    expected_counts = sorted(
        (junction["donor"], junction["acceptor"], row["barcode"], row["umis"])
        for junction in legacy["junctions"]
        for row in junction.get("cell_counts", [])
    )
    actual_counts = sorted(
        (row["donor"], row["acceptor"], row["entity"], row["umis"])
        for row in table(uniform, "counts")
    )
    return expected == actual and expected_counts == actual_counts


def semantic_jset(legacy: dict, uniform: dict) -> bool:
    requests = []
    for side, name in (("include", "inclusion_junctions"), ("exclude", "exclusion_junctions")):
        requests.extend(
            (side, row["locus"], row["present"], row.get("supporting_children"), row.get("posting_chunks"))
            for row in legacy[name]
        )
    actual_requests = [
        (
            row["side"], row["locus"], row["present"],
            row["archive_supporting_children"], row["archive_posting_chunks"],
        )
        for row in table(uniform, "junctions")
    ]
    expected_counts = sorted(
        (row["barcode"], row["include_only"], row["exclude_only"], row["both"])
        for row in legacy["cell_rows"]
    )
    actual_counts = sorted(
        (row["entity"], row["include_only"], row["exclude_only"], row["both"])
        for row in table(uniform, "counts")
    )
    return (
        legacy["totals"] == uniform["data"]["summary"]["totals"]
        and requests == actual_requests
        and expected_counts == actual_counts
    )


def semantic_events(legacy: dict, uniform: dict) -> bool:
    expected = sorted(
        (
            event["id"], event["event_type"], event["totals"]["include_only"],
            event["totals"]["exclude_only"], event["totals"]["both"],
            event["totals"]["informative_umis"],
        )
        for event in legacy["events"]
    )
    actual = sorted(
        (
            row["event_id"], row["event_type"], row["include_only"],
            row["exclude_only"], row["both"], row["informative_umis"],
        )
        for row in table(uniform, "events")
    )
    expected_counts = sorted(
        (event["id"], row["barcode"], row["include_only"], row["exclude_only"], row["both"])
        for event in legacy["events"]
        for row in event["cell_rows"]
    )
    actual_counts = sorted(
        (row["event_id"], row["entity"], row["include_only"], row["exclude_only"], row["both"])
        for row in table(uniform, "counts")
    )
    return expected == actual and expected_counts == actual_counts


def semantic_graph(legacy: dict, uniform: dict) -> bool:
    edges = sorted(
        (row["strand"], row["donor"], row["acceptor"], row["umis"], row["cells"])
        for row in legacy["edges"]
    )
    uniform_edges = sorted(
        (row["strand"], row["donor"], row["acceptor"], row["umis"], row["cells"])
        for row in table(uniform, "edges")
    )
    paths = sorted(
        (row["strand"], json.dumps(row["junctions"], sort_keys=True), row["umis"], row["cells"])
        for row in legacy["paths"]
    )
    uniform_paths = sorted(
        (row["strand"], json.dumps(row["junctions"], sort_keys=True), row["umis"], row["cells"])
        for row in table(uniform, "paths")
    )
    return legacy["totals"] == uniform["data"]["summary"]["totals"] and edges == uniform_edges and paths == uniform_paths


def semantic_federate(legacy_text: str, uniform: dict) -> bool:
    observed = []
    for line in legacy_text.splitlines():
        match = re.match(r"^.+: (\d+) UMIs / (\d+) cells ", line)
        if match:
            observed.append((int(match.group(1)), int(match.group(2))))
    actual = [(row["umis"], row["cells"]) for row in table(uniform, "archives")]
    summary = uniform["data"]["summary"]["totals"]
    return observed == actual and summary == {
        "umis": sum(row[0] for row in observed),
        "cells": sum(row[1] for row in observed),
    }


def semantic(case: str, legacy_bytes: bytes, uniform_bytes: bytes) -> bool:
    uniform = json.loads(uniform_bytes)
    if case == "federate":
        return semantic_federate(legacy_bytes.decode(), uniform)
    legacy = json.loads(legacy_bytes)
    return {
        "batch": semantic_batch,
        "junctions": semantic_junctions,
        "jset": semantic_jset,
        "events": semantic_events,
        "splice-graph": semantic_graph,
    }[case](legacy, uniform)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--junction-window", default="chr11:35138870-35232402")
    parser.add_argument("--junction-a", default="chr11:35139244-35139275")
    parser.add_argument("--junction-b", default="chr11:35139731-35140738")
    parser.add_argument("--dense-junction", default="chr16:89562391-89562883")
    parser.add_argument("--iterations", type=int, default=11)
    parser.add_argument("--warmups", type=int, default=2)
    parser.add_argument("--rss-runs", type=int, default=9)
    parser.add_argument("--scratch-root", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    args.scratch_root.mkdir(parents=True, exist_ok=True)

    archive = str(args.archive.resolve())
    cases = [
        Case("batch", ("query", archive, "batch", "--plan", str(args.plan), "--top", "0"),
             ("query", archive, "batch", "--plan", str(args.plan), "--top", "0", "--format", "json"), 1.10),
        Case("junctions", ("query", archive, "junctions", args.junction_window, "--with-cells", "--json"),
             ("query", archive, "junctions", args.junction_window, "--with-cells", "--format", "json"), 1.03),
        Case("jset", ("query", archive, "jset", "--include", args.junction_a, "--exclude", args.junction_b, "--top", "0", "--json"),
             ("query", archive, "jset", "--include", args.junction_a, "--exclude", args.junction_b, "--top", "0", "--format", "json"), 1.03),
        Case("events", ("query", archive, "events", args.junction_window, "--min-support", "1", "--min-informative", "1", "--top", "0", "--json"),
             ("query", archive, "events", args.junction_window, "--min-support", "1", "--min-informative", "1", "--top", "0", "--format", "json"), 1.03),
        Case("splice-graph", ("query", archive, "splice-graph", args.junction_window, "--min-support", "1", "--min-path-umis", "1", "--json"),
             ("query", archive, "splice-graph", args.junction_window, "--min-support", "1", "--min-path-umis", "1", "--format", "json"), 1.03),
        Case("federate", ("federate", archive, archive, args.dense_junction, "--top", "20"),
             ("federate", archive, archive, args.dense_junction, "--top", "20", "--format", "json"), 1.03),
    ]

    results = []
    for case in cases:
        baseline_cmd = [str(args.baseline.resolve()), *case.legacy]
        legacy_cmd = [str(args.candidate.resolve()), *case.legacy]
        uniform_cmd = [str(args.candidate.resolve()), *case.uniform]
        baseline_stdout = run(baseline_cmd, True).stdout
        legacy_stdout = run(legacy_cmd, True).stdout
        uniform_stdout = run(uniform_cmd, True).stdout
        legacy_exact = baseline_stdout == legacy_stdout
        semantic_equal = semantic(case.name, legacy_stdout, uniform_stdout)
        baseline_times, candidate_times = paired_times(
            baseline_cmd, legacy_cmd, args.iterations, args.warmups
        )
        legacy_times, uniform_times = paired_times(
            legacy_cmd, uniform_cmd, args.iterations, args.warmups
        )
        legacy_rss, uniform_rss = paired_rss(
            legacy_cmd, uniform_cmd, args.rss_runs, args.scratch_root, case.name
        )
        legacy_time = statistics.median(legacy_times)
        uniform_time = statistics.median(uniform_times)
        baseline_time = statistics.median(baseline_times)
        candidate_time = statistics.median(candidate_times)
        legacy_memory = statistics.median(legacy_rss)
        uniform_memory = statistics.median(uniform_rss)
        result = {
            "case": case.name,
            "legacy_stdout_exact": legacy_exact,
            "legacy_stdout_sha256": digest(legacy_stdout),
            "semantic_equivalence": semantic_equal,
            "legacy_stdout_bytes": len(legacy_stdout),
            "uniform_stdout_bytes": len(uniform_stdout),
            "legacy_compatibility": {
                "baseline_seconds": sample_summary(baseline_times),
                "candidate_seconds": sample_summary(candidate_times),
                "time_ratio": candidate_time / baseline_time,
                "gate": candidate_time <= baseline_time * 1.03,
            },
            "legacy_seconds": sample_summary(legacy_times),
            "uniform_seconds": sample_summary(uniform_times),
            "time_ratio": uniform_time / legacy_time,
            "time_limit": case.time_limit,
            "time_gate": uniform_time <= legacy_time * case.time_limit,
            "legacy_rss_kib": sample_summary(legacy_rss),
            "uniform_rss_kib": sample_summary(uniform_rss),
            "rss_ratio": uniform_memory / legacy_memory,
            "bounded_memory_gate": uniform_memory <= legacy_memory * 1.10,
        }
        result["gate"] = all(
            [legacy_exact, semantic_equal, result["legacy_compatibility"]["gate"],
             result["time_gate"], result["bounded_memory_gate"]]
        )
        results.append(result)

    report = {
        "schema": "gravlax.benchmark.uniform-io-query-family.v1",
        "candidate": historical_path_identity(args.candidate),
        "baseline": path_identity(args.baseline),
        "archive": path_identity(args.archive),
        "plan": path_identity(args.plan),
        "harness": historical_path_identity(Path(__file__)),
        "fixtures": {
            "junction_window": args.junction_window,
            "junction_a": args.junction_a,
            "junction_b": args.junction_b,
            "dense_junction": args.dense_junction,
        },
        "schedule": "paired alternating AB/BA",
        "iterations": args.iterations,
        "warmups": args.warmups,
        "rss_runs": args.rss_runs,
        "environment": {
            "rayon_num_threads": os.environ.get("RAYON_NUM_THREADS"),
            "cpu_affinity": sorted(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else None,
        },
        "cases": results,
        "passed": all(result["gate"] for result in results),
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"passed": report["passed"], "cases": [
        {"case": result["case"], "gate": result["gate"], "time_ratio": result["time_ratio"],
         "rss_ratio": result["rss_ratio"]} for result in results
    ]}, indent=2))
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
