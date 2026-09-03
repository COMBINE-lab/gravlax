#!/usr/bin/env python3
"""Reproducible semantic/performance gate for the region/junction I/O pilot."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import statistics
import subprocess
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Case:
    name: str
    kind: str
    locus: str
    extra: tuple[str, ...]
    output_class: str
    output_format: str = "json"


def command(binary: Path, archive: Path, case: Case, uniform: bool) -> list[str]:
    args = [str(binary), "query", str(archive), case.kind, case.locus, *case.extra]
    if uniform:
        args += ["--format", case.output_format]
    else:
        args += ["--json" if case.output_format == "json" else "--tsv"]
    return args


def run(command_line: list[str], capture: bool = True) -> subprocess.CompletedProcess[bytes]:
    return subprocess.run(
        command_line,
        check=True,
        stdout=subprocess.PIPE if capture else subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )


def sample_summary(samples: list[float]) -> dict[str, float]:
    median = statistics.median(samples)
    return {
        "median": median,
        "mad": statistics.median(abs(sample - median) for sample in samples),
        "min": min(samples),
        "max": max(samples),
    }


def path_identity(path: Path) -> dict[str, str | int]:
    digest = hashlib.sha256()
    if path.is_file():
        size = path.stat().st_size
        with path.open("rb") as reader:
            while chunk := reader.read(1024 * 1024):
                digest.update(chunk)
        return {
            "path": str(path),
            "kind": "file",
            "size_bytes": size,
            "sha256": digest.hexdigest(),
        }
    if not path.is_dir():
        raise FileNotFoundError(path)
    digest.update(b"gravlax-benchmark-tree-v1\0")
    size = 0
    files = 0
    for entry in sorted(candidate for candidate in path.rglob("*") if candidate.is_file()):
        relative = entry.relative_to(path).as_posix().encode()
        entry_size = entry.stat().st_size
        digest.update(len(relative).to_bytes(8, "little"))
        digest.update(relative)
        digest.update(entry_size.to_bytes(8, "little"))
        with entry.open("rb") as reader:
            while chunk := reader.read(1024 * 1024):
                digest.update(chunk)
        size += entry_size
        files += 1
    return {
        "path": str(path),
        "kind": "directory-tree-v1",
        "files": files,
        "size_bytes": size,
        "sha256": digest.hexdigest(),
    }


def historical_path_identity(path: Path) -> dict[str, object]:
    """Record a mutable artifact's measured location without claiming it is current."""
    result: dict[str, object] = path_identity(path)
    result["path_at_measurement"] = result.pop("path")
    result["path_semantics"] = (
        "Historical path at measurement time; size_bytes and sha256 are authoritative."
    )
    return result


def generated_fixture_identity(path: Path, argument_token: str) -> dict[str, object]:
    """Identify generated fixture bytes without serializing their temporary path."""
    result: dict[str, object] = path_identity(path)
    result.pop("path", None)
    result["kind"] = "generated-file"
    result["argument_token"] = argument_token
    return result


def affinity_ranges() -> str | None:
    if not hasattr(os, "sched_getaffinity"):
        return None
    cpus = sorted(os.sched_getaffinity(0))
    ranges: list[str] = []
    start = previous = cpus[0]
    for cpu in cpus[1:]:
        if cpu == previous + 1:
            previous = cpu
            continue
        ranges.append(str(start) if start == previous else f"{start}-{previous}")
        start = previous = cpu
    ranges.append(str(start) if start == previous else f"{start}-{previous}")
    return ",".join(ranges)


def paired_times(
    command_a: list[str], command_b: list[str], iterations: int, warmups: int
) -> tuple[list[float], list[float]]:
    """Measure a pair with alternating AB/BA order to limit warm-cache bias."""
    samples_a: list[float] = []
    samples_b: list[float] = []
    for index in range(warmups + iterations):
        order = (("a", command_a), ("b", command_b))
        if index % 2:
            order = tuple(reversed(order))
        for label, command_line in order:
            started = time.perf_counter()
            run(command_line, capture=False)
            elapsed = time.perf_counter() - started
            if index >= warmups:
                (samples_a if label == "a" else samples_b).append(elapsed)
    return samples_a, samples_b


def max_rss_kib(command_line: list[str], result: Path) -> int:
    subprocess.run(
        ["/usr/bin/time", "-f", "%M", "-o", str(result), *command_line],
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    return int(result.read_text().strip())


def paired_rss(
    command_a: list[str], command_b: list[str], runs: int, root: Path, prefix: str
) -> tuple[list[float], list[float]]:
    samples_a: list[float] = []
    samples_b: list[float] = []
    result = root / f"{prefix}-rss.txt"
    for index in range(runs):
        order = (("a", command_a), ("b", command_b))
        if index % 2:
            order = tuple(reversed(order))
        for label, command_line in order:
            rss = float(max_rss_kib(command_line, result))
            (samples_a if label == "a" else samples_b).append(rss)
    return samples_a, samples_b


def uniform_table(value: dict) -> tuple[dict, list[dict], dict]:
    table = value["data"]["tables"][0]
    fields = [field["name"] for field in table["schema"]["fields"]]
    rows = [dict(zip(fields, row, strict=True)) for row in table["rows"]]
    ordered_by = table["schema"].get("semantics", {}).get("ordered_by")
    if ordered_by not in (None, []):
        raise AssertionError("pilot count-table presentation order must remain unspecified")
    return value["data"]["summary"], rows, table["selection"]


def comparable_legacy(value: dict) -> tuple[dict, list[dict]]:
    if "start" in value:
        names = ("chrom", "start", "end", "molecules", "umis", "cells", "chunks_decoded")
        summary = {name: value[name] for name in names}
    else:
        names = ("chrom", "donor", "acceptor", "umis", "cells")
        summary = {name: value[name] for name in names}
        summary["archive_supporting_children"] = value["supporting_children"]
        summary["archive_posting_chunks"] = value["posting_chunks"]
    if "cell_rows" in value:
        rows = [
            {
                "aggregation": "cell",
                "entity": row["barcode"],
                "umis": row["umis"],
                "cells": None,
                "selected_cells": None,
            }
            for row in value["cell_rows"]
        ]
    elif "group_rows" in value:
        rows = [
            {
                "aggregation": "group",
                "entity": row["group"],
                "umis": row["umis"],
                "cells": row["cells"],
                "selected_cells": row["selected_cells"],
            }
            for row in value["group_rows"]
        ]
    else:
        rows = [
            {
                "aggregation": "bulk",
                "entity": "bulk",
                "umis": value["bulk"]["umis"],
                "cells": value["bulk"]["cells"],
                "selected_cells": value["scope"]["selected_cells"],
            }
        ]
    return summary, rows


def comparable_uniform(value: dict) -> tuple[dict, list[dict], dict]:
    summary, rows, selection = uniform_table(value)
    if "archive_supporting_children" in summary:
        names = (
            "chrom",
            "donor",
            "acceptor",
            "archive_supporting_children",
            "archive_posting_chunks",
            "umis",
            "cells",
        )
    else:
        names = ("chrom", "start", "end", "molecules", "umis", "cells", "chunks_decoded")
    return {name: summary[name] for name in names}, rows, selection


def with_all_cell_rows(case: Case) -> Case:
    extra = list(case.extra)
    top = "999999999" if case.kind == "region" else "0"
    if "--top" in extra:
        extra[extra.index("--top") + 1] = top
    else:
        extra += ["--top", top]
    return Case(case.name, case.kind, case.locus, tuple(extra), case.output_class, "json")


def verify_json_semantics(
    binary: Path, archive: Path, case: Case
) -> tuple[bytes, bytes, dict]:
    legacy_bytes = run(command(binary, archive, case, False)).stdout
    uniform_bytes = run(command(binary, archive, case, True)).stdout
    legacy = json.loads(legacy_bytes)
    uniform = json.loads(uniform_bytes)
    legacy_summary, legacy_rows = comparable_legacy(legacy)
    uniform_summary, uniform_rows, selection = comparable_uniform(uniform)
    if legacy_summary != uniform_summary:
        raise AssertionError(f"{case.name}: scientific summary mismatch")

    aggregation = uniform_rows[0]["aggregation"] if uniform_rows else "cell"
    if aggregation == "cell":
        full_value = json.loads(
            run(command(binary, archive, with_all_cell_rows(case), False)).stdout
        )
        _, full_rows = comparable_legacy(full_value)
        expected_rows = sorted(
            full_rows, key=lambda row: (-row["umis"], row["entity"])
        )[: selection["emitted_rows"]]
        if sorted(uniform_rows, key=repr) != sorted(expected_rows, key=repr):
            raise AssertionError(
                f"{case.name}: top-N subset does not follow (-umis, barcode)"
            )
        available = len(full_rows)
    else:
        if sorted(legacy_rows, key=repr) != sorted(uniform_rows, key=repr):
            raise AssertionError(f"{case.name}: complete row set mismatch")
        available = len(legacy_rows)

    expected_selection = {
        "available_rows": available,
        "emitted_rows": len(uniform_rows),
        "truncated": available > len(uniform_rows),
    }
    if selection != expected_selection:
        raise AssertionError(
            f"{case.name}: selection {selection!r} != {expected_selection!r}"
        )
    return legacy_bytes, uniform_bytes, selection


def parse_legacy_tsv(encoded: bytes) -> list[tuple[str, int]]:
    lines = encoded.decode().splitlines()
    if not lines or lines[0] != "barcode\tumis":
        raise AssertionError("legacy TSV header changed")
    return [(fields[0], int(fields[1])) for fields in (line.split("\t") for line in lines[1:])]


def parse_uniform_tsv(encoded: bytes) -> tuple[list[tuple[str, int]], dict]:
    lines = encoded.decode().splitlines()
    metadata: dict[str, int | bool] = {}
    for line in lines:
        if line.startswith("# available_rows="):
            metadata["available_rows"] = int(line.split("=", 1)[1])
        elif line.startswith("# emitted_rows="):
            metadata["emitted_rows"] = int(line.split("=", 1)[1])
        elif line.startswith("# truncated="):
            metadata["truncated"] = line.split("=", 1)[1] == "true"
    try:
        header_index = lines.index("aggregation\tentity\tumis\tcells\tselected_cells")
    except ValueError as error:
        raise AssertionError("uniform TSV header missing") from error
    rows = []
    for line in lines[header_index + 1 :]:
        if line == "# end_table":
            break
        fields = line.split("\t")
        if len(fields) != 5 or fields[0] != "cell":
            raise AssertionError(f"unexpected uniform TSV row: {line!r}")
        rows.append((fields[1], int(fields[2])))
    return rows, metadata


def verify_tsv_semantics(
    binary: Path, archive: Path, case: Case
) -> tuple[bytes, bytes, dict]:
    legacy_bytes = run(command(binary, archive, case, False)).stdout
    uniform_bytes = run(command(binary, archive, case, True)).stdout
    legacy_rows = parse_legacy_tsv(legacy_bytes)
    uniform_rows, selection = parse_uniform_tsv(uniform_bytes)
    if sorted(legacy_rows) != sorted(uniform_rows):
        raise AssertionError(f"{case.name}: TSV projected row set mismatch")
    expected_selection = {
        "available_rows": len(legacy_rows),
        "emitted_rows": len(uniform_rows),
        "truncated": False,
    }
    if selection != expected_selection:
        raise AssertionError(
            f"{case.name}: TSV selection {selection!r} != {expected_selection!r}"
        )
    return legacy_bytes, uniform_bytes, selection


def timed_file_publication(
    legacy_command: list[str],
    uniform_command: list[str],
    iterations: int,
    warmups: int,
    root: Path,
) -> tuple[list[float], list[float]]:
    legacy_samples: list[float] = []
    uniform_samples: list[float] = []
    for index in range(warmups + iterations):
        order = ("legacy", "uniform") if index % 2 == 0 else ("uniform", "legacy")
        for label in order:
            output = root / f"file-{index}-{label}.json"
            started = time.perf_counter()
            if label == "legacy":
                with output.open("xb") as writer:
                    subprocess.run(
                        legacy_command,
                        check=True,
                        stdout=writer,
                        stderr=subprocess.DEVNULL,
                    )
            else:
                run([*uniform_command, "--output", str(output)], capture=False)
            elapsed = time.perf_counter() - started
            if index >= warmups:
                (legacy_samples if label == "legacy" else uniform_samples).append(elapsed)
    return legacy_samples, uniform_samples


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--aie", type=Path, required=True, help="candidate release binary")
    parser.add_argument("--baseline-aie", type=Path, help="pre-migration release binary")
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--region", default="chr16:89550000-89575000")
    parser.add_argument("--junction", default="chr16:89562391-89562883")
    parser.add_argument("--iterations", type=int, default=15)
    parser.add_argument("--warmups", type=int, default=3)
    parser.add_argument("--rss-runs", type=int, default=9)
    parser.add_argument(
        "--scratch-root",
        type=Path,
        help="parent for temporary fixtures (default: target/uniform-io-pilot-scratch)",
    )
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()

    candidate = args.aie.resolve()
    baseline = args.baseline_aie.resolve() if args.baseline_aie else None
    archive = args.archive.resolve()
    if args.iterations < 1 or args.warmups < 0 or args.rss_runs < 1:
        parser.error("--iterations and --rss-runs must be positive; --warmups cannot be negative")
    scratch_root = (args.scratch_root or Path("target/uniform-io-pilot-scratch")).resolve()
    scratch_root.mkdir(parents=True, exist_ok=True)
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="run-", dir=scratch_root) as temporary:
        root = Path(temporary)
        seed_case = Case("seed", "junction", args.junction, ("--top", "0"), "output")
        seed = json.loads(run(command(candidate, archive, seed_case, False)).stdout)
        barcodes = [row["barcode"] for row in seed["cell_rows"]]
        if not barcodes:
            raise SystemExit("junction fixture has no cell rows")
        identity_case = Case("identity", "junction", args.junction, ("--top", "1"), "scan")
        identity_result = json.loads(run(command(candidate, archive, identity_case, True)).stdout)
        declared_archive_identity = identity_result["provenance"]["archives"][0]
        cells = root / "cells.txt"
        groups = root / "groups.tsv"
        cells.write_text("".join(f"{barcode}\n" for barcode in barcodes[:32]))
        groups.write_text(
            "".join(f"{barcode}\tgroup-{index % 4}\n" for index, barcode in enumerate(barcodes))
        )
        cases = [
            Case("junction-cell-small-json", "junction", args.junction, ("--top", "20"), "scan"),
            Case("junction-cell-large-json", "junction", args.junction, ("--top", "0"), "output"),
            Case(
                "junction-cell-large-tsv",
                "junction",
                args.junction,
                ("--top", "0"),
                "output",
                "tsv",
            ),
            Case(
                "junction-scoped-cell-json",
                "junction",
                args.junction,
                ("--cells", str(cells), "--agg", "cell", "--top", "0"),
                "scan",
            ),
            Case(
                "junction-group-json",
                "junction",
                args.junction,
                ("--groups", str(groups), "--agg", "group", "--top", "0"),
                "scan",
            ),
            Case(
                "junction-bulk-json",
                "junction",
                args.junction,
                ("--agg", "bulk", "--top", "0"),
                "scan",
            ),
            Case("region-cell-small-json", "region", args.region, ("--top", "20"), "scan"),
            Case(
                "region-cell-large-json",
                "region",
                args.region,
                ("--top", "999999999"),
                "output",
            ),
            Case(
                "region-bulk-json",
                "region",
                args.region,
                ("--agg", "bulk", "--top", "0"),
                "scan",
            ),
        ]

        records = []
        legacy_hashes = {}
        for case in cases:
            verifier = verify_json_semantics if case.output_format == "json" else verify_tsv_semantics
            legacy_bytes, uniform_bytes, selection = verifier(candidate, archive, case)
            legacy_command = command(candidate, archive, case, False)
            uniform_command = command(candidate, archive, case, True)
            legacy_hashes[case.name] = hashlib.sha256(legacy_bytes).hexdigest()

            legacy_times, uniform_times = paired_times(
                legacy_command, uniform_command, args.iterations, args.warmups
            )
            legacy_rss, uniform_rss = paired_rss(
                legacy_command, uniform_command, args.rss_runs, root, case.name
            )
            legacy_seconds = sample_summary(legacy_times)
            uniform_seconds = sample_summary(uniform_times)
            legacy_rss_kib = sample_summary(legacy_rss)
            uniform_rss_kib = sample_summary(uniform_rss)
            time_ratio = uniform_seconds["median"] / legacy_seconds["median"]
            rss_ratio = uniform_rss_kib["median"] / legacy_rss_kib["median"]
            limit = 1.03 if case.output_class == "scan" else 1.10

            compatibility = None
            if baseline:
                baseline_command = command(baseline, archive, case, False)
                baseline_bytes = run(baseline_command).stdout
                if baseline_bytes != legacy_bytes:
                    raise AssertionError(f"{case.name}: candidate legacy stdout changed from baseline")
                baseline_times, candidate_legacy_times = paired_times(
                    baseline_command, legacy_command, args.iterations, args.warmups
                )
                baseline_rss, candidate_legacy_rss = paired_rss(
                    baseline_command,
                    legacy_command,
                    args.rss_runs,
                    root,
                    f"{case.name}-baseline",
                )
                baseline_seconds = sample_summary(baseline_times)
                candidate_legacy_seconds = sample_summary(candidate_legacy_times)
                baseline_rss_kib = sample_summary(baseline_rss)
                candidate_legacy_rss_kib = sample_summary(candidate_legacy_rss)
                compatibility = {
                    "stdout_exact": True,
                    "baseline_seconds": baseline_seconds,
                    "candidate_seconds": candidate_legacy_seconds,
                    "time_ratio": candidate_legacy_seconds["median"]
                    / baseline_seconds["median"],
                    "baseline_rss_kib": baseline_rss_kib,
                    "candidate_rss_kib": candidate_legacy_rss_kib,
                    "rss_ratio": candidate_legacy_rss_kib["median"]
                    / baseline_rss_kib["median"],
                }
                compatibility["gate"] = (
                    compatibility["time_ratio"] <= 1.03
                    and compatibility["rss_ratio"] <= 1.10
                )

            record = {
                "case": case.name,
                "format": case.output_format,
                "output_class": case.output_class,
                "rows": selection["emitted_rows"],
                "legacy_stdout_bytes": len(legacy_bytes),
                "uniform_stdout_bytes": len(uniform_bytes),
                "legacy_seconds": legacy_seconds,
                "uniform_seconds": uniform_seconds,
                "time_ratio": time_ratio,
                "legacy_rss_kib": legacy_rss_kib,
                "uniform_rss_kib": uniform_rss_kib,
                "rss_ratio": rss_ratio,
                "semantic_equivalence": True,
                "time_limit": limit,
                "time_gate": time_ratio <= limit,
                "bounded_memory_gate": rss_ratio <= 1.10,
                "legacy_compatibility": compatibility,
            }
            record["gate"] = (
                record["time_gate"]
                and record["bounded_memory_gate"]
                and (compatibility is None or compatibility["gate"])
            )
            records.append(record)

        publication_case = next(
            case for case in cases if case.name == "junction-cell-large-json"
        )
        legacy_file_times, atomic_file_times = timed_file_publication(
            command(candidate, archive, publication_case, False),
            command(candidate, archive, publication_case, True),
            args.iterations,
            args.warmups,
            root,
        )
        legacy_file_seconds = sample_summary(legacy_file_times)
        atomic_file_seconds = sample_summary(atomic_file_times)
        file_ratio = atomic_file_seconds["median"] / legacy_file_seconds["median"]

        by_name = {record["case"]: record for record in records}
        bounded_checks = []
        for prefix in ("junction", "region"):
            small = by_name[f"{prefix}-cell-small-json"]
            large = by_name[f"{prefix}-cell-large-json"]
            ratio = large["uniform_rss_kib"]["median"] / small["uniform_rss_kib"]["median"]
            bounded_checks.append(
                {
                    "comparison": f"{prefix}-large-vs-small-json",
                    "emitted_row_ratio": large["rows"] / small["rows"],
                    "rss_ratio": ratio,
                    "limit": 1.10,
                    "gate": ratio <= 1.10,
                }
            )

        normalized_invocation = [
            str(Path(__file__).resolve()),
            "--aie",
            str(candidate),
            "--archive",
            str(archive),
            "--region",
            args.region,
            "--junction",
            args.junction,
            "--iterations",
            str(args.iterations),
            "--warmups",
            str(args.warmups),
            "--rss-runs",
            str(args.rss_runs),
            "--scratch-root",
            str(scratch_root),
        ]
        if baseline:
            normalized_invocation[3:3] = ["--baseline-aie", str(baseline)]
        if args.report:
            normalized_invocation += ["--report", str(args.report.resolve())]

        report = {
            "schema": "gravlax.benchmark.uniform-io-pilot.v4",
            "harness": historical_path_identity(Path(__file__).resolve()),
            "normalized_invocation": normalized_invocation,
            "archive": {
                **path_identity(archive),
                "declared_identity": declared_archive_identity,
            },
            "candidate": historical_path_identity(candidate),
            "baseline": path_identity(baseline) if baseline else None,
            "environment": {
                "rayon_num_threads": os.environ.get("RAYON_NUM_THREADS"),
                "cpu_affinity": affinity_ranges(),
                "logical_cpus_visible": len(os.sched_getaffinity(0))
                if hasattr(os, "sched_getaffinity")
                else os.cpu_count(),
            },
            "fixtures": {
                "region_locus": args.region,
                "junction_locus": args.junction,
                "cells": {
                    **generated_fixture_identity(cells, "<generated-cells-fixture>"),
                    "construction": "first 32 rows from the full legacy junction result",
                },
                "groups": {
                    **generated_fixture_identity(groups, "<generated-groups-fixture>"),
                    "construction": "all full legacy junction rows assigned round-robin to group-0..3",
                },
                "cases": [
                    {
                        "name": case.name,
                        "kind": case.kind,
                        "locus": case.locus,
                        "arguments": [
                            "<generated-cells-fixture>"
                            if value == str(cells)
                            else "<generated-groups-fixture>"
                            if value == str(groups)
                            else value
                            for value in case.extra
                        ],
                        "format": case.output_format,
                        "output_class": case.output_class,
                    }
                    for case in cases
                ],
            },
            "schedule": "paired alternating AB/BA",
            "iterations": args.iterations,
            "warmups": args.warmups,
            "rss_runs": args.rss_runs,
            "scratch_root": str(scratch_root),
            "cases": records,
            "legacy_stdout_sha256": legacy_hashes,
            "file_publication": {
                "case": publication_case.name,
                "legacy_redirect_seconds": legacy_file_seconds,
                "uniform_atomic_flush_seconds": atomic_file_seconds,
                "time_ratio": file_ratio,
                "limit": 1.10,
                "gate": file_ratio <= 1.10,
            },
            "bounded_output_memory": bounded_checks,
            "passed": (
                all(record["gate"] for record in records)
                and file_ratio <= 1.10
                and all(check["gate"] for check in bounded_checks)
            ),
            "notes": [
                "scan cases use a <=3% end-to-end median-time gate",
                "large-output cases use a <=10% end-to-end proxy for the output-only gate; this conservatively includes archive work",
                "uniform RSS must be <=10% above paired legacy RSS, and large-row RSS is also compared with small-row RSS",
                "RSS gates use alternating repeated runs and median/MAD because individual peak-RSS readings are noisy",
                "atomic publication uses Flush durability and is timed separately from stdout to /dev/null",
                "legacy compatibility requires byte-exact stdout, <=3% median time, and <=10% RSS when --baseline-aie is supplied",
            ],
        }
        encoded = json.dumps(report, indent=2, sort_keys=True) + "\n"
        if args.report:
            args.report.write_text(encoded)
        print(encoded, end="")
        if not report["passed"]:
            raise SystemExit(1)


if __name__ == "__main__":
    main()
