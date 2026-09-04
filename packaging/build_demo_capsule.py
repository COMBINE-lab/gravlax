#!/usr/bin/env python3
"""Build the path-independent data portion of a Gravlax demo capsule."""

from __future__ import annotations

import argparse
from contextlib import contextmanager
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import subprocess
import tempfile
from typing import Any, Iterator


BUILD_SCHEMA = "gravlax.demo-capsule-build.v1"
RECORD_SCHEMA = "gravlax.demo-capsule-build-record.v1"
HEX64 = re.compile(r"^[0-9a-f]{64}$")
IDENTIFIER = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")
BARCODE = re.compile(r"^[ACGT]{16}$")
COLLECTION_JUNCTION = re.compile(
    r"^([A-Za-z0-9][A-Za-z0-9_.-]*):([0-9]+)-([0-9]+)$"
)
PRIVATE_TEXT = ("/scratch", "/nfshomes", "/home/", "/Users/")
ABSOLUTE_PATH = re.compile(
    r"(?:^|[\s=:(])/(?!/)[A-Za-z0-9._~-]+(?:/[^\s]*)?|(?:^|[\s=:(])[A-Za-z]:[\\/]"
)


class CapsuleError(RuntimeError):
    """A capsule input or generated artifact violated the build contract."""


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(8 << 20):
            digest.update(block)
    return digest.hexdigest()


def _safe_filename(value: str) -> str:
    if not isinstance(value, str) or not IDENTIFIER.fullmatch(value):
        raise CapsuleError(f"unsafe capsule filename: {value!r}")
    if PurePosixPath(value).name != value:
        raise CapsuleError(f"capsule filename must be a basename: {value!r}")
    return value


def _safe_identifier(value: Any, label: str) -> str:
    if not isinstance(value, str) or not IDENTIFIER.fullmatch(value):
        raise CapsuleError(f"{label} must contain only letters, digits, '.', '_', and '-'")
    return value


def _public_text(value: Any, label: str) -> str:
    if not isinstance(value, str) or not value.strip() or "\t" in value or "\n" in value:
        raise CapsuleError(f"{label} must be nonempty, single-line text")
    if any(marker.lower() in value.lower() for marker in PRIVATE_TEXT) or ABSOLUTE_PATH.search(
        value
    ):
        raise CapsuleError(f"{label} contains a private filesystem locator")
    return value


def _https(value: Any, label: str) -> str:
    value = _public_text(value, label)
    if not value.startswith("https://"):
        raise CapsuleError(f"{label} must be an HTTPS URL")
    return value


def _expect_mapping(value: Any, label: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise CapsuleError(f"{label} must be an object")
    return value


def _expect_list(value: Any, label: str) -> list[Any]:
    if not isinstance(value, list) or not value:
        raise CapsuleError(f"{label} must be a nonempty array")
    return value


def _validate_input_declaration(value: Any, label: str) -> dict[str, Any]:
    declaration = _expect_mapping(value, label)
    relative = declaration.get("path")
    if not isinstance(relative, str) or not relative:
        raise CapsuleError(f"{label}.path must be nonempty")
    relative_path = Path(relative)
    if relative_path.is_absolute() or ".." in relative_path.parts:
        raise CapsuleError(f"{label}.path must stay below --source-root")
    expected = declaration.get("sha256")
    if not isinstance(expected, str) or not HEX64.fullmatch(expected):
        raise CapsuleError(f"{label}.sha256 must be 64 lowercase hexadecimal characters")
    expected_bytes = declaration.get("bytes")
    if type(expected_bytes) is not int or expected_bytes < 0:
        raise CapsuleError(f"{label}.bytes must be a nonnegative integer")
    provenance = _expect_mapping(declaration.get("provenance"), f"{label}.provenance")
    if not provenance:
        raise CapsuleError(f"{label}.provenance must not be empty")
    _public_text(provenance.get("accession"), f"{label}.provenance.accession")
    _https(provenance.get("source_url"), f"{label}.provenance.source_url")
    _assert_public_strings(provenance, f"{label}.provenance")
    return declaration


def _resolve_input(source_root: Path, declaration: dict[str, Any], label: str) -> Path:
    declaration = _validate_input_declaration(declaration, label)
    relative = declaration["path"]
    candidate = (source_root / relative).resolve(strict=True)
    if not candidate.is_relative_to(source_root):
        raise CapsuleError(f"{label}.path escapes --source-root")
    if not candidate.is_file():
        raise CapsuleError(f"{label}.path is not a regular file")
    expected = declaration["sha256"]
    observed = _sha256(candidate)
    if observed != expected:
        raise CapsuleError(f"SHA-256 mismatch for {label}: {observed} != {expected}")
    expected_bytes = declaration["bytes"]
    if expected_bytes != candidate.stat().st_size:
        raise CapsuleError(f"byte-size mismatch for {label}")
    return candidate


def _assert_public_strings(value: Any, label: str) -> None:
    if isinstance(value, dict):
        for key, item in value.items():
            _assert_public_strings(key, label)
            _assert_public_strings(item, label)
    elif isinstance(value, list):
        for item in value:
            _assert_public_strings(item, label)
    elif isinstance(value, str):
        if any(
            marker.lower() in value.lower() for marker in PRIVATE_TEXT
        ) or ABSOLUTE_PATH.search(value):
            raise CapsuleError(f"{label} contains a private filesystem locator")


def _public_input_record(declaration: dict[str, Any]) -> dict[str, Any]:
    return {
        "bytes": declaration["bytes"],
        "sha256": declaration["sha256"],
        "provenance": declaration["provenance"],
        "source_url_relationship": (
            "provenance-only; the exact local bytes are identified by sha256 and may be a "
            "derived or transformed object rather than bytes served directly by source_url"
        ),
    }


def _write_bytes_exclusive(path: Path, value: bytes) -> None:
    with path.open("xb") as handle:
        handle.write(value)


def _write_json_exclusive(path: Path, value: Any) -> None:
    _write_bytes_exclusive(
        path,
        (json.dumps(value, indent=2, sort_keys=True, ensure_ascii=True) + "\n").encode(),
    )


@contextmanager
def _new_output_directory(destination: Path) -> Iterator[Path]:
    destination = destination.resolve()
    if destination == Path(destination.anchor) or destination == Path.cwd().resolve():
        raise CapsuleError("refusing a broad output directory")
    destination.parent.mkdir(parents=True, exist_ok=True)
    try:
        destination.mkdir()
    except FileExistsError as error:
        raise FileExistsError(f"refusing to replace capsule directory: {destination}") from error
    try:
        yield destination
    except BaseException:
        shutil.rmtree(destination)
        raise


def _run(
    command: list[str],
    *,
    stdout: Any = subprocess.PIPE,
    cwd: Path | None = None,
) -> subprocess.CompletedProcess[bytes]:
    try:
        return subprocess.run(
            command,
            check=True,
            stdout=stdout,
            stderr=subprocess.PIPE,
            cwd=cwd,
        )
    except OSError as error:
        raise CapsuleError(f"cannot execute command: {command[0]}") from error
    except subprocess.CalledProcessError as error:
        stderr = error.stderr.decode("utf-8", "replace")[-4000:] if error.stderr else ""
        raise CapsuleError(f"command failed ({command[0]}): {stderr}") from error


def _tool_version(executable: Path, arguments: list[str]) -> str:
    result = _run([str(executable), *arguments])
    return result.stdout.decode("utf-8", "strict").strip()


def _stage_input(source: Path, destination: Path, expected_sha256: str) -> None:
    """Copy one verified input without sharing or mutating the source inode."""
    if destination.exists():
        raise FileExistsError(f"refusing to replace staged input: {destination}")
    with source.open("rb") as reader, destination.open("xb") as writer:
        shutil.copyfileobj(reader, writer, 8 << 20)
    if _sha256(destination) != expected_sha256:
        raise CapsuleError(f"staged input changed while copying: {source.name}")


def _stage_working_link(source: Path, destination: Path, expected_sha256: str) -> None:
    """Link a private per-build copy for sequential use within that same build."""
    if destination.exists():
        raise FileExistsError(f"refusing to replace staged input: {destination}")
    try:
        os.link(source, destination)
    except OSError:
        with source.open("rb") as reader, destination.open("xb") as writer:
            shutil.copyfileobj(reader, writer, 8 << 20)
    if _sha256(destination) != expected_sha256:
        raise CapsuleError(f"staged input changed while linking: {source.name}")


def _stage_executable(source: Path, destination: Path, expected_sha256: str) -> None:
    _stage_input(source, destination, expected_sha256)
    if not os.access(destination, os.X_OK):
        destination.chmod(destination.stat().st_mode | 0o100)


def _windows(spec: dict[str, Any]) -> list[tuple[str, int, int]]:
    parsed: list[tuple[str, int, int]] = []
    for index, raw in enumerate(_expect_list(spec.get("windows"), "windows")):
        item = _expect_mapping(raw, f"windows[{index}]")
        chrom = _safe_identifier(item.get("chrom"), f"windows[{index}].chrom")
        start, end = item.get("start"), item.get("end")
        if not isinstance(start, int) or not isinstance(end, int) or start < 0 or end <= start:
            raise CapsuleError(f"windows[{index}] must be a nonempty 0-based half-open interval")
        parsed.append((chrom, start, end))
    parsed.sort()
    for previous, current in zip(parsed, parsed[1:]):
        if previous[0] == current[0] and current[1] < previous[2]:
            raise CapsuleError("windows must not overlap; samtools -L could duplicate records")
    return parsed


GENE_ID = re.compile(r'(?:^|;\s*)gene_id "([^"]+)"')


def _subset_gtf(source: Path, windows: list[tuple[str, int, int]], destination: Path) -> int:
    by_chrom: dict[str, list[tuple[int, int]]] = {}
    for chrom, start, end in windows:
        by_chrom.setdefault(chrom, []).append((start, end))
    genes: set[str] = set()
    headers: list[str] = []
    with source.open("rt", encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, 1):
            if line.startswith("#"):
                headers.append(line)
                continue
            fields = line.rstrip("\n").split("\t", 8)
            if len(fields) != 9 or fields[0] not in by_chrom:
                continue
            try:
                one_based_start, inclusive_end = int(fields[3]), int(fields[4])
            except ValueError as error:
                raise CapsuleError(f"invalid GTF coordinates at {source}:{line_no}") from error
            overlaps = any(
                one_based_start <= window_end and inclusive_end >= window_start + 1
                for window_start, window_end in by_chrom[fields[0]]
            )
            if overlaps and (match := GENE_ID.search(fields[8])):
                genes.add(match.group(1))
    if not genes:
        raise CapsuleError(f"annotation {source.name} has no genes overlapping the demo windows")
    with destination.open("x", encoding="utf-8", newline="\n") as sink:
        sink.write("# gravlax_demo_subset: complete records for genes overlapping ")
        sink.write(",".join(f"{chrom}:{start}-{end}" for chrom, start, end in windows))
        sink.write("\n")
        sink.writelines(headers)
        with source.open("rt", encoding="utf-8") as handle:
            for line in handle:
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t", 8)
                if len(fields) == 9 and (match := GENE_ID.search(fields[8])):
                    if match.group(1) in genes:
                        sink.write(line if line.endswith("\n") else line + "\n")
    return len(genes)


def _sanitized_header(
    original: str,
    sample: str,
    accessions: list[str],
    windows: list[tuple[str, int, int]],
) -> bytes:
    header_lines = original.splitlines()
    hd_lines = [line for line in header_lines if line.startswith("@HD\t")]
    if len(hd_lines) != 1 or "SO:coordinate" not in hd_lines[0].split("\t")[1:]:
        raise CapsuleError("sorted BAM header must declare @HD SO:coordinate")
    retained = [hd_lines[0], *[line for line in header_lines if line.startswith("@SQ\t")]]
    if not any(line.startswith("@SQ\t") for line in retained):
        raise CapsuleError("source BAM header has no @SQ dictionary")
    for accession in accessions:
        _safe_identifier(accession, "donor accessions")
    region_text = ",".join(f"{chrom}:{start}-{end}" for chrom, start, end in windows)
    retained.append(
        "@CO\tPublic demo subset "
        f"{sample}; accessions={','.join(accessions)}; 0-based-half-open-windows={region_text}"
    )
    header = ("\n".join(retained) + "\n").encode("utf-8")
    _assert_public_strings(header.decode(), "sanitized BAM header")
    return header


def _read_group_map(path: Path, label: str) -> list[tuple[str, str]]:
    rows: list[tuple[str, str]] = []
    seen: set[str] = set()
    for line_no, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        fields = raw.split("\t")
        if len(fields) != 2 or not BARCODE.fullmatch(fields[0]):
            raise CapsuleError(f"{label}:{line_no} must be barcode<TAB>group")
        group = _safe_identifier(fields[1], f"{label}:{line_no} group")
        if fields[0] in seen:
            raise CapsuleError(f"{label} repeats barcode {fields[0]}")
        seen.add(fields[0])
        rows.append((fields[0], group))
    if not rows:
        raise CapsuleError(f"{label} is empty")
    return sorted(rows)


def _validate_star_log(
    path: Path, alignment: dict[str, Any], label: str
) -> dict[str, str]:
    try:
        contents = path.read_text(encoding="utf-8")
    except UnicodeDecodeError as error:
        raise CapsuleError(f"{label} is not a UTF-8 STAR Log.out") from error
    program = _safe_identifier(alignment.get("program"), "alignment.program")
    version = _public_text(alignment.get("version"), "alignment.version")
    command_line = _public_text(alignment.get("command_line"), "alignment.command_line")
    if program != "STAR" or not re.search(
        rf"(?m)^STAR version={re.escape(version)}\s*$", contents
    ):
        raise CapsuleError(f"{label} does not confirm STAR {version}")
    if "--twopassMode Basic" not in command_line or not re.search(
        r"(?m)^twopassMode\s+Basic(?:\s|$)", contents
    ):
        raise CapsuleError(f"{label} does not confirm STAR --twopassMode Basic")
    if "##### Final effective command line:" not in contents:
        raise CapsuleError(f"{label} lacks STAR's final effective command line")
    annotation_parameters = {
        "sjdbGTFfile": "-",
        "sjdbFileChrStartEnd": "-",
    }
    for parameter, expected in annotation_parameters.items():
        if not re.search(
            rf"(?m)^{re.escape(parameter)}\s+{re.escape(expected)}(?:\s|$)", contents
        ):
            raise CapsuleError(
                f"{label} does not confirm annotation-free alignment: "
                f"expected {parameter} {expected}"
            )
    return {
        "status": "validated-from-exact-Log.out-bytes",
        "program": program,
        "version": version,
        "twopass_mode": "Basic",
        "alignment_annotation_status": "absent-in-exact-Log.out",
        "sjdb_gtf_file": "-",
        "sjdb_file_chr_start_end": "-",
        "normalized_command_status": "caller-declared-public-summary",
    }


def _asset(path: Path, archive_root: str | None = None) -> dict[str, Any]:
    result: dict[str, Any] = {
        "filename": path.name,
        "sha256": _sha256(path),
        "bytes": path.stat().st_size,
    }
    if archive_root is not None:
        result["archive_root"] = archive_root
    return result


def _verified_archive_inspection(
    inspected: Any, label: str
) -> tuple[str, dict[str, Any], int]:
    try:
        identity = inspected["native_identity"]
        verification = inspected["verification"]
        molecular = inspected["molecular_evidence"]
        provenance = molecular["alignment_provenance"]
        alignment = provenance["alignment"]
        alignment_log = alignment["alignment_log"]
        catalogue = alignment["junction_catalogue"]
        tail = molecular["terminal_tail"]
        root = f"{identity['scheme']}:{identity['blake3']}"
    except (KeyError, TypeError) as error:
        raise CapsuleError(f"{label} has an incomplete archive inspection") from error
    if not re.fullmatch(r"aie-directory-root-v2:[0-9a-f]{64}", root):
        raise CapsuleError(f"{label} does not have a v2 directory root")
    if verification.get("directory_and_root") is not True or verification.get(
        "all_payloads"
    ) is not True:
        raise CapsuleError(f"{label} did not receive complete content verification")
    if molecular.get("schema") != "gravlax.molecular-evidence.v2":
        raise CapsuleError(f"{label} lacks logical molecular-evidence v2")
    if molecular.get("alignment_provenance_status") != "available" or (
        provenance.get("schema") != "gravlax.alignment-provenance.v1"
        or provenance.get("molecular_evidence_schema")
        != "gravlax.molecular-evidence.v2"
    ):
        raise CapsuleError(f"{label} lacks root-bound alignment provenance")
    if alignment.get("junction_discovery") != "per-library-two-pass" or (
        catalogue.get("role") != "per-library-pass1"
        or catalogue.get("section") != "alignment.junction-catalogue"
    ):
        raise CapsuleError(f"{label} lacks its per-library two-pass catalogue binding")
    if alignment.get("programs") != []:
        raise CapsuleError(
            f"{label} must not present a normalized command as verified BAM-header metadata"
        )
    if alignment_log.get("locator") != "STAR-Log.out" or not isinstance(
        alignment_log.get("identity"), dict
    ):
        raise CapsuleError(f"{label} lacks its path-independent STAR Log.out identity")
    if molecular.get("genome_reference_binding_status") != "available":
        raise CapsuleError(f"{label} lacks a genome reference binding")
    if molecular.get("terminal_tail_status") != "available":
        raise CapsuleError(f"{label} lacks terminal-tail capability")
    events = tail.get("events")
    if type(events) is not int or events < 0:
        raise CapsuleError(f"{label} has an invalid terminal-tail event count")
    _assert_public_strings(molecular, f"{label} provenance")
    return root, molecular, events


def _inspect_archive(aie: Path, archive: Path) -> tuple[str, int]:
    result = _run([str(aie), "inspect-archive", str(archive), "--verify-content", "--json"])
    try:
        inspected = json.loads(result.stdout)
    except json.JSONDecodeError as error:
        raise CapsuleError("aie emitted an invalid archive inspection") from error
    root, _, tail_events = _verified_archive_inspection(
        inspected, "generated archive"
    )
    return root, tail_events


def _stage_donor(
    *,
    aie: Path,
    samtools: Path,
    source_root: Path,
    output: Path,
    work: Path,
    windows_bed: Path,
    windows: list[tuple[str, int, int]],
    genome: Path,
    genome_sha256: str,
    alignment: dict[str, Any],
    donor_spec: dict[str, Any],
) -> tuple[str, dict[str, Any], list[tuple[str, str]], int, dict[str, Any]]:
    sample = _safe_identifier(donor_spec.get("sample"), "donor.sample")
    donor = _safe_identifier(donor_spec.get("donor"), f"{sample}.donor")
    resource = _safe_identifier(donor_spec.get("archive_resource"), f"{sample}.archive_resource")
    archive_filename = _safe_filename(donor_spec.get("archive_filename"))
    if not archive_filename.endswith(".aie"):
        raise CapsuleError(f"{sample}.archive_filename must end in .aie")
    accessions = donor_spec.get("accessions")
    if not isinstance(accessions, list) or not accessions:
        raise CapsuleError(f"{sample}.accessions must be a nonempty array")
    accessions = [_safe_identifier(value, f"{sample}.accessions") for value in accessions]

    declarations = {}
    paths = {}
    for key in ("bam", "whitelist", "groups", "junction_catalogue", "alignment_log"):
        declaration = _expect_mapping(donor_spec.get(key), f"{sample}.{key}")
        declarations[key] = declaration
        paths[key] = _resolve_input(source_root, declaration, f"{sample}.{key}")
    alignment_validation = _validate_star_log(
        paths["alignment_log"], alignment, f"{sample}.alignment_log"
    )

    donor_work = work / sample
    donor_work.mkdir()
    staged_inputs = {
        "whitelist": donor_work / "barcodes.txt",
        "groups": donor_work / "groups.tsv",
        "junction_catalogue": donor_work / "STAR-pass1-SJ.out.tab",
        "alignment_log": donor_work / "STAR-Log.out",
        "genome": donor_work
        / ("reference.fa.gz" if genome.name.endswith(".gz") else "reference.fa"),
    }
    for key in ("whitelist", "groups", "junction_catalogue", "alignment_log"):
        _stage_input(paths[key], staged_inputs[key], declarations[key]["sha256"])
    _stage_working_link(genome, staged_inputs["genome"], genome_sha256)
    group_rows = _read_group_map(staged_inputs["groups"], f"{sample}.groups")
    filtered = donor_work / "filtered.bam"
    _run(
        [
            str(samtools),
            "view",
            "--no-PG",
            "-b",
            "-L",
            str(windows_bed),
            "-o",
            str(filtered),
            str(paths["bam"]),
        ]
    )
    sorted_bam = donor_work / "sorted.bam"
    _run(
        [
            str(samtools),
            "sort",
            "--no-PG",
            "-o",
            str(sorted_bam),
            str(filtered),
        ]
    )
    original_header = _run(
        [str(samtools), "view", "--no-PG", "-H", str(sorted_bam)]
    ).stdout.decode("utf-8", "strict")
    header_path = donor_work / "public.header.sam"
    _write_bytes_exclusive(
        header_path,
        _sanitized_header(original_header, sample, accessions, windows),
    )
    sanitized = donor_work / "sanitized.bam"
    with sanitized.open("xb") as sink:
        _run(
            [str(samtools), "reheader", "--no-PG", str(header_path), str(sorted_bam)],
            stdout=sink,
        )
    _run([str(samtools), "quickcheck", "-v", str(sanitized)])
    records = int(
        _run([str(samtools), "view", "-c", str(sanitized)])
        .stdout.decode("ascii", "strict")
        .strip()
    )
    if records == 0:
        raise CapsuleError(f"{sample} has no alignments in the configured windows")

    staged_archive = donor_work / archive_filename
    command = [
        str(aie),
        "ingest-archive",
        sanitized.name,
        "--whitelist",
        staged_inputs["whitelist"].name,
        "--out",
        staged_archive.name,
        "--genome",
        staged_inputs["genome"].name,
        "--terminal-tails",
        "--junction-discovery",
        "per-library-two-pass",
        "--junction-catalogue",
        staged_inputs["junction_catalogue"].name,
        "--alignment-log",
        staged_inputs["alignment_log"].name,
        "--alignment-index-identity",
        _public_text(alignment.get("index_identity"), "alignment.index_identity"),
        "--alignment-chemistry",
        _public_text(alignment.get("chemistry"), "alignment.chemistry"),
    ]
    _run(command, cwd=donor_work)
    if _sha256(paths["bam"]) != declarations["bam"]["sha256"]:
        raise CapsuleError(f"{sample}.bam changed while filtering")
    for key in ("whitelist", "groups", "junction_catalogue", "alignment_log"):
        if _sha256(staged_inputs[key]) != declarations[key]["sha256"]:
            raise CapsuleError(f"{sample}.{key} changed while ingesting")
    if _sha256(staged_inputs["genome"]) != genome_sha256:
        raise CapsuleError("genome changed while ingesting")
    archive = output / archive_filename
    if archive.exists():
        raise FileExistsError(f"refusing to replace capsule file: {archive}")
    staged_archive.replace(archive)
    root, tail_events = _inspect_archive(aie, archive)
    provenance = {
        "sample": sample,
        "donor": donor,
        "accessions": accessions,
        "selected_alignment_records": records,
        "sanitized_bam": {
            "bytes": sanitized.stat().st_size,
            "sha256": _sha256(sanitized),
            "header_policy": (
                "coordinate-sort; retain-HD-SQ; remove-PG-RG-CO; add-public-scope-CO"
            ),
        },
        "inputs": {
            key: _public_input_record(declarations[key])
            for key in sorted(declarations)
        },
        "archive_root": root,
        "terminal_tail_events": tail_events,
        "alignment_log_validation": alignment_validation,
    }
    return resource, _asset(archive, root), group_rows, tail_events, provenance


ANNOTATION_STORY_FIELDS = {
    "archive",
    "annotation_before",
    "annotation_after",
    "assembly",
    "before_label",
    "after_label",
    "max_molecule_witnesses",
    "expected_gene_id",
    "expected_min_signed_delta",
    "expected_changed_molecule_records",
}
EVENT_STORY_FIELDS = {
    "archives",
    "shape_routes",
    "allow_unstamped",
    "collection_groups",
    "design",
    "annotation",
    "assembly",
    "annotation_label",
    "annotation_digest",
    "comparison_annotation",
    "comparison_annotation_label",
    "kinds",
    "require_groups",
    "min_group_umi_classes",
    "min_donors",
    "min_samples",
    "min_umi_classes",
    "min_side_umi_classes",
    "min_support",
    "terminal_cluster_bp",
    "max_terminal_events",
    "novel_only",
    "solo_strand",
    "max_candidates",
    "max_candidates_considered",
    "max_routed_entries",
    "max_exact_match_attempts",
    "max_annotation_comparisons",
    "expected_entity_id",
    "expected_min_exact_umi_classes",
    "expected_min_exact_donors",
    "expected_gap_primary_class",
    "expected_annotation_incompatible",
    "expected_rank",
    "expected_comparison_compatible_transcripts",
}
DRILLDOWN_STORY_FIELDS = {
    "archives",
    "archive_sample",
    "shape_routes",
    "allow_unstamped",
    "drilldown_groups",
    "junction",
    "predicates",
    "expression",
    "universe",
    "unit",
    "region_match",
    "placements",
    "allow_full_scan",
    "aggregation",
    "emit_membership",
    "max_memberships",
    "max_pattern_rows",
    "max_chunks",
    "max_evidence_records",
    "max_terminal_events",
    "expected_min_selected_units",
    "required_true_predicates",
    "story_note",
}
GROUP_MAP_SCOPE_FIELDS = {
    "label_derivation",
    "included_cells",
    "event_semantics",
    "drilldown_semantics",
}


def _story_fields(
    story: dict[str, Any], label: str, required: set[str], allowed: set[str]
) -> None:
    missing = sorted(required.difference(story))
    unknown = sorted(set(story).difference(allowed))
    if missing or unknown:
        raise CapsuleError(f"{label} has missing fields {missing} or unknown fields {unknown}")


def _archive_map(story: dict[str, Any], label: str, resources: set[str]) -> dict[str, str]:
    archives = _expect_mapping(story.get("archives"), f"{label}.archives")
    if not archives:
        raise CapsuleError(f"{label}.archives must not be empty")
    for sample, resource in archives.items():
        _safe_identifier(sample, f"{label}.archives sample")
        _safe_identifier(resource, f"{label}.archives resource")
        if resource not in resources:
            raise CapsuleError(f"{label}.archives references unknown resource {resource}")
    return archives


def _story_integer(story: dict[str, Any], field: str, minimum: int, label: str) -> None:
    if field in story and (
        type(story[field]) is not int or story[field] < minimum  # bool is not an integer here
    ):
        raise CapsuleError(f"{label}.{field} must be an integer of at least {minimum}")


def _collection_junction(value: Any, label: str) -> str:
    value = _public_text(value, label)
    match = COLLECTION_JUNCTION.fullmatch(value)
    if match is None or int(match.group(3)) <= int(match.group(2)):
        raise CapsuleError(
            f"{label} must be chrom:donor-acceptor without a strand suffix"
        )
    return value


def _validate_story_references(stories: dict[str, Any], resources: set[str]) -> None:
    required = {"annotation_reinterpretation", "event_discovery", "junction_drilldown"}
    if set(stories) != required:
        raise CapsuleError(f"stories must contain exactly {sorted(required)}")
    annotation = _expect_mapping(stories["annotation_reinterpretation"], "annotation story")
    _story_fields(
        annotation,
        "annotation_reinterpretation",
        {
            "archive",
            "annotation_before",
            "annotation_after",
            "assembly",
            "before_label",
            "after_label",
            "expected_gene_id",
            "expected_min_signed_delta",
        },
        ANNOTATION_STORY_FIELDS,
    )
    for field in ("archive", "annotation_before", "annotation_after"):
        resource = _safe_identifier(
            annotation.get(field), f"annotation_reinterpretation.{field}"
        )
        if resource not in resources:
            raise CapsuleError(
                f"annotation story references unknown resource {resource!r}"
            )
    for field in ("assembly", "before_label", "after_label"):
        _public_text(annotation[field], f"annotation_reinterpretation.{field}")
    _story_integer(annotation, "max_molecule_witnesses", 1, "annotation_reinterpretation")
    _safe_identifier(
        annotation["expected_gene_id"], "annotation_reinterpretation.expected_gene_id"
    )
    _story_integer(
        annotation,
        "expected_min_signed_delta",
        1,
        "annotation_reinterpretation",
    )
    _story_integer(
        annotation,
        "expected_changed_molecule_records",
        0,
        "annotation_reinterpretation",
    )

    event = _expect_mapping(stories["event_discovery"], "event_discovery")
    drilldown = _expect_mapping(stories["junction_drilldown"], "junction_drilldown")
    _story_fields(
        event,
        "event_discovery",
        {
            "archives",
            "collection_groups",
            "design",
            "annotation",
            "assembly",
            "annotation_label",
            "expected_entity_id",
            "expected_min_exact_umi_classes",
            "expected_min_exact_donors",
            "expected_gap_primary_class",
            "expected_annotation_incompatible",
            "expected_rank",
            "comparison_annotation",
            "comparison_annotation_label",
            "expected_comparison_compatible_transcripts",
        },
        EVENT_STORY_FIELDS,
    )
    _story_fields(
        drilldown,
        "junction_drilldown",
        {
            "archives",
            "archive_sample",
            "junction",
            "predicates",
            "expression",
            "universe",
            "expected_min_selected_units",
            "required_true_predicates",
            "story_note",
            "emit_membership",
        },
        DRILLDOWN_STORY_FIELDS,
    )
    _archive_map(event, "event_discovery", resources)
    drilldown_archives = _archive_map(drilldown, "junction_drilldown", resources)
    for field in ("collection_groups", "design"):
        resource = _safe_identifier(event.get(field), f"event_discovery.{field}")
        if resource not in resources:
            raise CapsuleError(f"event_discovery.{field} references an unknown resource")
    for field in ("annotation", "comparison_annotation"):
        annotation_resource = _safe_identifier(event[field], f"event_discovery.{field}")
        if annotation_resource not in resources:
            raise CapsuleError(f"event_discovery.{field} references an unknown resource")
    if "drilldown_groups" in drilldown:
        group_resource = _safe_identifier(
            drilldown["drilldown_groups"], "junction_drilldown.drilldown_groups"
        )
        if group_resource not in resources:
            raise CapsuleError(
                "junction_drilldown.drilldown_groups references an unknown resource"
            )
    archive_sample = _safe_identifier(
        drilldown.get("archive_sample"), "junction_drilldown.archive_sample"
    )
    if archive_sample not in drilldown_archives:
        raise CapsuleError("junction_drilldown.archive_sample is absent from archives")

    for field in ("shape_routes", "allow_unstamped", "novel_only"):
        if field in event and not isinstance(event[field], bool):
            raise CapsuleError(f"event_discovery.{field} must be boolean")
    for field in (
        "min_group_umi_classes",
        "min_donors",
        "min_samples",
        "min_umi_classes",
        "min_side_umi_classes",
        "max_candidates",
        "max_candidates_considered",
        "max_routed_entries",
        "max_exact_match_attempts",
        "max_annotation_comparisons",
    ):
        _story_integer(event, field, 1, "event_discovery")
    for field in ("min_support", "terminal_cluster_bp", "max_terminal_events"):
        _story_integer(event, field, 0, "event_discovery")
    for field in (
        "expected_min_exact_umi_classes",
        "expected_min_exact_donors",
        "expected_rank",
        "expected_comparison_compatible_transcripts",
    ):
        _story_integer(event, field, 1, "event_discovery")
    _public_text(event["expected_entity_id"], "event_discovery.expected_entity_id")
    if event["expected_gap_primary_class"] not in {
        "missing_junction",
        "boundary",
        "strand",
        "overlap",
    }:
        raise CapsuleError("event_discovery.expected_gap_primary_class is invalid")
    if not isinstance(event["expected_annotation_incompatible"], bool):
        raise CapsuleError("event_discovery.expected_annotation_incompatible must be boolean")
    kinds = event.get("kinds", [])
    allowed_kinds = {"junction", "alt-acceptor", "alt-donor", "cassette", "terminal-tail"}
    if (
        not isinstance(kinds, list)
        or any(not isinstance(kind, str) for kind in kinds)
        or len(kinds) != len(set(kinds))
        or any(kind not in allowed_kinds for kind in kinds)
    ):
        raise CapsuleError("event_discovery.kinds contains an invalid or duplicate kind")
    required_groups = event.get("require_groups", [])
    if (
        not isinstance(required_groups, list)
        or any(not isinstance(group, str) for group in required_groups)
        or len(required_groups) != len(set(required_groups))
    ):
        raise CapsuleError("event_discovery.require_groups must contain unique group names")
    for group in required_groups:
        _safe_identifier(group, "event_discovery.require_groups")
    if "solo_strand" in event and event["solo_strand"] not in {
        "forward",
        "reverse",
        "unstranded",
    }:
        raise CapsuleError("event_discovery.solo_strand is invalid")
    for field in ("assembly", "annotation_label", "comparison_annotation_label"):
        _public_text(event[field], f"event_discovery.{field}")
    if "annotation_digest" in event and not re.fullmatch(
        r"blake3:[0-9a-f]{64}", str(event["annotation_digest"])
    ):
        raise CapsuleError("event_discovery.annotation_digest is invalid")

    for field in ("shape_routes", "allow_unstamped", "allow_full_scan", "emit_membership"):
        if field in drilldown and not isinstance(drilldown[field], bool):
            raise CapsuleError(f"junction_drilldown.{field} must be boolean")
    if drilldown["emit_membership"] is not True:
        raise CapsuleError(
            "junction_drilldown.emit_membership must be true for the frozen query contract"
        )
    for field in ("max_memberships", "max_pattern_rows", "max_chunks", "max_evidence_records"):
        _story_integer(drilldown, field, 1, "junction_drilldown")
    _story_integer(drilldown, "max_terminal_events", 0, "junction_drilldown")
    _story_integer(drilldown, "expected_min_selected_units", 1, "junction_drilldown")
    _public_text(drilldown["story_note"], "junction_drilldown.story_note")
    _collection_junction(drilldown["junction"], "junction_drilldown.junction")
    for field in ("expression", "universe"):
        _public_text(drilldown[field], f"junction_drilldown.{field}")
    predicates = _expect_mapping(drilldown["predicates"], "junction_drilldown.predicates")
    if not predicates:
        raise CapsuleError("junction_drilldown.predicates must not be empty")
    for name, predicate in predicates.items():
        _safe_identifier(name, "junction_drilldown predicate name")
        _public_text(predicate, f"junction_drilldown predicate {name}")
    required_true = drilldown["required_true_predicates"]
    if (
        not isinstance(required_true, list)
        or not required_true
        or any(not isinstance(name, str) for name in required_true)
        or len(required_true) != len(set(required_true))
        or any(name not in predicates for name in required_true)
    ):
        raise CapsuleError(
            "junction_drilldown.required_true_predicates must be unique predicate names"
        )
    enums = {
        "unit": {"molecule-record", "umi-class"},
        "region_match": {"anchor", "aligned-block"},
        "placements": {"unique", "direct", "all"},
        "aggregation": {"auto", "cell", "group", "bulk"},
    }
    for field, values in enums.items():
        if field in drilldown and drilldown[field] not in values:
            raise CapsuleError(f"junction_drilldown.{field} is invalid")


def _notices(value: Any) -> tuple[str, list[dict[str, str]]]:
    rows: list[dict[str, str]] = []
    lines = ["# Third-party data notices", ""]
    seen_names: set[str] = set()
    for index, raw in enumerate(_expect_list(value, "third_party_notices")):
        item = _expect_mapping(raw, f"third_party_notices[{index}]")
        row = {
            "name": _public_text(item.get("name"), f"notice {index} name"),
            "source_url": _https(item.get("source_url"), f"notice {index} source_url"),
            "terms": _public_text(item.get("terms"), f"notice {index} terms"),
            "terms_url": _https(item.get("terms_url"), f"notice {index} terms_url"),
            "citation": _public_text(item.get("citation"), f"notice {index} citation"),
        }
        if row["name"] in seen_names:
            raise CapsuleError(f"third_party_notices repeats {row['name']}")
        seen_names.add(row["name"])
        rows.append(row)
        lines.extend(
            [
                f"## {row['name']}",
                "",
                f"- Source: {row['source_url']}",
                f"- Terms: {row['terms']} ({row['terms_url']})",
                f"- Citation: {row['citation']}",
                "",
            ]
        )
    return "\n".join(lines), rows


def _group_map_scope(value: Any) -> dict[str, str]:
    scope = _expect_mapping(value, "group_map_scope")
    if set(scope) != GROUP_MAP_SCOPE_FIELDS:
        raise CapsuleError(
            f"group_map_scope must contain exactly {sorted(GROUP_MAP_SCOPE_FIELDS)}"
        )
    return {
        field: _public_text(scope[field], f"group_map_scope.{field}")
        for field in sorted(GROUP_MAP_SCOPE_FIELDS)
    }


def _preflight_story_resources(spec: dict[str, Any]) -> None:
    _validate_input_declaration(spec.get("genome"), "genome")
    resource_names = {"collection_groups", "design", "drilldown_groups"}
    output_filenames = {
        "BUILD-RECORD.json",
        "THIRD-PARTY-NOTICES.md",
        "collection-groups.tsv",
        "donors.tsv",
        "drilldown-groups.tsv",
    }
    archive_map: dict[str, str] = {}
    for raw in _expect_list(spec.get("donors"), "donors"):
        donor = _expect_mapping(raw, "donor")
        sample = _safe_identifier(donor.get("sample"), "donor.sample")
        _safe_identifier(donor.get("donor"), f"{sample}.donor")
        accessions = donor.get("accessions")
        if not isinstance(accessions, list) or not accessions:
            raise CapsuleError(f"{sample}.accessions must be a nonempty array")
        for accession in accessions:
            _safe_identifier(accession, f"{sample}.accessions")
        for field in (
            "bam",
            "whitelist",
            "groups",
            "junction_catalogue",
            "alignment_log",
        ):
            _validate_input_declaration(donor.get(field), f"{sample}.{field}")
        resource = _safe_identifier(
            donor.get("archive_resource"), f"{sample}.archive_resource"
        )
        filename = _safe_filename(donor.get("archive_filename"))
        if not filename.endswith(".aie"):
            raise CapsuleError(f"{sample}.archive_filename must end in .aie")
        if sample in archive_map:
            raise CapsuleError(f"duplicate donor sample {sample}")
        if resource in resource_names:
            raise CapsuleError(f"duplicate or reserved resource name {resource}")
        if filename in output_filenames:
            raise CapsuleError(f"duplicate capsule filename {filename}")
        archive_map[sample] = resource
        resource_names.add(resource)
        output_filenames.add(filename)

    annotations = _expect_mapping(spec.get("annotations"), "annotations")
    if set(annotations) != {"before", "after"}:
        raise CapsuleError("annotations must contain exactly before and after")
    for side in ("before", "after"):
        declaration = _expect_mapping(annotations[side], f"annotations.{side}")
        _validate_input_declaration(declaration, f"annotations.{side}")
        resource = _safe_identifier(
            declaration.get("resource"), f"annotations.{side}.resource"
        )
        compiled = _safe_filename(declaration.get("filename"))
        subset = _safe_filename(declaration.get("source_filename"))
        if not compiled.endswith(".aic") or not subset.endswith(".gtf"):
            raise CapsuleError("annotation filenames must end in .aic and .gtf")
        if resource in resource_names or f"{resource}_gtf" in resource_names:
            raise CapsuleError(f"duplicate annotation resource {resource}")
        if compiled == subset or {compiled, subset}.intersection(output_filenames):
            raise CapsuleError(f"annotation {side} uses a duplicate capsule filename")
        resource_names.update((resource, f"{resource}_gtf"))
        output_filenames.update((compiled, subset))

    stories = _expect_mapping(spec.get("stories"), "stories")
    _validate_story_references(stories, resource_names)
    event = _expect_mapping(stories["event_discovery"], "event_discovery")
    drilldown = _expect_mapping(stories["junction_drilldown"], "junction_drilldown")
    if event["archives"] != archive_map:
        raise CapsuleError(
            "event_discovery.archives must map every staged sample to its archive resource"
        )
    if drilldown["archives"] != archive_map:
        raise CapsuleError(
            "junction_drilldown.archives must map every staged sample to its archive resource"
        )
    if event["collection_groups"] != "collection_groups" or event["design"] != "design":
        raise CapsuleError(
            "event_discovery must reference the generated collection_groups and design"
        )
    if drilldown.get("drilldown_groups") != "drilldown_groups":
        raise CapsuleError(
            "junction_drilldown.drilldown_groups must reference drilldown_groups"
        )


def build(spec_path: Path, source_root: Path, aie: Path, samtools: Path, output_dir: Path) -> None:
    spec = json.loads(spec_path.read_text(encoding="utf-8"))
    spec = _expect_mapping(spec, "build specification")
    if spec.get("schema") != BUILD_SCHEMA:
        raise CapsuleError(f"build specification schema must be {BUILD_SCHEMA}")
    capsule_id = _safe_identifier(spec.get("capsule_id"), "capsule_id")
    summary = _public_text(spec.get("summary"), "summary")
    software = _expect_mapping(spec.get("software"), "software")
    version = _public_text(software.get("version"), "software.version")
    revision = software.get("revision")
    if not isinstance(revision, str) or not re.fullmatch(r"[0-9a-f]{40}", revision):
        raise CapsuleError("software.revision must be a full 40-character Git commit")
    windows = _windows(spec)
    alignment = _expect_mapping(spec.get("alignment"), "alignment")
    if alignment.get("junction_discovery") != "per-library-two-pass":
        raise CapsuleError("the public demo builder requires per-library-two-pass provenance")
    _safe_identifier(alignment.get("program"), "alignment.program")
    _public_text(alignment.get("version"), "alignment.version")
    _public_text(alignment.get("command_line"), "alignment.command_line")
    _public_text(alignment.get("index_identity"), "alignment.index_identity")
    _public_text(alignment.get("chemistry"), "alignment.chemistry")
    _assert_public_strings(alignment, "alignment")
    public_alignment = {
        **alignment,
        "runtime_annotation_status": (
            "validated per donor from exact STAR Log.out: sjdbGTFfile=- and "
            "sjdbFileChrStartEnd=-"
        ),
        "index_manifest_status": (
            "caller-declared; the capsule builder does not consume the STAR index files or "
            "verify the listed component hashes"
        ),
    }
    assembly = _public_text(spec.get("assembly"), "assembly")
    group_map_scope = _group_map_scope(spec.get("group_map_scope"))
    _preflight_story_resources(spec)
    notice_text, notices = _notices(spec.get("third_party_notices"))
    source_root = source_root.resolve(strict=True)
    for raw_donor in _expect_list(spec.get("donors"), "donors"):
        donor = _expect_mapping(raw_donor, "donor")
        sample = _safe_identifier(donor.get("sample"), "donor.sample")
        log_declaration = _expect_mapping(
            donor.get("alignment_log"), f"{sample}.alignment_log"
        )
        log_path = _resolve_input(
            source_root, log_declaration, f"{sample}.alignment_log"
        )
        _validate_star_log(log_path, alignment, f"{sample}.alignment_log")
    aie = aie.resolve(strict=True)
    samtools = samtools.resolve(strict=True)
    aie_sha256 = _sha256(aie)
    samtools_sha256 = _sha256(samtools)
    if _tool_version(aie, ["--version"]) != f"aie {version}":
        raise CapsuleError("--aie version does not match software.version")
    samtools_version = _tool_version(samtools, ["--version"]).splitlines()[0]
    if _sha256(aie) != aie_sha256 or _sha256(samtools) != samtools_sha256:
        raise CapsuleError("a build executable changed during preflight")
    genome_decl = _expect_mapping(spec.get("genome"), "genome")
    genome = _resolve_input(source_root, genome_decl, "genome")

    with _new_output_directory(output_dir) as output:
        output_filenames = {
            "BUILD-RECORD.json",
            "THIRD-PARTY-NOTICES.md",
            "collection-groups.tsv",
            "donors.tsv",
            "drilldown-groups.tsv",
        }
        work = Path(tempfile.mkdtemp(prefix=".work-", dir=output))
        shared_work = work / "shared"
        shared_work.mkdir()
        staged_aie = shared_work / "aie"
        _stage_executable(aie, staged_aie, aie_sha256)
        staged_genome = shared_work / (
            "reference.fa.gz" if genome.name.endswith(".gz") else "reference.fa"
        )
        _stage_input(genome, staged_genome, genome_decl["sha256"])
        bed = work / "selection-windows.bed"
        _write_bytes_exclusive(
            bed,
            "".join(f"{chrom}\t{start}\t{end}\n" for chrom, start, end in windows).encode(),
        )
        resources: dict[str, dict[str, Any]] = {}
        public_inputs: dict[str, Any] = {
            "genome": _public_input_record(genome_decl),
            "donors": [],
            "annotations": {},
        }
        global_groups: list[tuple[str, str, str]] = []
        groups_by_sample: dict[str, list[tuple[str, str]]] = {}
        design: list[tuple[str, str]] = []
        archive_resources_by_sample: dict[str, str] = {}
        terminal_tail_events = 0
        donors = _expect_list(spec.get("donors"), "donors")
        samples_seen: set[str] = set()
        for raw in donors:
            donor_spec = _expect_mapping(raw, "donor")
            sample = _safe_identifier(donor_spec.get("sample"), "donor.sample")
            if sample in samples_seen:
                raise CapsuleError(f"duplicate donor sample {sample}")
            samples_seen.add(sample)
            archive_resource_declared = _safe_identifier(
                donor_spec.get("archive_resource"), f"{sample}.archive_resource"
            )
            if archive_resource_declared in {
                "collection_groups",
                "design",
                "drilldown_groups",
            }:
                raise CapsuleError(
                    f"archive resource name {archive_resource_declared} is reserved"
                )
            if archive_resource_declared in resources:
                raise CapsuleError(f"duplicate resource name {archive_resource_declared}")
            archive_filename_declared = _safe_filename(donor_spec.get("archive_filename"))
            if archive_filename_declared in output_filenames:
                raise CapsuleError(f"duplicate capsule filename {archive_filename_declared}")
            output_filenames.add(archive_filename_declared)
            (
                archive_resource,
                archive_asset,
                group_rows,
                tail_events,
                donor_provenance,
            ) = _stage_donor(
                aie=staged_aie,
                samtools=samtools,
                source_root=source_root,
                output=output,
                work=work,
                windows_bed=bed,
                windows=windows,
                genome=staged_genome,
                genome_sha256=genome_decl["sha256"],
                alignment=alignment,
                donor_spec=donor_spec,
            )
            if archive_resource != archive_resource_declared:
                raise AssertionError("validated archive resource changed while staging")
            resources[archive_resource] = archive_asset
            archive_resources_by_sample[sample] = archive_resource
            groups_by_sample[sample] = group_rows
            donor_name = donor_provenance["donor"]
            design.append((sample, donor_name))
            global_groups.extend((sample, barcode, group) for barcode, group in group_rows)
            terminal_tail_events += tail_events
            public_inputs["donors"].append(donor_provenance)
        if terminal_tail_events == 0:
            raise CapsuleError("demo windows yielded no terminal-tail events in any donor")

        design_path = output / "donors.tsv"
        _write_bytes_exclusive(
            design_path,
            ("sample\tdonor\n" + "".join(f"{s}\t{d}\n" for s, d in sorted(design))).encode(),
        )
        groups_path = output / "collection-groups.tsv"
        _write_bytes_exclusive(
            groups_path,
            (
                "sample\tbarcode\tgroup\n"
                + "".join(f"{s}\t{b}\t{g}\n" for s, b, g in sorted(global_groups))
            ).encode(),
        )
        resources["design"] = _asset(design_path)
        resources["collection_groups"] = _asset(groups_path)

        stories = _expect_mapping(spec.get("stories"), "stories")
        drilldown = _expect_mapping(stories.get("junction_drilldown"), "junction_drilldown")
        drilldown_sample = _safe_identifier(
            drilldown.get("archive_sample"), "junction_drilldown.archive_sample"
        )
        if drilldown_sample not in groups_by_sample:
            raise CapsuleError("junction_drilldown.archive_sample is not a staged donor")
        drilldown_groups_path = output / "drilldown-groups.tsv"
        _write_bytes_exclusive(
            drilldown_groups_path,
            "".join(
                f"{barcode}\t{group}\n"
                for barcode, group in groups_by_sample[drilldown_sample]
            ).encode(),
        )
        resources["drilldown_groups"] = _asset(drilldown_groups_path)

        annotations = _expect_mapping(spec.get("annotations"), "annotations")
        if set(annotations) != {"before", "after"}:
            raise CapsuleError("annotations must contain exactly before and after")
        annotation_work = work / "annotations"
        annotation_work.mkdir()
        for side in ("before", "after"):
            declaration = _expect_mapping(annotations[side], f"annotations.{side}")
            source = _resolve_input(source_root, declaration, f"annotations.{side}")
            staged_source = annotation_work / f"source-{side}.gtf"
            _stage_input(source, staged_source, declaration["sha256"])
            resource = _safe_identifier(declaration.get("resource"), f"annotations.{side}.resource")
            filename = _safe_filename(declaration.get("filename"))
            source_filename = _safe_filename(declaration.get("source_filename"))
            if not filename.endswith(".aic") or not source_filename.endswith(".gtf"):
                raise CapsuleError("annotation filenames must end in .aic and .gtf")
            if filename == source_filename or {filename, source_filename}.intersection(
                output_filenames
            ):
                raise CapsuleError(
                    f"annotation {side} uses a duplicate capsule filename"
                )
            output_filenames.update((filename, source_filename))
            if resource in resources or f"{resource}_gtf" in resources:
                raise CapsuleError(f"duplicate annotation resource {resource}")
            subset = output / source_filename
            genes = _subset_gtf(staged_source, windows, subset)
            if _sha256(staged_source) != declaration["sha256"]:
                raise CapsuleError(f"annotations.{side} changed while subsetting")
            compiled = output / filename
            _run(
                [str(staged_aie), "compile-annotation", subset.name, "--out", compiled.name],
                cwd=output,
            )
            if not compiled.is_file() or compiled.stat().st_size == 0:
                raise CapsuleError(f"compile-annotation did not create {filename}")
            resources[resource] = _asset(compiled)
            resources[f"{resource}_gtf"] = _asset(subset)
            public_inputs["annotations"][side] = {
                **_public_input_record(declaration),
                "complete_genes_selected": genes,
                "label": _public_text(declaration.get("label"), f"annotations.{side}.label"),
                "compiled_resource": resource,
                "subset_gtf_resource": f"{resource}_gtf",
            }

        event_archives = _expect_mapping(
            _expect_mapping(stories.get("event_discovery"), "event_discovery").get("archives"),
            "event_discovery.archives",
        )
        drilldown_archives = _expect_mapping(
            drilldown.get("archives"), "junction_drilldown.archives"
        )
        if event_archives != archive_resources_by_sample:
            raise CapsuleError(
                "event_discovery.archives must map every staged sample to its archive resource"
            )
        if drilldown_archives != archive_resources_by_sample:
            raise CapsuleError(
                "junction_drilldown.archives must map every staged sample to its archive resource"
            )
        if stories["event_discovery"].get("collection_groups") != "collection_groups":
            raise CapsuleError(
                "event_discovery.collection_groups must reference collection_groups"
            )
        if stories["event_discovery"].get("design") != "design":
            raise CapsuleError("event_discovery.design must reference design")
        if drilldown.get("drilldown_groups") != "drilldown_groups":
            raise CapsuleError(
                "junction_drilldown.drilldown_groups must reference drilldown_groups"
            )
        _validate_story_references(stories, set(resources))
        _assert_public_strings(stories, "stories")
        if _sha256(staged_aie) != aie_sha256 or _sha256(samtools) != samtools_sha256:
            raise CapsuleError("a build executable changed during capsule construction")
        record = {
            "schema": RECORD_SCHEMA,
            "capsule_id": capsule_id,
            "summary": summary,
            "assembly": assembly,
            "software": {
                "version": version,
                "revision": revision,
                "aie_sha256": aie_sha256,
                "samtools": samtools_version,
                "samtools_sha256": samtools_sha256,
            },
            "selection": {
                "coordinate_system": "0-based-half-open",
                "windows": [
                    {"chrom": chrom, "start": start, "end": end}
                    for chrom, start, end in windows
                ],
                "annotation_subset": (
                    "all records for every gene with any record overlapping a selection window"
                ),
                "bam_header": (
                    "coordinate-sort the selected records; retain the sorted @HD/@SQ with "
                    "SO:coordinate; remove source @PG/@RG/@CO; add a path-free public scope @CO; "
                    "do not relabel a normalized caller declaration as verified BAM-header "
                    "alignment metadata"
                ),
                "region_filtering": (
                    "retain complete BAM molecule/read-derived alignment records that overlap "
                    "at least one selected window; omit such records wholly outside every window"
                ),
                "junction_catalogue_scope": (
                    "embed each donor's complete genome-wide aggregate STAR pass-1 junction "
                    "catalogue as exact root-bound alignment provenance; it contains junction "
                    "coordinates and aggregate support columns but no sequences, cell barcodes, "
                    "or UMI identities"
                ),
                "multimapper_scope": (
                    "locus subsets omit alignment alternatives outside the selected windows, "
                    "so placement multiplicity and completeness are capsule-local and must not "
                    "be interpreted genome-wide"
                ),
            },
            "alignment": public_alignment,
            "group_map_scope": group_map_scope,
            "inputs": public_inputs,
            "resources": dict(sorted(resources.items())),
            "stories": stories,
            "terminal_tail_events": terminal_tail_events,
            "third_party_notices": notices,
        }
        _assert_public_strings(record, "build record")
        _write_json_exclusive(output / "BUILD-RECORD.json", record)
        _write_bytes_exclusive(output / "THIRD-PARTY-NOTICES.md", notice_text.encode())
        shutil.rmtree(work)
    print(output_dir.resolve())


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spec", required=True, type=Path)
    parser.add_argument("--source-root", required=True, type=Path)
    parser.add_argument("--aie", required=True, type=Path)
    parser.add_argument("--samtools", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    build(args.spec, args.source_root, args.aie, args.samtools, args.output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
