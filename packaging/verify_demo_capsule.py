#!/usr/bin/env python3
"""Verify a finalized demo capsule and execute its three scientific stories."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import subprocess
import tempfile
from typing import Any

from build_demo_capsule import (
    ABSOLUTE_PATH,
    CapsuleError,
    HEX64,
    PRIVATE_TEXT,
    RECORD_SCHEMA,
    _assert_public_strings,
    _notices,
    _safe_filename,
    _safe_identifier,
    _sha256,
    _validate_story_references,
    _verified_archive_inspection,
)
from finalize_demo_capsule import (
    MANIFEST_SCHEMA,
    FINALIZATION_SCHEMA,
    _binary_identity,
    _documentation,
    _github_release_tag,
    _immutable_url,
    _member_name,
    _wheel_identity,
)


EXPECTED_TABLES = {
    "find-events": {
        "capabilities",
        "entities",
        "components",
        "counts",
        "terminal_anchors",
        "terminal_counts",
    },
    "cooccur": {"predicates", "patterns", "memberships"},
}

MEMBERSHIP_COLUMNS = {
    "cell_id",
    "unit_id",
    "barcode",
    "umi_class",
    "chunk",
    "local_record",
    "global_record",
    "contributing_records",
    "pattern_mask",
    "completeness_mask",
    "matched_predicates",
    "selection_state",
    "selected",
}


def _run_json(command: list[str]) -> dict[str, Any]:
    try:
        result = subprocess.run(
            command,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
    except OSError as error:
        raise CapsuleError(
            f"cannot execute scientific verification command: {command[0]}"
        ) from error
    except subprocess.CalledProcessError as error:
        raise CapsuleError(
            f"scientific verification command failed ({' '.join(command[:3])}): "
            f"{error.stderr[-4000:]}"
        ) from error
    try:
        value = json.loads(result.stdout)
    except json.JSONDecodeError as error:
        raise CapsuleError(f"command emitted invalid JSON: {' '.join(command[:3])}") from error
    if not isinstance(value, dict):
        raise CapsuleError("command result must be a JSON object")
    return value


def _table_records(table: Any, label: str) -> list[dict[str, Any]]:
    try:
        fields = [field["name"] for field in table["schema"]["fields"]]
        rows = table["rows"]
    except (KeyError, TypeError) as error:
        raise CapsuleError(f"malformed typed result table {label}") from error
    if (
        not isinstance(rows, list)
        or any(not isinstance(field, str) for field in fields)
        or len(fields) != len(set(fields))
    ):
        raise CapsuleError(f"invalid typed result table {label}")
    records = []
    for row in rows:
        if not isinstance(row, list) or len(row) != len(fields):
            raise CapsuleError(f"table {label} contains a malformed row")
        records.append(dict(zip(fields, row)))
    return records


def _uniform_tables(bundle: dict[str, Any]) -> dict[str, list[dict[str, Any]]]:
    try:
        tables = bundle["data"]["tables"]
    except (KeyError, TypeError) as error:
        raise CapsuleError("uniform result lacks data.tables") from error
    if not isinstance(tables, list):
        raise CapsuleError("uniform result tables must be an array")
    result: dict[str, list[dict[str, Any]]] = {}
    for table in tables:
        try:
            name = table["name"]
        except (KeyError, TypeError) as error:
            raise CapsuleError("malformed typed result table") from error
        if not isinstance(name, str) or name in result:
            raise CapsuleError("invalid or duplicate typed table name")
        result[name] = _table_records(table, name)
    return result


def _named_table(bundle: dict[str, Any], name: str) -> list[dict[str, Any]]:
    try:
        table = bundle["data"][name]
    except (KeyError, TypeError) as error:
        raise CapsuleError(f"result lacks data.{name}") from error
    return _table_records(table, name)


def _summary(bundle: dict[str, Any]) -> dict[str, Any]:
    try:
        summary = bundle["data"]["summary"]
    except (KeyError, TypeError) as error:
        raise CapsuleError("uniform result lacks data.summary") from error
    if not isinstance(summary, dict):
        raise CapsuleError("uniform result summary must be an object")
    return summary


def _validate_asset(name: str, asset: Any, *, archive: bool = False) -> dict[str, Any]:
    if not isinstance(asset, dict):
        raise CapsuleError(f"asset {name} must be an object")
    allowed = {"url", "sha256", "filename", "member", "archive_root"}
    if set(asset).difference(allowed) or not {"url", "sha256", "filename"}.issubset(asset):
        raise CapsuleError(f"asset {name} has missing or unknown fields")
    filename = _safe_filename(asset["filename"])
    _immutable_url(
        asset["url"],
        f"asset {name} URL",
        filename=filename,
        kind="pypi-file" if name == "software.python_wheel" else "github-release",
    )
    if not isinstance(asset["sha256"], str) or not HEX64.fullmatch(asset["sha256"]):
        raise CapsuleError(f"asset {name} has an invalid SHA-256")
    root = asset.get("archive_root")
    if archive and not isinstance(root, str):
        raise CapsuleError(f"archive asset {name} lacks archive_root")
    if root is not None and not re.fullmatch(r"aie-directory-root-v2:[0-9a-f]{64}", root):
        raise CapsuleError(f"asset {name} has an invalid archive_root")
    return asset


def _validate_manifest(manifest: Any) -> dict[str, Any]:
    if not isinstance(manifest, dict) or manifest.get("schema") != MANIFEST_SCHEMA:
        raise CapsuleError(f"manifest schema must be {MANIFEST_SCHEMA}")
    if set(manifest) != {"schema", "software", "resources", "stories"}:
        raise CapsuleError("manifest has missing or unknown top-level fields")
    software = manifest["software"]
    if not isinstance(software, dict) or set(software) != {"version", "aie", "python_wheel"}:
        raise CapsuleError("manifest software object is malformed")
    if not isinstance(software["version"], str) or not software["version"]:
        raise CapsuleError("manifest software.version must be nonempty")
    _validate_asset("software.aie", software["aie"])
    _validate_asset("software.python_wheel", software["python_wheel"])
    if _github_release_tag(software["aie"]["url"]) != f"v{software['version']}":
        raise CapsuleError("software.aie URL tag must equal v<software.version>")
    if "member" not in software["aie"]:
        raise CapsuleError("software.aie must name its archive member")
    _member_name(software["aie"]["member"])
    resources = manifest["resources"]
    if not isinstance(resources, dict) or not resources:
        raise CapsuleError("manifest resources must be a nonempty object")
    filenames: set[str] = set()
    for name, asset in resources.items():
        _safe_identifier(name, "manifest resource name")
        _validate_asset(name, asset, archive=asset.get("filename", "").endswith(".aie"))
        if asset["filename"].endswith(".aicollection"):
            raise CapsuleError("path-bound .aicollection files must not be published")
        if asset["filename"] in filenames:
            raise CapsuleError(f"resources repeat filename {asset['filename']}")
        filenames.add(asset["filename"])
    stories = manifest["stories"]
    required_stories = {"annotation_reinterpretation", "event_discovery", "junction_drilldown"}
    if not isinstance(stories, dict) or set(stories) != required_stories:
        raise CapsuleError(f"manifest stories must contain exactly {sorted(required_stories)}")
    _validate_story_references(stories, set(resources))
    _assert_public_strings(manifest, "manifest")
    return manifest


def _verify_checksums(directory: Path) -> None:
    checksum_file = directory / "SHA256SUMS"
    if not checksum_file.is_file() or checksum_file.is_symlink():
        raise CapsuleError("capsule lacks a regular SHA256SUMS file")
    expected: dict[str, str] = {}
    for line_no, line in enumerate(checksum_file.read_text(encoding="ascii").splitlines(), 1):
        fields = line.split("  ")
        if len(fields) != 2 or not HEX64.fullmatch(fields[0]):
            raise CapsuleError(f"malformed SHA256SUMS line {line_no}")
        filename = _safe_filename(fields[1])
        if filename == "SHA256SUMS" or filename in expected:
            raise CapsuleError(f"invalid or duplicate SHA256SUMS filename {filename}")
        expected[filename] = fields[0]
    entries = list(directory.iterdir())
    unexpected = [path.name for path in entries if not path.is_file() or path.is_symlink()]
    if unexpected:
        raise CapsuleError(f"capsule contains a non-regular entry: {sorted(unexpected)}")
    actual = {path.name for path in entries if path.name != "SHA256SUMS"}
    if actual != set(expected):
        raise CapsuleError("SHA256SUMS does not cover exactly the finalized capsule files")
    for filename, digest in expected.items():
        if _sha256(directory / filename) != digest:
            raise CapsuleError(f"SHA-256 mismatch for {filename}")


def _verify_public_bytes(directory: Path) -> None:
    for path in directory.iterdir():
        if not path.is_file() or path.name == "SHA256SUMS":
            continue
        markers = tuple(marker.lower().encode() for marker in PRIVATE_TEXT)
        overlap = max(map(len, markers)) - 1
        carry = b""
        with path.open("rb") as handle:
            while block := handle.read(8 << 20):
                contents = (carry + block).lower()
                if any(marker in contents for marker in markers):
                    raise CapsuleError(f"{path.name} contains a private filesystem locator")
                carry = contents[-overlap:]
        if path.suffix in {".json", ".md", ".tsv", ".gtf", ".txt"}:
            try:
                text = path.read_text(encoding="utf-8")
            except UnicodeDecodeError as error:
                raise CapsuleError(f"text capsule asset {path.name} is not UTF-8") from error
            if ABSOLUTE_PATH.search(text):
                raise CapsuleError(f"{path.name} contains an absolute filesystem locator")


def _resource_paths(directory: Path, manifest: dict[str, Any]) -> dict[str, Path]:
    result = {}
    for name, asset in manifest["resources"].items():
        path = directory / asset["filename"]
        if not path.is_file() or path.is_symlink() or _sha256(path) != asset["sha256"]:
            raise CapsuleError(f"local resource {name} is missing or changed")
        result[name] = path
    return result


def _verify_capsule_records(directory: Path, manifest: dict[str, Any]) -> dict[str, Any]:
    required_files = {
        "BUILD-RECORD.json",
        "FINALIZATION-RECORD.json",
        "README.md",
        "RELEASE-NOTES.md",
        "THIRD-PARTY-NOTICES.md",
    }
    missing = sorted(name for name in required_files if not (directory / name).is_file())
    if missing:
        raise CapsuleError(f"capsule lacks required self-description files: {missing}")
    try:
        build_record = json.loads((directory / "BUILD-RECORD.json").read_text(encoding="utf-8"))
        final_record = json.loads(
            (directory / "FINALIZATION-RECORD.json").read_text(encoding="utf-8")
        )
    except json.JSONDecodeError as error:
        raise CapsuleError("capsule records must contain valid JSON") from error
    if not isinstance(build_record, dict) or build_record.get("schema") != RECORD_SCHEMA:
        raise CapsuleError(f"build record schema must be {RECORD_SCHEMA}")
    _assert_public_strings(build_record, "build record")
    build_software = build_record.get("software")
    if not isinstance(build_software, dict):
        raise CapsuleError("build record software object is malformed")
    if build_software.get("version") != manifest["software"]["version"]:
        raise CapsuleError("build-record software version differs from the manifest")
    revision = build_software.get("revision")
    if not isinstance(revision, str) or not re.fullmatch(r"[0-9a-f]{40}", revision):
        raise CapsuleError("build-record software revision must be full lowercase hexadecimal")
    if build_record.get("stories") != manifest["stories"]:
        raise CapsuleError("build-record stories differ from the manifest")
    expected_notices, _ = _notices(build_record.get("third_party_notices"))
    if (directory / "THIRD-PARTY-NOTICES.md").read_text(
        encoding="utf-8"
    ) != expected_notices:
        raise CapsuleError("third-party notices differ from the build record")
    built_resources = build_record.get("resources")
    if not isinstance(built_resources, dict) or set(built_resources) != set(manifest["resources"]):
        raise CapsuleError("build-record resources differ from the manifest")
    for name, asset in manifest["resources"].items():
        built = built_resources[name]
        if not isinstance(built, dict):
            raise CapsuleError(f"build-record resource {name} is malformed")
        for field in ("filename", "sha256", "archive_root"):
            if built.get(field) != asset.get(field):
                raise CapsuleError(f"build-record resource {name} differs in {field}")
        if built.get("bytes") != (directory / asset["filename"]).stat().st_size:
            raise CapsuleError(f"build-record resource {name} differs in byte size")

    expected_final_fields = {
        "schema",
        "capsule_id",
        "manifest",
        "software",
        "source",
        "data_base_url",
    }
    if not isinstance(final_record, dict) or set(final_record) != expected_final_fields:
        raise CapsuleError("finalization record has missing or unknown fields")
    if final_record.get("schema") != FINALIZATION_SCHEMA:
        raise CapsuleError(f"finalization record schema must be {FINALIZATION_SCHEMA}")
    if final_record.get("capsule_id") != build_record.get("capsule_id"):
        raise CapsuleError("finalization and build records name different capsules")
    if final_record.get("software") != manifest["software"]:
        raise CapsuleError("finalization-record software differs from the manifest")
    expected_source = {
        "revision": build_software.get("revision"),
        "software_release_tag": _github_release_tag(manifest["software"]["aie"]["url"]),
    }
    if final_record.get("source") != expected_source:
        raise CapsuleError("finalization-record source tag/revision binding is inconsistent")
    data_base_url = _immutable_url(final_record.get("data_base_url"), "data_base_url").rstrip(
        "/"
    )
    for name, asset in manifest["resources"].items():
        if asset["url"] != f"{data_base_url}/{asset['filename']}":
            raise CapsuleError(
                f"resource {name} URL is outside the finalized capsule data release"
            )
    expected_manifest = {
        "filename": "demo-manifest.json",
        "sha256": _sha256(directory / "demo-manifest.json"),
        "url": f"{data_base_url}/demo-manifest.json",
    }
    if final_record.get("manifest") != expected_manifest:
        raise CapsuleError("finalization-record manifest identity is inconsistent")
    expected_readme, expected_release_notes = _documentation(build_record, manifest)
    if (directory / "README.md").read_text(encoding="utf-8") != expected_readme:
        raise CapsuleError("README.md differs from the build record and manifest")
    if (
        (directory / "RELEASE-NOTES.md").read_text(encoding="utf-8")
        != expected_release_notes
    ):
        raise CapsuleError("RELEASE-NOTES.md differs from the build record and manifest")
    return build_record


def _verify_archives(
    aie: Path, manifest: dict[str, Any], resources: dict[str, Path]
) -> dict[str, Any]:
    inspected: dict[str, Any] = {}
    for name, asset in manifest["resources"].items():
        if not asset["filename"].endswith(".aie"):
            continue
        result = _run_json(
            [str(aie), "inspect-archive", str(resources[name]), "--verify-content", "--json"]
        )
        observed, _, _ = _verified_archive_inspection(result, f"archive {name}")
        if observed != asset["archive_root"]:
            raise CapsuleError(f"archive root mismatch for {name}")
        inspected[name] = result
    if not inspected:
        raise CapsuleError("manifest contains no archive resources")
    if (
        sum(
            item["molecular_evidence"]["terminal_tail"]["events"]
            for item in inspected.values()
        )
        == 0
    ):
        raise CapsuleError("capsule archives contain no terminal-tail events")
    return inspected


def _build_collection(
    aie: Path,
    destination: Path,
    archive_map: dict[str, str],
    resources: dict[str, Path],
    story: dict[str, Any],
) -> None:
    command = [str(aie), "collection", "build"]
    for sample, resource in sorted(archive_map.items()):
        command.append(f"--sample={sample}={resources[resource]}")
    if story.get("shape_routes", True):
        command.append("--shape-routes")
    if story.get("allow_unstamped", False):
        command.append("--allow-unstamped")
    command.extend([f"--out={destination}", "--format=json"])
    _run_json(command)


def _verify_annotation_story(
    aie: Path, story: dict[str, Any], resources: dict[str, Path]
) -> dict[str, Any]:
    command = [
        str(aie),
        "compare-annotations",
        str(resources[story["archive"]]),
        "--annotation-a",
        str(resources[story["annotation_before"]]),
        "--annotation-b",
        str(resources[story["annotation_after"]]),
        "--assembly",
        story["assembly"],
        "--annotation-a-label",
        story["before_label"],
        "--annotation-b-label",
        story["after_label"],
        "--max-molecule-witnesses",
        str(story.get("max_molecule_witnesses", 10000)),
        "--format=json",
    ]
    bundle = _run_json(command)
    delta_rows = _named_table(bundle, "count_deltas")
    if not delta_rows:
        raise CapsuleError("annotation story produces no count deltas")
    gene_id = story["expected_gene_id"]
    gene_rows = [row for row in delta_rows if row.get("comparison_gene_id") == gene_id]
    if not gene_rows:
        raise CapsuleError(f"annotation story has no count delta for locked gene {gene_id}")
    try:
        signed_delta = sum(row["signed_delta_b_minus_a"] for row in gene_rows)
    except (KeyError, TypeError) as error:
        raise CapsuleError(
            "annotation count_deltas lacks the locked signed-delta fields"
        ) from error
    if signed_delta < story["expected_min_signed_delta"]:
        raise CapsuleError(
            f"annotation story {gene_id} signed UMI delta is {signed_delta}, below "
            f"{story['expected_min_signed_delta']}"
        )
    if "expected_changed_molecule_records" in story:
        transitions = _named_table(bundle, "class_transitions")
        try:
            changed_records = sum(
                row["molecule_records"]
                for row in transitions
                if gene_id
                in {
                    row.get("annotation_a_selected_comparison_gene_id"),
                    row.get("annotation_b_selected_comparison_gene_id"),
                }
            )
        except (KeyError, TypeError) as error:
            raise CapsuleError(
                "annotation class_transitions has an invalid molecule count"
            ) from error
        if changed_records != story["expected_changed_molecule_records"]:
            raise CapsuleError(
                f"annotation story {gene_id} changed molecule records are {changed_records}, "
                f"expected {story['expected_changed_molecule_records']}"
            )
    return _summary(bundle)


def _verify_event_story(
    aie: Path,
    story: dict[str, Any],
    resources: dict[str, Path],
    temporary: Path,
) -> dict[str, Any]:
    collection = temporary / "event-discovery.aicollection"
    _build_collection(aie, collection, story["archives"], resources, story)

    def run_events(
        annotation_resource: str,
        annotation_label: str,
        *,
        novel_only: bool,
        annotation_digest: str | None = None,
    ) -> dict[str, Any]:
        command = [str(aie), "collection", "find-events", str(collection)]
        for kind in story.get("kinds", []):
            command.extend(["--kind", kind])
        command.extend(
            [
                "--design",
                str(resources[story["design"]]),
                "--groups",
                str(resources[story["collection_groups"]]),
                "--annotation",
                str(resources[annotation_resource]),
            ]
        )
        for group in story.get("require_groups", []):
            command.extend(["--require-group", group])
        for option, key, default in (
            ("--min-group-umi-classes", "min_group_umi_classes", 1),
            ("--min-donors", "min_donors", 1),
            ("--min-samples", "min_samples", 1),
            ("--min-umi-classes", "min_umi_classes", 1),
            ("--min-side-umi-classes", "min_side_umi_classes", 1),
            ("--min-support", "min_support", 2),
            ("--terminal-cluster-bp", "terminal_cluster_bp", 25),
            ("--max-terminal-events", "max_terminal_events", 10000000),
            ("--max-candidates", "max_candidates", 100000),
            ("--max-candidates-considered", "max_candidates_considered", 1000000),
            ("--max-routed-entries", "max_routed_entries", 10000000),
            ("--max-exact-match-attempts", "max_exact_match_attempts", 25000000),
            ("--max-annotation-comparisons", "max_annotation_comparisons", 10000000),
        ):
            command.extend([option, str(story.get(key, default))])
        command.extend(
            [
                "--assembly",
                story["assembly"],
                "--annotation-label",
                annotation_label,
            ]
        )
        if annotation_digest is not None:
            command.extend(["--annotation-digest", annotation_digest])
        if "solo_strand" in story:
            command.extend(["--solo-strand", str(story["solo_strand"])])
        if novel_only:
            command.append("--novel-only")
        command.append("--format=json")
        return _run_json(command)

    bundle = run_events(
        story["annotation"],
        story["annotation_label"],
        novel_only=story.get("novel_only", False),
        annotation_digest=story.get("annotation_digest"),
    )
    tables = _uniform_tables(bundle)
    if not EXPECTED_TABLES["find-events"].issubset(tables):
        raise CapsuleError("event story is missing a frozen typed table")
    entities = tables["entities"]
    if not entities:
        raise CapsuleError("event story retains no entity")
    if "terminal-tail" in story.get("kinds", []) and not any(
        row.get("kind") == "terminal_tail" for row in tables["entities"]
    ):
        raise CapsuleError("event story requests terminal tails but retains none")
    expected_id = story["expected_entity_id"]
    matching = [row for row in entities if row.get("entity_id") == expected_id]
    if len(matching) != 1:
        raise CapsuleError(f"event story does not contain exactly one locked entity {expected_id}")
    entity = matching[0]
    rank = entities.index(entity) + 1
    if rank != story["expected_rank"]:
        raise CapsuleError(
            f"locked event {expected_id} has rank {rank}, "
            f"expected {story['expected_rank']}"
        )
    for field, expectation in (
        ("exact_umi_classes", "expected_min_exact_umi_classes"),
        ("exact_donors", "expected_min_exact_donors"),
    ):
        value = entity.get(field)
        if type(value) is not int or value < story[expectation]:
            raise CapsuleError(
                f"locked event {expected_id} {field} is {value}, below {story[expectation]}"
            )
    if entity.get("gap_primary_class") != story["expected_gap_primary_class"]:
        raise CapsuleError(f"locked event {expected_id} has a different gap primary class")
    if entity.get("annotation_incompatible") is not story["expected_annotation_incompatible"]:
        raise CapsuleError(f"locked event {expected_id} changed annotation compatibility")

    comparison = run_events(
        story["comparison_annotation"],
        story["comparison_annotation_label"],
        novel_only=False,
    )
    comparison_tables = _uniform_tables(comparison)
    if not EXPECTED_TABLES["find-events"].issubset(comparison_tables):
        raise CapsuleError("comparison-annotation event story is missing a frozen typed table")
    comparison_matches = [
        row for row in comparison_tables["entities"] if row.get("entity_id") == expected_id
    ]
    if len(comparison_matches) != 1:
        raise CapsuleError(
            f"comparison annotation does not contain exactly one locked event {expected_id}"
        )
    comparison_entity = comparison_matches[0]
    if comparison_entity.get("annotation_incompatible") is not False or (
        comparison_entity.get("compatible_transcripts")
        != story["expected_comparison_compatible_transcripts"]
    ):
        raise CapsuleError(
            f"locked event {expected_id} changed comparison-annotation compatibility"
        )
    return {
        "discovery": _summary(bundle),
        "comparison_annotation": _summary(comparison),
    }


def _verify_drilldown_story(
    aie: Path,
    story: dict[str, Any],
    resources: dict[str, Path],
    temporary: Path,
) -> dict[str, Any]:
    collection = temporary / "junction-drilldown.aicollection"
    _build_collection(aie, collection, story["archives"], resources, story)
    federated = _run_json(
        [
            str(aie),
            "collection",
            "junction",
            str(collection),
            story["junction"],
            "--top=0",
            "--format=json",
        ]
    )
    federated_tables = _uniform_tables(federated)
    if not any(
        row.get("present") and row.get("umis", 0) > 0
        for row in federated_tables.get("samples", [])
    ):
        raise CapsuleError("junction drilldown has no supporting source archive")
    archive_resource = story["archives"][story["archive_sample"]]
    command = [str(aie), "query", str(resources[archive_resource]), "cooccur"]
    for name, predicate in sorted(story["predicates"].items()):
        command.extend(["--predicate", f"{name}={predicate}"])
    command.extend(["--where", story["expression"], "--universe", story["universe"]])
    for option, key, default in (
        ("--unit", "unit", "molecule-record"),
        ("--region-match", "region_match", "anchor"),
        ("--placements", "placements", "unique"),
        ("--agg", "aggregation", "cell"),
        ("--max-memberships", "max_memberships", 1000000),
        ("--max-pattern-rows", "max_pattern_rows", 1000000),
        ("--max-chunks", "max_chunks", 100000),
        ("--max-evidence-records", "max_evidence_records", 10000000),
        ("--max-terminal-events", "max_terminal_events", 10000000),
    ):
        command.extend([option, str(story.get(key, default))])
    if "drilldown_groups" in story:
        command.extend(["--groups", str(resources[story["drilldown_groups"]])])
    if story.get("allow_full_scan", False):
        command.append("--allow-full-scan")
    if story.get("emit_membership", False):
        command.append("--emit-membership")
    command.append("--format=json")
    cooccur = _run_json(command)
    tables = _uniform_tables(cooccur)
    if not EXPECTED_TABLES["cooccur"].issubset(tables):
        raise CapsuleError("co-occurrence story is missing a frozen typed table")
    cooccur_summary = _summary(cooccur)
    selected_units = cooccur_summary.get("selected_units")
    if (
        type(selected_units) is not int
        or selected_units < story["expected_min_selected_units"]
    ):
        raise CapsuleError(
            "co-occurrence story selected-unit count is below its locked minimum"
        )
    required_true = set(story["required_true_predicates"])
    required_pattern = False
    for row in tables["patterns"]:
        matched = row.get("matched_predicates")
        if (
            row.get("selection_state") == "true"
            and isinstance(matched, list)
            and all(isinstance(name, str) for name in matched)
            and required_true.issubset(matched)
            and isinstance(row.get("evidence_units"), int)
            and row["evidence_units"] > 0
        ):
            required_pattern = True
            break
    if not required_pattern:
        raise CapsuleError("co-occurrence story lacks its required true predicate pattern")
    memberships = tables["memberships"]
    if not memberships:
        raise CapsuleError("co-occurrence story emits no molecule membership witnesses")
    if any(not MEMBERSHIP_COLUMNS.issubset(row) for row in memberships):
        raise CapsuleError("co-occurrence memberships lack a frozen witness column")
    unit_ids = [row["unit_id"] for row in memberships]
    if any(not isinstance(unit_id, str) or not unit_id for unit_id in unit_ids) or len(
        unit_ids
    ) != len(set(unit_ids)):
        raise CapsuleError("co-occurrence memberships contain invalid or duplicate unit IDs")
    selected_witnesses = []
    for row in memberships:
        matched = row["matched_predicates"]
        if not isinstance(matched, list) or any(not isinstance(name, str) for name in matched):
            raise CapsuleError("co-occurrence membership predicates must be a string array")
        if not isinstance(row["pattern_mask"], str) or not isinstance(
            row["completeness_mask"], str
        ):
            raise CapsuleError("co-occurrence membership masks must be strings")
        if row["selection_state"] == "true" and row["selected"] is True:
            selected_witnesses.append(row)
    if len(selected_witnesses) != selected_units:
        raise CapsuleError(
            "co-occurrence selected-unit summary differs from membership witnesses"
        )
    if not any(
        required_true.issubset(row["matched_predicates"]) for row in selected_witnesses
    ):
        raise CapsuleError(
            "co-occurrence memberships lack a selected required-predicate witness"
        )
    return {
        "junction": _summary(federated),
        "cooccurrence": cooccur_summary,
    }


def verify(
    directory: Path,
    aie: Path,
    aie_asset: Path,
    python_wheel: Path,
    *,
    execute_stories: bool = True,
) -> dict[str, Any]:
    directory = directory.resolve(strict=True)
    if any(path.suffix == ".aicollection" for path in directory.iterdir()):
        raise CapsuleError("final capsule contains a path-bound .aicollection")
    _verify_checksums(directory)
    _verify_public_bytes(directory)
    manifest = _validate_manifest(json.loads((directory / "demo-manifest.json").read_text()))
    resources = _resource_paths(directory, manifest)
    build_record = _verify_capsule_records(directory, manifest)
    aie = aie.resolve(strict=True)
    aie_asset = aie_asset.resolve(strict=True)
    python_wheel = python_wheel.resolve(strict=True)
    if _sha256(aie_asset) != manifest["software"]["aie"]["sha256"]:
        raise CapsuleError("local aie release asset differs from the manifest")
    if _sha256(python_wheel) != manifest["software"]["python_wheel"]["sha256"]:
        raise CapsuleError("local Python wheel differs from the manifest")
    version = manifest["software"]["version"]
    released_version, released_binary_sha256 = _binary_identity(
        aie_asset, manifest["software"]["aie"]["member"]
    )
    if released_version != f"aie {version}" or _sha256(aie) != released_binary_sha256:
        raise CapsuleError("--aie is not the executable from the pinned release asset")
    if build_record["software"].get("aie_sha256") != released_binary_sha256:
        raise CapsuleError("build record is not bound to the pinned release executable")
    try:
        observed_version = subprocess.run(
            [str(aie), "--version"], check=True, capture_output=True, text=True
        ).stdout.strip()
    except (OSError, subprocess.CalledProcessError) as error:
        raise CapsuleError("cannot execute --aie") from error
    if observed_version != f"aie {version}":
        raise CapsuleError("--aie version differs from the manifest")
    if _wheel_identity(python_wheel) != ("gravlax-client", version):
        raise CapsuleError("Python wheel identity differs from the manifest")
    archive_inspections = _verify_archives(aie, manifest, resources)
    results: dict[str, Any] = {
        "schema": "gravlax.demo-capsule-verification.v1",
        "manifest_sha256": _sha256(directory / "demo-manifest.json"),
        "archives_verified": len(archive_inspections),
        "stories_executed": execute_stories,
    }
    if execute_stories:
        with tempfile.TemporaryDirectory(prefix="gravlax-demo-verify-") as temporary_name:
            temporary = Path(temporary_name)
            stories = manifest["stories"]
            results["annotation_reinterpretation"] = _verify_annotation_story(
                aie, stories["annotation_reinterpretation"], resources
            )
            results["event_discovery"] = _verify_event_story(
                aie, stories["event_discovery"], resources, temporary
            )
            results["junction_drilldown"] = _verify_drilldown_story(
                aie, stories["junction_drilldown"], resources, temporary
            )
    return results


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--aie", required=True, type=Path)
    parser.add_argument("--aie-asset", required=True, type=Path)
    parser.add_argument("--python-wheel", required=True, type=Path)
    parser.add_argument(
        "--structure-only",
        action="store_true",
        help="verify bytes, roots, and software identities without executing the three stories",
    )
    args = parser.parse_args()
    result = verify(
        args.directory,
        args.aie,
        args.aie_asset,
        args.python_wheel,
        execute_stories=not args.structure_only,
    )
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
