#!/usr/bin/env python3
"""Bind staged demo data to immutable released software and HTTPS asset URLs."""

from __future__ import annotations

import argparse
from email.parser import Parser
import json
from pathlib import Path, PurePosixPath
import re
import shutil
import stat
import subprocess
import tarfile
import tempfile
from typing import Any
from urllib.parse import urlsplit
import zipfile

from build_demo_capsule import (
    CapsuleError,
    HEX64,
    IDENTIFIER,
    RECORD_SCHEMA,
    _assert_public_strings,
    _expect_mapping,
    _group_map_scope,
    _https,
    _new_output_directory,
    _notices,
    _safe_filename,
    _safe_identifier,
    _sha256,
    _validate_story_references,
    _write_bytes_exclusive,
    _write_json_exclusive,
)


MANIFEST_SCHEMA = "gravlax.demo-capsule.v1"
FINALIZATION_SCHEMA = "gravlax.demo-capsule-finalization.v1"


def _member_name(value: str) -> str:
    if not isinstance(value, str) or not value:
        raise CapsuleError("--aie-member must be nonempty")
    member = PurePosixPath(value)
    if member.is_absolute() or ".." in member.parts or member.name != "aie":
        raise CapsuleError("--aie-member must be a relative archive path ending in /aie")
    return value


def _extract_binary(bundle: Path, member: str, destination: Path) -> None:
    payload: bytes
    if zipfile.is_zipfile(bundle):
        with zipfile.ZipFile(bundle) as archive:
            try:
                info = archive.getinfo(member)
            except KeyError as error:
                raise CapsuleError(f"aie member is absent from {bundle.name}: {member}") from error
            member_mode = info.external_attr >> 16
            if info.is_dir() or stat.S_ISLNK(member_mode):
                raise CapsuleError("aie member is not a regular archive file")
            payload = archive.read(info)
    else:
        try:
            with tarfile.open(bundle, "r:*") as archive:
                info = archive.getmember(member)
                if not info.isfile() or info.issym() or info.islnk():
                    raise CapsuleError("aie member is not a regular archive file")
                handle = archive.extractfile(info)
                if handle is None:
                    raise CapsuleError("cannot read aie member")
                payload = handle.read()
        except tarfile.TarError as error:
            raise CapsuleError(f"unsupported aie release archive: {bundle.name}") from error
    _write_bytes_exclusive(destination, payload)
    destination.chmod(destination.stat().st_mode | stat.S_IXUSR)


def _binary_identity(bundle: Path, member: str) -> tuple[str, str]:
    with tempfile.TemporaryDirectory(prefix="gravlax-demo-aie-") as temporary:
        binary = Path(temporary) / "aie"
        _extract_binary(bundle, member, binary)
        try:
            result = subprocess.run(
                [str(binary), "--version"],
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
        except (OSError, subprocess.CalledProcessError) as error:
            raise CapsuleError("released aie member cannot execute on this host") from error
        return result.stdout.strip(), _sha256(binary)


def _wheel_identity(wheel: Path) -> tuple[str, str]:
    if not zipfile.is_zipfile(wheel):
        raise CapsuleError("--python-wheel is not a wheel ZIP archive")
    with zipfile.ZipFile(wheel) as archive:
        metadata_names = [
            name for name in archive.namelist() if name.endswith(".dist-info/METADATA")
        ]
        if len(metadata_names) != 1:
            raise CapsuleError("Python wheel must contain exactly one METADATA file")
        metadata = Parser().parsestr(archive.read(metadata_names[0]).decode("utf-8", "strict"))
    name, version = metadata.get("Name"), metadata.get("Version")
    if not name or not version:
        raise CapsuleError("Python wheel METADATA lacks Name or Version")
    return name, version


def _copy_no_replace(source: Path, destination: Path) -> None:
    if destination.exists():
        raise FileExistsError(f"refusing to replace capsule file: {destination}")
    with source.open("rb") as reader, destination.open("xb") as writer:
        shutil.copyfileobj(reader, writer, 8 << 20)


def _immutable_url(
    value: Any,
    label: str,
    *,
    filename: str | None = None,
    kind: str = "github-release",
) -> str:
    value = _https(value, label)
    parsed = urlsplit(value)
    error = f"{label} must be an immutable, versioned {kind} HTTPS URL"
    if (
        parsed.scheme != "https"
        or parsed.username is not None
        or parsed.password is not None
        or parsed.port is not None
        or parsed.query
        or parsed.fragment
    ):
        raise CapsuleError(error)
    parts = tuple(part for part in PurePosixPath(parsed.path).parts if part != "/")
    if kind == "github-release":
        if (
            parsed.netloc.casefold() != "github.com"
            or parts[:4] != ("COMBINE-lab", "gravlax", "releases", "download")
            or len(parts) != (6 if filename is not None else 5)
            or not IDENTIFIER.fullmatch(parts[4])
            or parts[4].casefold() == "latest"
        ):
            raise CapsuleError(error)
    elif kind == "pypi-file":
        if (
            parsed.netloc.casefold() != "files.pythonhosted.org"
            or filename is None
            or len(parts) < 5
            or parts[0] != "packages"
        ):
            raise CapsuleError(error)
    else:
        raise AssertionError(f"unknown release URL kind: {kind}")
    if filename is not None and parts[-1] != filename:
        raise CapsuleError(f"{label} basename must be {filename}")
    return value


def _github_release_tag(url: str) -> str:
    """Return the fixed GitHub release tag from a URL already validated above."""
    return tuple(part for part in PurePosixPath(urlsplit(url).path).parts if part != "/")[4]


def _verify_software_tag(
    source_repository: Path, release_tag: str, expected_revision: str
) -> None:
    source_repository = source_repository.resolve(strict=True)
    try:
        result = subprocess.run(
            [
                "git",
                "-C",
                str(source_repository),
                "rev-parse",
                "--verify",
                f"refs/tags/{release_tag}^{{commit}}",
            ],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
    except (OSError, subprocess.CalledProcessError) as error:
        raise CapsuleError(
            f"cannot resolve software release tag {release_tag} in --source-repository"
        ) from error
    observed = result.stdout.strip()
    if observed != expected_revision:
        raise CapsuleError(
            f"software release tag {release_tag} peels to {observed}, not {expected_revision}"
        )


def _load_record(build_dir: Path) -> dict[str, Any]:
    path = build_dir / "BUILD-RECORD.json"
    try:
        record = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise CapsuleError("build directory lacks a valid BUILD-RECORD.json") from error
    if not isinstance(record, dict) or record.get("schema") != RECORD_SCHEMA:
        raise CapsuleError(f"build record schema must be {RECORD_SCHEMA}")
    _assert_public_strings(record, "build record")
    return record


def _validate_local_resources(build_dir: Path, record: dict[str, Any]) -> dict[str, dict[str, Any]]:
    resources = record.get("resources")
    if not isinstance(resources, dict) or not resources:
        raise CapsuleError("build record resources must be a nonempty object")
    seen: set[str] = set()
    validated: dict[str, dict[str, Any]] = {}
    for name, raw in resources.items():
        if not isinstance(name, str) or not isinstance(raw, dict):
            raise CapsuleError("invalid build-record resource")
        _safe_identifier(name, "build-record resource name")
        if set(raw).difference({"filename", "sha256", "bytes", "archive_root"}):
            raise CapsuleError(f"resource {name} has an unknown field")
        filename = _safe_filename(raw.get("filename"))
        if filename.endswith(".aicollection"):
            raise CapsuleError("path-bound .aicollection files must not be finalized")
        if filename in seen:
            raise CapsuleError(f"two resources use filename {filename}")
        seen.add(filename)
        digest = raw.get("sha256")
        if not isinstance(digest, str) or not HEX64.fullmatch(digest):
            raise CapsuleError(f"resource {name} has an invalid SHA-256")
        path = build_dir / filename
        if not path.is_file() or path.is_symlink():
            raise CapsuleError(f"resource {name} is missing or not a regular file")
        if path.stat().st_size != raw.get("bytes") or _sha256(path) != digest:
            raise CapsuleError(f"resource {name} changed after staging")
        root = raw.get("archive_root")
        if filename.endswith(".aie") and root is None:
            raise CapsuleError(f"archive resource {name} lacks an archive root")
        if root is not None and not re.fullmatch(r"aie-directory-root-v2:[0-9a-f]{64}", root):
            raise CapsuleError(f"resource {name} has an invalid archive root")
        validated[name] = raw
    return validated


def _manifest_resources(
    resources: dict[str, dict[str, Any]], data_base_url: str
) -> dict[str, dict[str, Any]]:
    result = {}
    for name, resource in sorted(resources.items()):
        asset = {
            "url": f"{data_base_url.rstrip('/')}/{resource['filename']}",
            "sha256": resource["sha256"],
            "filename": resource["filename"],
        }
        if "archive_root" in resource:
            asset["archive_root"] = resource["archive_root"]
        result[name] = asset
    return result


def _documentation(record: dict[str, Any], manifest: dict[str, Any]) -> tuple[str, str]:
    capsule_id = record.get("capsule_id")
    summary = record.get("summary")
    assembly = record.get("assembly")
    selection = record.get("selection")
    inputs = record.get("inputs")
    if not isinstance(selection, dict) or not isinstance(inputs, dict):
        raise CapsuleError("build record lacks selection or input provenance")
    windows = selection.get("windows")
    donors = inputs.get("donors")
    if not isinstance(capsule_id, str) or not isinstance(summary, str):
        raise CapsuleError("build record lacks capsule identity or summary")
    if not isinstance(assembly, str) or not isinstance(windows, list) or not windows:
        raise CapsuleError("build record lacks assembly or selection windows")
    if not isinstance(donors, list) or not donors:
        raise CapsuleError("build record lacks donor provenance")
    group_map_scope = _group_map_scope(record.get("group_map_scope"))
    window_lines = []
    for window in windows:
        if not isinstance(window, dict):
            raise CapsuleError("build record contains an invalid selection window")
        try:
            window_lines.append(f"- `{window['chrom']}:{window['start']}-{window['end']}`")
        except KeyError as error:
            raise CapsuleError("build record contains an incomplete selection window") from error
    donor_lines = []
    donor_ids: set[str] = set()
    for donor in donors:
        if not isinstance(donor, dict):
            raise CapsuleError("build record contains invalid donor provenance")
        try:
            accessions = ", ".join(donor["accessions"])
            donor_ids.add(donor["donor"])
            donor_lines.append(f"| {donor['sample']} | {donor['donor']} | {accessions} |")
        except (KeyError, TypeError) as error:
            raise CapsuleError("build record contains incomplete donor provenance") from error

    version = manifest["software"]["version"]
    aie_filename = manifest["software"]["aie"]["filename"]
    wheel_filename = manifest["software"]["python_wheel"]["filename"]
    annotation_story = manifest["stories"]["annotation_reinterpretation"]
    event_story = manifest["stories"]["event_discovery"]
    drilldown_note = manifest["stories"]["junction_drilldown"]["story_note"]
    readme = "\n".join(
        [
            f"# Gravlax demo capsule: {capsule_id}",
            "",
            summary,
            "",
            f"This immutable capsule was built with Gravlax {version} for `{assembly}`.",
            "Its `.aie` archives contain molecule-level genome-coordinate evidence and sparse ",
            "terminal-tail events. They contain no read sequences, base qualities, or read names.",
            "Selected alignment records were coordinate-sorted before archive ingest.",
            "Molecule/read-derived evidence is deliberately restricted to the curated windows ",
            "below; placement multiplicity and completeness outside those windows are not ",
            "represented. Each archive also embeds its donor's complete genome-wide aggregate ",
            "STAR pass-1 junction catalogue as root-bound alignment provenance. That catalogue ",
            "contains junction coordinates and aggregate support columns, but no sequences, cell ",
            "barcodes, or UMI identities.",
            "Each donor's exact STAR `Log.out` confirms two-pass mapping with ",
            "`sjdbGTFfile=-` and `sjdbFileChrStartEnd=-`; no external GTF or caller-supplied ",
            "junction list was provided, while STAR Basic two-pass used its own ",
            "pass-1 catalogue. The STAR index identity and component ",
            "hash list are caller-declared provenance because the builder does not consume the ",
            "original index files.",
            "They are curated event demonstrations, not complete donor archives or a genome-wide ",
            "atlas.",
            "",
            "## Curated windows",
            "",
            *window_lines,
            "",
            "## Public donor inputs",
            "",
            f"The capsule has {len(donors)} sample archives representing "
            f"{len(donor_ids)} distinct donor identifiers.",
            "",
            "| Sample | Donor | Accessions |",
            "| --- | --- | --- |",
            *donor_lines,
            "",
            "## Cell-group-map scope",
            "",
            group_map_scope["label_derivation"],
            group_map_scope["included_cells"],
            group_map_scope["event_semantics"],
            group_map_scope["drilldown_semantics"],
            "",
            "## Demonstrations",
            "",
            "1. Reinterpret one molecule-evidence archive using two annotation releases.",
            f"   Locked check: `{annotation_story['expected_gene_id']}` gains at least "
            f"{annotation_story['expected_min_signed_delta']:,} signed UMIs.",
            "2. Find a recurrent cassette-splicing event across donors and cell groups, then "
            "test its compatibility with the later annotation.",
            f"   Locked event: `{event_story['expected_entity_id']}` has at least "
            f"{event_story['expected_min_exact_umi_classes']} exact UMI classes across "
            f"{event_story['expected_min_exact_donors']} donors, and exactly "
            f"{event_story['expected_comparison_compatible_transcripts']} transcripts in "
            f"{event_story['comparison_annotation_label']} are compatible.",
            "3. Query a junction across the collection, then test biological predicates on the "
            "same molecule record.",
            "",
            f"**Drilldown interpretation:** {drilldown_note}",
            "",
            "Collections are intentionally not distributed because they contain local source "
            "paths. ",
            "The notebooks and verifier rebuild each `.aicollection` after download.",
            "",
            "## Verify a download",
            "",
            "Download this release's complete data assets. Separately download the official ",
            f"release archive `{aie_filename}` and wheel `{wheel_filename}` from the exact URLs ",
            "in `demo-manifest.json`, then extract the manifest's `aie` member. Run:",
            "",
            "```sh",
            "python <gravlax-checkout>/packaging/verify_demo_capsule.py . \\",
            "  --aie <extracted-release>/aie \\",
            f"  --aie-asset <release-downloads>/{aie_filename} \\",
            f"  --python-wheel <release-downloads>/{wheel_filename}",
            "```",
            "",
            "The verifier checks every SHA-256, each archive content root and provenance ",
            "manifest, both software identities, and live outputs from all three stories.",
            "`SHA256SUMS` covers every file in this flat capsule except `SHA256SUMS` itself. ",
            "`BUILD-RECORD.json` records source identities and selection details; ",
            "`THIRD-PARTY-NOTICES.md` records source terms and citations.",
            "",
        ]
    )
    release_notes = "\n".join(
        [
            f"# {capsule_id}",
            "",
            summary,
            "",
            f"Built for `{assembly}` with Gravlax {version}. The capsule includes {len(donors)} ",
            f"locus-restricted donor archives over {len(windows)} curated windows, annotation ",
            "subsets, cell-group and donor designs, recorded hash-pinned provenance, and ",
            "transport checksums.",
            group_map_scope["label_derivation"],
            group_map_scope["included_cells"],
            group_map_scope["event_semantics"],
            group_map_scope["drilldown_semantics"],
            "It exercises annotation reinterpretation, multi-donor cassette-event discovery, and ",
            "same-molecule predicate queries. See `README.md` for scope and verification.",
            f"Drilldown interpretation: {drilldown_note}",
            "",
        ]
    )
    _assert_public_strings(readme, "generated README")
    _assert_public_strings(release_notes, "generated release notes")
    return readme, release_notes


def finalize(
    *,
    build_dir: Path,
    output_dir: Path,
    data_base_url: str,
    aie_asset: Path,
    aie_url: str,
    aie_member: str,
    python_wheel: Path,
    python_wheel_url: str,
    source_repository: Path,
) -> str:
    build_dir = build_dir.resolve(strict=True)
    record = _load_record(build_dir)
    expected_notices, _ = _notices(record.get("third_party_notices"))
    notices_path = build_dir / "THIRD-PARTY-NOTICES.md"
    if not notices_path.is_file() or notices_path.read_text(encoding="utf-8") != expected_notices:
        raise CapsuleError("THIRD-PARTY-NOTICES.md differs from the build record")
    resources = _validate_local_resources(build_dir, record)
    stories = _expect_mapping(record.get("stories"), "build record stories")
    _validate_story_references(stories, set(resources))
    build_software = _expect_mapping(record.get("software"), "build record software")
    version = build_software.get("version")
    if not isinstance(version, str) or not version:
        raise CapsuleError("build record has no software version")
    revision = build_software.get("revision")
    if not isinstance(revision, str) or not re.fullmatch(r"[0-9a-f]{40}", revision):
        raise CapsuleError("build record has no full software source revision")
    data_base_url = _immutable_url(data_base_url, "--data-base-url").rstrip("/")
    aie_member = _member_name(aie_member)
    aie_asset = aie_asset.resolve(strict=True)
    python_wheel = python_wheel.resolve(strict=True)
    aie_filename = _safe_filename(aie_asset.name)
    wheel_filename = _safe_filename(python_wheel.name)
    aie_url = _immutable_url(aie_url, "--aie-url", filename=aie_filename)
    software_release_tag = _github_release_tag(aie_url)
    if software_release_tag != f"v{version}":
        raise CapsuleError("--aie-url release tag must equal v<build-record version>")
    _verify_software_tag(source_repository, software_release_tag, revision)
    python_wheel_url = _immutable_url(
        python_wheel_url,
        "--python-wheel-url",
        filename=wheel_filename,
        kind="pypi-file",
    )
    binary_version, binary_sha256 = _binary_identity(aie_asset, aie_member)
    if binary_version != f"aie {version}":
        raise CapsuleError("released aie asset version does not match the build record")
    if binary_sha256 != build_software.get("aie_sha256"):
        raise CapsuleError("demo data was not built with the released aie executable bytes")
    wheel_name, wheel_version = _wheel_identity(python_wheel)
    if wheel_name != "gravlax-client" or wheel_version != version:
        raise CapsuleError("released Python wheel identity does not match gravlax-client/version")

    with _new_output_directory(output_dir) as output:
        copied = sorted(
            path for path in build_dir.iterdir() if path.is_file() and not path.is_symlink()
        )
        if {path.name for path in copied} != {
            "BUILD-RECORD.json",
            "THIRD-PARTY-NOTICES.md",
            *(resource["filename"] for resource in resources.values()),
        }:
            raise CapsuleError("build directory contains missing or unexpected top-level files")
        for source in copied:
            _copy_no_replace(source, output / _safe_filename(source.name))

        manifest = {
            "schema": MANIFEST_SCHEMA,
            "software": {
                "version": version,
                "aie": {
                    "url": aie_url,
                    "sha256": _sha256(aie_asset),
                    "filename": aie_filename,
                    "member": aie_member,
                },
                "python_wheel": {
                    "url": python_wheel_url,
                    "sha256": _sha256(python_wheel),
                    "filename": wheel_filename,
                },
            },
            "resources": _manifest_resources(resources, data_base_url),
            "stories": stories,
        }
        _assert_public_strings(manifest, "demo manifest")
        manifest_path = output / "demo-manifest.json"
        _write_json_exclusive(manifest_path, manifest)
        manifest_sha256 = _sha256(manifest_path)
        finalization = {
            "schema": FINALIZATION_SCHEMA,
            "capsule_id": record["capsule_id"],
            "manifest": {
                "filename": manifest_path.name,
                "sha256": manifest_sha256,
                "url": f"{data_base_url}/{manifest_path.name}",
            },
            "software": manifest["software"],
            "source": {
                "revision": revision,
                "software_release_tag": software_release_tag,
            },
            "data_base_url": data_base_url,
        }
        _write_json_exclusive(output / "FINALIZATION-RECORD.json", finalization)
        readme, release_notes = _documentation(record, manifest)
        _write_bytes_exclusive(output / "README.md", readme.encode("utf-8"))
        _write_bytes_exclusive(output / "RELEASE-NOTES.md", release_notes.encode("utf-8"))
        checksum_rows = []
        for path in sorted(output.iterdir(), key=lambda item: item.name):
            if not path.is_file() or path.is_symlink():
                raise CapsuleError(f"unexpected non-file in finalized capsule: {path.name}")
            checksum_rows.append(f"{_sha256(path)}  {path.name}\n")
        _write_bytes_exclusive(output / "SHA256SUMS", "".join(checksum_rows).encode("ascii"))
    print(json.dumps({"directory": str(output_dir.resolve()), "manifest_sha256": manifest_sha256}))
    return manifest_sha256


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--data-base-url", required=True)
    parser.add_argument("--aie-asset", required=True, type=Path)
    parser.add_argument("--aie-url", required=True)
    parser.add_argument("--aie-member", required=True)
    parser.add_argument("--python-wheel", required=True, type=Path)
    parser.add_argument("--python-wheel-url", required=True)
    parser.add_argument("--source-repository", required=True, type=Path)
    args = parser.parse_args()
    finalize(
        build_dir=args.build_dir,
        output_dir=args.output_dir,
        data_base_url=args.data_base_url,
        aie_asset=args.aie_asset,
        aie_url=args.aie_url,
        aie_member=args.aie_member,
        python_wheel=args.python_wheel,
        python_wheel_url=args.python_wheel_url,
        source_repository=args.source_repository,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
