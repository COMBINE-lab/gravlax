#!/usr/bin/env python3
"""Verify the exact artifact set before an immutable GitHub release is created."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re


TARGETS = (
    ("aarch64-apple-darwin", ".tar.gz"),
    ("x86_64-apple-darwin", ".tar.gz"),
    ("x86_64-pc-windows-msvc", ".zip"),
    ("x86_64-unknown-linux-gnu", ".tar.gz"),
    ("x86_64-unknown-linux-musl", ".tar.gz"),
)
CHECKSUM_LINE = re.compile(r"([0-9a-f]{64}) [ *](\S+)")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def _read_checksum_manifest(path: Path) -> dict[str, str]:
    entries: dict[str, str] = {}
    for line in path.read_text(encoding="ascii").splitlines():
        match = CHECKSUM_LINE.fullmatch(line)
        if match is None:
            raise ValueError(f"invalid checksum line in {path.name}: {line!r}")
        digest, name = match.groups()
        if name in entries:
            raise ValueError(f"duplicate checksum entry in {path.name}: {name}")
        entries[name] = digest
    if not entries:
        raise ValueError(f"empty checksum manifest: {path.name}")
    return entries


def _verify_checksums(
    directory: Path,
    manifest_name: str,
    expected_names: set[str],
) -> None:
    entries = _read_checksum_manifest(directory / manifest_name)
    if set(entries) != expected_names:
        raise ValueError(
            f"{manifest_name} entries differ; expected={sorted(expected_names)}, "
            f"found={sorted(entries)}"
        )
    for name, expected in entries.items():
        actual = _sha256(directory / name)
        if actual != expected:
            raise ValueError(
                f"checksum mismatch for {name}: expected {expected}, found {actual}"
            )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--version", required=True)
    args = parser.parse_args()

    directory = args.directory.resolve(strict=True)
    archives = {
        f"gravlax-{target}{suffix}" for target, suffix in TARGETS
    }
    archive_checksums = {f"{name}.sha256" for name in archives}
    extras = {
        f"gravlax-{args.version}-source.tar.gz",
        f"gravlax-{args.version}.spdx.json",
        f"gravlax_client-{args.version}-py3-none-any.whl",
        f"gravlax_client-{args.version}.tar.gz",
    }
    required = (
        archives
        | archive_checksums
        | extras
        | {
            "gravlax-installer.ps1",
            "gravlax-installer.sh",
            "sha256.sum",
            "SHA256SUMS",
            "dist-manifest.json",
        }
    )
    present = {path.name for path in directory.iterdir() if path.is_file()}
    missing = required - present
    unexpected = present - required
    if missing or unexpected:
        raise ValueError(
            f"complete release set mismatch; missing={sorted(missing)}, "
            f"unexpected={sorted(unexpected)}"
        )

    _verify_checksums(directory, "sha256.sum", archives)
    _verify_checksums(directory, "SHA256SUMS", extras)
    for archive in archives:
        _verify_checksums(directory, f"{archive}.sha256", {archive})

    manifest = json.loads((directory / "dist-manifest.json").read_text(encoding="utf-8"))
    if manifest.get("announcement_tag") != f"v{args.version}":
        raise ValueError("cargo-dist manifest tag does not match the release version")
    releases = manifest.get("releases")
    identities = {
        (release.get("app_name"), release.get("app_version"))
        for release in releases
        if isinstance(release, dict)
    } if isinstance(releases, list) else set()
    if identities != {("gravlax", args.version)}:
        raise ValueError(f"cargo-dist release identity mismatch: {sorted(identities)}")

    print(f"verified the complete {len(required)}-file release set")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
