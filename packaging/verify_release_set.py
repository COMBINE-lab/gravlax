#!/usr/bin/env python3
"""Check custom release extras and write their non-overwriting SHA-256 manifest."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--version", required=True)
    args = parser.parse_args()

    directory = args.directory.resolve(strict=True)
    required = {
        f"gravlax-{args.version}-source.tar.gz",
        f"gravlax-{args.version}.spdx.json",
        f"gravlax_client-{args.version}-py3-none-any.whl",
        f"gravlax_client-{args.version}.tar.gz",
    }
    present = {path.name for path in directory.iterdir() if path.is_file()}
    missing = required - present
    unexpected = present - required
    if missing or unexpected:
        raise ValueError(
            f"release set mismatch; missing={sorted(missing)}, unexpected={sorted(unexpected)}"
        )

    checksum_path = directory / "SHA256SUMS"
    if checksum_path.exists():
        raise FileExistsError(f"refusing to overwrite checksum manifest: {checksum_path}")
    lines = [f"{_sha256(directory / name)}  {name}\n" for name in sorted(required)]
    checksum_path.write_text("".join(lines), encoding="ascii", newline="\n")
    print(f"verified {len(required)} release artifacts and wrote {checksum_path.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
