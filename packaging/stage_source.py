#!/usr/bin/env python3
"""Create an offline-buildable, deterministic Gravlax source release asset."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import subprocess
import tarfile
import tempfile

from stage_binary import _safe_version, _tar_gz


def _export_revision(repository: Path, revision: str, destination: Path) -> None:
    archive_path = destination.parent / "tracked-source.tar"
    subprocess.run(
        [
            "git",
            "-C",
            str(repository),
            "archive",
            "--format=tar",
            f"--output={archive_path}",
            revision,
        ],
        check=True,
    )
    destination.mkdir()
    with tarfile.open(archive_path, "r:") as archive:
        archive.extractall(destination, filter="data")


def _vendor_dependencies(source: Path) -> None:
    cargo_dir = source / ".cargo"
    cargo_dir.mkdir(exist_ok=True)
    config_path = cargo_dir / "config.toml"
    if config_path.exists() or (source / "vendor").exists():
        raise RuntimeError(
            "refusing to replace a tracked .cargo/config.toml or vendor directory"
        )
    result = subprocess.run(
        ["cargo", "vendor", "--locked", "vendor"],
        cwd=source,
        check=True,
        stdout=subprocess.PIPE,
        stderr=None,
        text=True,
    )
    config = result.stdout.rstrip() + "\n\n[net]\noffline = true\n"
    config_path.write_text(config, encoding="utf-8", newline="\n")
    # Resolve every Cargo.lock entry from the staged tree with networking
    # disabled before packaging it.
    subprocess.run(
        ["cargo", "metadata", "--locked", "--offline", "--format-version=1"],
        cwd=source,
        check=True,
        stdout=subprocess.DEVNULL,
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository", default=Path.cwd(), type=Path)
    parser.add_argument("--revision", required=True)
    parser.add_argument("--version", required=True, type=_safe_version)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    repository = args.repository.resolve(strict=True)
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    destination = output_dir / f"gravlax-{args.version}-source.tar.gz"
    if destination.exists():
        raise FileExistsError(f"refusing to overwrite source archive: {destination}")

    epoch = max(int(os.environ.get("SOURCE_DATE_EPOCH", "315532800")), 315532800)
    with tempfile.TemporaryDirectory(prefix="gravlax-source-") as temporary:
        root = Path(temporary) / f"gravlax-{args.version}"
        _export_revision(repository, args.revision, root)
        _vendor_dependencies(root)
        _tar_gz(root, destination, epoch)
    print(destination)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
