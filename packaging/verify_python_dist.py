#!/usr/bin/env python3
"""Verify that the Python wheel and sdist contain the public package contract."""

from __future__ import annotations

import argparse
from email.parser import BytesParser
from pathlib import Path, PurePosixPath
import re
import tarfile
import zipfile


REQUIRED_PACKAGE_FILES = {
    "gravlax/__init__.py",
    "gravlax/client.py",
    "gravlax/models.py",
    "gravlax/results.py",
    "gravlax/mex.py",
    "gravlax/py.typed",
}


def _metadata_from_wheel(path: Path) -> tuple[set[str], bytes]:
    with zipfile.ZipFile(path) as archive:
        names = set(archive.namelist())
        metadata = [name for name in names if name.endswith(".dist-info/METADATA")]
        if len(metadata) != 1:
            raise ValueError(f"{path}: expected one wheel METADATA file, found {len(metadata)}")
        return names, archive.read(metadata[0])


def _verify_wheel(path: Path, version: str) -> None:
    names, raw_metadata = _metadata_from_wheel(path)
    missing = REQUIRED_PACKAGE_FILES - names
    if missing:
        raise ValueError(f"{path}: missing package files: {sorted(missing)}")
    if not any(".dist-info/licenses/LICENSE" in name for name in names):
        raise ValueError(f"{path}: BSD license is not included in the wheel")
    metadata = BytesParser().parsebytes(raw_metadata)
    if metadata["Name"] != "gravlax-client" or metadata["Version"] != version:
        raise ValueError(f"{path}: package name/version do not match gravlax-client {version}")
    if metadata["License-Expression"] != "BSD-3-Clause":
        raise ValueError(f"{path}: expected a BSD-3-Clause License-Expression")


def _verify_sdist(path: Path, version: str) -> None:
    expected_root = f"gravlax_client-{version}"
    with tarfile.open(path, "r:gz") as archive:
        names = {PurePosixPath(name) for name in archive.getnames()}
    required = {
        PurePosixPath(expected_root, "LICENSE"),
        PurePosixPath(expected_root, "README.md"),
        PurePosixPath(expected_root, "pyproject.toml"),
        PurePosixPath(expected_root, "src/gravlax/py.typed"),
    }
    missing = required - names
    if missing:
        raise ValueError(f"{path}: missing sdist files: {sorted(map(str, missing))}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--version", required=True)
    args = parser.parse_args()

    wheel_pattern = re.compile(rf"gravlax_client-{re.escape(args.version)}-py3-none-any[.]whl$")
    sdist_pattern = re.compile(rf"gravlax_client-{re.escape(args.version)}[.]tar[.]gz$")
    files = [path for path in args.directory.iterdir() if path.is_file()]
    wheels = [path for path in files if wheel_pattern.fullmatch(path.name)]
    sdists = [path for path in files if sdist_pattern.fullmatch(path.name)]
    if len(wheels) != 1 or len(sdists) != 1:
        raise ValueError(
            f"expected one wheel and one sdist for {args.version}; "
            f"found wheels={len(wheels)}, sdists={len(sdists)}"
        )
    _verify_wheel(wheels[0], args.version)
    _verify_sdist(sdists[0], args.version)
    print(f"verified {wheels[0].name} and {sdists[0].name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
