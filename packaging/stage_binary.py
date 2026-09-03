#!/usr/bin/env python3
"""Create a deterministic, self-contained Gravlax binary release archive."""

from __future__ import annotations

import argparse
import gzip
import os
from pathlib import Path
import shutil
import stat
import subprocess
import tarfile
import tempfile
import zipfile


SHELLS = ("bash", "zsh", "fish")
TEXT_FILES = ("LICENSE", "README.md")


def _safe_version(value: str) -> str:
    if not value or any(character not in "0123456789." for character in value):
        raise argparse.ArgumentTypeError("version must contain only digits and dots")
    if len(value.split(".")) != 3:
        raise argparse.ArgumentTypeError("version must have major.minor.patch form")
    return value


def _run_completion(binary: Path, shell: str, destination: Path) -> None:
    result = subprocess.run(
        [str(binary), "completions", shell],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    if not result.stdout:
        raise RuntimeError(f"aie emitted an empty {shell} completion script")
    destination.write_bytes(result.stdout)


def _stage(binary: Path, repository: Path, root: Path) -> None:
    root.mkdir()
    executable = root / binary.name
    shutil.copyfile(binary, executable)
    executable.chmod(executable.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

    for name in TEXT_FILES:
        source = repository / name
        if not source.is_file():
            raise FileNotFoundError(f"required release file is missing: {source}")
        shutil.copyfile(source, root / name)

    completions = root / "completions"
    completions.mkdir()
    for shell in SHELLS:
        _run_completion(binary, shell, completions / f"aie.{shell}")


def _install_no_replace(temporary: Path, destination: Path) -> None:
    """Atomically install a same-filesystem file only if its name is unused."""
    try:
        os.link(temporary, destination)
    except FileExistsError as error:
        raise FileExistsError(f"refusing to overwrite release archive: {destination}") from error
    temporary.unlink()


def _tar_gz(source: Path, destination: Path, epoch: int) -> None:
    temporary = destination.with_suffix(destination.suffix + ".tmp")
    with temporary.open("xb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=epoch) as compressed:
            with tarfile.open(fileobj=compressed, mode="w", format=tarfile.PAX_FORMAT) as archive:
                for path in sorted(source.rglob("*"), key=lambda item: item.as_posix()):
                    relative = path.relative_to(source.parent)
                    info = archive.gettarinfo(str(path), arcname=relative.as_posix())
                    info.uid = 0
                    info.gid = 0
                    info.uname = "root"
                    info.gname = "root"
                    info.mtime = epoch
                    if path.is_file():
                        with path.open("rb") as handle:
                            archive.addfile(info, handle)
                    else:
                        archive.addfile(info)
    _install_no_replace(temporary, destination)


def _zip(source: Path, destination: Path) -> None:
    temporary = destination.with_suffix(destination.suffix + ".tmp")
    with zipfile.ZipFile(temporary, "x", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        for path in sorted(source.rglob("*"), key=lambda item: item.as_posix()):
            if path.is_dir():
                continue
            relative = path.relative_to(source.parent).as_posix()
            info = zipfile.ZipInfo(relative, date_time=(1980, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = (path.stat().st_mode & 0xFFFF) << 16
            with path.open("rb") as handle:
                archive.writestr(info, handle.read(), compress_type=zipfile.ZIP_DEFLATED, compresslevel=9)
    _install_no_replace(temporary, destination)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--target", required=True)
    parser.add_argument("--version", required=True, type=_safe_version)
    parser.add_argument("--format", required=True, choices=("tar.gz", "zip"))
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--repository", default=Path.cwd(), type=Path)
    args = parser.parse_args()

    binary = args.binary.resolve(strict=True)
    repository = args.repository.resolve(strict=True)
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    stem = f"gravlax-{args.version}-{args.target}"
    destination = output_dir / f"{stem}.{args.format}"
    if destination.exists():
        raise FileExistsError(f"refusing to overwrite release archive: {destination}")

    epoch = int(os.environ.get("SOURCE_DATE_EPOCH", "315532800"))
    # Tar and ZIP formats cannot represent dates before 1980 portably.
    epoch = max(epoch, 315532800)
    with tempfile.TemporaryDirectory(prefix="gravlax-release-") as temporary:
        root = Path(temporary) / stem
        _stage(binary, repository, root)
        if args.format == "tar.gz":
            _tar_gz(root, destination, epoch)
        else:
            _zip(root, destination)

    print(destination)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
