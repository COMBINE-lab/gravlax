#!/usr/bin/env python3
"""Render a Bioconda recipe from a published Gravlax source release asset."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path
import re


VERSION = re.compile(
    r"(?:0|[1-9][0-9]*)[.](?:0|[1-9][0-9]*)[.](?:0|[1-9][0-9]*)"
)
TOKENS = ("@VERSION@", "@SOURCE_SHA256@")


def _digest(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", required=True)
    parser.add_argument("--source-archive", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    if not VERSION.fullmatch(args.version):
        parser.error(
            "--version must be a stable major.minor.patch version without a leading v"
        )
    archive = args.source_archive.resolve(strict=True)
    expected_name = f"gravlax-{args.version}-source.tar.gz"
    if archive.name != expected_name:
        parser.error(
            f"--source-archive must be the release asset named {expected_name}; "
            f"received {archive.name}"
        )
    output = args.output.resolve()
    if output.name != "meta.yaml":
        parser.error("--output must name the Bioconda recipe file meta.yaml")
    if output.exists():
        raise FileExistsError(f"refusing to overwrite rendered recipe: {output}")

    template_path = Path(__file__).with_name("meta.yaml.in")
    template = template_path.read_text(encoding="utf-8")
    for token in TOKENS:
        if template.count(token) != 1:
            raise RuntimeError(
                f"Bioconda template must contain {token} exactly once"
            )
    rendered = template.replace("@VERSION@", args.version).replace(
        "@SOURCE_SHA256@", _digest(archive)
    )
    if any(token in rendered for token in TOKENS):
        raise RuntimeError("unresolved Bioconda template tokens remain")
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x", encoding="utf-8", newline="\n") as handle:
        handle.write(rendered)
    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
