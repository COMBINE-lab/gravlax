#!/usr/bin/env python3
"""Validate the single release version shared by every Gravlax package."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import tomllib


SEMVER = re.compile(r"[0-9]+[.][0-9]+[.][0-9]+")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tag", required=True)
    parser.add_argument("--repository", type=Path, default=Path.cwd())
    args = parser.parse_args()

    repository = args.repository.resolve(strict=True)
    cargo = tomllib.loads((repository / "Cargo.toml").read_text(encoding="utf-8"))
    python = tomllib.loads((repository / "python/pyproject.toml").read_text(encoding="utf-8"))
    python_init = (repository / "python/src/gravlax/__init__.py").read_text(encoding="utf-8")
    python_init_match = re.search(r'^__version__\s*=\s*"([^"\n]+)"\s*$', python_init, re.MULTILINE)
    if python_init_match is None:
        raise ValueError("could not find the Python runtime version")
    docs = json.loads((repository / "docs/package.json").read_text(encoding="utf-8"))
    docs_lock = json.loads((repository / "docs/package-lock.json").read_text(encoding="utf-8"))
    local_recipe = (repository / "packaging/conda/local/meta.yaml").read_text(encoding="utf-8")
    recipe_match = re.search(
        r'environ[.]get[(]"GRAVLAX_VERSION", "([^"\n]+)"[)]', local_recipe
    )
    if recipe_match is None:
        raise ValueError("could not find the local conda fallback version")
    versions = {
        "Rust workspace": cargo["workspace"]["package"]["version"],
        "Python package": python["project"]["version"],
        "Python runtime": python_init_match.group(1),
        "documentation package": docs["version"],
        "documentation lockfile": docs_lock["version"],
        "documentation lockfile root": docs_lock["packages"][""]["version"],
        "local conda recipe": recipe_match.group(1),
    }
    unique = set(versions.values())
    if len(unique) != 1:
        raise ValueError(f"release versions disagree: {versions}")
    version = unique.pop()
    if not SEMVER.fullmatch(version):
        raise ValueError(f"release version is not stable major.minor.patch: {version}")
    expected_tag = f"v{version}"
    if args.tag != expected_tag:
        raise ValueError(
            f"release tag {args.tag!r} does not match package version; expected {expected_tag!r}"
        )
    print(version)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
