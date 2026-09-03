from __future__ import annotations

import importlib.util
from importlib.machinery import SourceFileLoader
import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts/bump-and-publish"
LOADER = SourceFileLoader("gravlax_release_tool", str(SCRIPT))
SPEC = importlib.util.spec_from_loader("gravlax_release_tool", LOADER)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError(f"cannot load {SCRIPT}")
release_tool = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(release_tool)


class ReleaseToolTests(unittest.TestCase):
    def test_script_cannot_match_its_own_private_marker(self):
        marker = b"re" + b"comb"
        self.assertNotIn(marker, SCRIPT.read_bytes().lower())

    def test_marker_scan_is_case_insensitive(self):
        old_root = release_tool.ROOT
        old_paths = release_tool._tracked_paths
        try:
            with tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                clean = root / "clean.txt"
                flagged = root / "flagged.txt"
                clean.write_bytes(b"public source\n")
                flagged.write_bytes(b"former " + b"Re" + b"Comb" + b" project\n")
                release_tool.ROOT = root
                release_tool._tracked_paths = lambda: [clean, flagged]
                self.assertEqual(release_tool.forbidden_tracked_references(), ["flagged.txt"])
        finally:
            release_tool.ROOT = old_root
            release_tool._tracked_paths = old_paths

    def test_version_update_changes_every_release_surface(self):
        old_root = release_tool.ROOT
        try:
            with tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                (root / "python/src/gravlax").mkdir(parents=True)
                (root / "docs").mkdir()
                (root / "packaging/conda/local").mkdir(parents=True)
                (root / "Cargo.toml").write_text(
                    """[workspace.package]
version = "0.1.0"

[workspace.dependencies]
evidence-io = { package = "gravlax-evidence-io", version = "=0.1.0", path = "crates/evidence-io" }
anno = { package = "gravlax-anno", version = "=0.1.0", path = "crates/anno" }
ingest = { package = "gravlax-ingest", version = "=0.1.0", path = "crates/ingest" }
gravlax-output = { version = "=0.1.0", path = "crates/gravlax-output" }
""",
                    encoding="utf-8",
                )
                (root / "python/pyproject.toml").write_text(
                    '[project]\nname = "gravlax-client"\nversion = "0.1.0"\n',
                    encoding="utf-8",
                )
                (root / "python/src/gravlax/__init__.py").write_text(
                    '__version__ = "0.1.0"\n', encoding="utf-8"
                )
                (root / "docs/package.json").write_text(
                    json.dumps({"name": "gravlax-docs", "version": "0.1.0"}) + "\n",
                    encoding="utf-8",
                )
                (root / "docs/package-lock.json").write_text(
                    json.dumps(
                        {
                            "name": "gravlax-docs",
                            "version": "0.1.0",
                            "lockfileVersion": 3,
                            "packages": {
                                "": {"name": "gravlax-docs", "version": "0.1.0"}
                            },
                        },
                        indent=2,
                    )
                    + "\n",
                    encoding="utf-8",
                )
                (root / "packaging/conda/local/meta.yaml").write_text(
                    '{% set version = environ.get("GRAVLAX_VERSION", "0.1.0") %}\n',
                    encoding="utf-8",
                )
                release_tool.ROOT = root
                release_tool.update_versions("0.1.0", "0.2.0")
                self.assertEqual(set(release_tool.current_versions().values()), {"0.2.0"})
                cargo = (root / "Cargo.toml").read_text(encoding="utf-8")
                self.assertEqual(cargo.count('version = "=0.2.0"'), 4)
                docs_lock = json.loads(
                    (root / "docs/package-lock.json").read_text(encoding="utf-8")
                )
                self.assertEqual(docs_lock["version"], "0.2.0")
                self.assertEqual(docs_lock["packages"][""]["version"], "0.2.0")
        finally:
            release_tool.ROOT = old_root

    def test_version_update_validates_every_surface_before_writing(self):
        old_root = release_tool.ROOT
        try:
            with tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                (root / "python/src/gravlax").mkdir(parents=True)
                (root / "docs").mkdir()
                (root / "packaging/conda/local").mkdir(parents=True)
                fixtures = {
                    "Cargo.toml": """[workspace.package]
version = "0.1.0"

[workspace.dependencies]
evidence-io = { package = "gravlax-evidence-io", version = "=0.1.0", path = "crates/evidence-io" }
anno = { package = "gravlax-anno", version = "=0.1.0", path = "crates/anno" }
ingest = { package = "gravlax-ingest", version = "=0.1.0", path = "crates/ingest" }
gravlax-output = { version = "=0.1.0", path = "crates/gravlax-output" }
""",
                    "python/pyproject.toml": '[project]\nversion = "0.1.0"\n',
                    "python/src/gravlax/__init__.py": '__version__ = "0.1.0"\n',
                    "docs/package.json": '{"version": "0.1.0"}\n',
                    "docs/package-lock.json": json.dumps(
                        {
                            "version": "0.1.0",
                            "lockfileVersion": 3,
                            "packages": {"": {"version": "0.1.0"}},
                        },
                        indent=2,
                    )
                    + "\n",
                    "packaging/conda/local/meta.yaml": (
                        '{% set version = environ.get("GRAVLAX_VERSION", "0.1.0") %}\n'
                    ),
                }
                for relative, contents in fixtures.items():
                    (root / relative).write_text(contents, encoding="utf-8")
                before = {
                    relative: (root / relative).read_bytes() for relative in fixtures
                }
                release_tool.ROOT = root
                with self.assertRaises(release_tool.ReleaseError):
                    release_tool.update_versions("0.1.0", "0.2.0")
                self.assertEqual(
                    {relative: (root / relative).read_bytes() for relative in fixtures},
                    before,
                )
        finally:
            release_tool.ROOT = old_root

    def test_irreversible_publish_requires_a_matching_confirmation(self):
        parsed = release_tool.parse_args(
            ["1.2.3", "--publish-crates", "--confirm-publish", "v1.2.3"]
        )
        self.assertTrue(parsed.publish_crates)
        with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
            release_tool.parse_args(["1.2.3", "--dry-run", "--confirm-publish", "v1.2.3"])
        dispatched = release_tool.parse_args(["1.2.3", "--dispatch-python"])
        self.assertTrue(dispatched.dispatch_python)

    def test_release_versions_must_increase_strictly(self):
        release_tool._require_newer_version("0.1.9", "0.2.0")
        with self.assertRaises(release_tool.ReleaseError):
            release_tool._require_newer_version("0.2.0", "0.2.0")
        with self.assertRaises(release_tool.ReleaseError):
            release_tool._require_newer_version("1.0.0", "0.99.0")

    def test_current_version_dry_run_requires_a_clean_tree(self):
        versions = {"Rust workspace": "0.1.0", "Python package": "0.1.0"}
        with mock.patch.object(
            release_tool, "current_versions", return_value=versions
        ), mock.patch.object(release_tool, "_require_release_notes"), mock.patch.object(
            release_tool, "_require_clean_tree"
        ) as clean, mock.patch.object(release_tool, "validate_repository"), mock.patch.object(
            release_tool, "package_dry_runs"
        ), contextlib.redirect_stdout(io.StringIO()):
            release_tool.dry_run("0.1.0", check_history=True)
        clean.assert_called_once_with()

    def test_registry_propagation_is_polled_and_bounded(self):
        with mock.patch.object(
            release_tool, "_crate_version_exists", side_effect=[False, False, True]
        ) as exists, mock.patch.object(release_tool.time, "sleep") as sleep:
            release_tool._wait_for_crate_version(
                "gravlax-output", "0.2.0", timeout=10, poll_interval=1
            )
        self.assertEqual(exists.call_count, 3)
        self.assertEqual(sleep.call_args_list, [mock.call(1), mock.call(2)])

        with mock.patch.object(release_tool, "_crate_version_exists", return_value=False):
            with self.assertRaises(release_tool.ReleaseError):
                release_tool._wait_for_crate_version(
                    "gravlax-output", "0.2.0", timeout=0, poll_interval=0
                )

    def test_release_push_accepts_only_the_canonical_github_origin(self):
        accepted = (
            "https://github.com/COMBINE-lab/gravlax",
            "https://github.com/COMBINE-lab/gravlax.git",
            "git@github.com:COMBINE-lab/gravlax.git",
            "ssh://git@github.com/COMBINE-lab/gravlax.git",
        )
        for origin in accepted:
            with self.subTest(origin=origin), mock.patch.object(
                release_tool, "_git", return_value=origin
            ):
                release_tool._require_github_origin()

        with mock.patch.object(
            release_tool, "_git", return_value="/scratch/local/gravlax"
        ), self.assertRaises(release_tool.ReleaseError):
            release_tool._require_github_origin()

    def test_release_operations_require_main(self):
        with mock.patch.object(release_tool, "_git", return_value="main"):
            release_tool._require_main_branch()
        with mock.patch.object(
            release_tool, "_git", return_value="release-v0.1.0"
        ), self.assertRaises(release_tool.ReleaseError):
            release_tool._require_main_branch()


if __name__ == "__main__":
    unittest.main()
