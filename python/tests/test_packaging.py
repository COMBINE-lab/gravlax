from pathlib import Path
import tomllib
import unittest


PYTHON_ROOT = Path(__file__).resolve().parents[1]
REPOSITORY = PYTHON_ROOT.parent


class PackagingMetadataTests(unittest.TestCase):
    def test_python_and_cli_versions_stay_in_sync(self):
        python = tomllib.loads((PYTHON_ROOT / "pyproject.toml").read_text(encoding="utf-8"))
        cargo = tomllib.loads((REPOSITORY / "Cargo.toml").read_text(encoding="utf-8"))
        package_version = python["project"]["version"]
        self.assertEqual(package_version, cargo["workspace"]["package"]["version"])

        init = (PYTHON_ROOT / "src/gravlax/__init__.py").read_text(encoding="utf-8")
        self.assertIn(f'__version__ = "{package_version}"', init)

    def test_distribution_declares_and_contains_its_license(self):
        python = tomllib.loads((PYTHON_ROOT / "pyproject.toml").read_text(encoding="utf-8"))
        self.assertEqual(python["project"]["license"], "BSD-3-Clause")
        self.assertEqual(python["project"]["license-files"], ["LICENSE"])
        self.assertEqual(
            (PYTHON_ROOT / "LICENSE").read_bytes(),
            (REPOSITORY / "LICENSE").read_bytes(),
        )


if __name__ == "__main__":
    unittest.main()
