from __future__ import annotations

import hashlib
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[2]
RENDERER = ROOT / "packaging/bioconda/render_recipe.py"


class BiocondaRecipeTests(unittest.TestCase):
    def render(self, directory: Path) -> tuple[Path, Path]:
        archive = directory / "gravlax-0.1.0-source.tar.gz"
        archive.write_bytes(b"published Gravlax source asset\n")
        recipe = directory / "recipe/meta.yaml"
        subprocess.run(
            [
                sys.executable,
                str(RENDERER),
                "--version",
                "0.1.0",
                "--source-archive",
                str(archive),
                "--output",
                str(recipe),
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        return archive, recipe

    def test_rendered_recipe_is_bound_and_submission_complete(self):
        with tempfile.TemporaryDirectory() as temporary:
            archive, recipe = self.render(Path(temporary))
            text = recipe.read_text(encoding="utf-8")
            digest = hashlib.sha256(archive.read_bytes()).hexdigest()

            self.assertIn(f'set version = "0.1.0"', text)
            self.assertIn(f'set sha256 = "{digest}"', text)
            self.assertIn("releases/download/v{{ version }}", text)
            self.assertIn("cargo install --locked --offline --no-track", text)
            self.assertIn("compiler('c')", text)
            self.assertIn("stdlib('c')", text)
            self.assertIn("rust >=1.98,<2", text)
            self.assertIn("cargo-bundle-licenses", text)
            self.assertIn("THIRDPARTY.yml", text)
            self.assertIn("should_use_compilers", text)
            self.assertIn("missing_run_exports", text)
            self.assertNotIn("run_exports:", text)
            self.assertNotIn("pin_subpackage", text)
            self.assertNotIn("@VERSION@", text)
            self.assertNotIn("@SOURCE_SHA256@", text)

    def test_renderer_refuses_to_replace_recipe_file(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            _, recipe = self.render(directory)
            original_recipe = recipe.read_bytes()
            repeated = subprocess.run(
                [
                    sys.executable,
                    str(RENDERER),
                    "--version",
                    "0.1.0",
                    "--source-archive",
                    str(directory / "gravlax-0.1.0-source.tar.gz"),
                    "--output",
                    str(recipe),
                ],
                capture_output=True,
                text=True,
            )
            self.assertNotEqual(repeated.returncode, 0)
            self.assertIn("refusing to overwrite", repeated.stderr)
            self.assertEqual(recipe.read_bytes(), original_recipe)

    def test_renderer_accepts_only_stable_versions_and_meta_yaml_output(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            archive = directory / "gravlax-0.1.0-rc.1-source.tar.gz"
            archive.write_bytes(b"not a stable release")
            prerelease = subprocess.run(
                [
                    sys.executable,
                    str(RENDERER),
                    "--version",
                    "0.1.0-rc.1",
                    "--source-archive",
                    str(archive),
                    "--output",
                    str(directory / "recipe/meta.yaml"),
                ],
                capture_output=True,
                text=True,
            )
            self.assertNotEqual(prerelease.returncode, 0)
            self.assertIn("stable major.minor.patch", prerelease.stderr)

            stable_archive = directory / "gravlax-0.1.0-source.tar.gz"
            stable_archive.write_bytes(b"stable release")
            wrong_output = subprocess.run(
                [
                    sys.executable,
                    str(RENDERER),
                    "--version",
                    "0.1.0",
                    "--source-archive",
                    str(stable_archive),
                    "--output",
                    str(directory / "recipe/gravlax.yaml"),
                ],
                capture_output=True,
                text=True,
            )
            self.assertNotEqual(wrong_output.returncode, 0)
            self.assertIn("must name the Bioconda recipe file meta.yaml", wrong_output.stderr)


if __name__ == "__main__":
    unittest.main()
