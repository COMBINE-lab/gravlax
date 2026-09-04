from __future__ import annotations

import json
import unittest
from pathlib import Path


REPOSITORY = Path(__file__).resolve().parents[2]
NOTEBOOKS = REPOSITORY / "notebooks"


def code_source(name: str) -> str:
    notebook = json.loads((NOTEBOOKS / name).read_text(encoding="utf-8"))
    return "\n".join(
        "".join(cell.get("source", ()))
        for cell in notebook["cells"]
        if cell.get("cell_type") == "code"
    )


class DemoNotebookTests(unittest.TestCase):
    def test_every_code_cell_compiles(self):
        for path in sorted(NOTEBOOKS.glob("*.ipynb")):
            notebook = json.loads(path.read_text(encoding="utf-8"))
            for index, cell in enumerate(notebook["cells"]):
                if cell.get("cell_type") == "code":
                    compile(
                        "".join(cell.get("source", ())),
                        f"{path.name}:cell-{index}",
                        "exec",
                    )

    def test_event_discovery_uses_the_frozen_table_names(self):
        source = code_source("02_multi_donor_event_discovery.ipynb")
        self.assertIn("result.table('entities')", source)
        self.assertNotIn("result.table('events')", source)
        for table in (
            "capabilities",
            "entities",
            "components",
            "counts",
            "terminal_anchors",
            "terminal_counts",
        ):
            self.assertIn(repr(table), source)

    def test_collection_notebooks_rebuild_from_downloaded_source_archives(self):
        for name in (
            "02_multi_donor_event_discovery.ipynb",
            "03_federated_junction_cooccurrence.ipynb",
        ):
            source = code_source(name)
            self.assertIn("['collection', 'build']", source)
            self.assertIn("--sample=", source)
            self.assertNotIn(".aicollection'", source.split("fetch(")[0])

    def test_manifest_setup_is_fail_closed(self):
        for name in (
            "01_annotation_reinterpretation.ipynb",
            "02_multi_donor_event_discovery.ipynb",
            "03_federated_junction_cooccurrence.ipynb",
        ):
            source = code_source(name)
            self.assertRegex(source, r"MANIFEST_URL = ['\"]{2}")
            self.assertRegex(source, r"MANIFEST_SHA256 = ['\"]{2}")
            self.assertIn("raise RuntimeError", source)
            self.assertIn("hashlib.sha256", source)
            self.assertIn("EXPECTED_VERSION", source)
            self.assertIn("CLI_VERSION != f'aie {EXPECTED_VERSION}'", source)
            self.assertIn("PYTHON_VERSION != EXPECTED_VERSION", source)
            self.assertIn("required_story_fields", source)

    def test_each_notebook_plots_only_live_typed_output(self):
        expected_columns = {
            "01_annotation_reinterpretation.ipynb": (
                "comparison.count_deltas",
                "comparison_gene_id",
                "signed_delta_b_minus_a",
            ),
            "02_multi_donor_event_discovery.ipynb": (
                "result.table('entities')",
                "result.table('counts')",
                "informative_umi_classes",
            ),
            "03_federated_junction_cooccurrence.ipynb": (
                "federated.table('samples')",
                "cooccurrence.table('patterns')",
                "evidence_units",
            ),
        }
        for name, markers in expected_columns.items():
            source = code_source(name)
            self.assertIn("require_table_columns", source)
            self.assertIn("display(SVG(data=", source)
            for marker in markers:
                self.assertIn(marker, source)

    def test_manifest_schema_constrains_each_story(self):
        schema = json.loads(
            (NOTEBOOKS / "demo-manifest.schema.json").read_text(encoding="utf-8")
        )
        story_properties = schema["properties"]["stories"]["properties"]
        self.assertEqual(
            story_properties["annotation_reinterpretation"]["$ref"],
            "#/$defs/annotation_story",
        )
        self.assertEqual(
            story_properties["event_discovery"]["$ref"],
            "#/$defs/event_story",
        )
        self.assertEqual(
            story_properties["junction_drilldown"]["$ref"],
            "#/$defs/drilldown_story",
        )
        for name in ("annotation_story", "event_story", "drilldown_story"):
            definition = schema["$defs"][name]
            self.assertFalse(definition["additionalProperties"])
            self.assertTrue(definition["required"])
        event_story = schema["$defs"]["event_story"]
        self.assertEqual(
            event_story["allOf"][0]["then"]["required"],
            ["assembly", "annotation_label"],
        )
        self.assertIn("max_candidates_considered", event_story["properties"])
        self.assertIn("max_routed_entries", event_story["properties"])
        self.assertIn("max_exact_match_attempts", event_story["properties"])
        self.assertIn("max_annotation_comparisons", event_story["properties"])


if __name__ == "__main__":
    unittest.main()
