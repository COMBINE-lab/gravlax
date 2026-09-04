from __future__ import annotations

import json
import re
import unittest
from pathlib import Path


REPOSITORY = Path(__file__).resolve().parents[2]
NOTEBOOKS = REPOSITORY / "notebooks"
PUBLISHED_MANIFEST_URL = (
    "https://github.com/COMBINE-lab/gravlax/releases/download/"
    "demo-data-v1/demo-manifest.json"
)
PUBLISHED_MANIFEST_SHA256 = (
    "82c34aad442d478f1cb1243a6ccfe8ad9f937b81d9e1f946a8eb2cfc498214fd"
)


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
        self.assertIn("comparison_annotation", source)
        self.assertIn("comparison_result", source)
        self.assertIn("expected_comparison_compatible_transcripts", source)

    def test_collection_notebooks_rebuild_from_downloaded_source_archives(self):
        for name in (
            "02_multi_donor_event_discovery.ipynb",
            "03_federated_junction_cooccurrence.ipynb",
        ):
            source = code_source(name)
            self.assertIn("['collection', 'build']", source)
            self.assertIn("--sample=", source)
            self.assertNotRegex(source, r"fetch\([^\n]*\.aicollection")

    def test_manifest_setup_is_pinned_and_fail_closed(self):
        self.assertRegex(PUBLISHED_MANIFEST_SHA256, r"\A[0-9a-f]{64}\Z")
        for name in (
            "01_annotation_reinterpretation.ipynb",
            "02_multi_donor_event_discovery.ipynb",
            "03_federated_junction_cooccurrence.ipynb",
        ):
            source = code_source(name)
            self.assertEqual(
                source.count(f'MANIFEST_URL = "{PUBLISHED_MANIFEST_URL}"'),
                1,
            )
            self.assertEqual(
                source.count(
                    f'MANIFEST_SHA256 = "{PUBLISHED_MANIFEST_SHA256}"'
                ),
                1,
            )
            self.assertIn("raise RuntimeError", source)
            self.assertIn("hashlib.sha256", source)
            self.assertIn("EXPECTED_VERSION", source)
            self.assertIn("CLI_VERSION != f'aie {EXPECTED_VERSION}'", source)
            self.assertIn("PYTHON_VERSION != EXPECTED_VERSION", source)
            self.assertIn("required_story_fields", source)
            self.assertIn("RESERVED_ASSET_FILENAMES", source)
            self.assertIn("duplicate asset filename", source)
            self.assertIn("/latest/ is not allowed", source)
            self.assertIn("'--force-reinstall', '--no-deps'", source)

    def test_collection_and_drilldown_group_files_have_distinct_contracts(self):
        event_source = code_source("02_multi_donor_event_discovery.ipynb")
        drilldown_source = code_source(
            "03_federated_junction_cooccurrence.ipynb"
        )
        self.assertIn("story['collection_groups']", event_source)
        self.assertNotIn("story['drilldown_groups']", event_source)
        self.assertIn("story['drilldown_groups']", drilldown_source)
        self.assertNotIn("story['collection_groups']", drilldown_source)
        self.assertIn("story['emit_membership']", drilldown_source)
        self.assertIn("cooccurrence.table('memberships')", drilldown_source)
        self.assertIn("selected-unit summary differs from membership witnesses", drilldown_source)

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
                "cooccurrence.table('memberships')",
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
        self.assertIn("collection_groups", event_story["required"])
        self.assertNotIn("groups", event_story["properties"])
        self.assertIn(
            "drilldown_groups",
            schema["$defs"]["drilldown_story"]["properties"],
        )
        self.assertIn("emit_membership", schema["$defs"]["drilldown_story"]["required"])
        self.assertIs(
            schema["$defs"]["drilldown_story"]["properties"]["emit_membership"]["const"],
            True,
        )
        asset = schema["$defs"]["asset"]
        self.assertRegex(
            "https://example.test/releases/latest/file",
            re.compile(asset["properties"]["url"]["not"]["pattern"]),
        )
        self.assertIn(
            "manifest.json",
            asset["properties"]["filename"]["not"]["enum"],
        )
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
