from __future__ import annotations

import ast
import json
import re
import tempfile
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


COLLECTION_NOTEBOOKS = (
    "02_multi_donor_event_discovery.ipynb",
    "03_federated_junction_cooccurrence.ipynb",
)
ARCHIVE_ROOTS = {
    "donor-a": f"aie-directory-root-v2:{'a' * 64}",
    "donor-b": f"aie-directory-root-v2:{'b' * 64}",
}
COLLECTION_ROOT = f"aicollection-directory-root-v1:{'c' * 64}"


def loader_source(name: str) -> str:
    """Return the notebook's collection loader without executing any other cell."""
    source = code_source(name)
    tree = ast.parse(source)
    selected = [
        node
        for node in tree.body
        if (isinstance(node, ast.FunctionDef) and node.name == "prebuilt_collection")
        or (
            isinstance(node, ast.Assign)
            and any(
                getattr(target, "id", None) == "COLLECTION_ROOT_PATTERN"
                for target in node.targets
            )
        )
    ]
    if len(selected) != 2:
        raise AssertionError(f"{name} does not define exactly one collection loader")
    return "\n".join(
        segment
        for node in sorted(selected, key=lambda node: node.lineno)
        if (segment := ast.get_source_segment(source, node)) is not None
    )


class CollectionLoaderFixture:
    """A synthetic capsule download directory for the notebook loader."""

    def __init__(self, directory: Path, *, declare_collection: bool = True):
        self.directory = directory
        self.resources = {
            "archive_a": {
                "url": "https://example.invalid/archive-a.aie",
                "sha256": "0" * 64,
                "filename": "archive-a.aie",
                "archive_root": ARCHIVE_ROOTS["donor-a"],
            },
            "archive_b": {
                "url": "https://example.invalid/archive-b.aie",
                "sha256": "1" * 64,
                "filename": "archive-b.aie",
                "archive_root": ARCHIVE_ROOTS["donor-b"],
            },
        }
        self.archives = {}
        for sample, resource in (("donor-a", "archive_a"), ("donor-b", "archive_b")):
            path = directory / self.resources[resource]["filename"]
            path.write_bytes(b"archive\n")
            self.archives[sample] = path
        self.story = {
            "archives": {"donor-a": "archive_a", "donor-b": "archive_b"},
            "shape_routes": True,
        }
        self.locations = {
            "schema_version": 1,
            "locations": [
                {"identity": ARCHIVE_ROOTS["donor-a"], "path": "archive-a.aie"},
                {"identity": ARCHIVE_ROOTS["donor-b"], "path": "archive-b.aie"},
            ],
        }
        self.manifest = {"schema": "gravlax.demo-capsule.v1", "resources": self.resources}
        if declare_collection:
            self.manifest["collection"] = {
                "archives": {"donor-a": "archive_a", "donor-b": "archive_b"},
                "shape_routes": True,
                "allow_unstamped": False,
                "collection_root": COLLECTION_ROOT,
                "asset": {
                    "url": "https://example.invalid/demo.aicollection",
                    "sha256": "2" * 64,
                    "filename": "demo.aicollection",
                },
                "locations": {
                    "url": "https://example.invalid/demo-locations.json",
                    "sha256": "3" * 64,
                    "filename": "demo-locations.json",
                },
            }
        self.inspection = {
            "layers": [{"root_digest": "c" * 64}],
            "archives": [
                {
                    "id": sample,
                    "native_identity": {
                        "scheme": "aie-directory-root-v2",
                        "blake3": root.split(":", 1)[1],
                    },
                }
                for sample, root in ARCHIVE_ROOTS.items()
            ],
            "index": {"shape_route_archives": 2},
            "guard": {
                "content_identity_verified": True,
                "shape_route_payloads_verified": True,
                "shape_route_reconstruction_verified": True,
            },
        }
        self.commands = []

    def fetch(self, spec):
        path = self.directory / spec["filename"]
        if spec["filename"].endswith(".json"):
            path.write_text(json.dumps(self.locations), encoding="utf-8")
        elif not path.exists():
            path.write_bytes(b"collection\n")
        return path

    def story_resource(self, name):
        spec = self.resources.get(name)
        if not isinstance(name, str) or not isinstance(spec, dict):
            raise RuntimeError(f"story resource is not declared: {name!r}")
        return spec

    def client(self):
        fixture = self

        class Client:
            def result_raw(self, args):
                fixture.commands.append([str(item) for item in args])
                return fixture.inspection

        return Client()

    def namespace(self, name: str) -> dict:
        namespace = {
            "re": re,
            "json": json,
            "Path": Path,
            "manifest": self.manifest,
            "archives": self.archives,
            "fetch": self.fetch,
            "story_resource": self.story_resource,
            "client": self.client(),
            "asset_filenames": {
                asset["filename"]: name for name, asset in self.resources.items()
            },
            "RESERVED_ASSET_FILENAMES": {
                ".",
                "..",
                "manifest.json",
                "event-discovery.aicollection",
                "junction-drilldown.aicollection",
            },
        }
        exec(compile(loader_source(name), name, "exec"), namespace)  # noqa: S102
        return namespace


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

    def test_collection_notebooks_keep_the_v1_local_rebuild_path(self):
        for name in COLLECTION_NOTEBOOKS:
            source = code_source(name)
            self.assertIn("['collection', 'build']", source)
            self.assertIn("--sample=", source)
            self.assertIn("collection, locations = prebuilt_collection(story)", source)
            self.assertIn("if collection is None:", source)
            self.assertNotRegex(source, r"fetch\([^\n]*\.aicollection")

    def test_collection_notebooks_share_one_loader_and_pass_it_to_every_query(self):
        sources = {name: loader_source(name) for name in COLLECTION_NOTEBOOKS}
        self.assertEqual(len(set(sources.values())), 1)
        event = code_source("02_multi_donor_event_discovery.ipynb")
        self.assertEqual(event.count("collection, locations=locations,"), 2)
        drilldown = code_source("03_federated_junction_cooccurrence.ipynb")
        self.assertIn(
            "*([f'--locations={locations}'] if locations else ())", drilldown
        )
        self.assertNotIn("prebuilt_collection", code_source(
            "01_annotation_reinterpretation.ipynb"
        ))

    def test_loader_rebuilds_locally_when_the_manifest_declares_no_collection(self):
        for name in COLLECTION_NOTEBOOKS:
            with tempfile.TemporaryDirectory() as temporary:
                fixture = CollectionLoaderFixture(
                    Path(temporary), declare_collection=False
                )
                namespace = fixture.namespace(name)
                self.assertEqual(
                    namespace["prebuilt_collection"](fixture.story), (None, None)
                )
                self.assertEqual(fixture.commands, [])

    def test_loader_resolves_the_published_collection_through_its_locations(self):
        for name in COLLECTION_NOTEBOOKS:
            with tempfile.TemporaryDirectory() as temporary:
                fixture = CollectionLoaderFixture(Path(temporary))
                namespace = fixture.namespace(name)
                collection, locations = namespace["prebuilt_collection"](fixture.story)
                self.assertEqual(collection, fixture.directory / "demo.aicollection")
                self.assertEqual(locations, fixture.directory / "demo-locations.json")
                self.assertEqual(len(fixture.commands), 1)
                command = fixture.commands[0]
                self.assertEqual(command[:3], ["collection", "inspect", str(collection)])
                self.assertIn(f"--locations={locations}", command)
                self.assertIn("--verify-routes", command)

    def test_loader_fails_closed_on_every_published_collection_mismatch(self):
        cases = {
            "root": (
                lambda fixture: fixture.manifest["collection"].update(
                    collection_root=f"aicollection-directory-root-v1:{'d' * 64}"
                ),
                "root does not match the manifest",
            ),
            "unknown field": (
                lambda fixture: fixture.manifest["collection"].update(extra=1),
                "missing or unknown fields",
            ),
            "other archives": (
                lambda fixture: fixture.manifest["collection"].update(
                    archives={"donor-a": "archive_a"}
                ),
                "does not commit exactly",
            ),
            "other options": (
                lambda fixture: fixture.manifest["collection"].update(
                    shape_routes=False
                ),
                "collection options",
            ),
            "reserved filename": (
                lambda fixture: fixture.manifest["collection"]["asset"].update(
                    filename="event-discovery.aicollection"
                ),
                "reserved, or duplicate filename",
            ),
            "duplicate filename": (
                lambda fixture: fixture.manifest["collection"]["asset"].update(
                    filename="archive-a.aie"
                ),
                "reserved, or duplicate filename",
            ),
            "short locations": (
                lambda fixture: fixture.locations["locations"].pop(),
                "does not resolve exactly",
            ),
            "repeated identity": (
                lambda fixture: fixture.locations["locations"].append(
                    dict(fixture.locations["locations"][0])
                ),
                "repeats an identity",
            ),
            "non-basename location": (
                lambda fixture: fixture.locations["locations"][0].update(
                    path="nested/archive-a.aie"
                ),
                "capsule basenames",
            ),
            "unsupported locations schema": (
                lambda fixture: fixture.locations.update(schema_version=2),
                "unsupported location manifest schema",
            ),
            "committed root drift": (
                lambda fixture: fixture.inspection["archives"][0][
                    "native_identity"
                ].update(blake3="e" * 64),
                "verified archive roots",
            ),
            "layered collection": (
                lambda fixture: fixture.inspection["layers"].append(
                    {"root_digest": "f" * 64}
                ),
                "single rooted layer",
            ),
            "unverified routes": (
                lambda fixture: fixture.inspection["guard"].update(
                    shape_route_reconstruction_verified=False
                ),
                "did not authenticate",
            ),
            "missing routes": (
                lambda fixture: fixture.inspection["index"].update(
                    shape_route_archives=1
                ),
                "shape route for every",
            ),
        }
        for label, (mutate, message) in cases.items():
            with self.subTest(case=label):
                with tempfile.TemporaryDirectory() as temporary:
                    fixture = CollectionLoaderFixture(Path(temporary))
                    mutate(fixture)
                    namespace = fixture.namespace(COLLECTION_NOTEBOOKS[0])
                    with self.assertRaisesRegex(RuntimeError, message):
                        namespace["prebuilt_collection"](fixture.story)

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
        self.assertNotIn("collection", schema["required"])
        collection = schema["$defs"]["collection"]
        self.assertFalse(collection["additionalProperties"])
        self.assertEqual(
            set(collection["required"]),
            {"archives", "collection_root", "asset", "locations"},
        )
        self.assertRegex(
            COLLECTION_ROOT, re.compile(collection["properties"]["collection_root"]["pattern"])
        )
        self.assertEqual(
            collection["properties"]["asset"]["allOf"][0]["$ref"], "#/$defs/asset"
        )
        self.assertEqual(
            collection["properties"]["asset"]["properties"]["filename"]["pattern"],
            r"\.aicollection$",
        )
        self.assertIn("max_candidates_considered", event_story["properties"])
        self.assertIn("max_routed_entries", event_story["properties"])
        self.assertIn("max_exact_match_attempts", event_story["properties"])
        self.assertIn("max_annotation_comparisons", event_story["properties"])


if __name__ == "__main__":
    unittest.main()
