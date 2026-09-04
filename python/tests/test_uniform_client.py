from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from gravlax import Client, ExactSelection, RegionQuerySummary

from test_uniform_results import junction_bundle, region_bundle


def completed(argv, stdout: str, *, returncode: int = 0, stderr: str = ""):
    return subprocess.CompletedProcess(argv, returncode, stdout, stderr)


class UniformClientTests(unittest.TestCase):
    @patch("gravlax.client.subprocess.run")
    def test_region_wrapper_requests_the_uniform_json_contract(self, run):
        run.return_value = completed([], json.dumps(region_bundle()))

        result = Client(binary=("wrapper", "aie")).query_region(
            "sample with spaces.aie",
            "chr16:89550000-89575000",
            top=1,
            cells="selected cells.txt",
            aggregation="cell",
        )

        self.assertIsInstance(result.summary, RegionQuerySummary)
        self.assertIsInstance(result.table("counts").selection, ExactSelection)
        self.assertEqual(
            run.call_args.args[0],
            (
                "wrapper",
                "aie",
                "query",
                "sample with spaces.aie",
                "region",
                "chr16:89550000-89575000",
                "--top=1",
                "--format=json",
                "--cells=selected cells.txt",
                "--agg=cell",
            ),
        )
        self.assertIs(run.call_args.kwargs["shell"], False)

    @patch("gravlax.client.subprocess.run")
    def test_junction_wrapper_uses_top_zero_as_all(self, run):
        run.return_value = completed([], json.dumps(junction_bundle()))

        result = Client().query_junction(
            "sample.aie", "chr16:89562391-89562883", top=0
        )

        self.assertEqual(result.result_schema, "gravlax.query.junction.result.v1")
        self.assertIn("--top=0", run.call_args.args[0])
        self.assertIn("--format=json", run.call_args.args[0])

    @patch("gravlax.client.subprocess.run")
    def test_file_wrapper_streams_without_capturing_stdout(self, run):
        encoded = json.dumps(region_bundle()).encode("utf-8")

        def invoke(argv, **kwargs):
            kwargs["stdout"].write(encoded)
            return completed(argv, "", stderr="region complete\n")

        run.side_effect = invoke
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "region.json"
            result = Client().query_region_to_file(
                "sample.aie",
                "chr16:89550000-89575000",
                output,
                top=0,
            )
            self.assertEqual(result.bytes, len(encoded))
            self.assertEqual(output.read_bytes(), encoded)
            parsed = Client().result_bundle_from_file(output)
            self.assertEqual(parsed.table_names, ("counts",))

        self.assertNotEqual(run.call_args.kwargs["stdout"], subprocess.PIPE)
        self.assertIn("--format=json", run.call_args.args[0])

    def test_scope_and_numeric_validation_happen_before_execution(self):
        client = Client(binary="definitely-not-run")
        with self.assertRaisesRegex(ValueError, "mutually exclusive"):
            client.query_region(
                "sample.aie",
                "chr1:0-1",
                cells="cells.txt",
                groups="groups.tsv",
            )
        with self.assertRaisesRegex(ValueError, "requires groups"):
            client.query_junction(
                "sample.aie",
                "chr1:1-2",
                aggregation="group",
            )
        with self.assertRaisesRegex(ValueError, "nonnegative"):
            client.query_region("sample.aie", "chr1:0-1", top=-1)

    @patch("gravlax.client.subprocess.run")
    def test_cooccurrence_wrapper_keeps_each_predicate_and_expression_one_token(self, run):
        run.return_value = completed([], json.dumps(region_bundle()))

        Client(binary=("wrapper", "aie")).query_cooccurrence(
            "sample with spaces.aie",
            {
                "splice": "junction:chr1:100-200:+",
                "tail": "terminal:chr1:300-325:+",
            },
            "splice & !tail",
            universe="splice",
            groups="cell groups.tsv",
            aggregation="group",
            emit_membership=True,
            max_evidence_records=1234,
            max_terminal_events=5678,
        )

        argv = run.call_args.args[0]
        self.assertIn("--predicate=splice=junction:chr1:100-200:+", argv)
        self.assertIn("--predicate=tail=terminal:chr1:300-325:+", argv)
        self.assertIn("--where=splice & !tail", argv)
        self.assertIn("--universe=splice", argv)
        self.assertIn("--groups=cell groups.tsv", argv)
        self.assertIn("--emit-membership", argv)
        self.assertIn("--max-evidence-records=1234", argv)
        self.assertIn("--max-terminal-events=5678", argv)
        self.assertIn("--placements=unique", argv)
        self.assertIs(run.call_args.kwargs["shell"], False)

    def test_cooccurrence_wrapper_rejects_unbounded_class_union(self):
        client = Client(binary="definitely-not-run")
        with self.assertRaisesRegex(ValueError, "requires allow_full_scan"):
            client.query_cooccurrence(
                "sample.aie",
                {"gene": "region:chr1:1-100"},
                "gene",
                universe="gene",
                unit="umi-class",
            )
        with self.assertRaisesRegex(ValueError, "must name one"):
            client.query_cooccurrence(
                "sample.aie",
                {"gene": "region:chr1:1-100"},
                "gene",
                universe="missing",
            )

    @patch("gravlax.client.subprocess.run")
    def test_collection_find_events_wrapper_emits_repeatable_filters(self, run):
        run.return_value = completed([], json.dumps(region_bundle()))

        Client().collection_find_events(
            "atlas.aicollection",
            kinds=("junction", "cassette", "terminal-tail"),
            design="donors.tsv",
            groups="groups.tsv",
            require_groups=("neuron", "astrocyte"),
            min_donors=3,
            min_umi_classes=10,
            terminal_cluster_bp=12,
            max_terminal_events=12345,
            annotation="gencode.gtf",
            assembly="GRCh38",
            annotation_label="GENCODE 49",
            annotation_digest="blake3:" + "ab" * 32,
            novel_only=True,
            max_candidates_considered=4321,
            max_routed_entries=8765,
            max_exact_match_attempts=9876,
            max_annotation_comparisons=7654,
        )

        argv = run.call_args.args[0]
        self.assertEqual(argv[:3], ("aie", "collection", "find-events"))
        self.assertIn("--kind=junction", argv)
        self.assertIn("--kind=cassette", argv)
        self.assertIn("--kind=terminal-tail", argv)
        self.assertIn("--require-group=neuron", argv)
        self.assertIn("--require-group=astrocyte", argv)
        self.assertIn("--min-donors=3", argv)
        self.assertIn("--min-umi-classes=10", argv)
        self.assertIn("--min-group-umi-classes=1", argv)
        self.assertIn("--terminal-cluster-bp=12", argv)
        self.assertIn("--max-terminal-events=12345", argv)
        self.assertIn("--assembly=GRCh38", argv)
        self.assertIn("--annotation-label=GENCODE 49", argv)
        self.assertIn("--annotation-digest=blake3:" + "ab" * 32, argv)
        self.assertIn("--max-candidates-considered=4321", argv)
        self.assertIn("--max-routed-entries=8765", argv)
        self.assertIn("--max-exact-match-attempts=9876", argv)
        self.assertIn("--max-annotation-comparisons=7654", argv)
        self.assertIn("--novel-only", argv)
        self.assertIn("--format=json", argv)

    def test_collection_find_events_wrapper_rejects_invalid_dependencies(self):
        client = Client(binary="definitely-not-run")
        with self.assertRaisesRegex(ValueError, "requires groups"):
            client.collection_find_events(
                "atlas.aicollection", require_groups=("neuron",)
            )
        with self.assertRaisesRegex(ValueError, "requires annotation"):
            client.collection_find_events("atlas.aicollection", novel_only=True)
        with self.assertRaisesRegex(ValueError, "requires a non-empty assembly"):
            client.collection_find_events(
                "atlas.aicollection",
                annotation="gencode.gtf",
                annotation_label="GENCODE 49",
            )
        with self.assertRaisesRegex(ValueError, "require annotation"):
            client.collection_find_events(
                "atlas.aicollection", assembly="GRCh38"
            )
        with self.assertRaisesRegex(ValueError, "nonnegative"):
            client.collection_find_events(
                "atlas.aicollection", terminal_cluster_bp=-1
            )
        with self.assertRaisesRegex(ValueError, "requires solo_strand='forward'"):
            client.collection_find_events(
                "atlas.aicollection",
                kinds=("terminal-tail",),
                solo_strand="reverse",
            )

    @patch("gravlax.client.subprocess.run")
    def test_collection_find_events_does_not_send_group_threshold_without_groups(
        self, run
    ):
        run.return_value = completed([], json.dumps(region_bundle()))
        Client().collection_find_events("atlas.aicollection")
        argv = run.call_args.args[0]
        self.assertNotIn("--min-group-umi-classes=1", argv)
        self.assertNotIn("--explain", argv)

    @patch("gravlax.client.subprocess.run")
    def test_collection_find_events_file_wrapper_streams_the_shared_argv(self, run):
        encoded = json.dumps(region_bundle()).encode("utf-8")

        def invoke(argv, **kwargs):
            kwargs["stdout"].write(encoded)
            return completed(argv, "", stderr="find-events complete\n")

        run.side_effect = invoke
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "events.json"
            result = Client().collection_find_events_to_file(
                "atlas.aicollection",
                output,
                kinds=("terminal-tail",),
                terminal_cluster_bp=7,
            )
            self.assertEqual(result.bytes, len(encoded))
            self.assertTrue(output.is_file())
        argv = run.call_args.args[0]
        self.assertIn("--kind=terminal-tail", argv)
        self.assertIn("--terminal-cluster-bp=7", argv)
        self.assertIn("--format=json", argv)
        self.assertNotEqual(run.call_args.kwargs["stdout"], subprocess.PIPE)


if __name__ == "__main__":
    unittest.main()
