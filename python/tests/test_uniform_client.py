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


if __name__ == "__main__":
    unittest.main()
