from __future__ import annotations

import json
import subprocess
import unittest
from unittest.mock import patch

from gravlax import Client


def envelope():
    return {
        "$schema": "gravlax.result-envelope.v1",
        "result_schema": "gravlax.annotation.resolve.v1",
        "producer": {"name": "aie", "version": "0.1.0"},
        "provenance": {
            "assembly": "GRCh38.p14",
            "annotation": "GENCODE 49",
            "annotation_digest": "blake3:" + "a" * 64,
            "parameters": {"identifiers": ["TP53"]},
        },
        "warnings": [],
        "data": {
            "schema": {
                "id": "gravlax.annotation.resolve.v1",
                "fields": [
                    {"name": "requested", "data_type": "string"},
                    {"name": "stable_id", "data_type": "string"},
                ],
            },
            "rows": [["TP53", "ENSG00000141510.19"]],
        },
    }


class ResolveClientTests(unittest.TestCase):
    @patch("gravlax.client.subprocess.run")
    def test_resolve_returns_typed_envelope_with_literal_arguments(self, run):
        run.return_value = subprocess.CompletedProcess(
            [], 0, json.dumps(envelope()), ""
        )

        result = Client(binary=("wrapper", "aie")).resolve(
            "annotations/gencode 49.aic",
            ["TP53"],
            assembly="GRCh38.p14",
            annotation="GENCODE 49",
        )

        self.assertEqual(result.result_schema, "gravlax.annotation.resolve.v1")
        self.assertEqual(result.table.records()[0]["requested"], "TP53")
        self.assertEqual(
            run.call_args.args[0],
            (
                "wrapper",
                "aie",
                "resolve",
                "annotations/gencode 49.aic",
                "--assembly=GRCh38.p14",
                "--annotation=GENCODE 49",
                "--format=json",
                "--",
                "TP53",
            ),
        )
        self.assertFalse(run.call_args.kwargs["shell"])

    def test_resolve_requires_a_nonempty_identifier_sequence(self):
        with self.assertRaises(ValueError):
            Client().resolve(
                "annotation.gtf",
                [],
                assembly="GRCh38",
                annotation="test",
            )


if __name__ == "__main__":
    unittest.main()
