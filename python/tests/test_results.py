from __future__ import annotations

import importlib
import json
import unittest
from unittest.mock import patch

from gravlax import OptionalDependencyError, ProtocolError, parse_result


def table_envelope():
    return {
        "$schema": "gravlax.result-envelope.v1",
        "result_schema": "gravlax.query.region.v1",
        "producer": {"name": "aie", "version": "0.1.0"},
        "provenance": {
            "archives": ["sample.aie"],
            "assembly": "GRCh38",
            "parameters": {"top": 20},
        },
        "warnings": [],
        "data": {
            "schema": {
                "id": "gravlax.query.region.v1",
                "fields": [
                    {"name": "cell", "data_type": "string", "nullable": False},
                    {"name": "umis", "data_type": "uint64", "nullable": False},
                    {"name": "score", "data_type": "float64", "nullable": True},
                    {"name": "details", "data_type": "json", "nullable": False},
                ],
            },
            "rows": [["AAAC-1", 17, None, {"source": "unique"}]],
        },
    }


class ResultTests(unittest.TestCase):
    def test_parses_and_exposes_dependency_free_records(self):
        result = parse_result(json.dumps(table_envelope()))

        self.assertEqual(result.result_schema, "gravlax.query.region.v1")
        self.assertEqual(result.provenance.assembly, "GRCh38")
        self.assertEqual(result.table.columns, ("cell", "umis", "score", "details"))
        self.assertEqual(result.table.records()[0]["umis"], 17)

    def test_parses_role_named_annotation_provenance(self):
        value = table_envelope()
        value["provenance"]["annotations"] = [
            {
                "role": "before",
                "assembly": "GRCh38",
                "annotation": "GENCODE 48",
                "digest": "blake3:" + "a" * 64,
            },
            {
                "role": "after",
                "assembly": "GRCh38",
                "annotation": "GENCODE 49",
                "digest": "blake3:" + "b" * 64,
            },
        ]
        result = parse_result(value)
        self.assertEqual(
            [annotation.role for annotation in result.provenance.annotations],
            ["before", "after"],
        )
        self.assertEqual(
            result.provenance.as_dict()["annotations"][1]["annotation"],
            "GENCODE 49",
        )

        value["provenance"]["annotations"][1]["role"] = "before"
        with self.assertRaisesRegex(ProtocolError, "duplicate roles"):
            parse_result(value)

    def test_rejects_schema_disagreement(self):
        value = table_envelope()
        value["data"]["schema"]["id"] = "gravlax.query.other.v1"
        result = parse_result(value)
        with self.assertRaises(ProtocolError):
            _ = result.table

    def test_rejects_noncanonical_schema_and_field_names(self):
        value = table_envelope()
        value["result_schema"] = "gravlax query region v1"
        with self.assertRaisesRegex(ProtocolError, "ASCII"):
            parse_result(value)

        value = table_envelope()
        value["data"]["schema"]["fields"][0]["name"] = "bad\tname"
        with self.assertRaisesRegex(ProtocolError, "tab or line break"):
            _ = parse_result(value).table

        value = table_envelope()
        value["data"]["schema"]["fields"][1]["data_type"] = "u_int64"
        with self.assertRaisesRegex(ProtocolError, "unknown value"):
            _ = parse_result(value).table

    def test_rejects_wrong_logical_type_and_nonnullable_null(self):
        wrong_type = table_envelope()
        wrong_type["data"]["rows"][0][1] = True
        with self.assertRaisesRegex(ProtocolError, "uint64"):
            _ = parse_result(wrong_type).table

        wrong_null = table_envelope()
        wrong_null["data"]["rows"][0][0] = None
        with self.assertRaisesRegex(ProtocolError, "not nullable"):
            _ = parse_result(wrong_null).table

    def test_rejects_out_of_range_uint64(self):
        value = table_envelope()
        value["data"]["rows"][0][1] = 2**64
        with self.assertRaises(ProtocolError):
            _ = parse_result(value).table

    def test_optional_converter_has_actionable_install_hint(self):
        result = parse_result(table_envelope())
        real_import = importlib.import_module

        def missing(name, package=None):
            if name == "pandas":
                raise ModuleNotFoundError("No module named 'pandas'", name="pandas")
            return real_import(name, package)

        with patch("gravlax.results.importlib.import_module", side_effect=missing):
            with self.assertRaisesRegex(OptionalDependencyError, r"gravlax-client\[pandas\]"):
                result.to_pandas()

    def test_rejects_nonfinite_json_number(self):
        source = json.dumps(table_envelope()).replace("17", "NaN", 1)
        with self.assertRaises(ProtocolError):
            parse_result(source)

        value = table_envelope()
        value["data"]["rows"][0][3] = {"bad": float("nan")}
        with self.assertRaises(ProtocolError):
            parse_result(value)

        value = table_envelope()
        value["provenance"]["parameters"] = {"nested": [1, float("inf")]}
        with self.assertRaises(ProtocolError):
            parse_result(value)

    def test_rejects_float_that_overflows_binary64(self):
        value = table_envelope()
        value["data"]["rows"][0][2] = 10**10000
        with self.assertRaises(ProtocolError):
            _ = parse_result(value).table

    def test_integer_json_number_is_valid_for_float64(self):
        value = table_envelope()
        value["data"]["rows"][0][2] = 1
        self.assertEqual(parse_result(value).table.rows[0][2], 1)

    def test_rejects_invalid_provenance_and_blank_producer(self):
        for archives in (["sample.aie", "sample.aie"], ["  "]):
            with self.subTest(archives=archives):
                value = table_envelope()
                value["provenance"]["archives"] = archives
                with self.assertRaises(ProtocolError):
                    parse_result(value)

        value = table_envelope()
        value["provenance"]["assembly"] = "\t"
        with self.assertRaises(ProtocolError):
            parse_result(value)

        value = table_envelope()
        value["provenance"]["parameters"] = {" ": 1}
        with self.assertRaises(ProtocolError):
            parse_result(value)

        value = table_envelope()
        value["producer"]["name"] = "  "
        with self.assertRaises(ProtocolError):
            parse_result(value)


if __name__ == "__main__":
    unittest.main()
