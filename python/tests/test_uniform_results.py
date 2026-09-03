from __future__ import annotations

import copy
import json
import unittest
from typing import Callable

from gravlax import (
    BundleSummary,
    DeferredSelection,
    ExactSelection,
    JunctionQuerySummary,
    ProtocolError,
    RegionQuerySummary,
    RowSemantics,
    SortDirection,
    parse_result,
    parse_uniform_bundle,
)


def region_bundle() -> dict:
    return {
        "$schema": "gravlax.result-envelope.v1",
        "result_schema": "gravlax.query.region.result.v1",
        "producer": {"name": "aie", "version": "0.1.0"},
        "provenance": {
            "archives": ["aie-directory-root-v2:abc"],
            "parameters": {
                "aggregation": "cell",
                "archive_access": "range-index-selected archive chunks",
                "selection_policy": {
                    "requested_top": 1,
                    "top_zero_means_all": True,
                    "comparator": "umis descending, entity ascending (barcode)",
                },
            },
        },
        "warnings": [],
        "data": {
            "summary": {
                "coordinates": "0-based half-open",
                "anchor_semantics": True,
                "chrom": "chr16",
                "start": 89550000,
                "end": 89575000,
                "molecules": 14,
                "umis": 12,
                "cells": 2,
                "chunks_decoded": 3,
            },
            "tables": [
                {
                    "name": "counts",
                    "schema": {
                        "id": "gravlax.query.region.counts.v1",
                        "fields": [
                            {
                                "name": "aggregation",
                                "data_type": "string",
                                "nullable": False,
                            },
                            {
                                "name": "entity",
                                "data_type": "string",
                                "nullable": False,
                            },
                            {
                                "name": "umis",
                                "data_type": "uint64",
                                "nullable": False,
                            },
                            {
                                "name": "cells",
                                "data_type": "uint64",
                                "nullable": True,
                            },
                            {
                                "name": "selected_cells",
                                "data_type": "uint64",
                                "nullable": True,
                            },
                        ],
                        "semantics": {
                            "row_semantics": "set",
                            "key": ["aggregation", "entity"],
                        },
                    },
                    "selection": {
                        "available_rows": 2,
                        "emitted_rows": 1,
                        "truncated": True,
                    },
                    "rows": [["cell", "AACCGGTT", 7, None, None]],
                }
            ],
        },
    }


def junction_bundle() -> dict:
    value = region_bundle()
    value["result_schema"] = "gravlax.query.junction.result.v1"
    value["data"]["summary"] = {
        "coordinates": "0-based half-open junction boundaries",
        "chrom": "chr16",
        "donor": 89562391,
        "acceptor": 89562883,
        "archive_supporting_children": 9,
        "archive_posting_chunks": 2,
        "umis": 12,
        "cells": 2,
    }
    value["data"]["tables"][0]["schema"]["id"] = (
        "gravlax.query.junction.counts.v1"
    )
    return value


def future_bundle() -> dict:
    return {
        "$schema": "gravlax.result-envelope.v1",
        "result_schema": "gravlax.query.future.result.v1",
        "producer": {"name": "aie", "version": "0.1.0"},
        "provenance": {"archives": [], "parameters": {}},
        "warnings": [],
        "data": {
            "summary": {"scanned": 100, "method": "one-pass"},
            "tables": [
                {
                    "name": "records",
                    "schema": {
                        "id": "gravlax.query.future.records.v1",
                        "fields": [
                            {
                                "name": "ordinal",
                                "data_type": "uint64",
                                "nullable": False,
                            },
                            {
                                "name": "label",
                                "data_type": "string",
                                "nullable": False,
                            },
                        ],
                        "semantics": {
                            "row_semantics": "sequence",
                            "ordered_by": [
                                {"field": "ordinal", "direction": "ascending"}
                            ],
                        },
                    },
                    "rows": [[0, "first"]],
                    "selection": {
                        "available_rows": None,
                        "emitted_rows": 1,
                        "truncated": None,
                    },
                }
            ],
        },
    }


class UniformResultTests(unittest.TestCase):
    def test_parses_known_region_bundle_with_typed_summary_and_exact_selection(self):
        result = parse_uniform_bundle(json.dumps(region_bundle()))

        self.assertIsInstance(result.summary, RegionQuerySummary)
        assert isinstance(result.summary, RegionQuerySummary)
        self.assertEqual(result.summary.umis, 12)
        self.assertEqual(result.table_names, ("counts",))
        counts = result.table("counts")
        self.assertEqual(counts.records()[0]["entity"], "AACCGGTT")
        self.assertEqual(counts.semantics.row_semantics, RowSemantics.SET)
        self.assertEqual(counts.semantics.key, ("aggregation", "entity"))
        self.assertIsInstance(counts.selection, ExactSelection)
        assert isinstance(counts.selection, ExactSelection)
        self.assertEqual(counts.selection.available_rows, 2)
        self.assertTrue(counts.selection.truncated)

        envelope = parse_result(region_bundle())
        self.assertEqual(envelope.bundle.table_names, ("counts",))

    def test_parses_known_junction_summary(self):
        result = parse_uniform_bundle(junction_bundle())
        self.assertIsInstance(result.summary, JunctionQuerySummary)
        assert isinstance(result.summary, JunctionQuerySummary)
        self.assertEqual((result.summary.donor, result.summary.acceptor), (89562391, 89562883))

    def test_parses_generic_summary_and_deferred_selection(self):
        result = parse_uniform_bundle(future_bundle())

        self.assertIsInstance(result.summary, BundleSummary)
        assert isinstance(result.summary, BundleSummary)
        self.assertEqual(result.summary["scanned"], 100)
        table = result.table("records")
        self.assertEqual(table.semantics.row_semantics, RowSemantics.SEQUENCE)
        self.assertEqual(
            table.semantics.ordered_by[0].direction,
            SortDirection.ASCENDING,
        )
        self.assertIsInstance(table.selection, DeferredSelection)
        assert isinstance(table.selection, DeferredSelection)
        self.assertIsNone(table.selection.available_rows)
        self.assertIsNone(table.selection.truncated)

    def test_rejects_invalid_or_duplicate_bundle_identifiers(self):
        cases: list[tuple[str, Callable[[], object]]] = []

        bad_outer = future_bundle()
        bad_outer["result_schema"] = "../../result"
        cases.append(("outer", lambda: parse_uniform_bundle(bad_outer)))

        bad_name = future_bundle()
        bad_name["data"]["tables"][0]["name"] = "bad name"
        cases.append(("name", lambda: parse_uniform_bundle(bad_name)))

        bad_schema = future_bundle()
        bad_schema["data"]["tables"][0]["schema"]["id"] = "bad/schema"
        cases.append(("schema", lambda: parse_uniform_bundle(bad_schema)))

        duplicate = future_bundle()
        duplicate["data"]["tables"].append(
            copy.deepcopy(duplicate["data"]["tables"][0])
        )
        cases.append(("duplicate", lambda: parse_uniform_bundle(duplicate)))

        for label, parse in cases:
            with self.subTest(label=label), self.assertRaises(ProtocolError):
                parse()

    def test_rejects_bad_semantic_references_and_values(self):
        mutations = [
            {"key": []},
            {"key": ["missing"]},
            {"key": ["ordinal", "ordinal"]},
            {"ordered_by": []},
            {
                "ordered_by": [
                    {"field": "ordinal", "direction": "sideways"}
                ]
            },
            {
                "ordered_by": [
                    {"field": "missing", "direction": "ascending"}
                ]
            },
        ]
        for mutation in mutations:
            with self.subTest(mutation=mutation):
                value = future_bundle()
                semantics = value["data"]["tables"][0]["schema"]["semantics"]
                semantics.update(mutation)
                with self.assertRaises(ProtocolError):
                    parse_uniform_bundle(value)

        value = future_bundle()
        value["data"]["tables"][0]["schema"]["semantics"] = None
        with self.assertRaisesRegex(ProtocolError, "semantics is required"):
            parse_uniform_bundle(value)

    def test_rejects_selection_incoherence_and_row_count_mismatch(self):
        mutations = [
            {"available_rows": 0, "emitted_rows": 1, "truncated": False},
            {"available_rows": 2, "emitted_rows": 1, "truncated": False},
            {"available_rows": None, "emitted_rows": 1, "truncated": False},
            {"available_rows": 2, "emitted_rows": 2, "truncated": False},
            {"available_rows": 2**64, "emitted_rows": 1, "truncated": True},
            {"available_rows": 2, "emitted_rows": True, "truncated": True},
        ]
        for mutation in mutations:
            with self.subTest(mutation=mutation):
                value = region_bundle()
                value["data"]["tables"][0]["selection"] = mutation
                with self.assertRaises(ProtocolError):
                    parse_uniform_bundle(value)

    def test_rejects_duplicate_set_keys_and_wrong_declared_cell_types(self):
        duplicate = region_bundle()
        table = duplicate["data"]["tables"][0]
        table["rows"].append(copy.deepcopy(table["rows"][0]))
        table["selection"] = {
            "available_rows": 2,
            "emitted_rows": 2,
            "truncated": False,
        }
        with self.assertRaisesRegex(ProtocolError, "duplicates"):
            parse_uniform_bundle(duplicate)

        wrong_type = region_bundle()
        wrong_type["data"]["tables"][0]["rows"][0][2] = "seven"
        with self.assertRaisesRegex(ProtocolError, "uint64"):
            parse_uniform_bundle(wrong_type)

    def test_known_schema_ids_bind_the_summary_and_table_shape(self):
        no_summary = region_bundle()
        no_summary["data"].pop("summary")
        with self.assertRaisesRegex(ProtocolError, "typed summary"):
            parse_uniform_bundle(no_summary)

        wrong_summary = region_bundle()
        wrong_summary["data"]["summary"]["cells"] = -1
        with self.assertRaises(ProtocolError):
            parse_uniform_bundle(wrong_summary)

        wrong_table_id = region_bundle()
        wrong_table_id["data"]["tables"][0]["schema"]["id"] = (
            "gravlax.query.junction.counts.v1"
        )
        with self.assertRaisesRegex(ProtocolError, "requires table schema"):
            parse_uniform_bundle(wrong_table_id)

        wrong_shape = region_bundle()
        wrong_shape["data"]["tables"][0]["schema"]["fields"][2][
            "data_type"
        ] = "int64"
        with self.assertRaisesRegex(ProtocolError, "field types"):
            parse_uniform_bundle(wrong_shape)

    def test_rejects_malicious_strict_json_documents(self):
        source = json.dumps(future_bundle())
        malicious = source.replace(
            '"result_schema": "gravlax.query.future.result.v1",',
            '"result_schema": "gravlax.query.future.result.v1", '
            '"result_schema": "gravlax.query.future.result.v1",',
            1,
        )
        with self.assertRaisesRegex(ProtocolError, "duplicate JSON key"):
            parse_uniform_bundle(malicious)

        nonfinite = source.replace('"scanned": 100', '"scanned": NaN', 1)
        with self.assertRaisesRegex(ProtocolError, "non-finite"):
            parse_uniform_bundle(nonfinite)


if __name__ == "__main__":
    unittest.main()
