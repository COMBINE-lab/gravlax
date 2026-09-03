from __future__ import annotations

import copy
import json
import subprocess
import unittest
from unittest.mock import patch

from gravlax import (
    ANNOTATION_CLASS_TRANSITIONS_SCHEMA,
    ANNOTATION_COMPARISON_SCHEMA,
    ANNOTATION_CONTRIBUTING_CAUSES_SCHEMA,
    ANNOTATION_COUNT_DELTAS_SCHEMA,
    ANNOTATION_WITNESSES_SCHEMA,
    Client,
    ProtocolError,
    ResolvedPlan,
    TRANSCRIPT_EQUIVALENCE_CATALOG_SCHEMA,
    TRANSCRIPT_EQUIVALENCE_COUNTS_SCHEMA,
    TRANSCRIPT_EQUIVALENCE_MEMBERSHIP_SCHEMA,
    TRANSCRIPT_EQUIVALENCE_SCHEMA,
    parse_annotation_comparison,
    parse_transcript_equivalence,
)


def _field(name: str, data_type: str = "string", nullable: bool = False):
    return {"name": name, "data_type": data_type, "nullable": nullable}


def _table(schema: str, fields):
    return {"schema": {"id": schema, "fields": [dict(field) for field in fields]}, "rows": []}


CATALOG_FIELDS = (
    _field("ec_id"),
    _field("transcript_ids", "json"),
    _field("gene_ids", "json"),
    _field("ambiguous", "boolean"),
    _field("archived_umi_class_count", "uint64"),
    _field("cell_count", "uint64"),
    _field("complete_umi_class_count", "uint64"),
)
COUNT_FIELDS = (
    _field("aggregation"),
    _field("key"),
    _field("cell_id", "uint64", True),
    _field("ec_id", "string", True),
    _field("archived_umi_class_count", "uint64"),
    _field("ambiguous_umi_class_count", "uint64"),
    _field("no_compatible_transcript_umi_class_count", "uint64"),
    _field("conflicting_umi_class_count", "uint64"),
    _field("complete_umi_class_count", "uint64"),
    _field("incomplete_umi_class_count", "uint64"),
    _field("retained_record_count", "uint64"),
    _field("represented_alignment_count", "uint64"),
)
MEMBERSHIP_FIELDS = (
    _field("umi_class", "uint64"),
    _field("cell_id", "uint64"),
    _field("barcode"),
    _field("aggregation"),
    _field("key"),
    _field("ec_id", "string", True),
    _field("retained_record_count", "uint64"),
    _field("represented_alignment_count", "uint64"),
    _field("compatible_record_count", "uint64"),
    _field("unmatched_record_count", "uint64"),
    _field("ambiguous", "boolean"),
    _field("no_compatible_transcript", "boolean"),
    _field("conflict", "boolean"),
    _field("complete_within_archive_quotient", "boolean"),
    _field("retained_representatives_complete", "boolean"),
)

COUNT_DELTA_FIELDS = (
    _field("cell", "uint64"),
    _field("cell_barcode"),
    _field("comparison_gene_id"),
    _field("annotation_a_gene_id", nullable=True),
    _field("annotation_b_gene_id", nullable=True),
    _field("annotation_a_count", "uint64"),
    _field("annotation_b_count", "uint64"),
    _field("signed_delta_b_minus_a", "int64"),
)
CLASS_TRANSITION_FIELDS = (
    _field("cell", "uint64"),
    _field("cell_barcode"),
    _field("umi_class", "uint64"),
    _field("transition_kind"),
    _field("molecule_records", "uint64"),
    _field("evidence_rows", "uint64"),
    _field("changed_evidence_rows", "uint64"),
    _field("annotation_a_selected_comparison_gene_id", nullable=True),
    _field("annotation_a_selected_gene_id", nullable=True),
    _field("annotation_a_selected_weight", "uint64"),
    _field("annotation_a_counted", "boolean"),
    _field("annotation_a_canonical_class", "uint64", True),
    _field("annotation_a_gene_support", "json"),
    _field("annotation_a_same_gene_neighbors", "json"),
    _field("annotation_b_selected_comparison_gene_id", nullable=True),
    _field("annotation_b_selected_gene_id", nullable=True),
    _field("annotation_b_selected_weight", "uint64"),
    _field("annotation_b_counted", "boolean"),
    _field("annotation_b_canonical_class", "uint64", True),
    _field("annotation_b_gene_support", "json"),
    _field("annotation_b_same_gene_neighbors", "json"),
    _field("contributing_cause_count", "uint64"),
    _field("molecule_witnesses", "uint64"),
    _field("omitted_molecule_witnesses", "uint64"),
    _field("changed_row_witnesses", "uint64"),
    _field("omitted_changed_row_witnesses", "uint64"),
)
CAUSE_FIELDS = (
    _field("cell", "uint64"),
    _field("cell_barcode"),
    _field("umi_class", "uint64"),
    _field("transition_kind"),
    _field("contributing_cause"),
    _field("nonexclusive", "boolean"),
    _field("additive_count_attribution", "boolean"),
)
WITNESS_FIELDS = (
    _field("archive_ordinal", "uint64"),
    _field("cell", "uint64"),
    _field("cell_barcode"),
    _field("umi_class", "uint64"),
    _field("chrom"),
    _field("anchor", "uint64"),
    _field("evidence_rows", "uint64"),
    _field("changed_rows_total", "uint64"),
    _field("changed_rows_omitted", "uint64"),
    _field("annotation_a_selected_comparison_gene_id", nullable=True),
    _field("annotation_a_selected_gene_id", nullable=True),
    _field("annotation_a_counted", "boolean"),
    _field("annotation_a_canonical_class", "uint64", True),
    _field("annotation_b_selected_comparison_gene_id", nullable=True),
    _field("annotation_b_selected_gene_id", nullable=True),
    _field("annotation_b_counted", "boolean"),
    _field("annotation_b_canonical_class", "uint64", True),
    _field("contributing_causes", "json"),
    _field("changed_row_witnesses", "json"),
)


def annotation_comparison_envelope():
    return {
        "$schema": "gravlax.result-envelope.v1",
        "result_schema": ANNOTATION_COMPARISON_SCHEMA,
        "producer": {"name": "aie", "version": "0.1.0"},
        "provenance": {
            "archives": ["aie-directory-root-v2:abc"],
            "assembly": "GRCh38",
            "annotations": [
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
            ],
            "parameters": {},
        },
        "warnings": [],
        "data": {
            "summary": {
                "count_delta_rows": 0,
                "class_transition_rows": 0,
                "molecule_witness_rows": 0,
            },
            "semantics": {
                "final_count_deltas_are_exact": True,
                "class_transition_ledger_is_complete": True,
                "contributing_causes_are_nonexclusive": True,
                "contributing_causes_are_additive_attributions": False,
                "annotation_order_tie_break_is_biological_change": False,
                "molecule_witnesses_are_bounded": True,
            },
            "count_deltas": _table(
                ANNOTATION_COUNT_DELTAS_SCHEMA, COUNT_DELTA_FIELDS
            ),
            "class_transitions": _table(
                ANNOTATION_CLASS_TRANSITIONS_SCHEMA, CLASS_TRANSITION_FIELDS
            ),
            "contributing_causes": _table(
                ANNOTATION_CONTRIBUTING_CAUSES_SCHEMA, CAUSE_FIELDS
            ),
            "witnesses": _table(ANNOTATION_WITNESSES_SCHEMA, WITNESS_FIELDS),
        },
    }


def transcript_ec_envelope(*, membership: bool = True):
    data = {
        "scope": {"selection": {"kind": "gene", "stable_id": "g1"}},
        "semantics": {
            "compatibility": {
                "compatibility": "exact annotated junctions",
                "alternative_placements": "union within retained record",
                "record_reduction": "intersection across global UMI class",
                "exact_scope": "retained archive quotient and supplied annotation",
                "abundance_inferred": False,
                "full_isoform_phasing_claimed": False,
            },
            "count_unit": "archived UMI class",
        },
        "summary": {
            "scoped_umi_classes": 0,
            "transcript_ecs": 0,
            "count_rows": 0,
            "membership_rows": 0,
            "assigned_umi_classes": 0,
            "unassigned_umi_classes": 0,
            "complete_umi_classes": 0,
            "incomplete_umi_classes": 0,
        },
        "catalog": _table(TRANSCRIPT_EQUIVALENCE_CATALOG_SCHEMA, CATALOG_FIELDS),
        "counts": _table(TRANSCRIPT_EQUIVALENCE_COUNTS_SCHEMA, COUNT_FIELDS),
    }
    if membership:
        data["membership"] = _table(
            TRANSCRIPT_EQUIVALENCE_MEMBERSHIP_SCHEMA, MEMBERSHIP_FIELDS
        )
    return {
        "$schema": "gravlax.result-envelope.v1",
        "result_schema": TRANSCRIPT_EQUIVALENCE_SCHEMA,
        "producer": {"name": "aie", "version": "0.1.0"},
        "provenance": {
            "archives": ["aie-directory-root-v2:abc"],
            "assembly": "GRCh38",
            "annotation": "GENCODE 49",
            "annotation_digest": "blake3:" + "a" * 64,
            "parameters": {
                "emit_membership": membership,
                "max_ecs": 100_000,
                "max_memberships": 1_000_000,
                "max_count_rows": 1_000_000,
            },
        },
        "warnings": [],
        "data": data,
    }


def _annotation_input(role: str, resource: str):
    return {
        "role": role,
        "resource": resource,
        "annotation_path": f"/work/{resource}.gtf",
        "source_path": f"/work/{resource}.gtf",
        "assembly": "GRCh38",
        "annotation": resource,
        "source_identity": {
            "scheme": "full-file-blake3-v1",
            "digest": resource + "-digest",
        },
        "expected_command_identity": {
            "scheme": "full-file-blake3-v1",
            "digest": resource + "-digest",
        },
        "compatibility": [
            {
                "resource": "archive",
                "kind": "archive",
                "status": "verified",
                "declared_assembly": "GRCh38",
                "chromosome_digest": None,
                "genome_digest": None,
                "note": "declared assembly matches annotation",
            }
        ],
    }


def resolved_v5_compare():
    return {
        "schema_version": 5,
        "plan_schema_version": 1,
        "producer": {
            "name": "aie",
            "version": "0.1.0",
            "plan_engine": "aie-declarative-plan-v1",
            "executable_identity": {
                "scheme": "full-file-blake3-v1",
                "digest": "executable",
            },
        },
        "name": "comparison",
        "source_path": "/work/plan.yaml",
        "source_digest": "source",
        "project_name": "project",
        "project_root": "/work",
        "project_manifest": "/work/aie-project.yaml",
        "project_manifest_digest": "manifest",
        "resources": {},
        "embedded_resources": {},
        "steps": [
            {
                "id": "compare",
                "kind": "compare-annotations",
                "args": ["compare-annotations", "/work/archive.aie"],
                "outputs": [],
                "annotation_inputs": [
                    _annotation_input("a", "before"),
                    _annotation_input("b", "after"),
                ],
                "annotation_comparison": {
                    "annotation_a_resource": "before",
                    "annotation_b_resource": "after",
                    "assembly": "GRCh38",
                    "gene_key": "unversioned",
                    "solo_strand": "forward",
                    "final_count_delta_semantics": "exact B-minus-A",
                    "transition_evidence_semantics": "non-additive evidence",
                },
                "output_schema_ids": [ANNOTATION_COMPARISON_SCHEMA],
                "io_estimate": {
                    "known_selected_input_bytes": 123,
                    "known_selected_input_files": 3,
                    "unknown_prior_step_outputs": 0,
                    "read_bytes_lower_bound": 0,
                    "read_bytes_upper_bound": None,
                    "bound": "known-inputs-only",
                    "note": "no execution upper bound",
                },
                "explanation": ["two complete reductions"],
            }
        ],
    }


def resolved_v5_transcript_ec():
    document = resolved_v5_compare()
    document["name"] = "transcript equivalence"
    document["steps"] = [
        {
            "id": "ecs",
            "kind": "query-transcript-ecs",
            "args": [
                "query",
                "/work/archive.aie",
                "transcript-ecs",
                "--locus",
                "chr1:0-100",
            ],
            "outputs": [],
            "annotation_inputs": [_annotation_input("query", "genes")],
            "output_schema_ids": [
                TRANSCRIPT_EQUIVALENCE_SCHEMA,
                TRANSCRIPT_EQUIVALENCE_CATALOG_SCHEMA,
                TRANSCRIPT_EQUIVALENCE_COUNTS_SCHEMA,
            ],
            "io_estimate": {
                "known_selected_input_bytes": 100,
                "known_selected_input_files": 2,
                "unknown_prior_step_outputs": 0,
                "read_bytes_lower_bound": 0,
                "read_bytes_upper_bound": None,
                "bound": "known-inputs-only",
                "note": "full scan with no execution upper bound",
            },
            "explanation": ["compatibility, not abundance"],
        }
    ]
    return document


class ScientificResultTests(unittest.TestCase):
    def test_annotation_comparison_parses_four_distinct_evidence_tables(self):
        result = parse_annotation_comparison(annotation_comparison_envelope())

        self.assertEqual(result.count_deltas.schema_id, ANNOTATION_COUNT_DELTAS_SCHEMA)
        self.assertEqual(
            result.class_transitions.schema_id, ANNOTATION_CLASS_TRANSITIONS_SCHEMA
        )
        self.assertEqual(
            result.contributing_causes.schema_id,
            ANNOTATION_CONTRIBUTING_CAUSES_SCHEMA,
        )
        self.assertEqual(result.witnesses.schema_id, ANNOTATION_WITNESSES_SCHEMA)

    def test_annotation_comparison_rejects_additive_transition_claim_or_wrong_type(self):
        additive = annotation_comparison_envelope()
        additive["data"]["semantics"][
            "contributing_causes_are_additive_attributions"
        ] = True
        with self.assertRaisesRegex(ProtocolError, "additive_attributions"):
            parse_annotation_comparison(additive)

        wrong_type = annotation_comparison_envelope()
        wrong_type["data"]["count_deltas"]["schema"]["fields"][0][
            "data_type"
        ] = "string"
        with self.assertRaisesRegex(ProtocolError, "field types"):
            parse_annotation_comparison(wrong_type)

        wrong_count = annotation_comparison_envelope()
        wrong_count["data"]["summary"]["count_delta_rows"] = 1
        with self.assertRaisesRegex(ProtocolError, "disagrees"):
            parse_annotation_comparison(wrong_count)

        wrong_role = annotation_comparison_envelope()
        wrong_role["provenance"]["annotations"][0]["role"] = "query"
        with self.assertRaisesRegex(ProtocolError, "before.*after"):
            parse_annotation_comparison(wrong_role)

    def test_annotation_comparison_checks_row_level_scientific_invariants(self):
        value = annotation_comparison_envelope()
        value["data"]["count_deltas"]["rows"] = [
            [0, "AAAC", "g1", "g1", "g1", 2, 3, 1]
        ]
        value["data"]["summary"]["count_delta_rows"] = 1
        value["data"]["contributing_causes"]["rows"] = [
            [0, "AAAC", 7, "reassigned_final_count", "candidate_set_changed", True, False]
        ]
        parse_annotation_comparison(value)

        wrong_delta = copy.deepcopy(value)
        wrong_delta["data"]["count_deltas"]["rows"][0][-1] = 2
        with self.assertRaisesRegex(ProtocolError, "B-minus-A"):
            parse_annotation_comparison(wrong_delta)

        additive = copy.deepcopy(value)
        additive["data"]["contributing_causes"]["rows"][0][-1] = True
        with self.assertRaisesRegex(ProtocolError, "non-additive"):
            parse_annotation_comparison(additive)

    def test_transcript_equivalence_parses_nested_typed_tables(self):
        result = parse_transcript_equivalence(transcript_ec_envelope())

        self.assertEqual(result.catalog.schema_id, TRANSCRIPT_EQUIVALENCE_CATALOG_SCHEMA)
        self.assertEqual(result.counts.schema_id, TRANSCRIPT_EQUIVALENCE_COUNTS_SCHEMA)
        self.assertEqual(
            result.membership.schema_id, TRANSCRIPT_EQUIVALENCE_MEMBERSHIP_SCHEMA
        )

    def test_transcript_equivalence_membership_is_optional(self):
        result = parse_transcript_equivalence(transcript_ec_envelope(membership=False))
        self.assertIsNone(result.membership)

    def test_transcript_equivalence_rejects_abundance_or_phasing_claims(self):
        for name in ("abundance_inferred", "full_isoform_phasing_claimed"):
            with self.subTest(name=name):
                value = transcript_ec_envelope()
                value["data"]["semantics"]["compatibility"][name] = True
                with self.assertRaisesRegex(ProtocolError, name):
                    parse_transcript_equivalence(value)

    def test_transcript_equivalence_rejects_wrong_nested_schema_or_columns(self):
        value = transcript_ec_envelope()
        value["data"]["catalog"]["schema"]["id"] = "gravlax.other.v1"
        with self.assertRaisesRegex(ProtocolError, "table schema"):
            parse_transcript_equivalence(value)

        value = transcript_ec_envelope()
        value["data"]["counts"]["schema"]["fields"].pop()
        with self.assertRaisesRegex(ProtocolError, "columns"):
            parse_transcript_equivalence(value)

        value = transcript_ec_envelope()
        value["data"]["summary"]["scoped_umi_classes"] = 1
        with self.assertRaisesRegex(ProtocolError, "conservation"):
            parse_transcript_equivalence(value)

        value = transcript_ec_envelope()
        value["provenance"]["parameters"]["emit_membership"] = False
        with self.assertRaisesRegex(ProtocolError, "disagrees"):
            parse_transcript_equivalence(value)


class ScientificClientTests(unittest.TestCase):
    @patch("gravlax.client.subprocess.run")
    def test_compare_annotations_accepts_zero_witness_caps_and_returns_bundle(self, run):
        run.return_value = subprocess.CompletedProcess(
            [], 0, json.dumps(annotation_comparison_envelope()), ""
        )

        result = Client(binary=("wrapper", "aie")).compare_annotations(
            "sample.aie",
            "before.gtf",
            "after.gtf",
            assembly="GRCh38",
            annotation_a_label="before",
            annotation_b_label="after",
            annotation_a_digest="blake3:" + "a" * 64,
            annotation_b_digest="blake3:" + "b" * 64,
            gene_key="exact",
            solo_strand="reverse",
            max_molecule_witnesses=0,
            max_row_transitions_per_molecule=0,
            allow_identical=True,
        )

        self.assertEqual(result.envelope.result_schema, ANNOTATION_COMPARISON_SCHEMA)
        self.assertEqual(
            run.call_args.args[0],
            (
                "wrapper",
                "aie",
                "compare-annotations",
                "sample.aie",
                "--annotation-a=before.gtf",
                "--annotation-b=after.gtf",
                "--assembly=GRCh38",
                "--annotation-a-label=before",
                "--annotation-b-label=after",
                "--gene-key=exact",
                "--solo-strand=reverse",
                "--max-molecule-witnesses=0",
                "--max-row-transitions-per-molecule=0",
                "--format=json",
                "--annotation-a-digest=blake3:" + "a" * 64,
                "--annotation-b-digest=blake3:" + "b" * 64,
                "--allow-identical",
            ),
        )
        self.assertFalse(run.call_args.kwargs["shell"])

    @patch("gravlax.client.subprocess.run")
    def test_transcript_ecs_uses_literal_argv_and_returns_typed_bundle(self, run):
        run.return_value = subprocess.CompletedProcess(
            [], 0, json.dumps(transcript_ec_envelope()), ""
        )

        result = Client(binary=("wrapper", "aie")).transcript_ecs(
            "sample archive.aie",
            "GENCODE 49.gtf",
            assembly="GRCh38",
            annotation_label="GENCODE 49",
            feature="gene:TP53",
            groups="groups.tsv",
            aggregation="group",
            emit_membership=True,
        )

        self.assertIsNotNone(result.membership)
        self.assertEqual(
            run.call_args.args[0],
            (
                "wrapper",
                "aie",
                "query",
                "sample archive.aie",
                "transcript-ecs",
                "--annotation-file=GENCODE 49.gtf",
                "--assembly=GRCh38",
                "--annotation-label=GENCODE 49",
                "--solo-strand=forward",
                "--max-ecs=100000",
                "--max-memberships=1000000",
                "--format=json",
                "--feature=gene:TP53",
                "--groups=groups.tsv",
                "--agg=group",
                "--emit-membership",
            ),
        )
        self.assertFalse(run.call_args.kwargs["shell"])

    @patch("gravlax.client.subprocess.run")
    def test_scientific_client_validation_fails_before_execution(self, run):
        client = Client()
        with self.assertRaisesRegex(ValueError, "exactly one"):
            client.transcript_ecs(
                "sample.aie",
                "genes.gtf",
                assembly="GRCh38",
                annotation_label="release",
            )
        with self.assertRaisesRegex(ValueError, "mutually exclusive"):
            client.transcript_ecs(
                "sample.aie",
                "genes.gtf",
                assembly="GRCh38",
                annotation_label="release",
                locus="chr1:0-10",
                cells="cells.tsv",
                groups="groups.tsv",
            )
        with self.assertRaisesRegex(ValueError, "nonnegative integer"):
            client.compare_annotations(
                "sample.aie",
                "a.gtf",
                "b.gtf",
                assembly="GRCh38",
                annotation_a_label="A",
                annotation_b_label="B",
                max_molecule_witnesses=-1,
            )
        run.assert_not_called()


class ResolvedPlanV5Tests(unittest.TestCase):
    def test_v5_parses_role_named_comparison_inputs(self):
        plan = ResolvedPlan.from_mapping(resolved_v5_compare())
        self.assertEqual([item.role for item in plan.steps[0].annotation_inputs], ["a", "b"])
        self.assertEqual(plan.steps[0].annotation_comparison.assembly, "GRCh38")
        self.assertIsNone(plan.steps[0].io_estimate.read_bytes_upper_bound)

    def test_v5_rejects_missing_duplicate_or_inconsistent_roles(self):
        missing = resolved_v5_compare()
        missing["steps"][0]["annotation_inputs"].pop()
        with self.assertRaisesRegex(ProtocolError, "roles 'a' and 'b'"):
            ResolvedPlan.from_mapping(missing)

        duplicate = resolved_v5_compare()
        duplicate["steps"][0]["annotation_inputs"][1]["role"] = "a"
        with self.assertRaisesRegex(ProtocolError, "duplicate roles"):
            ResolvedPlan.from_mapping(duplicate)

        mismatch = resolved_v5_compare()
        mismatch["steps"][0]["annotation_comparison"]["assembly"] = "GRCh37"
        with self.assertRaisesRegex(ProtocolError, "disagrees"):
            ResolvedPlan.from_mapping(mismatch)

    def test_v4_remains_compatible_without_v5_fields(self):
        legacy = copy.deepcopy(resolved_v5_compare())
        legacy["schema_version"] = 4
        legacy["steps"][0].pop("annotation_inputs")
        legacy["steps"][0].pop("annotation_comparison")
        plan = ResolvedPlan.from_mapping(legacy)
        self.assertEqual(plan.steps[0].annotation_inputs, ())
        self.assertIsNone(plan.steps[0].annotation_comparison)

    def test_v5_parses_raw_locus_transcript_equivalence_contract(self):
        plan = ResolvedPlan.from_mapping(resolved_v5_transcript_ec())
        self.assertEqual(plan.steps[0].annotation_inputs[0].role, "query")
        self.assertIsNone(plan.steps[0].biological_intent)
        self.assertEqual(
            plan.steps[0].output_schema_ids[0], TRANSCRIPT_EQUIVALENCE_SCHEMA
        )

    def test_v5_rejects_unbound_scientific_inputs_and_foreign_schemas(self):
        unbound = resolved_v5_transcript_ec()
        unbound["steps"][0]["annotation_inputs"][0]["compatibility"] = []
        with self.assertRaisesRegex(ProtocolError, "must disclose"):
            ResolvedPlan.from_mapping(unbound)

        wrong_identity = resolved_v5_transcript_ec()
        wrong_identity["steps"][0]["annotation_inputs"][0]["source_identity"][
            "scheme"
        ] = "tree-blake3-v1"
        with self.assertRaisesRegex(ProtocolError, "full-file-blake3-v1"):
            ResolvedPlan.from_mapping(wrong_identity)

        foreign = resolved_v5_transcript_ec()
        foreign["steps"][0]["output_schema_ids"] = ["gravlax.abundance.v1"]
        with self.assertRaisesRegex(ProtocolError, "transcript-equivalence"):
            ResolvedPlan.from_mapping(foreign)


if __name__ == "__main__":
    unittest.main()
