from __future__ import annotations

import errno
import json
import os
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import gravlax.client as client_module
from gravlax import Client, CommandError, ProtocolError, ResolvedPlan


def completed(argv, stdout: str, *, returncode: int = 0, stderr: str = ""):
    return subprocess.CompletedProcess(argv, returncode, stdout, stderr)


class ClientTests(unittest.TestCase):
    @patch("gravlax.client.subprocess.run")
    def test_project_show_uses_literal_argv_without_a_shell(self, run):
        response = {
            "schema_version": 1,
            "name": "pbmc",
            "root": "/work/project",
            "manifest": "/work/project/aie-project.yaml",
            "resources": [
                {
                    "name": "sample",
                    "kind": "archive",
                    "path": "data/sample.aie",
                    "external": False,
                    "resolved_path": "/work/project/data/sample.aie",
                    "status": "ok",
                }
            ],
        }
        run.return_value = completed([], json.dumps(response))
        suspicious = Path("project with spaces; touch never")

        project = Client(binary=("wrapper", "aie")).project_show(project=suspicious)

        self.assertEqual(project.resources[0].name, "sample")
        argv = run.call_args.args[0]
        self.assertEqual(
            argv,
            (
                "wrapper",
                "aie",
                "project",
                "show",
                "--json",
                "--project=project with spaces; touch never",
            ),
        )
        self.assertIs(run.call_args.kwargs["shell"], False)

    @patch("gravlax.client.subprocess.run")
    def test_plan_check_returns_typed_plan_and_explanation(self, run):
        response = {
            "schema_version": 4,
            "plan_schema_version": 1,
            "producer": {
                "name": "aie",
                "version": "0.1.0",
                "plan_engine": "aie-declarative-plan-v1",
                "executable_identity": {
                    "scheme": "full-file-blake3-v1",
                    "digest": "executable-digest",
                },
            },
            "name": "replay",
            "source_path": "/work/plans/replay.yaml",
            "source_digest": "abc123",
            "project_name": "pbmc",
            "project_root": "/work",
            "project_manifest": "/work/aie-project.yaml",
            "project_manifest_digest": "def456",
            "resources": {
                "sample": {
                    "kind": "archive",
                    "path": "/work/data/sample.aie",
                    "external": False,
                    "bytes": 42,
                    "identity": {"scheme": "blake3", "digest": "aabbcc"},
                    "assembly": "GRCh38.p14",
                }
            },
            "embedded_resources": {
                "design/sample/archive": {
                    "owner_resource": "design",
                    "sample": "sample",
                    "role": "archive",
                    "declared_path": "../data/sample.aie",
                    "kind": "archive",
                    "path": "/work/data/sample.aie",
                    "external": False,
                    "bytes": 42,
                    "identity": {"scheme": "blake3", "digest": "aabbcc"},
                }
            },
            "steps": [
                {
                    "id": "inspect",
                    "kind": "inspect-archive",
                    "args": ["inspect-archive", "/work/data/sample.aie", "--json"],
                    "stdout": "/work/results/inspect.json",
                    "outputs": [
                        {
                            "name": "output",
                            "path": "/work/results/inspect.json",
                            "kind": "file",
                            "resource_kind": "file",
                            "staging_path": "/work/results/.inspect.aie-stage-inspect.json",
                        }
                    ],
                    "input_resources": ["sample"],
                    "step_inputs": [
                        {
                            "reference": "step:prepare:archive",
                            "producer_step": "prepare",
                            "output_name": "archive",
                            "path": "/work/results/prepared.aie",
                            "kind": "file",
                            "resource_kind": "archive",
                        }
                    ],
                    "embedded_resources": ["design/sample/archive"],
                    "prepared_inputs": [
                        {
                            "path": "/work/.aie/resolved-inputs/design.tsv",
                            "bytes": 8,
                            "identity": {
                                "scheme": "full-file-blake3-v1",
                                "digest": "prepared-digest",
                            },
                            "content": "design\n",
                        }
                    ],
                    "output_schema_ids": ["gravlax.archive.identity.v1"],
                    "io_estimate": {
                        "known_selected_input_bytes": 50,
                        "known_selected_input_files": 2,
                        "unknown_prior_step_outputs": 1,
                        "read_bytes_lower_bound": 0,
                        "read_bytes_upper_bound": None,
                        "bound": "known-inputs-only",
                        "note": "known selected file sizes only",
                    },
                    "explanation": ["archive sample -> /work/data/sample.aie"],
                },
                {
                    "id": "tp53",
                    "kind": "query-region",
                    "args": [
                        "query",
                        "/work/data/sample.aie",
                        "region",
                        "chr17:7661778-7687550",
                    ],
                    "stdout": "/work/results/tp53.json",
                    "outputs": [
                        {
                            "name": "output",
                            "path": "/work/results/tp53.json",
                            "kind": "file",
                            "resource_kind": "file",
                            "staging_path": "/work/results/.tp53.aie-stage-tp53.json",
                        }
                    ],
                    "input_resources": ["sample", "genes"],
                    "step_inputs": [],
                    "embedded_resources": [],
                    "prepared_inputs": [],
                    "biological_intent": {
                        "requested": "TP53",
                        "resolved_kind": "gene",
                        "stable_id": "ENSG00000141510.18",
                        "display_name": "TP53",
                        "matched_by": "gene_symbol",
                        "gene_ids": ["ENSG00000141510.18"],
                        "transcript_ids": [],
                        "annotation_resource": "genes",
                        "annotation_path": "/work/refs/gencode.gtf",
                        "assembly": "GRCh38.p14",
                        "annotation": "GENCODE 48",
                        "annotation_digest": "blake3:annotation-digest",
                        "contig": "chr17",
                        "start": 7661778,
                        "end": 7687550,
                        "strand": "reverse",
                        "locus": "chr17:7661778-7687550",
                        "compatibility": [
                            {
                                "resource": "sample",
                                "kind": "archive",
                                "status": "verified",
                                "declared_assembly": "GRCh38.p14",
                                "chromosome_digest": "chromosome-digest",
                                "genome_digest": None,
                                "note": "declared assembly matches annotation",
                            }
                        ],
                    },
                    "output_schema_ids": ["gravlax.query.region.v1"],
                    "io_estimate": {
                        "known_selected_input_bytes": 42,
                        "known_selected_input_files": 1,
                        "unknown_prior_step_outputs": 0,
                        "read_bytes_lower_bound": 0,
                        "read_bytes_upper_bound": 42,
                        "bound": "whole-selected-files",
                        "note": "whole selected archive bounds bytes read",
                    },
                    "explanation": ["feature TP53 -> chr17:7661778-7687550"],
                },
            ],
        }
        run.return_value = completed([], json.dumps(response), stderr="Plan replay\n")

        plan = Client().plan_check("plans/replay.yaml", explain=True)

        self.assertEqual(plan.resources["sample"].bytes, 42)
        self.assertEqual(plan.producer.plan_engine, "aie-declarative-plan-v1")
        self.assertEqual(
            plan.embedded_resources["design/sample/archive"].declared_path,
            "../data/sample.aie",
        )
        self.assertEqual(plan.steps[0].args[0], "inspect-archive")
        self.assertEqual(plan.steps[0].input_resources, ("sample",))
        self.assertEqual(plan.steps[0].step_inputs[0].producer_step, "prepare")
        self.assertEqual(plan.steps[0].prepared_inputs[0].content, "design\n")
        self.assertEqual(plan.resources["sample"].assembly, "GRCh38.p14")
        self.assertEqual(
            plan.steps[0].output_schema_ids, ("gravlax.archive.identity.v1",)
        )
        self.assertEqual(plan.steps[0].io_estimate.unknown_prior_step_outputs, 1)
        self.assertEqual(plan.steps[1].biological_intent.stable_id, "ENSG00000141510.18")
        self.assertEqual(
            plan.steps[1].biological_intent.compatibility[0].status, "verified"
        )
        self.assertEqual(
            plan.steps[0].outputs[0].staging_path,
            Path("/work/results/.inspect.aie-stage-inspect.json"),
        )
        self.assertEqual(plan.explanation_text, "Plan replay\n")
        self.assertEqual(
            run.call_args.args[0][-3:], ("--explain", "--", "plans/replay.yaml")
        )

        legacy = json.loads(json.dumps(response))
        legacy["schema_version"] = 3
        legacy["steps"][0].pop("output_schema_ids")
        legacy["steps"][0].pop("io_estimate")
        legacy_plan = ResolvedPlan.from_mapping(legacy)
        self.assertEqual(legacy_plan.schema_version, 3)
        self.assertEqual(legacy_plan.steps[0].output_schema_ids, ())
        self.assertIsNone(legacy_plan.steps[0].io_estimate)

    @patch("gravlax.client.subprocess.run")
    def test_project_add_passes_annotation_identity_as_literal_flags(self, run):
        run.return_value = completed([], "registered genes\n")

        Client(binary=("wrapper", "aie")).project_add(
            "genes",
            "refs/genes.gtf",
            kind="annotation",
            project="project with spaces",
            assembly="GRCh38.p14",
            annotation_label="GENCODE 48",
        )

        self.assertEqual(
            run.call_args.args[0],
            (
                "wrapper",
                "aie",
                "project",
                "add",
                "--kind=annotation",
                "--project=project with spaces",
                "--assembly=GRCh38.p14",
                "--annotation-label=GENCODE 48",
                "--",
                "genes",
                "refs/genes.gtf",
            ),
        )
        self.assertIs(run.call_args.kwargs["shell"], False)

    @patch("gravlax.client.subprocess.run")
    def test_doctor_returns_json_report_on_expected_nonzero_status(self, run):
        response = {
            "schema": "gravlax.doctor.v1",
            "version": "0.1.0",
            "ok": False,
            "strict": False,
            "summary": {"passed": 0, "warnings": 0, "failures": 1},
            "checks": [
                {
                    "id": "archive:missing.aie",
                    "status": "fail",
                    "summary": "missing",
                    "detail": None,
                    "remedy": "correct the path",
                    "data": None,
                }
            ],
        }
        run.return_value = completed(
            [], json.dumps(response), returncode=1, stderr="doctor found 1 failure\n"
        )

        report = Client().doctor(["missing.aie"])

        self.assertFalse(report.ok)
        self.assertEqual(report.exit_code, 1)
        self.assertEqual(report.checks[0].remedy, "correct the path")

    @patch("gravlax.client.subprocess.run")
    def test_non_json_failed_doctor_is_a_command_error(self, run):
        run.return_value = completed([], "", returncode=2, stderr="bad option")
        with self.assertRaises(CommandError):
            Client().doctor()

    @patch("gravlax.client.subprocess.run")
    def test_doctor_rejects_incoherent_ok_flag(self, run):
        response = {
            "schema": "gravlax.doctor.v1",
            "version": "0.1.0",
            "ok": True,
            "strict": False,
            "summary": {"passed": 0, "warnings": 0, "failures": 1},
            "checks": [
                {
                    "id": "archive:broken.aie",
                    "status": "fail",
                    "summary": "broken",
                    "detail": None,
                    "remedy": None,
                    "data": None,
                }
            ],
        }
        run.return_value = completed([], json.dumps(response))

        with self.assertRaisesRegex(ProtocolError, "ok flag"):
            Client().doctor()

    @patch("gravlax.client.subprocess.run")
    def test_doctor_rejects_success_document_from_failed_process(self, run):
        response = {
            "schema": "gravlax.doctor.v1",
            "version": "0.1.0",
            "ok": True,
            "strict": False,
            "summary": {"passed": 1, "warnings": 0, "failures": 0},
            "checks": [
                {
                    "id": "installation",
                    "status": "pass",
                    "summary": "ready",
                    "detail": None,
                    "remedy": None,
                    "data": None,
                }
            ],
        }
        run.return_value = completed([], json.dumps(response), returncode=1)

        with self.assertRaises(CommandError):
            Client().doctor()

    @patch("gravlax.client.subprocess.run")
    def test_doctor_strict_warning_is_not_ok(self, run):
        response = {
            "schema": "gravlax.doctor.v1",
            "version": "0.1.0",
            "ok": True,
            "strict": True,
            "summary": {"passed": 0, "warnings": 1, "failures": 0},
            "checks": [
                {
                    "id": "optional",
                    "status": "warn",
                    "summary": "not configured",
                    "detail": None,
                    "remedy": None,
                    "data": None,
                }
            ],
        }
        run.return_value = completed([], json.dumps(response))

        with self.assertRaisesRegex(ProtocolError, "ok flag"):
            Client().doctor()

    @patch("gravlax.client.subprocess.run")
    def test_run_json_rejects_duplicate_keys(self, run):
        run.return_value = completed([], '{"schema_version":1,"schema_version":1}')
        with self.assertRaises(ProtocolError):
            Client().run_json(["project", "show", "--json"])

    def test_command_strings_are_not_accepted(self):
        with self.assertRaises(TypeError):
            Client().run("doctor --json")  # type: ignore[arg-type]

    def test_reads_versioned_step_completion(self):
        response = {
            "schema_version": 2,
            "resolved_plan_digest": "plan-digest",
            "step_id": "quantify",
            "step_digest": "step-digest",
            "inputs": [
                {
                    "producer_step": "ingest",
                    "output_name": "archive",
                    "path": "/work/results/sample.aie",
                    "kind": "file",
                    "bytes": 99,
                    "identity": {
                        "scheme": "full-file-blake3-v1",
                        "digest": "input-digest",
                    },
                }
            ],
            "outputs": [
                {
                    "path": "/work/results/counts",
                    "kind": "directory",
                    "bytes": 123,
                    "identity": {
                        "scheme": "tree-blake3-v1",
                        "digest": "output-digest",
                    },
                }
            ],
        }
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "quantify.json"
            path.write_text(json.dumps(response), encoding="utf-8")
            completion = Client().step_completion(path)

        self.assertEqual(completion.step_id, "quantify")
        self.assertEqual(completion.inputs[0].producer_step, "ingest")
        self.assertEqual(completion.outputs[0].identity.scheme, "tree-blake3-v1")

    @patch("gravlax.client.subprocess.run")
    def test_plan_resume_is_a_literal_flag(self, run):
        run.return_value = completed([], "completed\n")
        Client().plan_run("plans/run.yaml", resume=True)
        self.assertEqual(
            run.call_args.args[0][-3:], ("--resume", "--", "plans/run.yaml")
        )

    @patch("gravlax.client.subprocess.run")
    def test_run_to_file_streams_and_installs_only_completed_output(self, run):
        def invoke(argv, **kwargs):
            kwargs["stdout"].write(b"one\ntwo\n")
            return completed(argv, "", stderr="note\n")

        run.side_effect = invoke
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "rows.tsv"
            result = Client().run_to_file(["query", "--tsv"], output)
            self.assertEqual(output.read_bytes(), b"one\ntwo\n")
            self.assertEqual(result.output_path, output)
            self.assertEqual(result.bytes, 8)
            self.assertEqual(list(output.parent.glob("*.gravlax-tmp")), [])
        self.assertIs(run.call_args.kwargs["shell"], False)

    @patch("gravlax.client.subprocess.run")
    def test_run_to_file_keeps_secure_temporary_descriptor(self, run):
        events: list[str] = []

        def invoke(argv, **kwargs):
            events.append("run")
            kwargs["stdout"].write(b"complete\n")
            return completed(argv, "")

        actual_close = os.close

        def record_close(descriptor):
            events.append("close")
            actual_close(descriptor)

        run.side_effect = invoke
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "rows.tsv"
            with patch("gravlax.client.os.close", side_effect=record_close):
                Client().run_to_file(["query", "--tsv"], output)
            self.assertEqual(output.read_bytes(), b"complete\n")
        # Some Python builds close FileIO descriptors below the os module. If os.close is
        # observable, it must happen only after the child has received the original descriptor.
        self.assertEqual(events[0], "run")

    @patch("gravlax.client.subprocess.run")
    def test_run_to_file_removes_partial_output_after_failure(self, run):
        def invoke(argv, **kwargs):
            kwargs["stdout"].write(b"partial")
            return completed(argv, "", returncode=1, stderr="failed")

        run.side_effect = invoke
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "rows.tsv"
            with self.assertRaises(CommandError):
                Client().run_to_file(["query", "--tsv"], output)
            self.assertFalse(output.exists())
            self.assertEqual(list(output.parent.glob("*.gravlax-tmp")), [])

    @patch("gravlax.client.subprocess.run")
    def test_run_to_file_refuses_existing_destination_before_execution(self, run):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "rows.tsv"
            output.write_text("keep", encoding="utf-8")
            with self.assertRaises(FileExistsError):
                Client().run_to_file(["query", "--tsv"], output)
            self.assertEqual(output.read_text(encoding="utf-8"), "keep")
        run.assert_not_called()

    @patch("gravlax.client.subprocess.run")
    def test_run_to_file_publishes_held_inode_after_staging_name_swap(self, run):
        with tempfile.TemporaryDirectory() as directory:
            parent = Path(directory)
            output = parent / "rows.tsv"

            def invoke(argv, **kwargs):
                kwargs["stdout"].write(b"trusted\n")
                kwargs["stdout"].flush()
                staging = next(parent.glob(".rows.tsv.*.gravlax-tmp"))
                displaced = parent / "displaced-original"
                os.replace(staging, displaced)
                staging.write_bytes(b"foreign\n")
                return completed(argv, "")

            run.side_effect = invoke
            result = Client().run_to_file(["query", "--tsv"], output)

            self.assertEqual(output.read_bytes(), b"trusted\n")
            self.assertEqual(result.bytes, len(b"trusted\n"))
            replacements = list(parent.glob(".rows.tsv.*.gravlax-tmp"))
            self.assertEqual(len(replacements), 1)
            self.assertEqual(replacements[0].read_bytes(), b"foreign\n")

    @patch("gravlax.client.subprocess.run")
    def test_run_to_file_replace_uses_held_inode_after_staging_swap(self, run):
        with tempfile.TemporaryDirectory() as directory:
            parent = Path(directory)
            output = parent / "rows.tsv"
            output.write_bytes(b"old\n")

            def invoke(argv, **kwargs):
                kwargs["stdout"].write(b"trusted replacement\n")
                kwargs["stdout"].flush()
                staging = next(parent.glob(".rows.tsv.*.gravlax-tmp"))
                displaced = parent / "displaced-original"
                os.replace(staging, displaced)
                staging.write_bytes(b"foreign replacement\n")
                return completed(argv, "")

            run.side_effect = invoke
            Client().run_to_file(["query", "--tsv"], output, replace=True)

            self.assertEqual(output.read_bytes(), b"trusted replacement\n")
            replacements = list(parent.glob(".rows.tsv.*.gravlax-tmp"))
            self.assertEqual(len(replacements), 1)
            self.assertEqual(replacements[0].read_bytes(), b"foreign replacement\n")

    @patch("gravlax.client.subprocess.run")
    def test_run_to_file_portable_fallback_fails_closed_on_staging_swap(self, run):
        with tempfile.TemporaryDirectory() as directory:
            parent = Path(directory)
            output = parent / "rows.tsv"

            def invoke(argv, **kwargs):
                kwargs["stdout"].write(b"trusted\n")
                kwargs["stdout"].flush()
                staging = next(parent.glob(".rows.tsv.*.gravlax-tmp"))
                displaced = parent / "displaced-original"
                os.replace(staging, displaced)
                staging.write_bytes(b"foreign\n")
                displaced.unlink()
                return completed(argv, "")

            run.side_effect = invoke
            unavailable = OSError(errno.ENOSYS, "descriptor link unavailable")
            with patch(
                "gravlax.client._link_descriptor_exact", side_effect=unavailable
            ):
                with self.assertRaisesRegex(OSError, "staging path changed"):
                    Client().run_to_file(["query", "--tsv"], output)

            self.assertFalse(output.exists())
            replacements = list(parent.glob(".rows.tsv.*.gravlax-tmp"))
            self.assertEqual(len(replacements), 1)
            self.assertEqual(replacements[0].read_bytes(), b"foreign\n")

    @patch("gravlax.client.subprocess.run")
    def test_run_to_file_detects_destination_swap_without_deleting_it(self, run):
        def invoke(argv, **kwargs):
            kwargs["stdout"].write(b"trusted\n")
            return completed(argv, "")

        run.side_effect = invoke
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "rows.tsv"
            real_link = client_module._link_held_descriptor

            def link_then_swap(descriptor, staging, destination):
                real_link(descriptor, staging, destination)
                if destination == output:
                    destination.unlink()
                    destination.write_bytes(b"foreign\n")

            with patch(
                "gravlax.client._link_held_descriptor",
                side_effect=link_then_swap,
            ):
                with self.assertRaisesRegex(OSError, "completed output file"):
                    Client().run_to_file(["query", "--tsv"], output)

            self.assertEqual(output.read_bytes(), b"foreign\n")

    @patch("gravlax.client.subprocess.run")
    def test_run_to_file_failure_does_not_unlink_swapped_staging_name(self, run):
        with tempfile.TemporaryDirectory() as directory:
            parent = Path(directory)
            output = parent / "rows.tsv"

            def invoke(argv, **kwargs):
                kwargs["stdout"].write(b"partial")
                kwargs["stdout"].flush()
                staging = next(parent.glob(".rows.tsv.*.gravlax-tmp"))
                displaced = parent / "displaced-original"
                os.replace(staging, displaced)
                staging.write_bytes(b"do not delete")
                displaced.unlink()
                return completed(argv, "", returncode=1, stderr="failed")

            run.side_effect = invoke
            with self.assertRaises(CommandError):
                Client().run_to_file(["query", "--tsv"], output)

            self.assertFalse(output.exists())
            replacements = list(parent.glob(".rows.tsv.*.gravlax-tmp"))
            self.assertEqual(len(replacements), 1)
            self.assertEqual(replacements[0].read_bytes(), b"do not delete")

    def test_resolved_plan_v6_parses_explicit_uniform_io(self):
        document = {
            "schema_version": 6,
            "plan_schema_version": 1,
            "producer": {
                "name": "aie",
                "version": "0.1.0",
                "plan_engine": "aie-declarative-plan-v2",
                "executable_identity": {
                    "scheme": "full-file-blake3-v1",
                    "digest": "executable",
                },
            },
            "name": "uniform",
            "source_path": "/work/plans/uniform.yaml",
            "source_digest": "source",
            "project_name": "project",
            "project_root": "/work",
            "project_manifest": "/work/aie-project.yaml",
            "project_manifest_digest": "manifest",
            "resources": {},
            "steps": [
                {
                    "id": "query",
                    "kind": "query-region",
                    "args": ["query", "sample.aie", "region", "chr1:0-1", "--format", "json"],
                    "outputs": [],
                    "uniform_io": {
                        "kind": "result",
                        "format": "json",
                        "publication": "stdout",
                    },
                    "output_schema_ids": [
                        "gravlax.query.region.result.v1",
                        "gravlax.query.region.counts.v1",
                    ],
                    "io_estimate": {
                        "known_selected_input_bytes": 0,
                        "known_selected_input_files": 0,
                        "unknown_prior_step_outputs": 0,
                        "read_bytes_lower_bound": 0,
                        "read_bytes_upper_bound": None,
                        "bound": "known-inputs-only",
                        "note": "no exact upper bound",
                    },
                    "explanation": ["uniform result -> stdout"],
                }
            ],
        }
        plan = ResolvedPlan.from_mapping(document)
        self.assertEqual(plan.schema_version, 6)
        self.assertEqual(plan.steps[0].uniform_io.kind, "result")
        self.assertEqual(plan.steps[0].uniform_io.format, "json")
        self.assertIsNone(plan.steps[0].uniform_io.output)

        invalid = json.loads(json.dumps(document))
        invalid["steps"][0]["uniform_io"]["output"] = "/work/results/query.json"
        with self.assertRaisesRegex(ProtocolError, "absent exactly when publication is stdout"):
            ResolvedPlan.from_mapping(invalid)


if __name__ == "__main__":
    unittest.main()
