from __future__ import annotations

from contextlib import redirect_stdout
import copy
import hashlib
from io import BytesIO, StringIO
import json
import os
from pathlib import Path
import stat
import subprocess
import sys
import tarfile
import tempfile
import unittest
from unittest.mock import patch
import zipfile


ROOT = Path(__file__).resolve().parents[2]
PACKAGING = ROOT / "packaging"
sys.path.insert(0, str(PACKAGING))

import build_demo_capsule as builder  # noqa: E402
import finalize_demo_capsule as finalizer  # noqa: E402
import verify_demo_capsule as verifier  # noqa: E402


VERSION = "9.9.9"
ARCHIVE_ROOT = f"aie-directory-root-v2:{'a' * 64}"
ARCHIVE_ROOT_B = f"aie-directory-root-v2:{'b' * 64}"
COLLECTION_ROOT = f"aicollection-directory-root-v1:{'c' * 64}"
GITHUB_RELEASES = "https://github.com/COMBINE-lab/gravlax/releases/download"
DATA_BASE_URL = f"{GITHUB_RELEASES}/demo-data-v1"
GROUP_MAP_SCOPE = {
    "label_derivation": "Cell-group labels were frozen without using FNBP1 counts.",
    "included_cells": (
        "Maps include only cells present in the prior FNBP1-targeted archive; donor A includes "
        "867 of 4,275 fully labeled cells."
    ),
    "event_semantics": (
        "The maps cover the FNBP1-support numerators but are not full cell-group denominators."
    ),
    "drilldown_semantics": (
        "The chr12 mechanics drilldown is evaluated only within this labeled subset."
    ),
}


FAKE_AIE = r'''#!/usr/bin/env python3
import json
from pathlib import Path
import sys


def table(name, fields, rows):
    return {
        "name": name,
        "schema": {"fields": [{"name": field} for field in fields]},
        "rows": rows,
    }


args = sys.argv[1:]
if args == ["--version"]:
    print("aie 9.9.9")
elif args and args[0] == "inspect-archive":
    print(json.dumps({
        "verification": {"directory_and_root": True, "all_payloads": True},
        "native_identity": {
            "scheme": "aie-directory-root-v2",
            "blake3": Path(args[1]).stem[-1] * 64,
        },
        "molecular_evidence": {
            "schema": "gravlax.molecular-evidence.v2",
            "alignment_provenance_status": "available",
            "alignment_provenance": {
                "schema": "gravlax.alignment-provenance.v1",
                "molecular_evidence_schema": "gravlax.molecular-evidence.v2",
                "alignment": {
                    "junction_discovery": "per-library-two-pass",
                    "programs": [],
                    "junction_catalogue": {
                        "role": "per-library-pass1",
                        "section": "alignment.junction-catalogue",
                    },
                    "alignment_log": {
                        "locator": "STAR-Log.out",
                        "identity": {"blake3": "b" * 64},
                    },
                },
            },
            "genome_reference_binding_status": "available",
            "genome_reference_binding": {"bound_by": "ingest-archive"},
            "terminal_tail_status": "available",
            "terminal_tail": {"events": 1},
        },
    }))
elif args[:2] == ["collection", "build"]:
    output = next(value.split("=", 1)[1] for value in args if value.startswith("--out="))
    samples = [
        value.split("=", 1)[1].split("=", 1)[0]
        for value in args
        if value.startswith("--sample=")
    ]
    Path(output).write_text(
        json.dumps({"samples": samples, "shape_routes": "--shape-routes" in args}),
        encoding="utf-8",
    )
    print(json.dumps({"data": {"tables": [], "summary": {"samples": len(samples)}}}))
elif args[:2] == ["collection", "inspect"]:
    payload = json.loads(Path(args[2]).read_text(encoding="utf-8"))
    locations = next(
        (value.split("=", 1)[1] for value in args if value.startswith("--locations=")), None
    )
    if locations is None:
        raise SystemExit("the demo collection must be opened through --locations")
    document = json.loads(Path(locations).read_text(encoding="utf-8"))
    if document.get("schema_version") != 1:
        raise SystemExit("unsupported location manifest")
    resolved = {
        entry["identity"]: Path(locations).parent / entry["path"]
        for entry in document["locations"]
    }
    if not resolved or not all(path.is_file() for path in resolved.values()):
        raise SystemExit("location manifest does not resolve every committed source")
    samples = payload["samples"]
    routed = payload.get("shape_routes", True)
    print(json.dumps({
        "layers": [{"root_digest": "c" * 64}],
        "archives": [
            {
                "id": sample,
                "native_identity": {
                    "scheme": "aie-directory-root-v2",
                    "blake3": sample[-1] * 64,
                },
            }
            for sample in sorted(samples)
        ],
        "index": {"shape_route_archives": len(samples) if routed else 0},
        "guard": {
            "content_identity_verified": True,
            "shape_route_payloads_verified": True,
            "shape_route_reconstruction_verified": (
                "--verify-routes" in args or "--verify-content" in args
            ),
        },
    }))
elif args and args[0] == "compare-annotations":
    print(json.dumps({"data": {
        "count_deltas": {
            "schema": {"fields": [
                {"name": "comparison_gene_id"},
                {"name": "signed_delta_b_minus_a"},
            ]},
            "rows": [["GENE1", 1]],
        },
        "summary": {"changed_features": 1},
    }}))
elif args[:2] == ["collection", "find-events"]:
    annotation = args[args.index("--annotation") + 1]
    comparison = annotation.endswith("annotation-after.aic")
    names = [
        "capabilities", "entities", "components", "counts",
        "terminal_anchors", "terminal_counts",
    ]
    tables = [table(name, [], []) for name in names]
    tables[1] = table(
        "entities",
        [
            "entity_id", "kind", "exact_umi_classes", "exact_donors",
            "gap_primary_class", "annotation_incompatible", "compatible_transcripts",
        ],
        [[
            "tail-1", "terminal_tail", 2, 1,
            (None if comparison else "missing_junction"), not comparison,
            (10 if comparison else 0),
        ]],
    )
    print(json.dumps({"data": {"tables": tables, "summary": {"entities": 1}}}))
elif args[:2] == ["collection", "junction"]:
    if args[3].endswith((":+", ":-")):
        raise SystemExit("collection junction locus must not have a strand suffix")
    print(json.dumps({"data": {
        "tables": [table("samples", ["sample", "present", "umis"], [["donor-a", True, 3]])],
        "summary": {"supporting_samples": 1},
    }}))
elif args and args[0] == "query" and "cooccur" in args:
    print(json.dumps({"data": {
        "tables": [
            table("predicates", ["name"], [["splice"], ["tail"]]),
            table(
                "patterns",
                ["matched_predicates", "selection_state", "evidence_units"],
                [[ ["splice", "tail"], "true", 1]],
            ),
            table(
                "memberships",
                [
                    "cell_id", "unit_id", "barcode", "umi_class", "chunk",
                    "local_record", "global_record", "contributing_records",
                    "pattern_mask", "completeness_mask", "matched_predicates",
                    "selection_state", "selected",
                ],
                [[
                    1, "record:1", "AAAAAAAAAAAAAAAA", 1, 0, 1, 1, 1,
                    "0x3", "0x3", ["splice", "tail"], "true", True,
                ]],
            ),
        ],
        "summary": {"selected_units": 1},
    }}))
else:
    print("unsupported fake aie invocation: " + repr(args), file=sys.stderr)
    raise SystemExit(2)
'''


FAKE_BUILDER_AIE = r'''#!/usr/bin/env python3
import json
from pathlib import Path
import sys


args = sys.argv[1:]
if args == ["--version"]:
    print("aie 9.9.9")
elif args and args[0] == "ingest-archive":
    names = {
        "bam": args[1],
        "whitelist": args[args.index("--whitelist") + 1],
        "genome": args[args.index("--genome") + 1],
        "junction_catalogue": args[args.index("--junction-catalogue") + 1],
        "alignment_log": args[args.index("--alignment-log") + 1],
    }
    if any(Path(value).is_absolute() or not Path(value).is_file() for value in names.values()):
        raise SystemExit("ingest inputs must be existing relative paths")
    Path(args[args.index("--out") + 1]).write_text(
        json.dumps({"locators": names}, sort_keys=True), encoding="utf-8"
    )
elif args and args[0] == "inspect-archive":
    payload = json.loads(Path(args[1]).read_text(encoding="utf-8"))
    print(json.dumps({
        "verification": {"directory_and_root": True, "all_payloads": True},
        "native_identity": {
            "scheme": "aie-directory-root-v2",
            "blake3": Path(args[1]).stem[-1] * 64,
        },
        "molecular_evidence": {
            "schema": "gravlax.molecular-evidence.v2",
            "alignment_provenance_status": "available",
            "alignment_provenance": {
                "schema": "gravlax.alignment-provenance.v1",
                "molecular_evidence_schema": "gravlax.molecular-evidence.v2",
                "alignment": {
                    "junction_discovery": "per-library-two-pass",
                    "programs": [],
                    "junction_catalogue": {
                        "role": "per-library-pass1",
                        "section": "alignment.junction-catalogue",
                    },
                    "alignment_log": {
                        "locator": payload["locators"]["alignment_log"],
                        "identity": {"blake3": "b" * 64},
                    },
                },
            },
            "genome_reference_binding_status": "available",
            "genome_reference_binding": {"bound_by": "ingest-archive"},
            "terminal_tail_status": "available",
            "terminal_tail": {"events": 1},
        },
    }))
elif args[:2] == ["collection", "build"]:
    samples = [
        value.split("=", 1)[1].split("=", 1)[0]
        for value in args
        if value.startswith("--sample=")
    ]
    paths = [value.split("=", 1)[1].split("=", 1)[1] for value in args if value.startswith("--sample=")]
    output = next(value.split("=", 1)[1] for value in args if value.startswith("--out="))
    if Path(output).is_absolute() or any(Path(value).is_absolute() for value in paths):
        raise SystemExit("collection build must use public relative basenames")
    if not all(Path(value).is_file() for value in paths):
        raise SystemExit("collection build sources must exist")
    Path(output).write_text(
        json.dumps({"samples": samples, "shape_routes": "--shape-routes" in args}),
        encoding="utf-8",
    )
    print(json.dumps({"data": {"tables": [], "summary": {"samples": len(samples)}}}))
elif args[:2] == ["collection", "inspect"]:
    payload = json.loads(Path(args[2]).read_text(encoding="utf-8"))
    locations = next(
        (value.split("=", 1)[1] for value in args if value.startswith("--locations=")), None
    )
    if locations is None:
        raise SystemExit("the demo collection must be opened through --locations")
    document = json.loads(Path(locations).read_text(encoding="utf-8"))
    if document.get("schema_version") != 1:
        raise SystemExit("unsupported location manifest")
    resolved = {
        entry["identity"]: Path(locations).parent / entry["path"]
        for entry in document["locations"]
    }
    if not resolved or not all(path.is_file() for path in resolved.values()):
        raise SystemExit("location manifest does not resolve every committed source")
    samples = payload["samples"]
    routed = payload.get("shape_routes", True)
    print(json.dumps({
        "layers": [{"root_digest": "c" * 64}],
        "archives": [
            {
                "id": sample,
                "native_identity": {
                    "scheme": "aie-directory-root-v2",
                    "blake3": sample[-1] * 64,
                },
            }
            for sample in sorted(samples)
        ],
        "index": {"shape_route_archives": len(samples) if routed else 0},
        "guard": {
            "content_identity_verified": True,
            "shape_route_payloads_verified": True,
            "shape_route_reconstruction_verified": (
                "--verify-routes" in args or "--verify-content" in args
            ),
        },
    }))
elif args and args[0] == "compile-annotation":
    if Path(args[1]).is_absolute() or Path(args[args.index("--out") + 1]).is_absolute():
        raise SystemExit("annotation paths must be relative")
    Path(args[args.index("--out") + 1]).write_bytes(b"compiled annotation\n")
else:
    raise SystemExit("unsupported builder aie invocation: " + repr(args))
'''


FAKE_SAMTOOLS = r'''#!/usr/bin/env python3
from pathlib import Path
import sys


args = sys.argv[1:]
if not Path(__file__).with_name("samtools-runtime.txt").is_file():
    raise SystemExit("samtools was relocated without its runtime dependency")
if args == ["--version"]:
    print("samtools 1.99")
elif args[:3] == ["view", "--no-PG", "-b"]:
    Path(args[args.index("-o") + 1]).write_bytes(b"filtered bam\n")
elif args[:3] == ["view", "--no-PG", "-H"]:
    if Path(args[3]).name != "sorted.bam":
        raise SystemExit("public header must be derived from sorted BAM")
    print("@HD\tVN:1.6\tSO:coordinate")
    print("@SQ\tSN:chr1\tLN:1000")
    print("@PG\tID:private\tCL:tool /scratch/private")
elif args and args[0] == "sort":
    Path(args[args.index("-o") + 1]).write_bytes(b"coordinate-sorted bam\n")
elif len(args) >= 4 and args[:2] == ["reheader", "--no-PG"]:
    if Path(args[3]).name != "sorted.bam" or not Path(args[3]).is_file():
        raise SystemExit("reheader input must be the sorted BAM")
    sys.stdout.buffer.write(b"sanitized bam\n")
elif args[:2] == ["quickcheck", "-v"]:
    pass
elif args[:2] == ["view", "-c"]:
    print("3")
else:
    raise SystemExit("unsupported fake samtools invocation: " + repr(args))
'''


class BuilderUnitTests(unittest.TestCase):
    def _source_declaration(self, root: Path, filename: str) -> dict[str, object]:
        path = root / filename
        return {
            "path": filename,
            "bytes": path.stat().st_size,
            "sha256": builder._sha256(path),
            "provenance": {
                "accession": f"fixture-{filename}",
                "source_url": f"https://example.org/source/{filename}",
            },
        }

    def test_builder_uses_relative_public_ingest_locators(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sources = root / "sources"
            sources.mkdir()
            payloads = {
                "reference.fa": b">chr1\nAAAAAAAAAAAAAAAAAAAA\n",
                "donor.bam": b"fixture bam\n",
                "whitelist.txt": b"AAAAAAAAAAAAAAAA\n",
                "groups.tsv": b"AAAAAAAAAAAAAAAA\timmune\n",
                "SJ.out.tab": b"chr1\t121\t180\t1\t1\t1\t1\t0\t20\n",
                "Log.out": (
                    b"STAR version=2.7.11b\n"
                    b"##### Final effective command line:\n"
                    b"STAR --twopassMode Basic\n"
                    b"twopassMode Basic\n"
                    b"sjdbFileChrStartEnd -\n"
                    b"sjdbGTFfile -\n"
                ),
                "before.gtf": (
                    b'chr1\ttest\tgene\t101\t200\t.\t+\t.\tgene_id "G1";\n'
                    b'chr1\ttest\texon\t101\t120\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n'
                ),
                "after.gtf": (
                    b'chr1\ttest\tgene\t101\t200\t.\t+\t.\tgene_id "G1";\n'
                    b'chr1\ttest\texon\t101\t180\t.\t+\t.\tgene_id "G1"; transcript_id "T2";\n'
                ),
            }
            for filename, payload in payloads.items():
                (sources / filename).write_bytes(payload)

            aie = root / "aie"
            aie.write_text(FAKE_BUILDER_AIE, encoding="utf-8")
            aie.chmod(aie.stat().st_mode | stat.S_IXUSR)
            samtools = root / "samtools"
            samtools.write_text(FAKE_SAMTOOLS, encoding="utf-8")
            samtools.chmod(samtools.stat().st_mode | stat.S_IXUSR)
            (root / "samtools-runtime.txt").write_text("required\n", encoding="utf-8")
            archive_map = {"donor-a": "archive_a"}
            spec = {
                "schema": builder.BUILD_SCHEMA,
                "capsule_id": "fixture",
                "summary": "A curated public fixture donor over one genomic window.",
                "software": {"version": VERSION, "revision": "b" * 40},
                "assembly": "GRCh38",
                "group_map_scope": GROUP_MAP_SCOPE,
                "windows": [{"chrom": "chr1", "start": 100, "end": 200}],
                "alignment": {
                    "program": "STAR",
                    "version": "2.7.11b",
                    "command_line": "STAR 2.7.11b --twopassMode Basic",
                    "junction_discovery": "per-library-two-pass",
                    "index_identity": "sha256:public-fixture-index",
                    "chemistry": "10x-3p-v3",
                },
                "genome": self._source_declaration(sources, "reference.fa"),
                "donors": [
                    {
                        "sample": "donor-a",
                        "donor": "person-a",
                        "archive_resource": "archive_a",
                        "archive_filename": "archive-a.aie",
                        "accessions": ["SRR123456"],
                        "bam": self._source_declaration(sources, "donor.bam"),
                        "whitelist": self._source_declaration(sources, "whitelist.txt"),
                        "groups": self._source_declaration(sources, "groups.tsv"),
                        "junction_catalogue": self._source_declaration(sources, "SJ.out.tab"),
                        "alignment_log": self._source_declaration(sources, "Log.out"),
                    }
                ],
                "annotations": {
                    "before": {
                        **self._source_declaration(sources, "before.gtf"),
                        "resource": "annotation_before",
                        "filename": "annotation-before.aic",
                        "source_filename": "annotation-before.gtf",
                        "label": "before",
                    },
                    "after": {
                        **self._source_declaration(sources, "after.gtf"),
                        "resource": "annotation_after",
                        "filename": "annotation-after.aic",
                        "source_filename": "annotation-after.gtf",
                        "label": "after",
                    },
                },
                "stories": {
                    "annotation_reinterpretation": {
                        "archive": "archive_a",
                        "annotation_before": "annotation_before",
                        "annotation_after": "annotation_after",
                        "assembly": "GRCh38",
                        "before_label": "before",
                        "after_label": "after",
                        "expected_gene_id": "GENE1",
                        "expected_min_signed_delta": 1,
                    },
                    "event_discovery": {
                        "archives": archive_map,
                        "collection_groups": "collection_groups",
                        "design": "design",
                        "annotation": "annotation_before",
                        "assembly": "GRCh38",
                        "annotation_label": "before",
                        "comparison_annotation": "annotation_after",
                        "comparison_annotation_label": "after",
                        "expected_entity_id": "tail-1",
                        "expected_min_exact_umi_classes": 1,
                        "expected_min_exact_donors": 1,
                        "expected_gap_primary_class": "missing_junction",
                        "expected_annotation_incompatible": True,
                        "expected_rank": 1,
                        "expected_comparison_compatible_transcripts": 10,
                    },
                    "junction_drilldown": {
                        "archives": archive_map,
                        "archive_sample": "donor-a",
                        "drilldown_groups": "drilldown_groups",
                        "junction": "chr1:120-180",
                        "predicates": {"splice": "junction:chr1:120-180:+"},
                        "expression": "splice",
                        "universe": "splice",
                        "emit_membership": True,
                        "expected_min_selected_units": 1,
                        "required_true_predicates": ["splice"],
                        "story_note": "A technical positive control for the query mechanics.",
                    },
                },
                "third_party_notices": [
                    {
                        "name": "Public fixture",
                        "source_url": "https://example.org/source",
                        "terms": "Fixture data terms",
                        "terms_url": "https://example.org/terms",
                        "citation": "Fixture citation",
                    }
                ],
            }
            spec_path = root / "spec.json"
            spec_path.write_text(json.dumps(spec), encoding="utf-8")
            output = root / "capsule-build"
            with redirect_stdout(StringIO()):
                builder.build(spec_path, sources, aie, samtools, output)

            archived = json.loads((output / "archive-a.aie").read_text(encoding="utf-8"))
            self.assertEqual(
                archived["locators"],
                {
                    "alignment_log": "STAR-Log.out",
                    "bam": "sanitized.bam",
                    "genome": "reference.fa",
                    "junction_catalogue": "STAR-pass1-SJ.out.tab",
                    "whitelist": "barcodes.txt",
                },
            )
            record_text = (output / "BUILD-RECORD.json").read_text(encoding="utf-8")
            self.assertNotIn(str(root), record_text)
            record = json.loads(record_text)
            self.assertEqual(record["group_map_scope"], GROUP_MAP_SCOPE)
            self.assertEqual(
                record["inputs"]["donors"][0]["alignment_log_validation"][
                    "alignment_annotation_status"
                ],
                "absent-in-exact-Log.out",
            )
            self.assertIn("caller-declared", record["alignment"]["index_manifest_status"])
            self.assertIn("multimapper_scope", record["selection"])
            self.assertIn("collection_groups", record["resources"])
            self.assertIn("drilldown_groups", record["resources"])
            self.assertEqual(
                (output / "drilldown-groups.tsv").read_text(encoding="utf-8"),
                "AAAAAAAAAAAAAAAA\timmune\n",
            )

            no_memberships = copy.deepcopy(spec)
            no_memberships["stories"]["junction_drilldown"]["emit_membership"] = False
            spec_path.write_text(json.dumps(no_memberships), encoding="utf-8")
            with self.assertRaisesRegex(builder.CapsuleError, "emit_membership must be true"):
                builder.build(
                    spec_path, sources, aie, samtools, root / "invalid-membership-build"
                )
            self.assertFalse((root / "invalid-membership-build").exists())

            malformed = copy.deepcopy(spec)
            del malformed["stories"]["event_discovery"]["comparison_annotation_label"]
            spec_path.write_text(json.dumps(malformed), encoding="utf-8")
            marker = root / "tool-was-invoked"
            sentinel = (
                "#!/usr/bin/env python3\n"
                "from pathlib import Path\n"
                f"Path({str(marker)!r}).write_text('called', encoding='utf-8')\n"
            )
            for executable in (aie, samtools):
                executable.write_text(sentinel, encoding="utf-8")
                executable.chmod(executable.stat().st_mode | stat.S_IXUSR)
            with self.assertRaisesRegex(builder.CapsuleError, "missing fields"):
                builder.build(spec_path, sources, aie, samtools, root / "invalid-build")
            self.assertFalse(marker.exists(), "malformed stories must fail before tools run")
            self.assertFalse((root / "invalid-build").exists())

    def test_star_log_must_confirm_no_runtime_annotation_inputs(self) -> None:
        alignment = {
            "program": "STAR",
            "version": "2.7.11b",
            "command_line": "STAR 2.7.11b --twopassMode Basic",
        }
        with tempfile.TemporaryDirectory() as temporary:
            log = Path(temporary) / "Log.out"
            base = (
                "STAR version=2.7.11b\n"
                "##### Final effective command line:\n"
                "twopassMode Basic\n"
                "sjdbFileChrStartEnd -\n"
            )
            log.write_text(base, encoding="utf-8")
            with self.assertRaisesRegex(builder.CapsuleError, "sjdbGTFfile -"):
                builder._validate_star_log(log, alignment, "fixture.alignment_log")

            log.write_text(base + "sjdbGTFfile annotations.gtf\n", encoding="utf-8")
            with self.assertRaisesRegex(builder.CapsuleError, "sjdbGTFfile -"):
                builder._validate_star_log(log, alignment, "fixture.alignment_log")

            log.write_text(base + "sjdbGTFfile -\n", encoding="utf-8")
            validation = builder._validate_star_log(log, alignment, "fixture.alignment_log")
            self.assertEqual(
                validation["alignment_annotation_status"],
                "absent-in-exact-Log.out",
            )

    def test_staging_never_shares_an_external_source_inode(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "source.bin"
            source.write_bytes(b"immutable source bytes\n")
            digest = builder._sha256(source)
            before = source.stat()
            first = root / "build-one.bin"
            second = root / "build-two.bin"

            builder._stage_input(source, first, digest)
            builder._stage_input(source, second, digest)

            after = source.stat()
            self.assertEqual(after.st_ino, before.st_ino)
            self.assertEqual(after.st_nlink, before.st_nlink)
            self.assertEqual(after.st_ctime_ns, before.st_ctime_ns)
            self.assertNotEqual(first.stat().st_ino, source.stat().st_ino)
            self.assertNotEqual(second.stat().st_ino, source.stat().st_ino)
            self.assertNotEqual(first.stat().st_ino, second.stat().st_ino)
            original = first.read_bytes()
            source.write_bytes(b"x" * len(original))
            self.assertEqual(source.stat().st_size, len(original))
            self.assertNotEqual(builder._sha256(source), digest)
            self.assertEqual(first.read_bytes(), original)
            self.assertEqual(second.read_bytes(), original)
            declaration = {
                "bytes": len(original),
                "sha256": digest,
                "provenance": {"accession": "fixture", "source_url": "https://example.org"},
            }
            source.write_bytes(b"larger mutable source bytes\n")
            self.assertEqual(builder._public_input_record(declaration)["bytes"], len(original))

    def test_windows_are_sorted_and_overlaps_are_rejected(self) -> None:
        parsed = builder._windows(
            {
                "windows": [
                    {"chrom": "chr2", "start": 20, "end": 30},
                    {"chrom": "chr1", "start": 10, "end": 20},
                    {"chrom": "chr1", "start": 20, "end": 25},
                ]
            }
        )
        self.assertEqual(parsed, [("chr1", 10, 20), ("chr1", 20, 25), ("chr2", 20, 30)])
        with self.assertRaisesRegex(builder.CapsuleError, "must not overlap"):
            builder._windows(
                {
                    "windows": [
                        {"chrom": "chr1", "start": 10, "end": 21},
                        {"chrom": "chr1", "start": 20, "end": 30},
                    ]
                }
            )

    def test_annotation_subset_retains_complete_overlapping_genes(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            source = directory / "source.gtf"
            destination = directory / "subset.gtf"
            source.write_text(
                "##description: fixture\n"
                'chr1\ttest\tgene\t101\t150\t.\t+\t.\tgene_id "G1";\n'
                'chr1\ttest\texon\t501\t550\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n'
                'chr1\ttest\tgene\t701\t750\t.\t+\t.\tgene_id "G2";\n',
                encoding="utf-8",
            )
            selected = builder._subset_gtf(source, [("chr1", 100, 151)], destination)
            text = destination.read_text(encoding="utf-8")
            self.assertEqual(selected, 1)
            self.assertIn("complete records for genes overlapping", text)
            self.assertIn('gene_id "G1"', text)
            self.assertIn("501\t550", text)
            self.assertNotIn('gene_id "G2"', text)

    def test_sanitized_header_keeps_dictionary_and_removes_private_history(self) -> None:
        original = (
            "@HD\tVN:1.6\tSO:coordinate\n"
            "@SQ\tSN:chr1\tLN:1000\n"
            "@RG\tID:old\n"
            "@PG\tID:STAR\tCL:STAR --genomeDir /scratch/private\n"
            "@CO\t/nfshomes/private/input.bam\n"
        )
        header = builder._sanitized_header(
            original,
            "donor-a",
            ["SRR123456"],
            [("chr1", 100, 200)],
        ).decode()
        self.assertEqual(sum(line.startswith("@HD\t") for line in header.splitlines()), 1)
        self.assertIn("SO:coordinate", header.splitlines()[0])
        self.assertEqual(sum(line.startswith("@SQ\t") for line in header.splitlines()), 1)
        self.assertNotIn("@PG\t", header)
        self.assertIn("accessions=SRR123456", header)
        self.assertNotIn("@RG\t", header)
        self.assertNotIn("/scratch", header)
        self.assertNotIn("/nfshomes", header)
        with self.assertRaisesRegex(builder.CapsuleError, "private filesystem locator"):
            builder._public_text("aligner --input /mnt/private/sample.bam", "command")
        with self.assertRaisesRegex(builder.CapsuleError, "SO:coordinate"):
            builder._sanitized_header(
                original.replace("\tSO:coordinate", ""),
                "donor-a",
                ["SRR123456"],
                [("chr1", 100, 200)],
            )

    def test_group_map_is_headerless_unique_and_deterministic(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            groups = Path(temporary) / "groups.tsv"
            groups.write_text(
                "TTTTTTTTTTTTTTTT\tT-cell\nAAAAAAAAAAAAAAAA\tB-cell\n",
                encoding="utf-8",
            )
            self.assertEqual(
                builder._read_group_map(groups, "groups"),
                [("AAAAAAAAAAAAAAAA", "B-cell"), ("TTTTTTTTTTTTTTTT", "T-cell")],
            )
            groups.write_text(
                "AAAAAAAAAAAAAAAA\tB-cell\nAAAAAAAAAAAAAAAA\tT-cell\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(builder.CapsuleError, "repeats barcode"):
                builder._read_group_map(groups, "groups")

    def _collection_output(self, root: Path) -> tuple[Path, Path, dict[str, dict[str, object]]]:
        aie = root / "aie"
        aie.write_text(FAKE_BUILDER_AIE, encoding="utf-8")
        aie.chmod(aie.stat().st_mode | stat.S_IXUSR)
        output = root / "capsule"
        output.mkdir()
        resources: dict[str, dict[str, object]] = {}
        for letter in ("a", "b"):
            path = output / f"archive-{letter}.aie"
            path.write_bytes(f"archive {letter}\n".encode())
            resources[f"archive_{letter}"] = builder._asset(
                path, f"aie-directory-root-v2:{letter * 64}"
            )
        return aie, output, resources

    def test_staged_collection_publishes_an_identity_keyed_location_manifest(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            aie, output, resources = self._collection_output(root)
            staging_root = root / "public-staging"
            staging_root.mkdir()
            with patch.object(builder, "PRIVATE_TEXT", ("/definitely-not-a-real-root",)):
                record = builder._stage_collection(
                    aie=aie,
                    output=output,
                    declaration={
                        "filename": "demo.aicollection",
                        "locations_filename": "demo-locations.json",
                    },
                    resources=resources,
                    archive_map={"donor-a": "archive_a", "donor-b": "archive_b"},
                    output_filenames={"BUILD-RECORD.json"},
                    staging_root=staging_root,
                )
            self.assertEqual(record["collection_root"], COLLECTION_ROOT)
            self.assertTrue(record["shape_routes"])
            self.assertFalse(record["allow_unstamped"])
            self.assertEqual(
                record["archives"], {"donor-a": "archive_a", "donor-b": "archive_b"}
            )
            collection = output / "demo.aicollection"
            locations = output / "demo-locations.json"
            self.assertTrue(collection.is_file() and locations.is_file())
            self.assertEqual(record["sha256"], builder._sha256(collection))
            self.assertEqual(record["locations"]["sha256"], builder._sha256(locations))
            document = json.loads(locations.read_text(encoding="utf-8"))
            self.assertEqual(
                document,
                {
                    "schema_version": 1,
                    "locations": [
                        {
                            "identity": f"aie-directory-root-v2:{'a' * 64}",
                            "path": "archive-a.aie",
                        },
                        {
                            "identity": f"aie-directory-root-v2:{'b' * 64}",
                            "path": "archive-b.aie",
                        },
                    ],
                },
            )
            self.assertFalse(list(staging_root.iterdir()))

    def test_staged_collection_refuses_a_private_staging_root(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            aie, output, resources = self._collection_output(root)
            staging_root = root / "private-staging"
            staging_root.mkdir()
            with patch.object(builder, "PRIVATE_TEXT", ("/private-staging",)):
                with self.assertRaisesRegex(
                    builder.CapsuleError, "public-neutral directory"
                ):
                    builder._stage_collection(
                        aie=aie,
                        output=output,
                        declaration={
                            "filename": "demo.aicollection",
                            "locations_filename": "demo-locations.json",
                        },
                        resources=resources,
                        archive_map={"donor-a": "archive_a"},
                        output_filenames=set(),
                        staging_root=staging_root,
                    )
            self.assertFalse((output / "demo.aicollection").exists())

    def test_collection_junction_rejects_a_strand_suffix(self) -> None:
        self.assertEqual(
            builder._collection_junction("chr12:18715181-18853140", "junction"),
            "chr12:18715181-18853140",
        )
        with self.assertRaisesRegex(builder.CapsuleError, "without a strand suffix"):
            builder._collection_junction("chr12:18715181-18853140:+", "junction")


class FinalizerVerifierTests(unittest.TestCase):
    def test_story_validator_matches_the_checked_in_manifest_schema(self) -> None:
        schema = json.loads(
            (ROOT / "notebooks/demo-manifest.schema.json").read_text(encoding="utf-8")
        )
        definitions = schema["$defs"]
        self.assertEqual(
            set(definitions["annotation_story"]["properties"]),
            builder.ANNOTATION_STORY_FIELDS,
        )
        self.assertEqual(
            set(definitions["event_story"]["properties"]), builder.EVENT_STORY_FIELDS
        )
        self.assertEqual(
            set(definitions["drilldown_story"]["properties"]),
            builder.DRILLDOWN_STORY_FIELDS,
        )
        self.assertIn("collection_groups", definitions["event_story"]["required"])
        self.assertNotIn("drilldown_groups", definitions["event_story"]["properties"])
        self.assertIn("drilldown_groups", definitions["drilldown_story"]["properties"])
        self.assertNotIn("collection_groups", definitions["drilldown_story"]["properties"])
        self.assertEqual(schema["properties"]["collection"]["$ref"], "#/$defs/collection")
        self.assertEqual(
            set(definitions["collection"]["properties"]),
            {
                "archives",
                "shape_routes",
                "allow_unstamped",
                "collection_root",
                "asset",
                "locations",
            },
        )
        self.assertFalse(definitions["collection"]["additionalProperties"])
        self.assertNotIn("collection", schema["required"])

    def _write_fake_aie(self, path: Path) -> bytes:
        payload = FAKE_AIE.encode()
        path.write_bytes(payload)
        path.chmod(path.stat().st_mode | stat.S_IXUSR)
        return payload

    def _write_release_bundle(self, path: Path, payload: bytes) -> str:
        member = "gravlax-x86_64-unknown-linux-gnu/aie"
        with tarfile.open(path, "w:gz") as archive:
            info = tarfile.TarInfo(member)
            info.mode = 0o755
            info.size = len(payload)
            archive.addfile(info, BytesIO(payload))
        return member

    def _write_wheel(self, path: Path) -> None:
        with zipfile.ZipFile(path, "w") as archive:
            archive.writestr(
                "gravlax_client-9.9.9.dist-info/METADATA",
                "Metadata-Version: 2.1\nName: gravlax-client\nVersion: 9.9.9\n",
            )

    def _write_build_directory(self, directory: Path, *, collection: bool = False) -> None:
        directory.mkdir()
        resource_payloads = {
            "archive-a.aie": b"archive a\n",
            "archive-b.aie": b"archive b\n",
            "annotation-before.aic": b"annotation before\n",
            "annotation-after.aic": b"annotation after\n",
            "collection-groups.tsv": (
                b"sample\tbarcode\tgroup\n"
                b"donor-a\tAAAAAAAAAAAAAAAA\timmune\n"
                b"donor-b\tCCCCCCCCCCCCCCCC\timmune\n"
            ),
            "drilldown-groups.tsv": b"AAAAAAAAAAAAAAAA\timmune\n",
            "donors.tsv": b"sample\tdonor\ndonor-a\tperson-a\ndonor-b\tperson-b\n",
        }
        for filename, payload in resource_payloads.items():
            (directory / filename).write_bytes(payload)
        resources = {
            "archive_a": builder._asset(directory / "archive-a.aie", ARCHIVE_ROOT),
            "archive_b": builder._asset(directory / "archive-b.aie", ARCHIVE_ROOT_B),
            "annotation_before": builder._asset(directory / "annotation-before.aic"),
            "annotation_after": builder._asset(directory / "annotation-after.aic"),
            "collection_groups": builder._asset(directory / "collection-groups.tsv"),
            "drilldown_groups": builder._asset(directory / "drilldown-groups.tsv"),
            "design": builder._asset(directory / "donors.tsv"),
        }
        archive_map = {"donor-a": "archive_a", "donor-b": "archive_b"}
        stories = {
            "annotation_reinterpretation": {
                "archive": "archive_a",
                "annotation_before": "annotation_before",
                "annotation_after": "annotation_after",
                "assembly": "GRCh38",
                "before_label": "before",
                "after_label": "after",
                "expected_gene_id": "GENE1",
                "expected_min_signed_delta": 1,
            },
            "event_discovery": {
                "archives": archive_map,
                "collection_groups": "collection_groups",
                "design": "design",
                "annotation": "annotation_before",
                "assembly": "GRCh38",
                "annotation_label": "before",
                "comparison_annotation": "annotation_after",
                "comparison_annotation_label": "after",
                "kinds": ["junction", "terminal-tail"],
                "expected_entity_id": "tail-1",
                "expected_min_exact_umi_classes": 2,
                "expected_min_exact_donors": 1,
                "expected_gap_primary_class": "missing_junction",
                "expected_annotation_incompatible": True,
                "expected_rank": 1,
                "expected_comparison_compatible_transcripts": 10,
            },
            "junction_drilldown": {
                "archives": archive_map,
                "archive_sample": "donor-a",
                "drilldown_groups": "drilldown_groups",
                "junction": "chr12:18715181-18853140",
                "predicates": {
                    "splice": "junction:chr12:18715181-18853140:+",
                    "tail": "terminal:chr12:18853181-18853182:+",
                },
                "expression": "splice & tail",
                "universe": "tail",
                "region_match": "aligned-block",
                "aggregation": "group",
                "emit_membership": True,
                "expected_min_selected_units": 1,
                "required_true_predicates": ["splice", "tail"],
                "story_note": "A technical positive control for the query mechanics.",
            },
        }
        collection_record = None
        if collection:
            collection_payload = json.dumps(
                {"samples": sorted(archive_map), "shape_routes": True}
            ).encode()
            (directory / "demo.aicollection").write_bytes(collection_payload)
            locations = builder._locations_document(
                {
                    ARCHIVE_ROOT: resources["archive_a"]["filename"],
                    ARCHIVE_ROOT_B: resources["archive_b"]["filename"],
                }
            )
            (directory / "demo-locations.json").write_bytes(
                (json.dumps(locations, indent=2, sort_keys=True, ensure_ascii=True) + "\n").encode()
            )
            collection_record = {
                **builder._asset(directory / "demo.aicollection"),
                "collection_root": COLLECTION_ROOT,
                "archives": archive_map,
                "shape_routes": True,
                "allow_unstamped": False,
                "locations": builder._asset(directory / "demo-locations.json"),
            }
        record = {
            "schema": builder.RECORD_SCHEMA,
            "capsule_id": "test-capsule",
            "summary": "A curated demonstration from two independent public fixture donors.",
            "assembly": "GRCh38",
            "group_map_scope": GROUP_MAP_SCOPE,
            "software": {
                "version": VERSION,
                "revision": "0" * 40,
                "aie_sha256": hashlib.sha256(FAKE_AIE.encode()).hexdigest(),
            },
            "selection": {
                "windows": [
                    {"chrom": "chr9", "start": 129907500, "end": 129926000},
                    {"chrom": "chr12", "start": 18714000, "end": 18854000},
                ]
            },
            "inputs": {
                "donors": [
                    {
                        "sample": "donor-a",
                        "donor": "person-a",
                        "accessions": ["SRR000001"],
                    },
                    {
                        "sample": "donor-b",
                        "donor": "person-b",
                        "accessions": ["SRR000002"],
                    },
                ]
            },
            "resources": resources,
            **({"collection": collection_record} if collection_record else {}),
            "stories": stories,
            "third_party_notices": [
                {
                    "name": "Fixture data",
                    "source_url": "https://example.org/source",
                    "terms": "Fixture data terms",
                    "terms_url": "https://example.org/terms",
                    "citation": "Fixture citation",
                }
            ],
        }
        (directory / "BUILD-RECORD.json").write_text(
            json.dumps(record, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        notice_text, _ = builder._notices(record["third_party_notices"])
        (directory / "THIRD-PARTY-NOTICES.md").write_text(notice_text, encoding="utf-8")

    def _source_repository(self, root: Path, build_dir: Path) -> Path:
        repository = root / "source-repository"
        repository.mkdir()
        subprocess.run(["git", "init", "--quiet", str(repository)], check=True)
        (repository / "fixture.txt").write_text("fixture\n", encoding="utf-8")
        subprocess.run(["git", "-C", str(repository), "add", "fixture.txt"], check=True)
        subprocess.run(
            [
                "git",
                "-C",
                str(repository),
                "-c",
                "user.name=Capsule Fixture",
                "-c",
                "user.email=capsule@example.org",
                "commit",
                "--quiet",
                "-m",
                "fixture",
            ],
            check=True,
        )
        revision = subprocess.run(
            ["git", "-C", str(repository), "rev-parse", "HEAD"],
            check=True,
            stdout=subprocess.PIPE,
            text=True,
        ).stdout.strip()
        subprocess.run(
            ["git", "-C", str(repository), "tag", f"v{VERSION}"], check=True
        )
        record_path = build_dir / "BUILD-RECORD.json"
        record = json.loads(record_path.read_text(encoding="utf-8"))
        record["software"]["revision"] = revision
        record_path.write_text(
            json.dumps(record, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        return repository

    def _finalize_fixture(
        self, root: Path, *, collection: bool = False
    ) -> tuple[Path, Path, Path, Path, Path]:
        build_dir = root / "build"
        output_dir = root / "final"
        self._write_build_directory(build_dir, collection=collection)
        source_repository = self._source_repository(root, build_dir)
        aie = root / "aie"
        payload = self._write_fake_aie(aie)
        release = root / "gravlax-x86_64-unknown-linux-gnu.tar.gz"
        member = self._write_release_bundle(release, payload)
        wheel = root / "gravlax_client-9.9.9-py3-none-any.whl"
        self._write_wheel(wheel)
        with redirect_stdout(StringIO()):
            finalizer.finalize(
                build_dir=build_dir,
                output_dir=output_dir,
                data_base_url=DATA_BASE_URL,
                aie_asset=release,
                aie_url=f"{GITHUB_RELEASES}/v{VERSION}/{release.name}",
                aie_member=member,
                python_wheel=wheel,
                python_wheel_url=(
                    f"https://files.pythonhosted.org/packages/aa/bb/fixture/{wheel.name}"
                ),
                source_repository=source_repository,
            )
        return output_dir, aie, release, wheel, source_repository

    def test_finalized_capsule_is_independent_non_overwriting_and_executable(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            output, aie, release, wheel, source_repository = self._finalize_fixture(root)
            build_archive = root / "build/archive-a.aie"
            final_archive = output / "archive-a.aie"
            self.assertNotEqual(os.stat(build_archive).st_ino, os.stat(final_archive).st_ino)
            self.assertTrue((output / "SHA256SUMS").is_file())
            readme = (output / "README.md").read_text(encoding="utf-8")
            self.assertIn("no read sequences, base qualities, or read names", readme)
            self.assertIn("three stories", readme)
            release_notes = (output / "RELEASE-NOTES.md").read_text(encoding="utf-8")
            for phrase in (
                "frozen without using FNBP1 counts",
                "867 of 4,275",
                "not full cell-group denominators",
                "chr12 mechanics drilldown",
            ):
                self.assertIn(phrase, readme)
                self.assertIn(phrase, release_notes)
            expected_checksums = {
                path.name
                for path in output.iterdir()
                if path.is_file() and path.name != "SHA256SUMS"
            }
            covered = {
                line.split("  ", 1)[1]
                for line in (output / "SHA256SUMS").read_text(encoding="ascii").splitlines()
            }
            self.assertEqual(covered, expected_checksums)

            result = verifier.verify(output, aie, release, wheel)
            self.assertEqual(result["archives_verified"], 2)
            self.assertTrue(result["stories_executed"])
            self.assertEqual(result["junction_drilldown"]["cooccurrence"]["selected_units"], 1)

            with self.assertRaisesRegex(FileExistsError, "refusing to replace capsule directory"):
                with redirect_stdout(StringIO()):
                    finalizer.finalize(
                        build_dir=root / "build",
                        output_dir=output,
                        data_base_url=DATA_BASE_URL,
                        aie_asset=release,
                        aie_url=f"{GITHUB_RELEASES}/v{VERSION}/{release.name}",
                        aie_member="gravlax-x86_64-unknown-linux-gnu/aie",
                        python_wheel=wheel,
                        python_wheel_url=(
                            "https://files.pythonhosted.org/packages/aa/bb/fixture/"
                            f"{wheel.name}"
                        ),
                        source_repository=source_repository,
                    )

            final_archive.write_bytes(final_archive.read_bytes() + b"changed\n")
            with self.assertRaisesRegex(builder.CapsuleError, "SHA-256 mismatch"):
                verifier._verify_checksums(output)

    def test_manifest_validation_rejects_mutable_or_path_bound_assets(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output, _, _, _, _ = self._finalize_fixture(Path(temporary))
            manifest = json.loads((output / "demo-manifest.json").read_text(encoding="utf-8"))
            verifier._validate_manifest(manifest)

            mutable = copy.deepcopy(manifest)
            mutable["resources"]["archive_a"]["url"] = (
                "https://github.com/COMBINE-lab/gravlax/releases/latest/archive-a.aie"
            )
            with self.assertRaisesRegex(builder.CapsuleError, "immutable"):
                verifier._validate_manifest(mutable)

            moving_branch = copy.deepcopy(manifest)
            moving_branch["resources"]["archive_a"]["url"] = (
                "https://raw.githubusercontent.com/COMBINE-lab/gravlax/main/archive-a.aie"
            )
            with self.assertRaisesRegex(builder.CapsuleError, "immutable"):
                verifier._validate_manifest(moving_branch)

            non_pypi_wheel = copy.deepcopy(manifest)
            wheel = non_pypi_wheel["software"]["python_wheel"]
            wheel["url"] = f"{GITHUB_RELEASES}/v{VERSION}/{wheel['filename']}"
            with self.assertRaisesRegex(builder.CapsuleError, "pypi-file"):
                verifier._validate_manifest(non_pypi_wheel)

            wrong_software_tag = copy.deepcopy(manifest)
            aie_asset = wrong_software_tag["software"]["aie"]
            aie_asset["url"] = f"{GITHUB_RELEASES}/v0.0.1/{aie_asset['filename']}"
            with self.assertRaisesRegex(builder.CapsuleError, "URL tag"):
                verifier._validate_manifest(wrong_software_tag)

            scattered = copy.deepcopy(manifest)
            archive = scattered["resources"]["archive_a"]
            archive["url"] = f"{GITHUB_RELEASES}/demo-data-v2/{archive['filename']}"
            with self.assertRaisesRegex(builder.CapsuleError, "outside the finalized"):
                verifier._verify_capsule_records(output, scattered)

            path_bound = copy.deepcopy(manifest)
            path_bound["resources"]["archive_a"]["filename"] = "atlas.aicollection"
            path_bound["resources"]["archive_a"]["url"] = (
                f"{DATA_BASE_URL}/atlas.aicollection"
            )
            with self.assertRaisesRegex(builder.CapsuleError, "path-bound"):
                verifier._validate_manifest(path_bound)

            build_record_path = output / "BUILD-RECORD.json"
            build_record = json.loads(build_record_path.read_text(encoding="utf-8"))
            build_record["software"]["revision"] = "not-a-git-revision"
            build_record_path.write_text(
                json.dumps(build_record, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(builder.CapsuleError, "software revision"):
                verifier._verify_capsule_records(output, manifest)

    def test_published_collection_is_resolved_through_its_location_manifest(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            output, aie, release, wheel, _ = self._finalize_fixture(root, collection=True)
            manifest = json.loads((output / "demo-manifest.json").read_text(encoding="utf-8"))
            collection = manifest["collection"]
            self.assertEqual(collection["collection_root"], COLLECTION_ROOT)
            self.assertEqual(collection["archives"], {"donor-a": "archive_a", "donor-b": "archive_b"})
            self.assertTrue(collection["shape_routes"])
            self.assertEqual(
                collection["asset"]["url"], f"{DATA_BASE_URL}/demo.aicollection"
            )
            self.assertEqual(
                collection["locations"]["url"], f"{DATA_BASE_URL}/demo-locations.json"
            )
            self.assertNotIn("demo.aicollection", {
                asset["filename"] for asset in manifest["resources"].values()
            })
            readme = (output / "README.md").read_text(encoding="utf-8")
            self.assertIn("demo-locations.json", readme)
            self.assertNotIn("rebuild each `.aicollection`", readme)
            covered = {
                line.split("  ", 1)[1]
                for line in (output / "SHA256SUMS").read_text(encoding="ascii").splitlines()
            }
            self.assertIn("demo.aicollection", covered)
            self.assertIn("demo-locations.json", covered)

            result = verifier.verify(output, aie, release, wheel)
            self.assertEqual(result["collection_root"], COLLECTION_ROOT)
            self.assertEqual(result["archives_verified"], 2)
            self.assertTrue(result["stories_executed"])

            document = json.loads(
                (output / "demo-locations.json").read_text(encoding="utf-8")
            )
            self.assertEqual(document["schema_version"], 1)
            self.assertEqual(
                {entry["path"] for entry in document["locations"]},
                {"archive-a.aie", "archive-b.aie"},
            )

    def test_published_collection_fails_closed_on_manifest_drift(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            output, aie, release, wheel, _ = self._finalize_fixture(root, collection=True)
            manifest = json.loads((output / "demo-manifest.json").read_text(encoding="utf-8"))
            verifier._validate_manifest(copy.deepcopy(manifest))

            wrong_root = copy.deepcopy(manifest)
            wrong_root["collection"]["collection_root"] = f"aicollection-directory-root-v1:{'d' * 64}"
            with self.assertRaisesRegex(builder.CapsuleError, "root differs from the manifest"):
                verifier._verify_collection(
                    aie,
                    wrong_root,
                    output / "demo.aicollection",
                    output / "demo-locations.json",
                )

            other_archives = copy.deepcopy(manifest)
            other_archives["collection"]["archives"] = {"donor-a": "archive_a"}
            with self.assertRaisesRegex(builder.CapsuleError, "differs from the published collection"):
                verifier._validate_manifest(other_archives)

            other_options = copy.deepcopy(manifest)
            other_options["collection"]["shape_routes"] = False
            with self.assertRaisesRegex(builder.CapsuleError, "options differ"):
                verifier._validate_manifest(other_options)

            scattered = copy.deepcopy(manifest)
            scattered["collection"]["asset"]["url"] = (
                f"{GITHUB_RELEASES}/demo-data-v3/demo.aicollection"
            )
            with self.assertRaisesRegex(builder.CapsuleError, "outside the finalized"):
                verifier._verify_capsule_records(output, scattered)

            undeclared = copy.deepcopy(manifest)
            del undeclared["collection"]
            with self.assertRaisesRegex(builder.CapsuleError, "disagree about a published"):
                verifier._verify_capsule_records(output, undeclared)

            (output / "demo-locations.json").write_text(
                json.dumps(
                    {
                        "schema_version": 1,
                        "locations": [
                            {"identity": ARCHIVE_ROOT, "path": "archive-a.aie"}
                        ],
                    }
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(builder.CapsuleError, "does not resolve exactly"):
                verifier._verify_collection(
                    aie,
                    manifest,
                    output / "demo.aicollection",
                    output / "demo-locations.json",
                )

    def test_v1_capsules_without_a_collection_stay_valid(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output, aie, release, wheel, _ = self._finalize_fixture(Path(temporary))
            manifest = json.loads((output / "demo-manifest.json").read_text(encoding="utf-8"))
            self.assertNotIn("collection", manifest)
            self.assertNotIn(
                "demo.aicollection", {path.name for path in output.iterdir()}
            )
            readme = (output / "README.md").read_text(encoding="utf-8")
            self.assertIn("rebuild each `.aicollection`", readme)
            result = verifier.verify(output, aie, release, wheel)
            self.assertNotIn("collection_root", result)
            self.assertTrue(result["stories_executed"])

    def test_drilldown_verifier_requires_membership_witness_table(self) -> None:
        def table(name: str, fields: list[str], rows: list[list[object]]) -> dict[str, object]:
            return {
                "name": name,
                "schema": {"fields": [{"name": field} for field in fields]},
                "rows": rows,
            }

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            build_dir = root / "build"
            self._write_build_directory(build_dir)
            record = json.loads((build_dir / "BUILD-RECORD.json").read_text(encoding="utf-8"))
            resources = {
                name: build_dir / asset["filename"]
                for name, asset in record["resources"].items()
            }
            junction = {
                "data": {
                    "tables": [
                        table(
                            "samples",
                            ["sample", "present", "umis"],
                            [["donor-a", True, 1]],
                        )
                    ],
                    "summary": {"supporting_samples": 1},
                }
            }
            missing_memberships = {
                "data": {
                    "tables": [
                        table("predicates", ["name"], [["splice"], ["tail"]]),
                        table(
                            "patterns",
                            ["matched_predicates", "selection_state", "evidence_units"],
                            [[ ["splice", "tail"], "true", 1]],
                        ),
                    ],
                    "summary": {"selected_units": 1},
                }
            }
            with patch.object(
                verifier,
                "_run_json",
                side_effect=[{}, junction, missing_memberships],
            ):
                with self.assertRaisesRegex(builder.CapsuleError, "frozen typed table"):
                    verifier._verify_drilldown_story(
                        root / "aie",
                        record["stories"]["junction_drilldown"],
                        resources,
                        root,
                    )


if __name__ == "__main__":
    unittest.main()
