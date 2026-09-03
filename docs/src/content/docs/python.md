---
title: Python and AnnData
description: Validate projects, run plans, and load typed Gravlax results from Python.
---

The `gravlax-client` package is a dependency-light companion to the `aie`
binary. It executes explicit argument arrays without a shell, so paths and
identifiers remain single literal arguments. Archive interpretation stays in
the Rust executable; Python consumes the same versioned JSON plans and result
schemas as the command line and Explorer.

## Install from a checkout

Install `aie` first, then choose only the notebook integrations you need:

```sh
python -m pip install ./python
python -m pip install './python[pandas,arrow]'
python -m pip install './python[anndata]'
```

The base package has no third-party runtime dependencies.

## Check and run a plan

```python
from gravlax import Client

aie = Client()
project = aie.project_show(project="analysis")
plan = aie.plan_check(
    "analysis/plans/replay.yaml",
    project="analysis",
    explain=True,
)

print(project.name, plan.source_digest)
for step in plan.steps:
    print(step.id, step.args)

report = aie.doctor(project="analysis")
if report.ok:
    aie.plan_run("analysis/plans/replay.yaml", project="analysis")
```

Resolved-plan v3, v4, v5, and v6 fields are typed rather than discarded. Version 4
adds biological intent and assembly-compatibility evidence, output-schema IDs,
exact known selected-input sizes, and explicitly available or unavailable
execution-read bounds. Version 5 adds paired annotation-comparison intent and
validates the annotation roles and typed output schemas for comparison and
transcript-equivalence steps. Version 6 adds the selected uniform I/O mode,
publication semantics, and report/result destinations. `plan.producer` binds the exact `aie` executable
identity, `plan.embedded_resources` records files named inside cohort designs,
and each step exposes its named resources, typed prior-step inputs, canonical
prepared inputs, semantic outputs, and no-clobber staging paths.

For an interrupted run, `plan_run(..., resume=True)` only skips steps whose
versioned completion records still match the resolved plan and output
identities. External project resources must be registered intentionally with
`project_add(..., external=True)` and remain absolute, read-only inputs.

Failed doctor checks are data, so `doctor()` returns the report with
`report.ok == False` and the command's exit code. Unexpected command failures
raise `CommandError` and retain their exact argument vector and diagnostics.

## Existing JSON and large output

```python
raw = aie.result_raw([
    "query", "sample.aie", "region", "chr1:1000000-2000000", "--json"
])
print(raw["schema"])

written = aie.run_to_file(
    ["query", "sample.aie", "junctions", "chr1:1-100000000", "--tsv"],
    "junctions.tsv",
)
```

Existing query JSON remains command-specific. `result_raw()` parses it without
claiming a shared schema. `run_to_file()` keeps large JSON, TSV, or binary
stdout out of Python memory, installs the destination only after success, and
refuses an existing file unless replacement was explicit.

## Shared typed result envelopes

The client implements `gravlax.result-envelope.v1` for producers that
explicitly advertise that contract. Identifier resolution is an envelope
producer and fails the whole batch on ambiguity:

```python
resolved = aie.resolve(
    "gencode.v49.annotation.aic",
    ["TP53", "transcript:ENST00000269305"],
    assembly="GRCh38.p14",
    annotation="GENCODE 49",
)

print(resolved.table.records())
print(resolved.provenance.annotation_digest)
```

The assembly, release, and observed annotation-content digest are carried in
the result provenance. For direct Python calls, `assembly=` is a caller
assertion; use registered resources in a checked project plan when assembly
compatibility must be verified rather than merely recorded.

Paired annotation comparison and experimental transcript equivalence also have
typed client methods:

```python
comparison = aie.compare_annotations(
    "sample.aie",
    "gencode.v44.gtf",
    "gencode.v49.aic",
    assembly="GRCh38.p14",
    annotation_a_label="GENCODE 44",
    annotation_b_label="GENCODE 49",
)
print(comparison.count_deltas.records())
print(comparison.contributing_causes.records())

ecs = aie.transcript_ecs(
    "sample.aie",
    "gencode.v49.aic",
    assembly="GRCh38.p14",
    annotation_label="GENCODE 49",
    feature="gene:ENSG00000141510",
    aggregation="bulk",
)
print(ecs.catalog.records())
print(ecs.counts.records())
```

Comparison deltas are exact only within the retained archive quotient under the
fixed alignment/barcode policy; its causes are non-additive explanations.
Transcript ECs are retained-evidence compatibility sets, not abundance,
full-read equivalence, isoform calls, or phasing. Transcript ECs can differ
from classes derived from every original read, and their
`no_compatible_transcript` and `conflict` flags are
non-exclusive.

Historical query `--json` responses remain command-specific for byte
compatibility. Supported result-streaming query families expose the typed
envelope explicitly through `--format=json`. Commands that advertise a
separate uniform operation report use `--report-format=json`; that is a
per-command interface, not a rule for every artifact producer. A compatible
saved envelope can also be loaded directly:

```python
from gravlax import ResultEnvelope

result = ResultEnvelope.from_file("typed-result.json")
rows = result.table.records()
frame = result.to_pandas()
arrow_table = result.to_arrow()
```

The envelope parser verifies the result schema, field widths, nullability, and
logical scalar types. Arrow output preserves the envelope metadata and marks
logical JSON columns explicitly.

`to_arrow()` constructs an in-memory Arrow table in Python; it does not imply
that CLI commands support `--format arrow` or emit Arrow IPC. Supported
scientific commands can emit shared text, TSV, and JSON results; native Arrow
IPC and an R client are not currently provided.

For row-oriented envelopes, `result.to_anndata(obs_names="cell")` puts the
table in `.obs` without inventing a numeric matrix. A MEX result with the
shared `metadata.json` completion marker has a matrix-aware path:

```python
from gravlax import read_mex

counts = read_mex("analysis/results/counts")
adata = counts.to_anndata()
```

This validates the completion marker, every declared coordinate and label,
bounds, duplicates, and truncation; transposes the stored feature-by-barcode
matrix to AnnData's cell-by-feature convention; fills `.obs` and `.var`; and
preserves provenance in `.uns["gravlax"]`.
The package's full API and dependency-free test command are in
[`python/README.md`](https://github.com/COMBINE-lab/gravlax/blob/main/python/README.md).
