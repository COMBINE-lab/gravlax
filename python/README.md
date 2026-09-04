# Gravlax Python client

The `gravlax-client` package is a small, dependency-free wrapper around the
`aie` executable. It uses argument arrays—never a shell—and parses the
versioned project, resolved-plan, doctor, result-envelope, and MEX contracts.
The Rust executable remains the only implementation of archive semantics.

## Install

Build or install `aie`, ensure it is on `PATH`, then install the client from a
checkout:

```sh
python -m pip install ./python
```

Install only the integrations a notebook needs:

```sh
python -m pip install './python[pandas,arrow]'
python -m pip install './python[anndata]'
```

## Projects, plans, and diagnostics

```python
from gravlax import Client

aie = Client()  # or Client(binary="/opt/gravlax/bin/aie")

project = aie.project_show(project="analysis")
resolved = aie.plan_check(
    "analysis/plans/replay.yaml",
    project="analysis",
    explain=True,
)
for step in resolved.steps:
    print(step.id, step.args)

report = aie.doctor(["analysis/data/sample.aie"], project="analysis")
for check in report.checks:
    print(check.status, check.summary, check.remedy or "")

if report.ok:
    aie.plan_run("analysis/plans/replay.yaml", project="analysis")
```

The client reads resolved-plan v3, v4, v5, and v6. Version 4 adds biological
intent, assembly-compatibility evidence, output-schema IDs, and conservative
I/O estimates. Version 5 adds paired annotation-comparison intent and validates
the annotation roles and typed output schemas for comparison and transcript-
equivalence steps. Version 6 exposes each explicit uniform result/report format,
stdout or atomic-file publication mode, and destination through
`step.uniform_io`. All preserve access to `resolved.producer`, embedded cohort-
design resources, typed prior-step inputs, and every step's named inputs,
canonical prepared inputs, final outputs, semantic output roles, and staging
paths. These provenance fields are validated rather than silently ignored by
the client. Plan fields accept `step:<id>` or
`step:<id>:<output-name>` to consume a compatible output from an earlier
declaration.

Use `project_add(..., external=True)` only for an intentional absolute,
read-only input outside the otherwise portable project. An interrupted plan
can be continued with `plan_run(..., resume=True)`; `aie` verifies the input,
step, and output identities in its versioned completion records before it
skips anything.

`doctor()` returns its complete report even when checks fail; inspect
`report.ok` or `report.exit_code`. Other unsuccessful commands raise
`CommandError`, whose `.result` retains stdout, stderr, return code, and the
exact argument vector.

## Command-specific JSON and large output

Commands with their own `--json` modes retain those command-specific schemas
for byte compatibility. Parse one without treating it as a shared envelope:

```python
raw = aie.result_raw([
    "query", "sample.aie", "region", "chr1:1000000-2000000", "--json"
])
print(raw["schema"])
```

For a large JSON, TSV, or binary response, stream stdout to a new file instead
of retaining it in Python memory:

```python
written = aie.run_to_file(
    ["query", "sample.aie", "junctions", "chr1:1-100000000", "--tsv"],
    "junctions.tsv",
)
print(written.bytes)
```

The destination is installed only after a successful command and is not
overwritten unless `replace=True` is explicit.

### Uniform named-table bundles

Region and exact-junction queries have dedicated Python convenience methods for
their opt-in uniform JSON interface. They return a strict named-table bundle;
other commands that emit this contract can use the generic bundle methods
shown below:

```python
region = aie.query_region(
    "sample.aie",
    "chr16:89550000-89575000",
    top=20,
)
print(region.summary.umis, region.summary.cells)

counts = region.table("counts")
print(counts.semantics.row_semantics, counts.semantics.key)
print(counts.selection.available_rows, counts.selection.truncated)
print(counts.records())
```

`query_junction()` returns the corresponding typed junction summary and the
same count-table shape. `top=0` means all rows in this uniform interface.
Physical row order remains distinct from logical set/multiset/sequence
semantics; the parser validates declared keys and ordering-field references but
does not invent an ordering.

Use the file variants when a result may be large. They keep subprocess stdout
out of Python memory and atomically install the completed file:

```python
written = aie.query_junction_to_file(
    "sample.aie",
    "chr16:89562391-89562883",
    "junction.json",
    top=0,
)
bundle = aie.result_bundle_from_file(written.output_path)
```

Boolean evidence-unit queries and atlas-wide event discovery have argument-safe
wrappers and matching bounded-memory file variants:

```python
cooccurrence = aie.query_cooccurrence(
    "sample.aie",
    {
        "locus": "region:chr1:155230000-155240000:+",
        "splice": "junction:chr1:155234452-155235327:+",
        "tail": "terminal:chr1:155239900-155240025:+",
    },
    "locus & splice & !tail",
    universe="locus",
)
for pattern in cooccurrence.table("patterns").records():
    print(pattern["pattern_mask"], pattern["selection_state"])

written = aie.collection_find_events_to_file(
    "atlas.aicollection",
    "events.json",
    kinds=("junction", "cassette", "terminal-tail"),
    design="donors.tsv",
    groups="groups.tsv",
    min_donors=3,
)
```

`selection_state="unknown"` preserves an unresolvable absence when a
two-representative chain omitted middle read placements; a positive predicate
always has a retained witness. `unit="umi-class"` requires
`allow_full_scan=True` and describes a barcode-corrected cell plus exact raw
UMI-value class; it does not collapse one-mismatch UMI edges and is not proof
of one physical molecule. The default `placements="unique"` excludes
multimappers; `direct` and `all` are explicitly diagnostic placement modes.
`collection_find_events_to_file()` rebuilds unique-chain, exact raw-UMI-value
class counts by sample, donor, and group from the collection's rooted source
archives and keeps the result stream out of Python memory.

For any command that emits the same JSON contract, `result_bundle(args)` and
`parse_uniform_bundle(document)` expose unique named tables, typed fields,
row semantics, and exact or deferred selection metadata. A deferred one-pass
selection represents unknown availability and truncation explicitly as null;
it is never silently treated as complete.

## Shared typed results

The package also implements the `gravlax.result-envelope.v1` contract for
producers that explicitly advertise it. It is deliberately separate from
`result_raw()` because each command's own `--json` output keeps its own schema:

```python
from gravlax import ResultEnvelope

resolved = aie.resolve(
    "gencode.v49.aic",
    ["TP53", "transcript:ENST00000269305"],
    assembly="GRCh38.p14",
    annotation="GENCODE 49",
)
print(resolved.table.records())
print(resolved.provenance.annotation_digest)

comparison = aie.compare_annotations(
    "sample.aie",
    "gencode.v44.gtf",
    "gencode.v49.aic",
    assembly="GRCh38.p14",
    annotation_a_label="GENCODE 44",
    annotation_b_label="GENCODE 49",
)
print(comparison.count_deltas.records())

ecs = aie.transcript_ecs(
    "sample.aie",
    "gencode.v49.aic",
    assembly="GRCh38.p14",
    annotation_label="GENCODE 49",
    feature="gene:ENSG00000141510",
    aggregation="bulk",
)
print(ecs.catalog.records())

result = ResultEnvelope.from_file("typed-result.json")
records = result.table.records()  # dependency-free list of dictionaries
frame = result.to_pandas()        # requires the pandas extra
arrow = result.to_arrow()         # requires the arrow extra
observations = result.to_anndata(obs_names="cell")
```

Annotation comparison is exact only within the retained archive quotient and
fixed alignment/barcode policy; its causes are non-additive explanations.
Transcript ECs are retained-evidence compatibility sets rather than abundance,
isoform calls, or phasing. They are derived from the representatives retained
in the archive rather than every source read. Their conflict and
no-compatible-transcript flags are non-exclusive.

Command-specific query `--json` responses retain their own schemas. Every
supported result-streaming query family exposes the typed contract explicitly
through `--format=json`. Commands that advertise a separate typed operation
report use `--report-format=json`; that is a per-command interface, not a rule
for every artifact producer. Native Arrow IPC is not yet a general CLI output,
and no R client is currently provided.

The AnnData conversion above preserves a row-oriented envelope in `.obs` and
does not guess a count matrix. For a MEX result carrying the shared
`metadata.json` completion marker, use the matrix-aware reader instead:

```python
from gravlax import read_mex

mex = read_mex("analysis/results/counts")
matrix = mex.to_scipy()       # exact feature-by-barcode orientation
adata = mex.to_anndata()      # conventional cell-by-feature orientation
```

MEX loading requires `metadata.json`, validates every coordinate, bounds,
duplicate, declared nonzero, label count, and file path, and carries the result
schema and provenance into `adata.uns["gravlax"]`.

## Test

The core test suite needs no optional scientific Python packages:

```sh
cd python
PYTHONPATH=src python -m unittest discover -s tests -v
```
