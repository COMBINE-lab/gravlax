---
title: Projects and plans
description: Name inputs once, validate an analysis, and preserve its exact resolved form.
---

Projects make repeated Gravlax work easier to read and reproduce. A project is
a portable directory with one `aie-project.yaml` manifest, source plans under
`plans/`, ordinary outputs under `results/`, and immutable resolved-plan
snapshots under `.aie/resolved-plans/`.

To try the complete interface without a large dataset, use the repository's
`examples/demo-project`. Its checked plan compiles a tiny demonstration GTF,
persists a resolved snapshot, and produces an artifact the Explorer can show.

## Create a project

```sh
aie project init experiments/pbmc --name pbmc-reanalysis
cd experiments/pbmc
```

For a fully movable workspace, keep inputs inside the project and give each one
a stable biological or workflow name:

```sh
aie project add pbmc-v1 data/pbmc-v1.aie --kind archive
aie project add gencode-v49 annotations/gencode.v49.aic --kind annotation \
  --assembly GRCh38.p14 --annotation-label "GENCODE 49"
aie project add filtered-cells metadata/barcodes.tsv --kind barcodes
aie project show
```

Large read-only inputs can remain in shared storage when copying them would be
wasteful. Make that portability tradeoff explicit:

```sh
aie project add pbmc-atlas /shared/archives/pbmc-atlas.aie \
  --kind archive --external
```

An external entry stores its canonical absolute path and `external: true`.
Ordinary entries remain relative and movable. Symlinks that escape the project
are still rejected; use `--external` so reviewers can see the dependency.

Resource kinds cover archives, collections, annotations, genomes, BAMs,
barcodes, whitelists, groups, cell lists, cohort designs, metadata, and a
generic file escape hatch. `--kind auto` recognizes ordinary `.aie`,
`.aicollection`, `.aic`, GTF/FASTA/BAM, barcode, whitelist, group, cell,
design, and metadata names; use an explicit kind when a generic filename would
be ambiguous.

Annotations used for identifier-based plans carry both an assembly and an
immutable annotation label. The resolved plan computes and records the exact
annotation-content digest separately, so an editable label never substitutes
for byte identity. Archives, collections, and genomes may also declare an
assembly:

```sh
aie project add pbmc-v1 data/pbmc-v1.aie --kind archive \
  --assembly GRCh38.p14 --replace
```

When a feature plan selects a coordinate resource with an assembly label, the
checker requires an exact match. Without such a label it reports compatibility
as `unverified`; a shared chromosome name alone is not treated as proof of an
assembly. For readable archives the checker also records chromosome and genome
signatures and rejects a resolved contig absent from the archive. This is
stronger than a direct path-based command, where `--assembly` is a caller
assertion recorded in provenance: project compatibility is what verifies that
assertion against registered coordinate-resource metadata when available.

Commands contain every output beneath the canonical project root. Inputs are
either contained there too or are visibly declared external; a symlink cannot
silently change that choice.

### Project command options

| Command | Argument or option | Default | Description |
|---|---|---|---|
| `project init` | `[DIRECTORY]` | `.` | Directory to initialize |
| `project init` | `--name <NAME>` | directory name | Human-readable project name |
| `project add` | `<NAME> <PATH>` | required | Stable resource name and existing input path |
| `project add` | `--kind <KIND>` | `auto` | Resource type: `auto`, `archive`, `collection`, `annotation`, `genome`, `bam`, `barcodes`, `whitelist`, `groups`, `cells`, `design`, `metadata`, or `file` |
| `project add` | `--external` | off | Register a canonical absolute read-only path outside the project |
| `project add` | `--project <PATH>` | search upward | Project directory or manifest |
| `project add` | `--replace` | off | Replace a resource with the same name |
| `project add` | `--assembly <NAME>` | — | Assembly identity for an annotation, archive, collection, or genome |
| `project add` | `--annotation-label <LABEL>` | — | Annotation label; requires `--assembly` and an annotation resource |
| `project show` | `--project <PATH>` | search upward | Project directory or manifest |
| `project show` | `--json` | off | Emit a versioned JSON document instead of the text view |

## Write a plan

Plans are strict, versioned YAML or JSON. This `plans/replay.yaml` example uses
the names registered above:

```yaml
schema_version: 1
name: replay-gencode-v49
steps:
  - id: quantify
    kind: replay-rows
    archive: pbmc-v1
    annotation: gencode-v49
    barcodes: filtered-cells
    out_dir: results/counts-v49
```

Each step delegates to an existing `aie` operation. The plan layer resolves
names, paths, output ownership, and exact argument lists; it does not contain a
second implementation of replay or query logic.

Plans currently cover archive inspection and ingest, annotation compilation
and extension, replay, paired annotation comparison, transcript-equivalence,
single-archive region/junction/junction-set/event/terminal queries, federation,
collection queries, and cohort event or splice-graph workflows. Query plans can
add a named cell/group scope without carrying raw paths from step to step:

```yaml
  - id: cell-type-events
    kind: query-events
    archive: pbmc-v1
    feature: gene:ENSG00000166349
    annotation: gencode-v49
    scope:
      groups: cell-types
      aggregation: group
    uniform_output:
      format: json
      output: results/cell-type-events.json
```

`uniform_output` explicitly opts a step into the common result envelope. Its
format is `text`, `tsv`, or `json`. Omit the nested `output` to stream the
schema-bearing result to standard output, or give a project-relative path for
atomic no-clobber publication. The older top-level `format` and `output` fields
retain their command-specific behavior, so existing source-plan v1 files keep
their prior command lines. A step cannot combine uniform output with a
non-default legacy format or a legacy output destination.

Each declared artifact is installed complete and without replacing an existing
path. Filesystems do not provide a portable transaction spanning several output
paths. If another process creates a later destination during a multi-output
install, the runner therefore preserves any earlier complete outputs and its
remaining staging artifacts, reports their paths, and stops. It never removes
an installed pathname as rollback, because that pathname could already have
been replaced by the other process.

Every region-like step accepts exactly one of an explicit 0-based half-open
`locus` or a biological `feature`. A compact feature string uses the registered
annotation identity. A qualified request can additionally pin the expectation
in the plan itself and fails if the project metadata differs:

```yaml
    feature:
      identifier: TP53
      assembly: GRCh38.p14
      annotation: GENCODE 49
```

Genes, transcripts, and exons are resolved through the same ambiguity-detecting
resolver as `aie resolve`. The resolved snapshot preserves the requested and
stable IDs, match basis, genes/transcripts, strand, explicit coordinates,
assembly, annotation label, content digest, and compatibility evidence.

Terminal-boundary exploration has a canonical plan step too:

```yaml
  - id: tp53-ends
    kind: query-apa
    archive: pbmc-v1
    feature: TP53
    annotation: gencode-v49
    strand: forward
    groups: cell-types
    site_gap: 24
    uniform_output:
      format: json
      output: results/tp53-ends.json
```

`query-apa` also supports a registered genome, internal-priming exclusion,
permutations and seed, and an SVG plot.

The two annotation-dependent capability queries use the same named-resource,
identity, and compatibility checks. After registering both annotation releases,
a comparison step is concise:

```yaml
  - id: compare-releases
    kind: compare-annotations
    archive: pbmc-v1
    annotation_a: gencode-v44
    annotation_b: gencode-v49
    format: json
    output: results/gencode-v44-v49.json
```

Its signed count deltas are exact within the retained archive quotient and fixed
alignment/barcode policy; transition causes are explanatory and non-additive.
See [annotation comparison](/gravlax/cli/compare-annotations/) for the four
typed tables.

The experimental transcript-equivalence query can use the same biological
resolver and group resource:

```yaml
  - id: tp53-transcript-ecs
    kind: query-transcript-ecs
    archive: pbmc-v1
    annotation: gencode-v49
    feature: gene:ENSG00000141510
    scope:
      groups: cell-types
      aggregation: group
    format: json
    output: results/tp53-transcript-ecs.json
```

This produces compatibility sets and archived UMI-class counts, not transcript
abundance or phasing, and it can differ from equivalence classes derived from
every original read. See the
[transcript-equivalence reference](/gravlax/cli/transcript-ecs/) for the exact
scope and output flags.

An input can also consume a typed output from an earlier declaration. Use
`step:<id>` when that step has one output, or `step:<id>:<output-name>` when it
has several. For example, one plan can ingest an archive, compile an
annotation, and replay both results:

```yaml
  - id: ingest
    kind: ingest-archive
    bam: reads
    whitelist: whitelist
    output: results/sample.aie
    uniform_report:
      format: json
      output: results/ingest-report.json
  - id: compile
    kind: compile-annotation
    annotation: source-gtf
    output: results/genes.aic
    uniform_report:
      format: json
      output: results/compile-report.json
  - id: replay
    kind: replay-rows
    archive: step:ingest:archive
    annotation: step:compile:annotation
    barcodes: filtered-cells
    out_dir: results/counts
    uniform_report:
      format: json
      output: results/replay-report.json
```

Artifact-producing steps use `uniform_report`: their archive, annotation, or
MEX directory remains the primary output while the typed report describes the
operation and its artifacts. `extend-annotation` supports the same report
object. Report files participate in staging, collision checks, completion
records, and resume identity checks.

References are declaration-ordered: forward references and cycles fail, as do
ambiguous output names or role mismatches such as using an annotation output
where an archive is required.

For cohort event plans, `samples` maps each stable sample ID to a named archive
resource and `groups` maps the same sample IDs to named group resources.
Replicate-aware `cohort-splice-graph` steps instead refer to one registered
`design` resource. The checker rejects missing resources and mismatched sample
maps before any evidence is decoded.

## Check before running

```sh
aie plan check plans/replay.yaml --explain
```

`--explain` shows every named-resource resolution, biological resolution,
assembly-compatibility status, output-schema identifier, exact sizes for known
selected input files, output, and delegated command. Execution read bounds are
reported as unavailable unless the checker can prove one; command read
multiplicity, route-dependent collection members, and not-yet-produced step
outputs are not presented as exact I/O predictions. It is safe to use while
reviewing a plan: nothing is written or run.
Use `--json` to obtain the complete resolved-plan document.

### Plan command options

Both commands take one `<PLAN>` YAML or JSON path.

| Command | Option | Default | Description |
|---|---|---|---|
| `plan check` / `plan run` | `--project <PATH>` | search upward from plan | Project directory or manifest |
| `plan check` / `plan run` | `--explain` | off | Show every resource, output, and delegated-command resolution |
| `plan check` | `--json` | off | Emit the complete resolved plan as JSON |
| `plan run` | `--dry-run` | off | Resolve and display the plan without writing or executing |
| `plan run` | `--resume` | off | Skip only steps whose recorded plan and output identities still match |

Resolved-plan schema v6 records an explicit `uniform_io` object for every
opted-in step: result versus report, selected format, stdout versus atomic
no-clobber file publication, and the final destination when present.
`output_schema_ids` names each primary artifact contract plus the outer
result/report schema and its emitted named tables. File arguments in the
command vector point at recorded staging paths;
the runner installs them only after successful completion. Artifact reports
still name the final logical destinations rather than those private staging
paths.

For a rehearsal of the run interface, use:

```sh
aie plan run plans/replay.yaml --dry-run
```

A real run resolves the plan again, writes its deterministic JSON snapshot
beneath `.aie/resolved-plans/`, and only then invokes the existing commands:

```sh
aie plan run plans/replay.yaml --explain
```

Every resource path, input size and content identity, project-manifest digest,
source-plan digest, output, and exact argument vector is preserved in the
snapshot. Rooted archives and collections use their cheap committed identity;
other inputs receive a full-file BLAKE3 identity. Child commands write to
recorded staging paths and final artifacts are installed without overwriting.

Resume an interrupted or already completed plan explicitly:

```sh
aie plan run plans/replay.yaml --resume
```

A step is skipped only when its resolved-plan digest, step digest, producer
executable, upstream step-output identities, and freshly observed full output
identities match the versioned record under `.aie/completions/`. Missing,
partial, changed, or unrecorded outputs stop resume with an error. Keep `.aie/`; removing it
discards both resolved snapshots and the evidence required for exact resume.

Use [`aie explore`](/gravlax/cli/explore/) for a local, read-only view of the
manifest, resource status, resolved snapshots, and results. Its scientific
builder also resolves registered biological identifiers and compiles a
single-step plan in memory, showing exact coordinate, identity, route, scope,
schema, and I/O information before you copy an export. Explorer never runs or
saves the generated plan.
