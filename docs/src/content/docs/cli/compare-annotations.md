---
title: aie compare-annotations
description: Measure and explain annotation changes exactly within one retained archive quotient.
---

`aie compare-annotations` replays one evidence archive against two annotations
under the same assignment policy. It reports the exact signed change in final
gene/UMI counts within that retained archive quotient and a structural ledger
explaining which UMI-class states changed. The upstream genome alignments and
barcode correction remain fixed. The archive is decoded once; each annotation
nevertheless receives its own class aggregation and final 1-mismatch UMI
collapse.

```sh
aie compare-annotations sample.aie \
  --annotation-a gencode.v44.annotation.gtf \
  --annotation-b gencode.v49.annotation.aic \
  --assembly GRCh38.p14 \
  --annotation-a-label "GENCODE 44" \
  --annotation-b-label "GENCODE 49"
```

Source GTF and compiled `.aic` inputs can be mixed. Each file is opened once,
hashed, verified when an expected digest is supplied, rewound, and parsed from
that same open file. This binds the model used for replay to the digest in the
result provenance.

For this direct path-based command, `--assembly` is a caller assertion recorded
in provenance. Neither matching labels nor content digests prove that the
annotation and archive use the same reference assembly. A
[project plan](/gravlax/cli/projects/) can additionally verify declared
archive/annotation assembly compatibility and records `unverified` when the
required coordinate-resource identity is unavailable.

## What is exact

The `count-deltas` table contains nonzero signed `B - A` changes after two
complete, independent final collapses. These rows are not obtained by simply
subtracting class assignments. A changed class can alter which neighboring
UMIs survive collapse, so each annotation side is collapsed globally before
counts are compared.

This exactness is scoped to the evidence retained in `.aie`, the fixed
alignment/barcode policy that produced it, and the selected strand, gene-key,
assignment, and collapse semantics. It is not a full-read counterfactual and
does not predict the result of fresh annotation-aware realignment or different
barcode correction.

The `class-transitions` table is a complete ledger of changed UMI-class states.
It includes each side's selected gene, support, same-gene neighbors, canonical
class, and final-count status. It intentionally has no class-level count-delta
column: class transitions are not additive contributions to gene-count
deltas.

The `contributing-causes` table names observed changes such as candidate-set,
class-winner, collapse-neighborhood, or final-contribution changes. Causes are
non-exclusive state descriptions, not unique counterfactual attributions, and
must not be summed. `annotation_order_tie_break_changed` is an explicit method
artifact: both sides have the same comparison-key support with a tied maximum,
but their local gene order selects different winners. It is reported so exact
replay parity is explainable; it is not a biological structural change.

The `witnesses` table contains deterministic examples of changed molecule
records and rows. Witnesses are bounded by
`--max-molecule-witnesses` (default 10,000) and
`--max-row-transitions-per-molecule` (default 32). Exact totals and omitted
counts remain in the summary and class ledger even when no witnesses are
requested.

## Gene identity across releases

`--gene-key unversioned` is the default. It removes exactly one terminal
`.digits` suffix, so a change such as `ENSG….12` to `ENSG….13` is not reported
as a lost and gained gene. The command rejects an annotation if this
normalization would merge two IDs within that annotation. Use
`--gene-key exact` when versioned IDs must remain distinct. Gene symbols are
never used as joins. Output retains the comparison key and the original ID
from each side.

Every table identifies a cell by its human-readable 16-base barcode. The
archive's dense numeric cell ID is retained in a separate `cell` column for
stable joins and diagnostics.

## Output

Text is the default and gives a concise scientific summary. JSON emits one
`gravlax.result-envelope.v1` object whose result schema is
`gravlax.annotation.compare.v1`. Its `data` contains a summary, an explicit
semantics block, and four typed tables with independent versioned schemas:

| JSON member / TSV selector | Table schema |
|---|---|
| `count_deltas` / `count-deltas` | `gravlax.annotation.compare.count-deltas.v1` |
| `class_transitions` / `class-transitions` | `gravlax.annotation.compare.class-transitions.v1` |
| `contributing_causes` / `contributing-causes` | `gravlax.annotation.compare.contributing-causes.v1` |
| `witnesses` / `witnesses` | `gravlax.annotation.compare.witnesses.v1` |

TSV emits one selected table and therefore requires `--table`. `--table` is
rejected with text or JSON output.

```sh
# Complete machine-readable comparison
aie compare-annotations sample.aie \
  --annotation-a old.gtf --annotation-b new.gtf \
  --assembly GRCh38.p14 \
  --annotation-a-label old --annotation-b-label new \
  --format json -o comparison.json

# One pipeline-safe table
aie compare-annotations sample.aie \
  --annotation-a old.gtf --annotation-b new.gtf \
  --assembly GRCh38.p14 \
  --annotation-a-label old --annotation-b-label new \
  --format tsv --table count-deltas > count-deltas.tsv
```

JSON and TSV carry paired annotation provenance with `before` and `after`
roles, observed annotation digests, replay parameters, and the archive's
authenticated v2 directory root. Legacy v1 archives have no rooted
commitment; the command reports that limitation and recommends sealing or
rewriting the archive as v2.

`-o/--output` installs a fully rendered result atomically at a new path.
Existing files are never replaced, and validation or digest failures leave no
partial output. If the two observed annotation digests are identical, the
command fails unless `--allow-identical` is given for an explicit A/A control.

## Options

| Option | Default | Description |
|---|---|---|
| `--annotation-a <A>` | required | Before annotation, as GTF or `.aic` |
| `--annotation-b <B>` | required | After annotation, as GTF or `.aic` |
| `--assembly <ASSEMBLY>` | required | Caller-asserted shared reference assembly identity; project plans can additionally verify compatibility |
| `--annotation-a-label <LABEL>` | required | Before annotation provenance label |
| `--annotation-b-label <LABEL>` | required | After annotation provenance label |
| `--annotation-a-digest <DIGEST>` | — | Expected `blake3:<64 lowercase hex>` digest |
| `--annotation-b-digest <DIGEST>` | — | Expected `blake3:<64 lowercase hex>` digest |
| `--gene-key <POLICY>` | `unversioned` | `unversioned` or `exact` |
| `--solo-strand <STRAND>` | `forward` | `forward`, `reverse`, or `unstranded`, shared by both sides |
| `--max-molecule-witnesses <N>` | `10000` | Global witness-row bound; zero retains none |
| `--max-row-transitions-per-molecule <N>` | `32` | Changed-row bound inside each witness; zero retains none |
| `--allow-identical` | off | Permit an explicit identical-digest A/A comparison |
| `--format <FORMAT>` | `text` | `text`, `json`, or `tsv` |
| `--table <TABLE>` | — | Required only for TSV |
| `-o, --output <FILE>` | stdout | Atomically create a new output file |
