---
title: Diagnostic commands
description: BAM comparison, replay inspection, archive accounting, and model-analysis commands under aie dev.
---

`aie dev` contains diagnostic and research-oriented commands. They expose
intermediate comparisons and measurements rather than the supported archive,
replay, and query result contracts. The command names `gate-c` and `gate-d`
are retained for compatibility; the descriptions below state the quantities
they compute.

The former top-level spellings are hidden aliases. Prefer `aie dev <command>`.
As with every command, `aie dev <command> --help` prints the option reference
for the installed version.

## `gate-c` — compare alignment evidence

Compares the annotation-independent evidence tuple for reads in two BAM files,
usually an annotation-free alignment and an annotation-aware alignment of the
same reads. The optional STARsolo BAM supplied with `--gene-assigned-from`
restricts the comparison to read names with a `GX` gene assignment.

```sh
aie dev gate-c annotation-free.bam annotation-aware.bam \
  --gene-assigned-from annotation-aware.bam --json-out alignment-diff.json
```

| Argument or option | Description |
|---|---|
| `<A>` | Baseline BAM, normally the annotation-free alignment |
| `<B>` | Comparison BAM |
| `--gene-assigned-from <BAM>` | Keep only reads with a STARsolo `GX` assignment in this BAM |
| `--json-out <PATH>` | Write per-category counts as JSON |

## `gate-d` — compare UMI grouping

Compares annotation-independent locus/UMI grouping with STARsolo per-cell,
per-gene UMI grouping. The input BAM must be coordinate sorted and carry
`CB`, `UB`, `UR`, and `GX` tags.

| Argument or option | Default | Description |
|---|---|---|
| `<BAM>` | required | STARsolo BAM used for the comparison |
| `--locus-gap <BP>` | `500` | Maximum gap between consecutive read starts assigned to one molecule |
| `--cells <PATH>` | all barcodes | Restrict summaries to barcodes in a STARsolo `filtered/barcodes.tsv` file |
| `--correct-umi` | off | Apply annotation-free one-mismatch UMI merging within each locus |
| `--json-out <PATH>` | — | Write the summary as JSON |

## `build` — measure evidence size

Extracts molecular evidence from a coordinate-sorted, annotation-free BAM and
reports the size and composition of alternative evidence representations. It
does not create the production archive used by `ingest-archive`.

| Argument or option | Default | Description |
|---|---|---|
| `<BAM>` | required | BAM carrying raw `CR` and `UR` tags |
| `--whitelist <PATH>` | required | One accepted 16-base barcode per line |
| `--locus-gap <BP>` | `50000` | Gap used to join reads into a locus |
| `--zstd-level <N>` | `19` | zstd level used in the size calculation |
| `--no-umi-collapse` | off | Keep exact UMI classes and defer one-mismatch merging |
| `--json-out <PATH>` | — | Write measurements as JSON |

## `umi-graph` — inspect UMI adjacency

Measures the cell-scoped one-mismatch UMI graph and compares graph-based
replay with grouping derived from `CB`, `UB`, and `GX` tags when those tags are
present.

| Argument or option | Default | Description |
|---|---|---|
| `<BAM>` | required | BAM with `CR`/`UR` and, for the comparison, `CB`/`UB`/`GX` |
| `--cells <PATH>` | all barcodes | Restrict per-cell summaries to the listed barcodes |
| `--locus-gap <BP>` | `2000` | Gap used by the base ingest-time molecule definition |
| `--store-window <BP>` | `3000000` | Maximum genomic span for retained graph edges |
| `--json-out <PATH>` | — | Write measurements as JSON |

## `replay` — replay directly from a BAM

Extracts rows from an annotation-free BAM and quantifies a GTF without first
creating an archive. It also exposes switches that isolate barcode correction,
assignment, multimapper handling, representative selection, and UMI collapse.
For ordinary archive quantification, use
[`replay-rows`](/gravlax/cli/replay-rows/).

| Argument or option | Default | Description |
|---|---|---|
| `<BAM>` | required | Coordinate-sorted annotation-free BAM with `CR` and `UR` tags |
| `--gtf <PATH>` | required | Annotation to quantify |
| `--whitelist <PATH>` | required | Barcode whitelist used to correct raw `CR` tags |
| `--barcodes <PATH>` | required | Barcode list and output-column order |
| `--out-dir <PATH>` | required | Directory for `matrix.mtx`, `features.tsv`, and `barcodes.tsv` |
| `--locus-gap <BP>` | `2000` | Gap used to join reads into a locus |
| `--use-cb-tag` | off | Use the BAM's corrected `CB` tag instead of correcting `CR` |
| `--simple-bc` | off | Use exact-or-unique-one-mismatch barcode correction |
| `--no-collapse` | off | Count distinct raw UMIs without one-mismatch collapse |
| `--no-multigene-filter` | off | Keep UMI classes compatible with more than one gene |
| `--read-level` | off | Classify every read rather than molecule representatives |
| `--chain-representative` | off | With `--read-level`, keep the most-contained unique read for each junction chain |
| `--two-reps` | off | With `--chain-representative`, also keep the most-extended read when distinct |
| `--no-multimappers` | off | Ignore reads with more than one genomic alignment |

## `assign-diff` — compare gene assignment

Compares the STARsolo `GX` assignment in a BAM with Gravlax's local
`alignToTranscript` assignment for the same alignments and GTF.

| Argument or option | Default | Description |
|---|---|---|
| `<BAM>` | required | STARsolo BAM containing `GX` tags |
| `--gtf <PATH>` | required | The annotation used to produce the STARsolo assignments |
| `--examples <N>` | `5` | Maximum example reads printed for each disagreement category |
| `--json-out <PATH>` | — | Write category counts as JSON |

## `sig-stats` — inspect signatures and coding

Measures evidence-signature multiplicity, alternative-placement sharing, site
offsets, and compressed stream sizes.

| Argument or option | Default | Description |
|---|---|---|
| `<BAM>` | required | Annotation-free BAM with raw `CR`/`UR` tags and secondary alignments |
| `--whitelist <PATH>` | required | Barcode whitelist |
| `--locus-gap <BP>` | `2000` | Gap used to join reads into a locus |
| `--site-gap <BP>` | `64` | Gap used to cluster molecule starts into a 3′ site |
| `--zstd-level <N>` | `19` | zstd level used in compressed-size measurements |
| `--json-out <PATH>` | — | Write measurements as JSON |

## `debug` — inspect an encoded archive

Reports per-section and per-stream byte accounting, order-0 entropy estimates,
and dictionary-sharing ratios for an existing archive.

| Argument or option | Description |
|---|---|
| `<ARCHIVE>` | Archive to inspect |
| `--dump-dir <PATH>` | Write concatenated raw chunk streams and dictionary/index sections as `.bin` files |

## `em` — analyze multimapper recovery

Compares per-cell, pooled, group-aware, and evidence-depth-dependent recovery
models on masked evidence. With `--mask 0`, `--emit`, and `--barcodes`, it can
write the additive pooled recovered-count layer. The other model-only modes
write evaluation results and cannot emit counts. See the complete
[`aie dev em` option reference](/gravlax/cli/em/).
