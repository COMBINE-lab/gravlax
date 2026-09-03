---
title: aie replay-rows
description: Quantify a compatible GTF from an .aie index using Gene or Velocyto semantics.
---

Quantify a compatible GTF from an `.aie` archive. Replay uses the fixed genome
alignments and barcode correction recorded at ingest; it does not rerun an
annotation-aware aligner.

## Usage

```sh
aie replay-rows [OPTIONS] --gtf <GTF> --barcodes <BARCODES> --out-dir <OUT_DIR> <INPUT>
```

```sh
aie replay-rows sample.aie \
  --gtf gencode.v49.annotation.gtf \
  --barcodes barcodes.tsv \
  --out-dir counts-v49/
```

## Arguments

| Argument | Description |
|---|---|
| `<INPUT>` | `.aie` archive, ingest BAM with `--from-bam`, or molecule BAM with `--from-molecule-bam` |

## Options

| Option | Default | Description |
|---|---|---|
| `--gtf <GTF/AIC>` | required | The GTF or compiled `.aic` annotation to replay; coordinates and contig names must match the archive's reference genome |
| `--barcodes <BARCODES>` | required | Barcode list defining output column order (e.g. a raw `barcodes.tsv`) |
| `--out-dir <OUT_DIR>` | required | Receives `matrix.mtx`, `features.tsv`, `barcodes.tsv` |
| `--velocity` | off | Emit STARsolo Velocyto semantics (spliced/unspliced/ambiguous matrices) instead of Gene |
| `--audit-multigene` | off | Print the multi-gene ambiguity audit (the EM upside bound) instead of emitting a matrix |
| `--solo-strand <STRAND>` | `forward` | STARsolo-compatible assignment strand: `forward`, `reverse`, or `unstranded`; 10x 5′ Gene expression uses `reverse` |
| `--from-bam` | off | Interpret `input` as an ingest BAM and extract rows on the fly |
| `--from-molecule-bam` | off | Read a sequence-free BAM produced by `export-molecule-bam`, including its opaque UMI-class IDs and one-mismatch edges |
| `--whitelist <WHITELIST>` | — | Barcode whitelist; required with `--from-bam` |
| `--locus-gap <GAP>` | `2000` | Locus grouping gap (bp); must match the ingest setting |
| `--eager` | off | Decode the entire archive before replay instead of streaming bounded chunk batches; results are identical |
| `--report-format <FORMAT>` | — | Opt in to a versioned `text`, `tsv`, or `json` operation report |
| `--report-output <PATH>` | stdout | Atomically publish the report without replacing an existing file; requires `--report-format` |

## Semantics

The replay implements the relevant STARsolo assignment and collapse rules:
junction concordance against exon structure (`alignToTranscript`), stranded
assignment, per-(cell, gene) aggregation, best-gene filtering, multi-gene UMI
filtering, and greedy 1-mismatch tie-merging over the stored UMI adjacency
graph. `--velocity` additionally ports the Velocyto transcript-set
intersection with its flank-tolerance classifier, applying the same UMI
correction map the Gene collapse produces.

For repeated replay, compile the GTF once with
[`aie compile-annotation`](/gravlax/cli/compile-annotation/) and pass the
resulting `.aic` to the unchanged `--gtf` option. This only removes repeated
GTF parsing and overlap-index construction; matrices are byte-identical.

## Equivalent archive and BAM input

`--from-bam` and archive replay share one row abstraction and produce
byte-identical matrices after the same molecule reduction.

This input equivalence is not a claim of byte identity with a fresh
STARsolo run. Fresh annotation-aware processing can change genome alignments;
the evaluated end-to-end deviation is 0.22–0.45% of UMI mass for the human
3′ datasets and 0.75% in one mouse 5′ v2 dataset for Gene replay, with separate,
larger component-wise values reported for velocity.

```sh
# Compare archive-sourced and BAM-sourced matrices:
aie replay-rows sample.aie      --gtf v49.gtf --barcodes bc.tsv --out-dir a/
aie replay-rows ingest.bam --from-bam --whitelist wl.txt \
                                --gtf v49.gtf --barcodes bc.tsv --out-dir b/
diff a/matrix.mtx b/matrix.mtx
```

## Uniform report and MEX metadata

When `--report-format` is present, replay emits
`gravlax.archive.replay-report.v1`. The typed summary binds the exact archive
(or BAM), the exact annotation snapshot parsed for assignment, the barcode
list, all replay parameters, molecule/assignment/UMI totals, and a manifest of
the MEX files. Matrix, feature, and barcode identities are computed during the
existing write pass; reporting neither clones nor resorts matrix entries.

The scientific `matrix.mtx`, `features.tsv`, and `barcodes.tsv` bytes are
identical with and without reporting. Reporting additionally publishes
`metadata.json` last in the MEX directory using the shared
`gravlax.result-envelope.v1` envelope and
`gravlax.replay.mex-artifact.v1` schema. Gene replay exposes the standard MEX
manifest fields directly, so the shared `read_mex()` path can validate and
open it. Velocity metadata names `spliced.mtx` as its canonical matrix view and
also lists all three component matrices plus `entries.complete.tsv`. The
legacy velocity files retain STARsolo's shared-coordinate layout, including
explicit zero values in an individual component, so the strict single-matrix
reader does not yet claim full velocity-bundle support. The multi-gene audit
also has a typed summary; its legacy human output remains unchanged when
reporting is omitted.

The existing MEX writer creates or updates component files directly in
`--out-dir`. Therefore `metadata.json` is a completion marker, not an atomic
transaction for the whole directory. Report and metadata destinations are
preflighted before replay; each metadata/report file is itself atomically
installed without replacement. For a fresh, unambiguous result, use a new
output directory.
