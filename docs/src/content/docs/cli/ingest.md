---
title: Ingest setup and preflight
description: Generate an annotation-free STAR recipe and validate BAM/whitelist inputs.
---

The `aie ingest` command prepares inputs for the one-time archive build. It
does not build an archive itself; the checked inputs are passed unchanged to
[`aie ingest-archive`](/gravlax/cli/ingest-archive/).

## Chemistry-specific STAR recipe

Print an explicit recipe instead of adapting a generic alignment command by
hand:

```sh
aie ingest recipe --chemistry 10x-3p-v3
```

Supported choices are `10x-3p-v2` (10 bp UMI) and `10x-3p-v3` (12 bp UMI,
including v3.1). Use `--genome-dir`, `--read1`, `--read2`, `--whitelist`,
`--out-prefix`, and `--threads` to substitute real paths. `--plain-fastq`
omits the decompression command.

The recipe is deliberately annotation-free: it uses the junction-only
STARsolo feature, retains secondary alignments, writes `CR`, `CY`, `UR`, and
`NH`, and never adds a GTF. The selected STAR genome directory must itself
have been built without a GTF or annotation-derived splice junctions; the
recipe states this because a command cannot prove it from the directory name.

| `ingest recipe` option | Default | Description |
|---|---|---|
| `--chemistry <NAME>` | required | `10x-3p-v2` (10-base UMI) or `10x-3p-v3` (12-base UMI, including v3.1) |
| `--genome-dir <DIR>` | `star-index-nogtf` | STAR genome directory built without annotation-derived junctions |
| `--read1 <FASTQ>` | `sample_R1.fastq.gz` | Barcode/UMI read |
| `--read2 <FASTQ>` | `sample_R2.fastq.gz` | cDNA read |
| `--whitelist <FILE>` | chemistry-specific filename | Barcode whitelist passed to STAR |
| `--out-prefix <PREFIX>` | `align/` | STAR output prefix |
| `--threads <N>` | `24` | STAR worker threads |
| `--plain-fastq` | off | Treat inputs as uncompressed FASTQ and omit `zcat` |

## Full preflight

Before a long ingest, scan the complete BAM and validate the exact whitelist:

```sh
aie ingest check align/Aligned.sortedByCoord.out.bam \
  --whitelist 3M-february-2018.txt \
  --chemistry 10x-3p-v3
```

The check reads every BAM record. It verifies the reference dictionary and
actual coordinate order; required raw `CR`, `UR`, and `NH` tags; barcode and
UMI lengths; supported sequence alphabets; barcode qualities; and secondary
alignment retention when multimappers are observed. It also rejects malformed
whitelist lines that the ingest loader would otherwise ignore. When a dataset
contains no multimappers, secondary retention cannot be demonstrated and is
reported as a warning rather than silently claimed.

Use `--format json` for a stable `gravlax.ingest.preflight.v1` report.
Failures return a nonzero status while leaving the JSON report parseable.
`--strict` also treats warnings as unsuccessful.

| `ingest check` argument or option | Default | Description |
|---|---|---|
| `<BAM>` | required | Coordinate-sorted annotation-free BAM intended for `ingest-archive` |
| `--whitelist <FILE>` | required | Exact barcode whitelist intended for ingest |
| `--chemistry <NAME>` | infer UMI length | Require the observed UMI length to match `10x-3p-v2` or `10x-3p-v3` |
| `--strict` | off | Return an unsuccessful status for warnings as well as errors |
| `--format <FORMAT>` | `text` | `text` or `json` report |
