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

The default recipe is annotation-free: it uses the junction-only
STARsolo feature, retains secondary alignments, writes `CR`, `CY`, `UR`, and
`NH`, and never adds a GTF. Junction seeding is available as an explicit
option (below) and is recorded in archive provenance. The selected STAR genome directory must itself
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
| `--junction-seed <FILE>` | off | Insert a fixed splice-junction list (STAR `--sjdbFileChrStartEnd` format) at mapping time; see below |
| `--sjdb-overhang <N>` | chemistry-derived (`90` for 10x 3' v3/v3.1, `97` for 10x 3' v2) | STAR `--sjdbOverhang` emitted with `--junction-seed`: the cDNA read length minus one. Override only for non-standard read lengths |
| `--one-pass` | off | Omit `--twopassMode Basic`, so junctions come only from the aligner's own detection and any seed |

## Optional junction seeding and one-pass alignment

Both options are off by default. The default recipe inserts no junctions and
runs per-library two-pass discovery.

Splice-junction sets change far less between annotation releases than
transcript sets do (across GENCODE v32 to v49, 99% of v32's junctions persist,
and junctions added since v32 carry well under 1% of observed junction reads).
Seeding the alignment with a junction list therefore leaves the archive
annotation-independent with respect to gene models while letting the aligner
place reads across known junctions it might otherwise miss. A seed file can be
derived from any GTF or compiled annotation:

```sh
aie ingest junctions --gtf gencode.v32.annotation.gtf --out v32.junctions.tab
aie ingest recipe --chemistry 10x-3p-v3 --junction-seed v32.junctions.tab
```

`--sjdbOverhang` follows STAR's rule of the cDNA read length minus one and is
derived from `--chemistry`: 90 for 10x 3' v3/v3.1 (91 bp cDNA reads) and 97 for
10x 3' v2 (98 bp cDNA reads). Pass `--sjdb-overhang <N>` to override it when a
library was sequenced at a non-standard read length.

The seed file has one junction per line: chromosome, 1-based inclusive intron
start and end, and strand. With `--one-pass`, only the seed and the aligner's
own detection supply junctions; the recipe then declares the seed as a
`frozen-catalogue` at ingest. Without `--one-pass`, the recipe declares the
pass-1 table as `per-library-two-pass` and records the seed as the alignment
annotation. In every case the printed `aie ingest-archive` line carries the
matching `--junction-discovery`, `--junction-catalogue`, and
`--alignment-annotation` flags, so the archive records exactly how junctions
were supplied, with the seed file's digest.

| `ingest junctions` option | Default | Description |
|---|---|---|
| `--gtf <PATH>` | required | Uncompressed GTF or compiled `.aic` annotation |
| `--out <FILE>` | required | Seed file to write; an existing file is not overwritten |

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
