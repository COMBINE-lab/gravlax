---
title: aie ingest-archive
description: Build the .aie index from the annotation-free ingest BAM.
---

Build the `.aie` archive from the annotation-free ingest BAM. This is the
one-time indexing step: after it completes, every replay and query runs from
the archive alone.

## Usage

```sh
aie ingest-archive [OPTIONS] --whitelist <WHITELIST> --out <OUT> <BAM>
```

```sh
aie ingest check align/Aligned.sortedByCoord.out.bam \
  --whitelist 3M-february-2018.txt \
  --chemistry 10x-3p-v3

aie ingest-archive align/Aligned.sortedByCoord.out.bam \
  --whitelist 3M-february-2018.txt \
  --out sample.aie
```

## Arguments

| Argument | Description |
|---|---|
| `<BAM>` | Annotation-free ingest BAM (`CR`/`UR`/`CY` tags, secondaries included, coordinate-sorted) |

## Options

| Option | Default | Description |
|---|---|---|
| `--whitelist <WHITELIST>` | required | 10x barcode whitelist (one 16 bp barcode per line) |
| `--out <OUT>` | required | Output `.aie` path |
| `--locus-gap <GAP>` | `2000` | Single-linkage gap (bp) for grouping reads into loci |
| `--zstd-level <LEVEL>` | `19` | zstd level for chunk streams (level 12 trades ~6% size for ~25% faster ingest; query latency is flat across levels) |
| `--chunk-mb <MB>` | `4` | Genomic chunk size in megabases — the size-vs-random-access granularity knob; the ingest report prints total size per setting |
| `--genome <FASTA>` | — | Stamp the alignment reference's normalized per-contig BLAKE3 signature during ingest |
| `--report-format <FORMAT>` | — | Opt in to a versioned `text`, `tsv`, or `json` operation report |
| `--report-output <PATH>` | stdout | Atomically publish the report without replacing an existing file; requires `--report-format` |

## Notes

- The input BAM must come from an **annotation-free** alignment (see the
  [quick start](/gravlax/quickstart/) or generate an explicit command with
  [`aie ingest recipe`](/gravlax/cli/ingest/)). Secondary alignments are
  evidence — do not filter them. `aie ingest check` scans the BAM and
  whitelist before this longer build.
- Barcode correction runs here, against the whitelist, using an
  annotation-independent corrector. It is fixed in the index;
  re-correction under a different whitelist is the one replay the index
  cannot perform.
- The ingest report prints per-section byte accounting, which is also
  embedded in the archive footer (`aie dev debug` reprints it later).

## Uniform operation report

Omitting `--report-format` preserves the established stdout report exactly.
With reporting enabled, stdout contains only the uniform report (or is empty
when `--report-output` is used). The `gravlax.archive.ingest-report.v1`
summary records molecule, class, dictionary, and byte totals; exact BLAKE3
identities for the BAM, whitelist, optional FASTA, and completed archive; and
all build parameters. Its `sections` table streams one row per archive section
without buffering a second copy of the section data.

The report path is checked before BAM extraction starts and is installed
atomically without replacing an existing path. It is deliberately separate
from `--out`: the archive remains the primary artifact.
