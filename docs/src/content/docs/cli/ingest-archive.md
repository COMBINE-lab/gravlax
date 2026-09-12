---
title: aie ingest-archive
description: Build the .aie index from a tagged, coordinate-sorted ingest BAM.
---

Build the `.aie` archive from a tagged, coordinate-sorted ingest BAM. This is the
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
  --junction-discovery per-library-two-pass \
  --junction-catalogue align/_STARpass1/SJ.out.tab \
  --alignment-chemistry 10x-3p-v3 \
  --terminal-tails \
  --out sample.aie
```

## Arguments

| Argument | Description |
|---|---|
| `<BAM>` | Tagged ingest BAM (`CR`/`UR`/`CY` tags, secondaries included, coordinate-sorted) |

## Options

| Option | Default | Description |
|---|---|---|
| `--whitelist <WHITELIST>` | required | 10x barcode whitelist (one 16 bp barcode per line) |
| `--out <OUT>` | required | Output `.aie` path |
| `--locus-gap <GAP>` | `2000` | Single-linkage gap (bp) for grouping reads into loci |
| `--zstd-level <LEVEL>` | `19` | zstd level for chunk streams (level 12 trades ~6% size for ~25% faster ingest; query latency is flat across levels) |
| `--chunk-mb <MB>` | `4` | Genomic chunk size in megabases — the size-vs-random-access granularity knob; the ingest report prints total size per setting |
| `--genome <FASTA>` | — | Bind a reference for sequence-consulting queries by recording its exact identity and normalized per-contig BLAKE3 signature; its relationship to the alignment is caller-declared, not inferred |
| `--terminal-tails` | off | Apply the frozen forward-stranded 10x 3′ cDNA rule (forward trailing A / reverse leading T) and retain its sparse, sequence-free terminal-tail observable; this is not a chemistry-generic detector. See [The `.aie` format](/gravlax/format/#sparse-terminal-tail-evidence) |
| `--junction-discovery <MODE>` | `unspecified` | Record the alignment's junction source: `one-pass`, `per-library-two-pass`, `frozen-catalogue`, or `unspecified`; this is an explicit declaration, not an inference from `@PG` |
| `--junction-catalogue <PATH>` | — | Exact STAR-style junction table supplied to pass 2; for `per-library-two-pass`, this is the pass-1 output (normally `_STARpass1/SJ.out.tab` for STAR Basic), not pass 2's final `SJ.out.tab`. Required by `per-library-two-pass` and `frozen-catalogue`, forbidden for the other modes, and embedded with its digest and parsed data-row count. Gravlax verifies its bytes; its role in alignment is caller-declared |
| `--alignment-annotation <PATH>` | — | Hash an annotation supplied to the aligner or its index and record its path and caller-declared role |
| `--alignment-index-identity <TEXT>` | — | Record a caller-supplied content identity or reproducible locator for the aligner index; directories are not hashed implicitly |
| `--alignment-input <PATH>` | — | Hash one source-read or other aligner input and retain its locator; repeat in the aligner's original input order |
| `--alignment-log <PATH>` | — | Hash an aligner log that records resolved defaults and retain its locator |
| `--alignment-chemistry <TEXT>` | — | Record the caller-declared library chemistry used for alignment and strand interpretation |
| `--geometry-fidelity` | off | Retain every distinct accepted unique-read geometry with its exact read multiplicity, instead of two coordinate-extreme representatives per junction chain. See [When to use geometry fidelity](#when-to-use-geometry-fidelity) |
| `--report-format <FORMAT>` | — | Opt in to a versioned `text`, `tsv`, or `json` operation report |
| `--report-output <PATH>` | stdout | Atomically publish the report without replacing an existing file; requires `--report-format` |

## Notes

- The archive representation leaves downstream annotation choice open, but
  alignment provenance still matters. For discovery-oriented use, an
  annotation-free alignment can learn junctions from the same library with a
  two-pass run; declare that run with `per-library-two-pass` and preserve its
  pass-1 junction table. If an annotation was supplied to alignment or index
  construction, record it with `--alignment-annotation` so later comparisons
  can distinguish retained evidence from evidence the aligner never emitted.
  Secondary alignments are evidence — do not filter them. `aie ingest check`
  scans the BAM and whitelist before this longer build.
- `per-library-two-pass` means the embedded catalogue came from pass 1 of this
  library. `frozen-catalogue` means an already fixed external catalogue was
  reused, for example a junction seed written by `aie ingest junctions` and
  aligned with `aie ingest recipe --junction-seed ... --one-pass`; the recipe
  prints the matching declaration. Gravlax verifies the supplied file's current bytes but relies on the
  caller for that relationship; it does not infer either mode from the BAM
  header.
- Barcode correction runs here, against the whitelist, using an
  annotation-independent corrector. It is fixed in the index;
  re-correction under a different whitelist is the one replay the index
  cannot perform.
- The ingest report prints per-section byte accounting, which is also
  embedded in the archive footer (`aie dev debug` reprints it later).
- New archives use logical schema `gravlax.molecular-evidence.v2` inside the
  rooted v2 container. The core molecule streams remain compatible with older
  readers. New readers distinguish an absent optional capability from a
  capability that was evaluated and produced zero events.

:::caution[Review provenance before sharing]
The archive records the path strings supplied to provenance options and the
BAM header's `@PG` command lines. They may contain user names, directories,
sample identifiers, or command arguments. The declared junction catalogue is
embedded exactly; source reads, annotations, and logs are represented by
identities and locators rather than copied into the archive.
:::

## When to use geometry fidelity

By default a UMI class stores each junction chain as a read count plus two
coordinate-extreme representatives. `--geometry-fidelity` (0.2.2) instead keeps
every distinct unique-read geometry once, with its exact multiplicity, inside
the same cell/UMI/locus record. Gene counting does not need it. Use it when a
downstream analysis depends on every read in a UMI class — the ambiguous
component of RNA velocity, or per-read geometry queries — and the storage
premium is acceptable.

Measured on 10x PBMC 5k (383.9M reads, GENCODE v49, 5,038 called cells,
gravlax 0.2.2; ingest took 7.4 minutes at 24 threads in both modes):

| Quantity | Default | `--geometry-fidelity` |
|---|---|---|
| Archive size | 544.0 MB (11.3 bits/read) | 692.7 MB (14.4 bits/read) |
| Gene count-matrix deviation from STARsolo | 0.307% | 0.253% |
| RNA-velocity spliced deviation | 2.67% | 2.35% |
| RNA-velocity unspliced deviation | 0.86% | 0.85% |
| RNA-velocity ambiguous deviation | 6.13% | 4.67% |

The archive grows by 27%. Gene counts barely move, so a workflow that only
replays Gene or GeneFull matrices gains little for that cost. The largest
change is the ambiguous velocity component, and even there representative
reduction accounts for only about a quarter of the 6.13% deviation; the rest
comes from alignment differences and no archive setting recovers it.

Older readers reject the `distinct-unique-geometries-v1` provenance rule, so an
archive built with this switch is not readable by pre-0.2.0 tools. The other
experimental ingest switches are described on
[The `.aie` format](/gravlax/format/#experimental-archive-improvements).

## Uniform operation report

Omitting `--report-format` preserves the established stdout report exactly.
With reporting enabled, stdout contains only the uniform report (or is empty
when `--report-output` is used). The `gravlax.archive.ingest-report.v1`
summary records molecule, class, dictionary, and byte totals; exact BLAKE3
identities for the BAM, whitelist, optional FASTA, and completed archive; and
all build parameters. It also reports the complete parsed alignment manifest,
terminal-tail availability and cardinalities, and every alignment/tail option.
Its `sections` table streams one row per archive section without buffering a
second copy of the section data.

The report path is checked before BAM extraction starts and is installed
atomically without replacing an existing path. It is deliberately separate
from `--out`: the archive remains the primary artifact.
