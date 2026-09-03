---
title: The .aie format
description: What the evidence index stores, how it is laid out on disk, and the conditions for replay.
---

An `.aie` file is a single seekable container holding the annotation-free
molecular evidence of one sample: every molecule's placement geometry, UMI
equivalence structure, and multimapping structure, plus the indexes that make
it queryable. This page describes the format at the level a user or tool
author needs.

## Design principles

1. **Retained evidence is explicit.** The molecule payload
   summarizes read signatures with up to two span-extreme representatives per
   junction chain. Archive and BAM replay are byte-identical after this shared
   reduction; comparison with fresh STARsolo measures both the reduction and
   changes caused by annotation-aware alignment.
2. **Compression comes from structure the data actually has.** Some of the
   same structures also support indexed queries:
   - geometry → a **shape dictionary** of block-length/junction-offset
     vectors;
   - junctions → a **junction catalogue** with inverted postings (the
     junction-query index);
   - UMI values → an **adjacency graph** of equivalence classes and
     1-mismatch edges (19× smaller than storing the values, and replay
     accuracy *improves*, because merge decisions defer to replay time when
     the gene context is known);
   - multimapping → a **paralog pattern dictionary** (anchor-relative offset
     vectors, interned at 7.9× sharing; the pattern is also the equivalence
     class EM consumes);
   - genomic chunks and cell postings support region and APA scans without a
     whole-archive decode.
3. **Access is genome-major.** Replay is a sequential genome-ordered scan;
   junction/region/discovery queries are range reads; per-cell access goes
   through postings, never a scan.
4. **Evolution is explicit.** A versioned container and named sections allow
   optional evidence to be added. The current core is intentionally lossy:
   barcode correction is fixed at ingest, sequence and qualities are absent, and most
   chains retain two representatives rather than all reads.

## What is stored per molecule

The core representation stores *quotients* rather than identities —
annotation-dependent processing never consumes raw values, only relations
among them:

- **Junction chains, two representatives each.** Unique-mapping reads group
  into loci (single-linkage, 2 kb gap) and, within a locus, into junction
  chains keyed by absolute junction coordinates. Each chain stores its two
  span-extreme reads — the most-contained and most-extended — plus the chain
  read count. The extremes bracket the containment behaviour of every read in
  between. Chains are position-sorted and the span-minimum representative
  comes first, so the anchor offset of the first representative is implied
  and never stored.
- **UMI classes and edges, never values.** UMIs are stored as global
  (cell, value) equivalence classes plus cell-scoped 1-mismatch adjacency
  edges. Every collapse policy — 1MM tie-merging, multi-gene UMI filtering,
  velocity's cell-scoped intersection — is a function of these relations,
  applied at replay with the target annotation in hand.
- **Cell-of-class, once per class.** The cell barcode is a pure function of
  the UMI class, so cell identity is stored once per class in a global
  section, not per molecule.
- **Paralog patterns, interned.** A multimapper's alternative placements are
  stored as an anchor-relative difference set — offsets, strand flips, and a
  same-shape flag — interned in a dictionary shared across molecules and
  cells.

## Layout

```text
header         magic "AIE0", container version
dictionaries   cells (packed corrected barcodes, frequency-ordered)
               shapes (block-length / junction-offset vectors)
               junction catalogue (chrom, donor, acceptor)
               paralog patterns (offset/flip vectors relative to anchor)
               rANS tables (global static tables for the memoryless streams)
cell-of-class  one cell id per UMI class, in 65,536-class blocks,
               each block its own compressed frame with per-block codec choice
chunks         per 4 Mb genomic bin: ten columnar streams, independently
               compressed (anchor, class, layout, weight, rep.pos, rep.shape,
               mm.pos, mm.shape, mm.pattern, mm.weight)
graph          cell-scoped 1-mismatch edges among global UMI classes
indexes        genomic range → chunk · junction id → chunk postings
               (with genome-wide support totals) · cell id → chunk postings
footer         section directory (name, offset, raw length, compressed length)
```

### Chunks are the access unit

Genomic 4 Mb chunks are self-contained given the dictionaries: every stream is
zstd-compressed per chunk (level 19 by default), with a static-table rANS
stage on the five streams whose values are memoryless — the fitted tables live
in a dictionary section so chunk decode stays independent. Gene-scale queries
decode one or two chunks.

### Lazy open

Opening an index reads only the header, chromosome table, and chunk directory;
dictionaries load on first use and cell-of-class blocks decompress
individually as classes are touched. Open cost is **~10 ms, independent of
archive size**. Full replays decode chunks and cell-of-class blocks in
parallel.

### Compatibility fields

`meta` records the stream layout (`chunk_streams`), cell-of-class block size
(`coc_block`), and codec generation (`codec`). Readers **must refuse** values
they do not understand rather than misread the file — the format versions
forward, never silently.

## Size and fidelity

Across four 10x datasets (blood, tumor, brain nuclei; 66.6M–383.9M reads),
the index costs **11–18 bits per input read**:

| artifact | typical size (66.6M-read PBMC) | vs index |
|---|---|---|
| FASTQ (R1+R2) | 5.0 GB | 32–49× |
| BAM (all tags) | 4.2 GB | 26–37× |
| CRAM 3.1 (archive mode, all tags) | 1.40 GB | 9.0–12.7× |
| **evidence index (`.aie`)** | **111.1 MB** | — |

Section sizes are available from the container directory, and `aie dev debug`
reconstructs detailed per-stream accounting for an archive.

## Observable format properties

- **Round-trip**: write → read → identical molecule set, bit for bit.
- **Equivalent replay inputs**: replay-from-archive and replay-from-BAM produce
  **byte-identical matrices** after the shared molecule reduction.
- **Equivalent execution modes**: bounded archive replay and eager archive replay
  produce the same result,
  including UMI classes whose evidence crosses physical chunk boundaries.
- **Fresh-pipeline comparison**: replay from fixed annotation-free alignments
  is distinct from running a fresh annotation-aware alignment. Dataset-specific
  differences between those workflows are reported separately.
