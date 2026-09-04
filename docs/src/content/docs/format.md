---
title: The .aie format
description: What the evidence index stores, how it is laid out on disk, and the conditions for replay.
---

An `.aie` file is a single seekable container holding molecular evidence from
one fixed tagged alignment: every molecule's placement geometry, UMI
equivalence structure, and multimapping structure, plus the indexes that make
it queryable. The recommended workflow uses annotation-free alignment;
root-bound provenance records a caller declaration when the supplied alignment
used an annotation. This page describes the format at the level a user or tool author needs.

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

## Container version and evidence schema

Two versions describe different layers of an archive:

- **Container v2** is the outer `AIE0` file layout. Its footer commits to the
  ordered section directory and every compressed section digest. Adding a named
  section therefore does not require a new container grammar.
- **`gravlax.molecular-evidence.v2`** is the logical schema written by current
  `ingest-archive`. It adds root-bound alignment provenance and can add the
  optional terminal-tail observable while leaving the ten core molecule streams
  unchanged.

A rooted archive can predate the logical v2 schema. `seal-archive` authenticates
the bytes of an older archive but does not invent provenance or evidence that was
not captured at ingest. Readers treat a missing capability as **unavailable**,
not as an observed count of zero. Unknown logical schemas and partial capability
section sets are errors.

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
provenance     canonical alignment manifest; optional exact junction catalogue
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
optional tail  sparse tail route index · one event list per selected core chunk
footer         section directory (name, offset, raw length, compressed length)
```

### Alignment provenance

Logical v2 archives contain `alignment.provenance`, canonical compact JSON with
schema `gravlax.alignment-provenance.v1`. It records:

- identities computed from the same held BAM and whitelist files that ingest
  consumed, plus the exact genome FASTA identity and normalized genome signature
  when `--genome` is supplied at ingest (the relationship to alignment remains a
  caller declaration);
- BAM `@PG` records in header order, marked as observations from the verified BAM;
- an explicit caller declaration of junction discovery as `one-pass`,
  `per-library-two-pass`, `frozen-catalogue`, or `unspecified`;
- the ordered identities and locators supplied with repeatable
  `--alignment-input`, and optional annotation, aligner log, index identity, and
  chemistry declarations;
- the Gravlax version and every effective molecule-reduction, chunking, codec,
  compression, and optional tail-extraction parameter.

The physical-layout values duplicated in `meta` (`chunk_bp`, `chunk_streams`,
and `codec`) must equal their provenance-manifest counterparts; opening or
inspecting a logical-v2 archive rejects any disagreement. `locus_gap` and
`zstd_level` are root-committed records of the construction settings, but the
finished archive does not independently reveal enough history to verify that
those declared settings were used.

Gravlax never infers junction-discovery mode from an `@PG` command line. A
`verified_from_consumed_bytes` status means Gravlax hashed stable bytes;
`verified_bam_header` means a value was read from that BAM header; and
`declared_by_caller` describes the asserted relationship between a verified file
and the earlier alignment run.

The currently bound reference is recorded separately in
`meta.genome_reference_binding`. Its `bound_by` field distinguishes an
ingest-time binding from a later `stamp-genome` operation, and it records the
verified FASTA identity, normalized signature, and contig-coverage check. A
later stamp never rewrites `alignment.provenance` and does not prove that its
FASTA was the reference used to create the BAM.

For two-pass modes, `alignment.junction-catalogue` stores the exact supplied
catalogue bytes. The manifest records its full-file BLAKE3 digest, byte length,
parsed nonblank/noncomment data-row count, and whether it is the same library's
pass-1 table or an externally frozen catalogue.

:::caution[Provenance can contain sensitive paths]
The manifest retains caller-supplied path strings and BAM `@PG` command lines.
Those values can include user names, directories, sample identifiers, or command
arguments. Inspect them before distributing an archive. Source-read, annotation,
and log contents are not embedded; their identities and locators are. A declared
junction catalogue is embedded exactly.
:::

### Sparse terminal-tail evidence

`--terminal-tails` evaluates the frozen rule
`forward-cdna-terminal-softclip-v1`. A qualifying uniquely mapped (`NH=1`),
primary, mapped, nonsupplementary record must have a corrected cell barcode and
UMI and:

- on the forward genomic strand, a trailing soft clip with at least 6 bases,
  at least 4/5 A, and a terminal A run of at least 4; or
- on the reverse genomic strand, a leading soft clip satisfying the same
  thresholds for T.

Hard clips and padding outside the terminal soft clip do not change which clip
is inspected. The cleavage anchor is the 0-based exclusive aligned end on `+`
and the 0-based inclusive aligned start on `-`. Qualification uses the exact
unsaturated clip counts with overflow-safe integer arithmetic. Duplicate
witnesses are ranked by exact unsaturated tail fraction, then terminal run,
then clip length.

Every qualifying read-level anchor survives the normal span-extreme molecule
reduction, including an anchor observed only on a middle read. Duplicate
witnesses are collapsed globally by `(chromosome, cleavage anchor, strand,
corrected cell, UMI class)`. Each resulting event is attached to a stable
serialized molecule-record ordinal, so terminal predicates can participate in
same-record Boolean queries. When duplicate-key witnesses originated in
different molecule records, the event remains attached to the record containing
the strongest read-level signal selected by the ranking above; exact-signal ties
choose the first serialized record. The signal and its origin record are never
selected independently.

Multi-mapping (`NH>1`) records are excluded because the v1 side section does
not bind a tail observation to one member of a placement pattern. Treating the
BAM-primary placement as an exact cleavage anchor would overstate what the
alignment establishes. Records without an explicit integer `NH=1` tag are also
excluded because unique placement cannot be established from this section.

The retained signal is deliberately bounded: clip length, matching-tail bases,
and terminal run each saturate at 31. Raw clipped sequence, base qualities,
and exact counts above 31 are not stored. “Lossless terminal-tail evidence” means
that every qualifying globally deduplicated cleavage-anchor key selected by this
versioned rule is retained. It does not mean preservation of exact raw counts,
clipped sequence, or base qualities.

The sparse encoding is:

```text
index.tail
  "TAILIDX1" · route_count(varint)
  repeated: chunk_delta · chrom · min_anchor · anchor_span
            · selected_molecules · event_count       (all varints)

tail.cN
  "TAILCHN1" · ordinary_chunk_molecules · selected_molecules · events
  repeated selected molecule:
    local_molecule_ordinal_delta · event_count        (varints)
    repeated event:
      cleavage_anchor_delta_from_molecule             (signed varint)
      signal                                           (little-endian u16)
```

Signal bit 15 is the genomic reverse-strand flag. Bits 10–14, 5–9, and
0–4 hold the three five-bit counts. Route, cardinality, coordinate-envelope,
chromosome, local molecule identity, and molecule/event strand agreement are
validated when the side section is read. `meta.terminal_tail` declares the rule,
unique-mapping alignment scope, thresholds, coordinate convention, total
selected molecules, events, and routed chunks. If that typed declaration is
absent, a tail query fails as unavailable;
only a present declaration with zero counts is a measured zero.

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
