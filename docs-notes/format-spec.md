# `.aie` molecular-evidence archive format

Gravlax writes a seekable, authenticated `.aie` v2 container. The archive stores corrected,
molecule-resolved genome-alignment evidence without gene, transcript, or exon annotations. A
compatible annotation can therefore be supplied later for quantification, discovery, or indexed
queries without repeating alignment and barcode/UMI correction.

This document distinguishes two versioned layers:

- **container version 2** defines section layout, content identity, and integrity; and
- the `meta` payload guards the current molecular codec (`rans2`, ten chunk streams, and blocked
  cell-of-class coding).

Readers also accept the earlier seekable v1 container. Writers emit v2.

## Evidence model

The archive records the post-correction molecular abstraction used by Gravlax replay:

- an ordered chromosome dictionary and optional digest of the aligned reference genome;
- corrected cell barcodes;
- global UMI classes and their one-mismatch adjacency edges;
- genome-ordered molecule records containing chromosome, strand, anchor, UMI class, and weights;
- one or two representative placements for each observed alignment chain, expressed as a genomic
  start and a reusable block shape; and
- multimapper placements expressed with a dictionary of relative alternative-placement patterns.

A shape is an ordered vector of aligned blocks `(offset, length)` relative to a placement start.
Adjacent blocks define exact junction boundaries. When members of one chain have different
footprints, the stored representatives are the most-contained and most-extended placements; their
chain weight records the number of supporting children. Multimapper patterns retain alternative
chromosome, displacement, orientation, and shape information relative to the anchor placement.

The archive does not store gene assignments or transcript identifiers. It also does not claim to
reproduce a new annotation-aware alignment: replay interprets the fixed archived alignments under a
new annotation. Archive-backed and BAM-backed Gravlax replay are exact with respect to the same
post-correction molecular evidence, while comparison with a fresh aligner is a separate fidelity
measurement.

## Authenticated v2 container

All integer fields in the outer container are little-endian. Each payload is an independent zstd
frame.

```text
fixed header
  magic                 4 bytes: "AIE0"
  container version     u32: 2

section area, repeated in physical order
  name length           u8, nonzero
  UTF-8 name             name-length bytes
  raw length            u64
  compressed length     u64
  compressed payload    compressed-length bytes

section terminator      u8: 0

authenticated directory
  section count         u32
  for each physical section:
    name length + name
    physical offset     u64
    raw length          u64
    compressed length   u64
    compressed BLAKE3   32 bytes

footer                  44 bytes
  directory offset      u64
  directory root        32 bytes
  magic                 4 bytes: "AIED"
```

The directory root is a domain-separated BLAKE3 commitment to `AIE0`, container version 2, the
directory offset, and the exact directory bytes. It therefore commits section names, ordering,
physical offsets, raw and compressed lengths, and every compressed-payload digest. This is a
content commitment, not a publisher signature.

Opening a v2 archive verifies the footer and root, requires a contiguous canonical section layout,
and checks every directory entry against its inline section header. It does not read unrelated
payloads. When a section is selected, the reader hashes its exact compressed bytes, compares that
digest with the authenticated directory entry, and then performs bounded decompression to exactly
the declared raw length. `aie inspect-archive --verify-content` selects and verifies every payload.

The native archive identity is the directory root under the scheme `aie-directory-root-v2`.
Gravlax also derives an `aie-encoded-sections-v1` identity from the ordered section names, lengths,
and compressed-payload digests. This second identity intentionally excludes container offsets and
footer representation, so byte-preserving v1-to-v2 sealing retains the same encoded-evidence
identity.

## Molecular sections

| Section | Contents |
|---|---|
| `meta` | JSON cardinalities and hard layout guards: molecule, edge, cell, shape, pattern, and class counts; chunk width; chunk-stream count; cell-of-class block width; codec name; optional reference-genome signature. |
| `chroms` | Newline-delimited chromosome names in archive identifier order. |
| `cells` | Corrected barcodes packed as little-endian `u32` values. |
| `shapes` | Varint-coded block vectors `(offset, length)`. |
| `patterns` | Varint-coded multimapper alternative-placement patterns; a flag elides shape IDs when the alternative uses the anchor shape. |
| `rans.tables` | Six archive-wide static rANS tables: five molecule value streams and the cell-of-class alternative codec. |
| `cN` | Molecule chunk `N`, normally a 4-Mb anchor-coordinate bin on one chromosome. |
| `coc.N` | Cell ID for each UMI class in a 65,536-class block. A one-byte tag selects delta-varints or rANS over absolute IDs for that block. |
| `index.chunks` | Chromosome, bin start, molecule count, class base, maximum anchor, and distinct-cell count for every molecule chunk. |
| `index.junctions` | Coordinate-sorted junction catalogue `(chromosome, donor, acceptor)`. |
| `index.jpost` | Per-junction supporting-child total and delta-coded molecule-chunk postings. |
| `index.cellpost` | Delta-coded molecule-chunk postings for each cell. |
| `edges` | Delta-coded one-mismatch edges between global UMI classes. |

Genomic coordinates are zero-based. Region intervals are half-open; junction donor and acceptor
values are the boundaries between adjacent aligned blocks.

## Molecule chunks

Molecules are sorted by chromosome and anchor, then partitioned by `meta.chunk_bp` (4,000,000 by
default). A `cN` payload contains ten length-prefixed column streams:

```text
anchor
class
layout
chain weight / representative count
representative position
representative shape
multimapper position
multimapper shape
multimapper pattern
multimapper weight
```

Anchor positions are nonnegative deltas within a chunk. A class token either introduces the next
global UMI class or refers backward to an existing class; the class's cell is stored once in the
appropriate `coc.N` block. `layout` records strand, chain count, and multimapper count. The first
representative of a molecule with a chain is its anchor and its position is therefore implicit.
Other positions are anchor-relative. Shape and pattern identifiers reference the global
dictionaries.

The class, chain-weight, representative-position, multimapper-position, and multimapper-weight
streams use the archive-wide static rANS tables. The remaining structured streams use varints.
The complete ten-stream chunk is then compressed as one zstd section, making a genomic chunk the
unit of random access and bounded decoding.

The `meta.codec`, `meta.chunk_streams`, and `meta.coc_block` fields are mandatory decoding guards.
A reader refuses an unfamiliar combination instead of interpreting bytes under a different
molecular layout.

## Indexed access

Ordinary replay walks `index.chunks` in genome order and releases bounded batches after global
UMI-class aggregation. A region query selects chunks by chromosome and anchor overlap. A junction
query locates the coordinate in `index.junctions`, uses `index.jpost` to find candidate chunks, and
checks placement geometry against `shapes`. Cell and group scopes are applied while reducing those
selected chunks, and only touched `coc.N` blocks are decompressed. Multimapper patterns and other
dictionaries are loaded lazily when the selected operation requires them. The archive also carries
`index.cellpost` for direct cell-to-chunk routing.

Indexes route work but do not store query answers. Counts are recomputed from molecule records,
UMI classes, and adjacency edges. The same principle applies to optional federated shape routes in
`.aicollection`: they narrow the shape candidates for exact point and junction-set predicates but
never supply a count.

## Legacy v1 archives

Seekable v1 has the same named compressed molecular sections and a tail directory, but its
directory has neither per-payload digests nor a root commitment. The reader validates directory
entries against inline headers, but persistent content identity requires a complete file scan and
ordinary selected reads cannot authenticate an uncommitted payload.

```sh
# Report native and cross-container identities.
aie inspect-archive legacy.aie --json

# Copy compressed sections exactly into a new rooted v2 container.
aie seal-archive legacy.aie --out rooted.aie --json

# Authenticate the directory and verify every compressed payload.
aie inspect-archive rooted.aie --verify-content --json
```

`seal-archive` accepts only v1 input and refuses to overwrite its destination. It reads and bounded-
decodes every source frame, copies the compressed bytes without recompression, writes a canonical
v2 directory/footer, verifies every output payload, and requires the source and output
encoded-section identities to match. The source remains unchanged.

For v2, ordinary `inspect-archive` verifies the directory root and reports zero payload bytes read
for content identity; selected payloads are verified when used. For v1, inspection performs a
complete scan to derive its full-file and encoded-section identities. `--verify-content` is the
explicit whole-payload audit for either format.

## Validation and limits

Container readers reject unknown versions, empty or duplicate section names, noncanonical or
overlapping physical layouts, directory/header disagreement, truncation, trailing layout bytes,
unsafe section sizes or compression ratios, arithmetic overflow, root mismatch, and selected
payload digest mismatch. The v2 implementation limits the directory to 1,000,000 sections, each
raw section to 512 MiB, and each compressed section to 528 MiB; decompression is capped at the
declared raw length plus one byte.

Molecular decoders independently validate dictionary identifiers, stream cardinalities, varint and
rANS termination, coordinate arithmetic, class back-references, cell-of-class block lengths,
junction and chunk postings, and the layout fields in `meta`. Integrity failures are errors; they do
not produce partial biological results.
