# Scaling UMI-class identifiers beyond 2^32

## Status

Design note. The companion change on this branch makes construction **fail
loudly** at the current limit. Widening the identifier is future work; this note
scopes it so the format change can be made deliberately.

## The current limit

A UMI class is a distinct `(corrected cell, raw UMI value)` pair within one
archive. Its identifier is a `u32` in every place it appears:

- `MolRec.umi_class: u32` and `Mol.class: u32` (the per-molecule field).
- `n_classes: u32` on the in-memory archive struct, and the `classes` value in
  archive metadata, read back through `required_meta_u32`.
- The construction counters: the `class_of` closure in `rows.rs`
  (`*n_classes += 1`) and the encoder counter in `archivecmd.rs`
  (`next_class += 1`).
- The on-disk per-molecule class stream, which stores `0` for a newly
  introduced class and `next_class - umi_class` for a back-reference, decoded
  as `u32`.
- The packed discovery reducer, `(event << 34) | (class << 2) | mask` in a
  `u64`: 30 event bits, 32 class bits, 2 side bits.

So an archive holds at most `2^32 - 1` UMI classes. That is about 4.29 billion
`(cell, UMI)` pairs, well beyond any single sample today. Before this change,
crossing it wrapped silently in release builds and produced an archive that
re-opened without error while aliasing distinct molecules onto shared ids. The
guard on this branch converts that into an error at ingest.

## What a wider identifier has to touch

Widening is not a one-line type change. The identifier is load-bearing in four
places, in rough order of difficulty.

### 1. The packed discovery reducer (the hard constraint)

`collection find-events` and `query events` reduce evidence by packing
`(event, class, mask)` into one integer and doing a sort plus bitwise-OR merge.
The word must hold `event_bits + class_bits + 2`. Today that is
`30 + 32 + 2 = 64`, exactly a `u64`.

The event count is not small. The cohort discovery run in the paper enumerated
about 5.4 million candidate definitions, which needs 23 event bits. A 40-bit
class with a 23-bit event space needs `40 + 23 + 2 = 65` bits, which no longer
fits a `u64`. So a wider class forces one of:

- **Runtime bit allocation with a `u128` fallback.** Compute
  `class_bits = bits(n_classes)` and `event_bits = bits(n_events)` per reduction
  pass. If `class_bits + event_bits + 2 <= 64`, keep the current `u64` fast
  path. Otherwise pack into a `u128`. The `u128` sort and OR-reduce are slower
  and use twice the memory per hit, but they are correct and only pay that cost
  on archives that actually exceed the `u64` budget. This keeps every realistic
  archive (fewer than `2^30` classes) on the fast path unchanged.
- **A struct-sort fallback.** Sort `(event, class)` tuples with a `u32`/`u64`
  mask alongside, instead of one packed word. Fully general, simplest to reason
  about, slowest. Reasonable as the >`u64` path if `u128` packing is not wanted.

Recommendation: runtime bit allocation, `u64` fast path preserved, `u128`
fallback. Gate it on a single computed `class_bits` so the two paths cannot
diverge.

### 2. On-disk encoding and metadata

- Promote the per-molecule class field to `u64` in memory (a `ClassId(u64)`
  newtype keeps the intent visible and prevents accidental narrowing).
- The class stream is already delta and back-reference coded, so it extends to
  wider values by widening the varint/rANS path to `u64` deltas. No structural
  change to the stream, only its element width.
- Write `classes` in metadata as a `u64` and read it with a `meta_u64` that
  enforces the new documented cap, not `required_meta_u32`.
- Bump the archive format version and add a capability flag so an old reader
  rejects a wide-class archive loudly instead of truncating its metadata. New
  readers keep reading old `u32` archives unchanged.

### 3. Dense per-class side arrays (the real memory limit)

Several places allocate a vector with one entry per class,
`vec![_; n_classes as usize]`, for example `cell_of_class`, and the direct and
outcome arrays in `annotationcompare.rs`, and the EM support arrays. At `2^40`
classes a dense per-class array is on the order of terabytes and is not viable.

`cell_of_class` already decompresses in fixed `COC_BLOCK` blocks on demand, so
its on-disk form scales; the dense in-memory `Vec` is the issue. The wide-class
path needs these structures to be block-lazy or hash-backed rather than dense.
Enumerating and converting every `n_classes as usize` allocation is the bulk of
the implementation work and should be a checklist item, since a missed one
reintroduces a hard `usize`-sized dense array keyed by class.

### 4. Choosing the cap

`u64` everywhere is more than needed and makes the dense-array problem worse
without bound. A documented cap of `2^40 - 1` (about 1.1 trillion classes) gives
roughly 256x the current headroom while keeping the identifier inside 40 bits,
which is what the packed reducer's bit budget reasons about. Store it as `u64`
on disk for alignment simplicity, but enforce the `2^40 - 1` limit at
construction with the same loud error this branch adds at `2^32 - 1`.

## Suggested phasing

1. **(this branch)** Fail loudly at `2^32 - 1` during construction. No format
   change. Old and new archives read identically.
2. Introduce `ClassId(u64)` and thread it through memory and the on-disk stream
   behind a new format version and capability flag, keeping the reader
   backward compatible. Raise the construction cap to `2^40 - 1`.
3. Make the packed reducer allocate bits at runtime with a `u128` fallback.
4. Convert dense per-class arrays to block-lazy or hash-backed forms.

Phases 2 through 4 can land independently. Only phase 3 is needed before an
archive with more than `2^30` classes can actually be queried, and only phase 4
before one can be queried without a very large machine.
