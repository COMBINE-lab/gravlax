# Structural compression and the cost of full geometry fidelity

These are read-only, reversible coding experiments on `astra-improvements`.
They do not change the ingest defaults or introduce a supported archive codec.

## Main result

The best tested factorization removes approximately a quarter of the additional
cost of fidelity. It does **not** make fidelity nearly free.

| Representation, same optional index and chunk boundaries | Compact evidence | Full fidelity | Additional fidelity bytes | Fidelity premium |
|---|---:|---:|---:|---:|
| Existing encoding | 384,486 | 466,431 | 81,945 | 21.31% |
| Position deltas + local rANS + separated layout | 380,094 | 441,674 | 61,580 | 16.20% |

New values are projected whole-file bytes from **actually compressed, verified
candidate frames**, retaining existing dictionaries, entropy tables, and indexes.
They charge one experiment-envelope byte in each transformed chunk. A production
codec/capability declaration has not been designed, so its metadata overhead is
not included. The outputs are measurements, not newly usable `.aie` archives.

The fidelity archive shrinks 24,757 bytes (5.31%). However, the compact archive
also shrinks 4,392 bytes; the fair reduction in the *additional fidelity cost* is
therefore 20,365 bytes (24.85%), not 24,757. Comparing the new fidelity encoding
only against the old compact encoding would suggest a 14.87% premium and would
overstate the matched-encoding improvement.

## Why this helps

The full-fidelity fixture has 161,818 retained unique representatives, compared
with 108,844 in the compact encoding. The position stream grows from 54,908 to
107,882 values because the first representative of each record is implicit.

The current position stream codes each representative's distance from its
molecule anchor. The candidate instead codes the signed difference from the
preceding representative, resetting at each molecule. This preserves ordinary
two-representative chains too, including negative steps between chains. Full
fidelity's position-canonical singleton entries make these increments especially
useful. A per-chunk rANS table models the transformed values; its complete stored
table cost is included. Strand, chain count and multimapper count are separated
into three layout columns. No counts, identities or geometry are discarded.

## What did not work

| Candidate | Compact projected bytes | Fidelity projected bytes |
|---|---:|---:|
| Unchanged | 384,486 | 466,431 |
| Class tokens as varints | 385,669 | 466,985 |
| Birth bitmap + class backward distances | 385,172 | 466,592 |
| Birth bitmap + repeated-class ID deltas | 384,340 | 465,825 |
| Local frequency-ranked shape IDs | 389,774 | 471,226 |
| Separated layout columns | 383,412 | 463,228 |
| Representative position deltas, varint/zstd | 388,119 | 453,703 |
| Position deltas + layout, varint/zstd | 387,219 | 449,953 |
| Position deltas + layout + local rANS | **380,094** | **441,674** |
| Per-molecule shared splice chains and terminal-length residuals | 395,386 | 487,266 |
| Position deltas + chain residuals + layout + local rANS | 390,901 | 462,887 |
| Shared complete molecule-geometry templates | 455,447 | 575,680 |

The chain-residual prototype shares exact absolute junction boundaries within a
molecule, references the previous representative on that chain, and stores only
the change in terminal-block length. Position changes determine the first-block
length and internal offsets. Unspliced reads share the empty junction chain.
This is exact, including nonoverlapping unspliced reads and changing first/last
exon lengths, but references/residuals cost more than the current shape-ID stream
on this fixture. This rejects this particular encoding, not every possible chain
factorization.

Whole-geometry templates share strand, chain/representative layout, relative
positions, shapes and multimapper patterns across records, with weights and UMI
classes kept separately. The complete template dictionary and ID stream outweigh
the eliminated repetition. The same observation applies to local shape-ID maps.

Earlier shape-*dictionary* factoring remains useful: its compressed dictionary
falls from 15,689 to 10,254 bytes. That is different from factoring per-record
fidelity geometry. It affects both evidence modes and does not explain away the
extra fidelity cost. Of the previous 5,684-byte whole-file compression-tuning
saving, 5,435 bytes came from this dictionary factorization, 267 from cell-map
selection, offset by 18 extra section-name bytes.

## Validation and limitations

- Every input payload is authenticated. Each candidate is compressed, decompressed,
  inversely transformed, and compared with the **exact original chunk bytes** and
  all decoded `MolRec` values, including representative ordering and weights.
- Seven repetitions per candidate, all chunks, on each of three real CD45 archives:
  default compact, indexed compact, and indexed full-fidelity. Same zstd level 19.
- Tests cover empty and repeated classes, cross-chunk class bases, malformed
  templates, nonoverlapping unspliced geometry, spliced terminal changes, and
  preservation of both unique and multimapper geometry columns.
- Full workspace/all-target/all-feature suite: **402 tests pass**, 23 targets.
  Strict clippy with `-D warnings` passes.
- Timings in the raw reports are research-harness timings. Candidate restoration
  reconstructs canonical streams, including rANS re-encoding, before invoking the
  current decoder. These are **not native candidate reader or query/replay speed
  measurements** and must not be used to claim a runtime regression or improvement.
- This is one narrow-locus corpus. There is no genome-wide result or proof that
  16.2% is a lower bound. No new production encoding has been enabled.

Follow-up: this experiment has now been performed; see `astra-sparse-fidelity.md`.
The tested sparse capsules reconstruct exactly but are larger than the existing
full-fidelity archive under both tested chunk layouts.

The next distinct experiment proposed here was a compact base plus sparse
geometry/multiplicity exceptions, with an explicit reconstruction contract.
Simply labeling the current full-fidelity encoding sparse does not provide that
contract. Such an exception codec must account for its routing, references and
count corrections, not merely measure omitted coordinates in isolation. Native
decoding and wider-corpus benchmarks are needed before promoting any winner.

## Reproduce

```sh
cargo build --release --locked
target/release/aie dev structural-bench INPUT.aie --level 19 --repeats 7
```

The command emits JSON to stdout and never writes the source archive. Reports
are checked in alongside this note as `astra-structural-{default,compact,fidelity}.json`.
Inputs are the `reader.aie`, `indexed.aie` and `fidelity.aie` outputs under
`/tmp/gravlax-astra-benchmark-v3`, reproducible with the earlier benchmark script.
Measured binary SHA256:
`99adf95145dff90817744283557b307722b997fc59d17412e50bf0e4b18182e7`.
