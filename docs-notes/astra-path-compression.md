# Path geometry, independent streams, and OpenZL experiments

## Outcome

Two candidates deserve further work, for different reasons:

- **Position-predicted shared geometry:** 5.05% smaller under fine indexed
  chunking, 7.68% smaller under default chunking, with essentially unchanged
  native row decoding time on this fixture.
- **Grouped path geometry with adaptive entropy columns:** 7.85% and 10.19%
  smaller respectively, but native row decoding takes approximately 1.9x and
  1.8x as long as the existing representation.

These percentages compare complete experimental files against the unchanged
encoding in the same experimental container. The files include all original
dictionaries, indexes, auxiliary evidence, new framing, manifests and integrity
metadata. The actual source archives are slightly smaller than these matched
baselines because they lack the experiment declaration and longer section names.

No production archive format, default, ingest path, query adapter, or dependency
has changed. Use `aie dev path-bench` to reproduce the research codecs. `.gexp`
files are experimental containers, not production `.aie` archives.

## Complete-file results

All variants are lossless relative to the full-fidelity input. Decode times are
medians of 15 warm runs over all chunks: decompression plus direct reconstruction
of complete `MolRec` rows. Dictionary setup and disk I/O are excluded. These are
not end-to-end query or replay timings. Encoding time is reported separately in
the JSON results and is substantially higher for trial-based entropy selection.

| Encoding | Indexed: bytes | Decode ms | Default: bytes | Decode ms |
|---|---:|---:|---:|---:|
| Existing full-fidelity encoding, matched wrapper | 466,776 | 6.66 | 419,578 | 6.59 |
| Existing columns, independently compressed | 473,074 | 6.68 | 418,655 | 6.60 |
| Grouped paths, one frame | 441,446 | 9.67 | 389,777 | 9.28 |
| Grouped paths, four frame groups | 439,616 | 9.71 | 387,679 | 9.32 |
| Grouped paths, independent columns | 440,126 | 9.61 | 384,609 | 9.21 |
| Grouped paths, adaptive local entropy | **430,154** | 12.58 | **376,819** | 12.08 |
| Grouped paths, constrained counts | 443,269 | 9.78 | 390,405 | 9.26 |
| Shared geometry, absolute references | 498,474 | 6.39 | 444,275 | 5.92 |
| Shared geometry, reference deltas | 487,173 | 6.35 | 430,426 | 5.86 |
| Shared geometry, position-predicted reference deltas | **443,216** | **6.54** | **387,364** | **6.44** |

Actual existing full-fidelity source sizes are 466,431 and 419,312 bytes.
The corresponding existing compact archives are 384,486 and 340,008 bytes.
Thus the smallest new experimental files have fidelity premiums of **11.88%**
and **10.83%**, down from **21.31%** and **23.32%**. The near-baseline-speed shared
candidate has premiums of 15.27% and 13.93%. New experimental overhead is included
in these premiums; no hypothetical removal of old dictionaries is credited.

## What was implemented, and why it matters

### 1. Observed paths, sorted starts, and aligned-length residuals

The fixture has 161,818 geometry entries in 72,600 molecule/path groups. Of these
entries, 154,608 (95.54%) sum to 91 reference-consuming aligned bases. This is a
property of this dataset, never a decoder assumption or an assertion about
sequenced read length.

Store absolute observed splice junction paths in a chunk-local dictionary. Group
each molecule's geometries by path, store sorted starts as gaps from its anchor
and then its preceding start, and store exact total-aligned-length differences
from a chunk-local mode. The path and start determine the first and internal
blocks; the remaining aligned length determines the last block. Unspliced
geometry needs only start and length. Preserve a positive multiplicity per
geometry and recover the original canonical order and shape IDs on decoding.

This is not a new annotation dependency: paths come entirely from retained
evidence. Insertions/deletions retain the existing block-geometry semantics;
neither read sequence nor additional CIGAR detail is inferred or added.

The constrained-count experiment stores a group's total excess above one count
per geometry. It omits all counts when this excess is zero and infers the final
count otherwise. That simple variant lost; it does not rule out a more elaborate
combinatorial count-allocation codec. Start gaps were tested; occupancy bitmaps
and enumerative ordered-set codes were not implemented in this pass.

### 2. Independently authenticated side-by-side streams

Instead of an outer frame that must all be read and decompressed, substreams are
adjacent independent sections. The existing container supplies offsets, lengths,
per-section compressed digests and a rooted directory. No cross-chunk decoding
state is introduced. The four-group layout separates identity/MM evidence,
path dictionary, group structure, and geometry/count details. The fully split
layout exposes each geometry column independently.

An authenticated identity-only scan reads **78,004 vs 300,974 payload bytes** for
the indexed grouped candidates, and **72,935 vs 292,689** for default chunking:
74.1% and 75.1% less. These counts exclude already-loaded directory/dictionaries;
they do not describe every query. Identity values are checked against the source.

Naively splitting the old ten columns saves payload bytes but increases the
indexed file by 1.35% after framing. The effect is smaller under default chunking.
This is why payload-only comparisons are insufficient.

Local entropy selection compares complete zstd frames of varints against a
local rANS table and coded values, charging every mode/count/table byte. It keeps
varints where smaller. There is no new compression dependency. Selection is
per-column, not a claim of globally optimal compression. This prototype still
materializes some intermediate buffers, including varints after rANS decoding.

### 3. OpenZL backend ablation

An optional C helper builds against official OpenZL **v0.2.0**, commit
`3dceb64867840201fb8f57a29d179995f700c9b8`. It is not linked into Gravlax.
Inputs are typed, pre-entropy u64 streams exported by `path-bench`; signed
residuals already have their reversible zigzag mapping. Ten OpenZL graph/level
configurations test generic numeric compression, delta transforms, field-LZ,
range packing, FSE, zstd, flatpack, tokenization and bitpacking. All use standard
nodes, so no custom OpenZL decoder is required.

| Backend frame totals, geometry columns only | Indexed | Default |
|---|---:|---:|
| Varints + zstd 19 control | 187,176 | 184,146 |
| OpenZL generic numeric, level 19 | 267,469 | 250,449 |
| Best tested OpenZL graph per stream | 192,968 | 194,125 |
| Hybrid: OpenZL or varint/zstd per stream | 180,551 | 178,753 |

OpenZL alone loses even with per-stream selection. Hybrid savings are 6,625 and
5,393 bytes, before archive integration and backend-choice metadata. The main
wins are length residuals and multiplicities. Its measured backend decoding is
also slower here; conversion between typed values and row geometry is excluded.
The native Rust local-entropy variant obtains useful gains without introducing
OpenZL. This is evidence against adoption for these streams today, not a verdict
against trained multi-input graphs, other datasets, or the framework generally.

These totals are actual persisted, reopened, exactly verified codec frames,
**not whole archive sizes**. Bare C zstd frame sizes differ slightly from the
Rust streaming wrapper; the C control is used consistently within this ablation.
No cross-field training or held-out graph optimization was performed. Exhaustive
per-stream selection charges the sum of candidate encoding costs, not just the
winning candidate. Hybrid selection uses different decoded intermediate types,
so it still needs an explicit format adapter before production use.

### 4. Shared geometry needs position-conditioned references

There are 35,037 distinct absolute geometries across the fixture, and 140,351
entries belong to geometries occurring more than once. This does not itself
guarantee savings. A chunk-local dictionary with one-off geometries inline lost,
even after differencing successive references across the chunk.

The successful refinement uses the molecule's already-known anchor to locate a
lower bound in the position-sorted shared dictionary. Encode the first reference
relative to that bound, then subsequent references as positive rank gaps within
the same molecule. This preserves molecule order and membership without a
permutation stream. It reduces indexed reference/group payload from 206,574 to
151,312 bytes, turning a losing dictionary into a useful speed/size candidate.

This implementation retains the original access index. It does not claim savings
from unifying evidence and index postings, or from sharing whole molecule sets.

## Validation and limitations

- Every source and output container is fully integrity-verified.
- Every complete decoded row is compared for exact equality before and after
  writing/reopening each candidate: cell/class/chromosome/strand, ordered unique
  geometries and multiplicities, plus every multimapper tuple.
- Molecule and chunk order are unchanged; all non-chunk sections are copied
  verbatim in their original compressed form, including optional evidence.
- Five new tests cover fixed examples, 24 deterministic randomized mixed-length
  datasets, one through seven blocks, same-start geometries, MM-only rows,
  shared geometries, maximum counts, truncation, malformed lengths/IDs, overflow,
  conservation failures, entropy envelopes and deterministic encoding.
- Workspace all-feature/all-target tests pass (411 tests); strict clippy passes.
  The optional C helper also compiles with `-Wall -Wextra -Werror`.
- The sandbox initially denied a pre-existing loopback listener test; the full
  suite passes with local socket permission. No unrelated server code changed.

Only one real narrow-locus CD45 dataset is measured, under two chunk layouts;
these are not independent biological datasets. Synthetic tests establish codec
correctness, not general compression performance or improved biological accuracy.
The experiment retains full input chunks in memory and has a 16-million-item
decoder budget: it is a research harness, not a whole-atlas streaming replacement.

The next production-oriented candidate is position-predicted shared geometry;
the path/entropy representation remains the size-oriented alternative. Both
need broader libraries/loci and end-to-end query/replay and memory measurements
before changing defaults. No universal compression ratio is claimed.

## Reproduction and artifacts

```sh
cargo build --release --locked -p gravlax --bin aie
target/release/aie dev path-bench \
  /tmp/gravlax-astra-benchmark-v3/fidelity.aie \
  --out /tmp/gravlax-path-new-indexed --repeats 15
target/release/aie dev path-bench \
  /tmp/gravlax-fidelity-default-for-sparse.aie \
  --out /tmp/gravlax-path-new-default --repeats 15
```

Output directories must not already exist. Commands never rewrite the sources.
Each output includes the experimental files, exported typed streams and report.
Final artifacts are in `/tmp/gravlax-path-indexed-complete` and
`/tmp/gravlax-path-default-complete`; these are temporary local artifacts, not a
published or immutable release bundle. Checked-in full native reports are
`astra-path-indexed-results.json` and `astra-path-default-results.json`.

Optional backend benchmark (source/dependency download requires network access):

```sh
git clone --depth 1 --branch v0.2.0 https://github.com/facebook/openzl.git /tmp/gravlax-openzl-new
make -C /tmp/gravlax-openzl-new -j8 libopenzl.a
cc -O3 -std=c11 -Wall -Wextra -Werror \
  -I/tmp/gravlax-openzl-new/include -I/tmp/gravlax-openzl-new/deps/zstd/lib \
  scripts/pathbench-openzl.c /tmp/gravlax-openzl-new/libopenzl.a \
  /tmp/gravlax-openzl-new/deps/zstd/lib/libzstd.a \
  /tmp/gravlax-openzl-new/deps/lz4/lib/liblz4.a \
  -lm -pthread -o /tmp/gravlax-path-openzl-new
python3 scripts/pathbench-openzl.py /tmp/gravlax-path-new-indexed \
  /tmp/gravlax-openzl-new-results --helper /tmp/gravlax-path-openzl-new
```

The OpenZL reports bind the source archive root, pin the release commit, and hash
each exported typed input. Complete frame artifacts and reports are under
`/tmp/gravlax-path-openzl-indexed-complete` and
`/tmp/gravlax-path-openzl-default-complete`. Compact backend summaries are checked
in as `astra-path-openzl-indexed.json` and `astra-path-openzl-default.json`.
