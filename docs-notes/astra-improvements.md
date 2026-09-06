# Astra archive improvements

Branch: `astra-improvements`, baseline `ef94f220` (v0.1.6).

Follow-up structural/fidelity-cost experiments are recorded in
`astra-structural-compression.md`; they are not enabled production codecs.

Implementation order:

1. Faster selective access: byte-compatible section lookup/read and decoder improvements;
   authenticated optional class/geometry postings; bounded access units; projected tail identities.
2. Sparse geometry fidelity: retain distinct omitted unique-read placements with explicit
   completeness, original record association, multiplicity and capability semantics.
3. Deeper compression: final-size codec selection, geometry factoring and class-to-cell coding experiments.

Validation requirements:

- Existing archives remain readable and existing default semantics remain unchanged.
- Derived indexes bind to the authenticated source root and cannot replace molecular evidence.
- Index absence has a correct fallback; partial/corrupt indexes fail closed.
- Additional geometry must be validated against direct BAM evidence, including adversarial
  non-nested spans and same-UMI records at separate loci.
- Compression experiments must report complete archive plus index size, construction time,
  query/replay time and peak memory. Do not promote an experiment solely on encoded size.
- Preserve full payload verification, deterministic construction, bounded decoding and explicit
  unsupported-codec rejection.

Available initial real-data fixture: public 10x CD45/PTPRC BAM and v2 archive under
`/tmp/gravlax-cd45-case-study`. This is a narrow-locus fixture, not evidence of genome-wide
performance. Whole-genome corpus access must be restored before making atlas-scale claims.

## Implemented

1. Section-name lookup is indexed; directory reads are buffered; inline header
   checks use one read per section. Integrity checks and old wire bytes remain
   unchanged. rANS expands values during symbol decoding, removing an intermediate
   vector. Optional `index.access` routes cross-record classes, aligned blocks and
   multimapper-alternative junctions. `--chunk-records` permits smaller access
   units without splitting equal anchors. Collection tail queries decode only
   the core identity columns needed to attach authenticated tail evidence.
2. `--geometry-fidelity` retains distinct unique geometry and exact multiplicity
   within the original molecule record. This is a sparse replacement of affected
   chains with singleton geometry entries, not a separate sidecar. The root-bound
   provenance rule explicitly changes to `distinct-unique-geometries-v1`.
   Direct-BAM replay has a corresponding reference switch. No sequences or
   qualities are retained and multimapper capture is unchanged.
3. `--compression-tuning` compares actual compressed frames for delta/rANS/run
   cell mappings and an optional factored splice-shape dictionary. It reuses the
   selected frame and keeps old encodings on ties. High-level compression uses
   at most four concurrent contexts for either fine chunks or codec tuning.

The format and CLI usage are documented in `docs/src/content/docs/format.md`.
The ten core streams still share one frame per chunk. A broader replay cache was
not retained: it needs explicit annotation/dictionary lifetime management and an
independent benchmark before promotion. Molecule-class stream redesign, separately
compressed core columns and atlas-scale partitioned postings remain future work;
they are not represented as completed optimizations here.

## Validation

- `cargo test --workspace --all-targets --all-features --locked`: **399 passed**,
  23 test targets, including existing rooted integrity/legacy archive tests. The
  existing Explorer loopback tests require execution outside the network sandbox.
- Strict all-target/all-feature clippy with `-D warnings`: passed.
- New regressions cover omitted middle geometry, duplicate reads, a same-UMI
  separate locus, exact BAM/archive Gene and velocity MEX parity, index fallback,
  all-alternative routing, truncated indexes/shapes, missing/extra routes, and
  bounded cell-run expansion.
- Every benchmark archive passes `doctor --verify-content`. The default candidate
  archive is byte-identical to baseline (SHA256
  `75b7ab2be2550cbc589636b8055bf778285d2daa777e9cfb3a3ad5acaec40486`).
- Lossless access/compression variants have identical cooccur tables and replay
  matrices. Fidelity and fidelity+compression likewise match each other.
- On the real fixture, direct-BAM fidelity replay and the all-options archive
  match byte-for-byte for Gene and spliced/unspliced/ambiguous velocity matrices.
  Gene totals change from 23,160 (compact) to 23,268 (fidelity). This is fidelity to
  accepted geometry, not independent STARsolo/biological accuracy validation.

## Measured results

Machine: 128 available CPUs, Rust 1.89 release build. Fixture: 54,173 molecule
records, 37,895 UMI classes, 10,845 cells, 4,567 shapes. Warm-process-launch query
medians use seven runs, replay medians three; ingest is one measured run per
configuration and its RSS is noisy. The query is the CD45 A-exon/junction
intersection, exact UMI-class unit. Replay uses the supplied CD45 annotation,
unstranded mode and the full barcode whitelist (including MEX output cost).

| Configuration | Whole archive bytes | Query median ms | Replay median ms | Ingest ms | Ingest peak MiB |
|---|---:|---:|---:|---:|---:|
| v0.1.6 baseline | 340,008 | 18.18 | 72.33 | 338.8 | 1,047.6 |
| New reader, default encoding | 340,008 | 18.15 | 75.57 | 338.5 | 941.6 |
| Access index + 4,096-record target | 384,486 | 15.95 | 71.08 | 347.8 | 365.2 |
| Compression tuning, level 19 | 334,324 | 18.24 | 74.28 | 333.6 | 366.3 |
| Index + fine chunks + fidelity | 466,431 | 16.88 | 73.10 | 357.7 | 370.5 |
| All options, level 19 | 460,747 | 17.28 | 73.65 | 376.5 | 370.0 |
| Compression tuning, level 3 | 377,509 | 18.86 | 76.27 | 259.6 | 42.3 |
| Compression tuning, level 9 | 361,678 | 18.65 | 75.88 | 259.5 | 108.2 |

Interpretation:

- Indexed query latency is 12.3% lower, at 13.1% additional total archive size.
  The index itself is 25,047 compressed bytes; remaining overhead comes from finer
  chunking and its effects on compression. This query reads 11 of 14 chunks.
- Compression tuning saves 5,684 whole-file bytes (1.67%). Shapes shrink from
  15,689 to 10,254 compressed bytes (34.6%); cell mapping saves another 267 bytes.
- Fidelity adds 81,945 bytes over the indexed compact archive (21.3%). The query
  still selects 3,956 UMI classes, but its retained unique-geometry negative scope
  becomes complete. This does not establish biological absence.
- Lower zstd levels trade about 28% less ingest time for 6–11% larger files versus
  the old level-19 baseline. These single-run ingest results are exploratory.
- There is **no demonstrated whole-archive replay speedup** on this fixture.
  The larger bounded classification cache was therefore not enabled.
- Synthetic 4,096-section measurements improve median open time from 4.457 to
  1.652 ms, 20,000 section reads from 43.960 to 5.666 ms, and one-million-value
  rANS decoding from 6.301 to 4.155 ms, with identical encoded rANS bytes. These
  isolate mechanisms and must not be reported as genome-wide query gains.

## Reproduction and artifacts

Checked-in results: `astra-benchmark-results.json`, `astra-access-baseline.json`
and `astra-access-final.json` in this directory. Raw per-run outputs, archives,
RSS and logs are in `/tmp/gravlax-astra-benchmark-v3`; these are temporary local
artifacts, not a durable published bundle. Binary/input SHA256 values are recorded
in the benchmark JSON. The original binary is `/tmp/gravlax-astra-baseline/aie`.

```sh
cargo build --release --locked
python benchmarks/astra_selective.py \
  --baseline /tmp/gravlax-astra-baseline/aie --candidate target/release/aie \
  --bam /tmp/gravlax-cd45-case-study/source/ptprc-grch38-complete-multiplicity.bam \
  --whitelist /tmp/gravlax-cd45-case-study/source/737K-august-2016.txt \
  --gtf /tmp/gravlax-cd45-case-study/source/CD45_exons_nochr.hg38.txt \
  --barcodes /tmp/gravlax-cd45-case-study/source/737K-august-2016.txt \
  --out /tmp/gravlax-astra-new-run
cargo run -p gravlax-evidence-io --example access_benchmark --release --locked
```

The output directory must not already exist. The filtered cell list is not a
valid replay output scope for this BAM because additional whitelist cells are
present. Full-atlas benchmarking and independent accuracy evaluation remain
necessary before changing defaults. No release or deployment is part of this branch.
