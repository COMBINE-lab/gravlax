# Changelog

This file records user-visible changes in each Gravlax release.

## Unreleased

### Demonstration capsules

- A demo capsule may now publish a prebuilt, rooted `.aicollection` with shape
  routes plus a `--locations` manifest keyed by the committed source identities.
  `packaging/build_demo_capsule.py` builds it from a public-neutral staging root
  (`--collection-staging-root`, defaulting to the system temporary directory) and
  records its `aicollection-directory-root-v1` root; `finalize_demo_capsule.py`
  binds both files to the immutable data-release URLs; `verify_demo_capsule.py`
  verifies that root through the published location manifest and runs both
  collection stories against the relocated bytes without rebuilding.
- `demo-manifest.schema.json` gains an optional top-level `collection` object.
  Existing `demo-data-v1` manifests, which declare none, remain valid.
- Notebooks 02 and 03 use a published collection when the manifest declares one
  and otherwise rebuild it locally as before. The published path is fail-closed:
  hash-pinned bytes, no fallback address, the declared collection must commit
  exactly that story's archives and build options, and the location manifest must
  resolve exactly the archives already verified by archive root. Notebook 01 is
  unchanged.

## [0.2.2] - 2026-09-10

### Compatibility and assignment statistics

Existing `.aie` archives, including legacy v1 archives, and compiled `.aic`
annotations remain usable without re-ingestion or migration. Archive and
annotation encodings and default Gene count matrices are unchanged.

**Reporting change:** `replay-rows` reports now define `assigned_molecules` as
molecule records with at least one uniquely assigned representative. Previously
this field counted assigned representative rows, so its value can decrease even
when the emitted count matrix is identical. It now aliases
`assignment_statistics.assigned_molecule_records`; use
`assignment_statistics.molecule_records` as its denominator. Consumers needing
the previous row count must use `assignment_statistics.assigned_representative_rows`,
with `assignment_statistics.representative_rows` as the denominator.

Assignment statistics separately report molecule records, representative rows
and raw UMI classes across the full consumed input, including barcodes outside
the called-nucleus set. These are not collapsed-UMI counts. `counted_umis` reports
the full-input post-collapse total; use the emitted matrix sum for a selected
barcode population. This reporting correction also applies to default Gene
replay; archive compatibility does not imply unchanged report-field semantics.

### GeneFull replay and EM

- Add `replay-rows --gene-full` for intron-inclusive, strand-aware gene-span
  assignment from aligned blocks, with the existing global UMI collapse.
  Support streaming/eager archives, BAM input, and compiled annotations.
  Skipped alignment gaps alone do not assign a gene.
- Add Python `Client.replay(..., gene_full=True, solo_strand=...)` and record
  the counting model in replay reports and MEX provenance. Reject incompatible
  velocity/audit flags.
- Add `dev em --gene-full` for recovery evaluation and pooled emission, with
  explicit masked-class coverage and counting-model metadata. Accuracy is
  conditional on evaluable raw UMI classes before one-mismatch collapse.
- Add `dev em --solo-strand forward|reverse|unstranded` to recovery EM and its
  eager reference, and support `--star --gene-full` under all three strand
  policies. Apply the same model and strand to unique and ambiguous evidence;
  record them and the output counting units in metadata. Gene/forward remains
  the command default.

### Execution controls and validation

- Add independent `dev em --eval-barcodes` and `--modes` selection, preserving
  full-input prior fitting and the existing default models. Selecting fewer
  models reduces work without changing retained models' results.
- Spill packed EM supports according to actual retained candidate volume;
  expose `--support-memory-mib` (512 MiB by default; zero forces spilling) and
  report storage scope and capacity. This budget is not a whole-process memory
  cap and does not establish a reduction in peak process memory.
- Query GeneFull directly from aligned blocks using separate strand indexes.
  Isolated forward/reverse lookups were approximately 1.77 times faster; whole
  replay and EM timings do not establish a comparable end-to-end speedup.
- Validate against an exhaustive overlap oracle, synthetic STARsolo comparisons
  for both counting models and all strands, streaming/eager and forced-spill
  equivalence, and matched brain-nucleus matrices and EM. Fixed-nucleus EM keeps
  the original 6,460 nuclei and preserves the manuscript's shared metrics;
  nucleus calling is assessed separately from this fixed-population comparison.
  Default Gene replay matrices and EM results are preserved. The existing Gene
  UMI tie-breaking difference from STARsolo remains; synthetic EM agreement is
  not a claim of universal STARsolo equivalence.

See [GeneFull replay validation](docs/genefull-validation.md),
[EM validation](docs/genefull-em-validation.md), and
[strand and optimization validation](docs/genefull-em-strands-optimization.md)
for counting units, denominators, biological comparisons and measured tradeoffs.

## [0.2.1] - 2026-09-07

- Make rooted collection sources and parent layers relocatable through an optional
  `--locations` manifest keyed by committed content identities. Existing archive and
  collection formats, index bytes and roots remain unchanged; relocation adds no
  archive scan or decoding beyond normal authenticated access.
- Stop treating historical inode, device, size and timestamps as source acceptance
  criteria or cross-layer duplicate identities. Keep content-based duplicate checks,
  same-open-file authentication, consumed-payload checks and explicit full audits.
- Add location manifests to Python collection event search and collection-backed
  analysis plans, with resolved paths and manifest digests in uniform provenance.
- Legacy unrooted archives now authenticate each source open by full-file hash;
  their relocation remains supported but cannot provide the rooted fast path.
- Test moved layered bundles, unchanged archive/index hashes and I/O, cross-file
  replacement, malformed mappings, wrong roots, legacy identity and corruption.
- Normalize canonical path aliases in relocation acceptance tests on macOS and
  Windows while preserving content, query-result and I/O comparisons.

## [0.2.0] - 2026-09-05

### Composable query language

- Reuse built-in shape kernels for compatible GQ predicates, with a fixed-size
  placement cache and allocation-free shape validation. Offer opt-in
  `--parallel-decode` with bounded, ordered windows and reuse empty/sample-only aggregation keys
  per member. Preserve reference results, logical work and uncertainty diagnostics.
- Add an opt-in `mimalloc` executable feature (mimalloc 0.1.52, v2 backend), with
  allocator identification in GQ profiles and an interleaved latency/RSS benchmark.
  Keep the system allocator by default: measured speed gains carry a peak-memory cost.
  Reusable library crates and native C-library allocation are not overridden.
- Add GQ v1 with `aie gq validate`, `explain`, and `run`, source-located diagnostics,
  typed expressions, bounded nonrecursive functions, project annotations, metadata,
  archive federations, and native evidence execution. Parsing uses winnow 1.0.4;
  the tested Rust minimum remains 1.89.
- Distinguish unique observations, diagnostic stored representatives, grouped
  multimapping signatures and their alternatives. Universal quantification is
  classical; nonempty requirements are explicit. Record, exact raw-UMI class and
  cell scopes cannot silently substitute for one another.
- Add conservative three-valued geometry proofs and predicate-specific read-support
  bounds. Compact representatives bracket starts, not generally ends; fixed splice
  gaps and internal blocks provide additional sound proofs without archive changes.
- Support explicit alignment/transcript strand frames, fixed paths, union overlap,
  record-scoped tails, accepted-observation totals, grouped Truth tallies, summaries,
  exact denominator reporting, projection, sorting and disclosed pagination.
- Add bounded unique-junction enumeration with exact record/class/cell support.
  Full-scan fallbacks, closure and denominator costs require explicit permission.
- Add automatic physical GQ planning: fused filter/Truth-derive/aggregate kernels,
  compact Truth slots, sparse typed accumulators, bounded hash-based distinct counts,
  direct record iteration, unused identity projection pruning and stable top-k.
  Preserve logical work budgets, unknown diagnostics, floating-point input order and
  complete-unit/denominator semantics. Expose `--engine reference` and `--profile`
  for differential evaluation and post-parse phase timing.
- Emit typed GQ tables through the existing result envelope, content-bound provenance
  and atomic no-clobber output. Add Python `gq_validate`, `gq_explain`, `gq_run`, and
  `gq_run_to_file`, plus executable examples and a language reference.

Existing command semantics remain unchanged. In particular, the legacy co-occurrence
command's non-vacuous universal mode is not redefined. Its specialized executor can
be faster than the new general GQ interpreter; GQ is an expressiveness addition, not
a blanket performance replacement. No mandatory archive-format migration is needed.

### Archive access and optional storage improvements

- Speed up authenticated section lookup, directory/header reads and rANS decoding
  without changing their wire encodings. Project terminal-query identities without
  decoding unused molecule columns.
- Add opt-in `--access-index` class/geometry routes and `--chunk-records` finer
  access units, with authenticated completeness checks and bounded index expansion.
- Add opt-in `--geometry-fidelity` to retain distinct accepted unique geometries
  and multiplicities within the original molecule records. This can change query,
  assignment and velocity results and declares a new provenance reduction rule.
- Add experimental `--compression-tuning`, selecting final compressed cell-map
  and factored-shape encodings losslessly. New selected codecs and the fidelity
  provenance rule require a compatible reader; older readers fail closed.
- Keep all four ingest switches disabled by default. Default evidence, chunking,
  codecs and container/schema versions are unchanged; producer-version provenance
  changes archive bytes and content roots across releases.
- Retain structural, sparse-correction and split/path/shared-geometry experiments
  under `aie dev`; their experimental files are not supported production archives.
  OpenZL is an external research helper, not a Gravlax dependency.

### Built-in queries and release validation

- Add co-occurrence any/all-placement scopes, overlap/start/end predicates, explicit
  junction tolerances, paths/subpaths and payload-free query explanations. Preserve
  record-level defaults and non-vacuous legacy universal semantics.
- Accelerate larger predicate panels with compiled matching and early cell
  filtering. Small panels retain scalar matching; no archive storage is added.
- Add deterministic parser mutation tests and offline executable acceptance of
  archive/federation GQ, serial/parallel/reference execution and no-clobber output.
- Exercise GQ on release-target executables and extracted Linux release archives.
  Keep Rust 1.89 as the minimum supported version.

## [0.1.6] - 2026-09-04

### Correctness and compatibility

- Encode an empty multimap-pattern stream with the canonical empty rANS table,
  allowing BAMs containing only unique placements to be ingested, read, and
  validated normally.
- Lower the declared minimum supported Rust version from 1.98 to the tested
  dependency floor, Rust 1.89. Pull requests and pushes to `main` now run the
  locked workspace tests and deny all Clippy warnings on that exact toolchain.
- Activate the checked-in Google Colab notebooks against the immutable,
  checksum-bound `demo-data-v1` release locators, preserving fail-closed
  downloads and executable demonstrations.
- Correct the Bioconda Intel and Apple Silicon linker selection so both macOS
  builds use the active Conda C compiler driver. With the lower MSRV, the recipe
  returns to the standard Rust compiler activation package and no longer needs
  a compiler-policy lint exception.

Archive schemas, archive commitments, query result schemas, and command
interfaces are unchanged from version 0.1.5.

## [0.1.5] - 2026-09-04

### Molecular evidence and provenance

- Add logical molecular-evidence schema v2 within the authenticated `.aie` v2
  container. New archives carry a root-bound alignment-provenance manifest
  recording exact consumed-input identities, construction parameters, genome
  binding, alignment metadata, and an optional exact two-pass junction
  catalogue. Legacy archives remain readable and report unavailable provenance
  rather than inventing it.
- Add optional sparse terminal-tail evidence for uniquely mapped 10x 3′ cDNA
  reads. The side section retains every qualifying globally deduplicated
  cleavage-anchor key without adding read sequence, qualities, or names, and
  remains outside the unchanged core evidence streams when disabled.

### Coordinate-free and same-molecule queries

- Add `aie collection find-events` for coordinate-free junction, alternative
  donor/acceptor, cassette, and terminal-tail discovery across samples, donors,
  and cell groups. Queries can require recurrent exact UMI-class support and
  classify evidence against a supplied annotation as missing-junction,
  boundary, strand, or overlapping-model gaps.
- Add `aie query … cooccur` for bounded Boolean region, junction, and terminal
  predicates on the same retained molecule record, with three-valued handling
  of incomplete evidence and an explicit diagnostic exact-raw-UMI-class union
  mode.
- Add typed Python wrappers for collection event discovery and co-occurrence,
  including streaming result files and the same explicit resource limits as the
  native CLI.

### Demonstrations and maintenance

- Add three fail-closed Google Colab demonstrations for annotation
  reinterpretation, multi-donor event discovery, and federated
  junction/co-occurrence analysis. Add non-overwriting build, finalization, and
  verification tools for an independently rooted, checksum-bound demo-data
  capsule; the notebooks require immutable published locators and never fall
  back to invented or cached results.
- Resolve the existing workspace Clippy diagnostics and require the complete
  workspace, all targets, and all features to pass with warnings denied.

## [0.1.4] - 2026-09-03

This is the first complete Gravlax distribution. Version 0.1.0 established the
Rust packages on crates.io, but its immutable GitHub release contains only the
distribution manifest: native archives, installers, the vendored source
archive, and Python packages were not published for that version. The 0.1.1
release attempt stopped before publication when the native builds exposed two
platform-specific defects. The 0.1.2 attempt also stopped before publication
when its publisher jobs did not select the supported Python runtime. The 0.1.3
attempt built every release artifact but stopped before publication
because its checksum verifier rejected the trailing blank line emitted by
cargo-dist 0.32. Install version 0.1.4 when using a packaged release.

### Distribution changes

- Publish native archives and installers for 64-bit GNU and musl Linux, Intel
  and Apple Silicon macOS, and 64-bit Windows, together with checksums, build
  attestations, an SPDX dependency inventory, and a vendored source archive.
- Publish the `gravlax-client` wheel and source distribution to PyPI from the
  exact files attached to the immutable GitHub release.
- Pin portable release code generation explicitly: `x86-64` for 64-bit Linux
  and Windows, `penryn` for Intel macOS, and `apple-m1` for Apple Silicon.
  Intel macOS artifacts target macOS 10.12 or newer; Apple Silicon artifacts
  target macOS 11.0 or newer.
- Require successful validation and artifact assembly before creating the
  GitHub release, and keep Python 3.10 packaging validation compatible with
  the client's supported Python versions.
- Support concurrent positioned archive reads on Windows and use the Linux
  system-call interface needed by fully static musl builds.
- Build and smoke-test the Windows and musl targets on ordinary changes, before
  a release tag is created.
- Select the supported Python runtime explicitly in the registry publisher
  jobs.
- Accept cargo-dist 0.32's trailing blank line in `sha256.sum` while continuing
  to reject interior blank lines and malformed checksum records.

There are no archive-format, result-format, or command-interface changes from
version 0.1.0.

## [0.1.3] - 2026-09-03

The annotated 0.1.3 tag was retained after every release artifact was built.
Its checksum verifier rejected cargo-dist 0.32's trailing blank line in
`sha256.sum`, so the workflow stopped before obtaining registry credentials,
publishing to either registry, or creating a GitHub release. No 0.1.3 Python
package, Rust crate, native archive, or installer was publicly published. The
corrected distribution is version 0.1.4.

## [0.1.2] - 2026-09-03

The annotated 0.1.2 tag was retained after its release workflow found that the
registry publisher jobs did not select the supported Python runtime. The
workflow stopped before registry publication or GitHub release creation: no
0.1.2 Rust crates, Python package, native archives, or installers were
published. The corrected distribution is version 0.1.4.

## [0.1.1] - 2026-09-03

The annotated 0.1.1 tag was retained after its release workflow found Windows
and static-musl portability defects. Validation stopped before registry
publication or GitHub release creation: no 0.1.1 Rust crates, Python package,
native archives, or installers were published. The corrected distribution is
version 0.1.4.

## [0.1.0] - 2026-09-03

This is the first public release of Gravlax. It introduces compact,
molecule-resolved evidence archives for single-cell RNA-seq and lets an
annotation be changed or compared without realigning the source reads.

### Highlights

- Build authenticated `.aie` v2 archives from annotation-free alignments;
  validate, inspect, seal, stamp, extend, replay, and combine those archives.
- Replay Gene, GeneFull, and Velocyto-style matrices against a GTF or compiled
  AIC annotation while retaining fixed alignment and barcode-correction
  decisions.
- Query regions, junctions, junction sets, splice events, splice graphs, 3′
  endpoints, and transcript-compatibility classes at cell, group, bulk, or
  multi-sample scope where supported.
- Compare two annotations on the same retained evidence and report signed
  count changes, class transitions, non-exclusive causes, and bounded
  witnesses.
- Define reusable projects and analysis plans with content-bound inputs,
  explicit biological intent, safe resume, and a loopback-only Explorer for
  building and inspecting plans.
- Request typed text, TSV, or JSON results from supported commands without
  changing their established default output; operation reports and directory
  bundles use the same schema and provenance model.
- Control the `aie` executable and validate typed results from Python with the
  dependency-light `gravlax-client` source package.
- Install the `aie` executable from crates.io. Native archives, installers,
  the vendored source archive, and the PyPI package first become available in
  version 0.1.4.

### Compatibility notes

- New archives use `.aie` v2. Seekable v1 archives remain readable and can be
  sealed into authenticated v2 containers without recompressing section
  payloads.
- Transcript-equivalence results describe compatibility with the retained
  archive representatives; they are not transcript-abundance estimates or
  full-isoform phasing.
- The Python distribution does not embed the Rust executable. Install `aie`
  separately and keep its version aligned with `gravlax-client`.
- Gravlax 0.1 is an initial public interface. Result and archive formats are
  versioned so readers can reject incompatible future changes explicitly.

[0.1.6]: https://github.com/COMBINE-lab/gravlax/releases/tag/v0.1.6
[0.1.5]: https://github.com/COMBINE-lab/gravlax/releases/tag/v0.1.5
[0.1.4]: https://github.com/COMBINE-lab/gravlax/releases/tag/v0.1.4
[0.1.3]: https://github.com/COMBINE-lab/gravlax/tree/v0.1.3
[0.1.2]: https://github.com/COMBINE-lab/gravlax/tree/v0.1.2
[0.1.1]: https://github.com/COMBINE-lab/gravlax/tree/v0.1.1
[0.1.0]: https://github.com/COMBINE-lab/gravlax/releases/tag/v0.1.0
