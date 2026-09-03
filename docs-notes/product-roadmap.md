# Gravlax product and capability roadmap

This note records the September 2026 ranking of opportunities to reduce archive
size, add scientifically useful queries, and improve usability. The ranking is
based on the current authenticated `.aie` format, the content-addressed
`.aicollection` implementation, the eight-donor SEZ archive audit, and the
implemented command surfaces. It is a decision record, not a claim that the
unimplemented size estimates have been experimentally established.

## Summary priorities

1. **User experience (delivered):** project workspaces and declarative,
   inspectable analysis plans, with the other interfaces built on the same
   resolved plan and result schemas.
2. **Capability (delivered):** paired annotation comparison and explanation,
   the most direct expression of the archive's purpose without new evidence
   fields.
3. **Space:** share canonical geometry and multimapper-pattern dictionaries in
   optional cohort packs. For standalone archives, investigate cell-local UMI
   edge coding. Remove the currently unused cell-posting index as an immediate
   small win.

## Implementation status

The first workflow release is implemented in this tree. It includes portable
projects and named resources; strict source plans, explained resolved snapshots,
and identity-checked resume; ambiguity-detecting gene/transcript/exon resolution;
the shared `gravlax.result-envelope.v1` contract and a typed subprocess-based
Python client; `doctor`, ingest recipe/preflight, shell completion, and the
public/development command split; and a loopback-only, read-only Explorer that
uses the same resolution and plan APIs. Paired annotation comparison is
implemented as a planned workflow and typed result. Annotation-conditional
transcript equivalence classes are also implemented, but remain experimental
under the validation gate described below.

Distribution machinery is implemented in-tree: native and Python release
assembly, integrity metadata, a Dockerfile and CI image gate, a local conda
recipe, and a Bioconda submission template. Static distribution tests pass, but
local container execution was blocked by the available rootless runtime rather
than claimed as tested. No GitHub release, PyPI package, registry image, or
Bioconda package has been published. Publication remains an explicit maintainer
action.

The uniform-output cascade is implemented across supported direct CLI query,
federation, cohort, collection, archive-lifecycle, compilation, export, and
extension surfaces, plus the corresponding operation subset represented in
source-plan v1; not every direct command has a project step kind. Python
validates both single-table envelopes and named-table bundles, converts
compatible tables to pandas or Arrow, and provides the existing row/MEX AnnData
paths. Historical command-specific JSON/TSV modes remain for byte compatibility
and are selected separately from uniform output. Native Arrow IPC is not
generally available from the CLI, and there is no R client; those remain
follow-up work rather than properties of the current release.

## 1. Archive-space opportunities

Across the audited eight SEZ archives, compressed section bytes are dominated
by molecule chunks (45.1%), cell-of-class mappings (29.7%), multimapper
patterns (17.5%), UMI adjacency edges (4.3%), and shapes (2.6%). Generic codec
tuning is therefore unlikely to yield a large improvement: the current stream
codes are already close to their measured memoryless bounds.

### Ranked opportunities

1. **Cohort-wide canonical shape and multimapper-pattern pool.** The eight
   archives contain 7.23 million local shapes that canonicalize to 4.77 million
   values and 17.27 million local patterns that canonicalize to 13.69 million.
   Shape and pattern sections occupy 265.4 MB in aggregate. An optional
   content-addressed cohort pack could plausibly reduce aggregate storage by
   roughly 2--4% after local-ID maps. Standalone archives should remain
   materializable so a missing shared object never makes evidence inaccessible.
2. **Cell-local UMI adjacency-edge coding.** Genome-ordered class identifiers
   make cell-local graph edges look nearly random. Encoding endpoints as ranks
   within a cell should preserve the exact graph while plausibly saving
   approximately 1--3.5% of total archive bytes, depending on edge density.
3. **Called-cell core plus a specialized ambient/rescue tier.** A separately
   encoded ambient tier may save roughly 4--10% on suitable datasets, but this
   estimate is less certain. The tier must remain lossless if later cell recall
   is supported; dropping it is a reduced-capability profile, not a transparent
   compression change.
4. **Remove or make optional `index.cellpost`.** The section is written but is
   not read by any implemented command. Removing it saves about 0.2--1.3%
   depending on the dataset. Either remove the corresponding direct-cell-route
   documentation or implement and evaluate a consumer.
5. **Explicit reduced-capability profiles.** One-representative or unique-only
   archives might save 10--18%, but they impair replay fidelity or multimapper
   inference. They must never replace the standard profile and must carry
   machine-readable capability flags that make unsupported queries fail.

### Low-priority codec experiments

A short block-local cell-of-class entropy audit is reasonable, but cell-ID
permutation cannot reduce rANS entropy and prior order experiments were
negative. Even an optimistic 6% cell-of-class gain is only 1.8% of aggregate
source bytes before local-table overhead. Promote a block-local model only if
a real re-encode, including table bytes and lazy-read costs, saves at least 2%
of total archive bytes. Wider chunks or stronger generic compression are not
priority projects.

Collection route bytes are accounted separately from `.aie` archive bytes.
Omitting local-shape routes makes a collection sidecar much smaller but loses
the corresponding junction-query acceleration and does not shrink the source
archives.

## 2. New query capabilities

### Ranked opportunities

1. **Counterfactual annotation comparison and explanation.** Scan the archive
   once, classify its evidence under annotations A and B, run both collapse
   reductions independently, and emit exact signed cell/gene deltas,
   class transitions, contributing causes, and bounded molecular witnesses.
   Final counts can be exact relative to the archived quotient and fixed
   alignment/barcode policy; this is not equivalence to full reads or to fresh
   annotation-aware alignment. Because collapse is nonlinear, explanations
   report contributing cause sets rather than claim a unique cause for every
   final count. This capability is implemented as `aie compare-annotations`
   and the `compare-annotations` plan step. Group deltas can be obtained by
   summing signed per-cell rows against an explicit group map; the current CLI
   does not accept a group scope or emit a group table directly.
2. **Lossless terminal-tail events and direct terminal-site queries.** Add an
   optional chunk-local, one-to-many event list attached to selected molecule
   ordinals, containing signed cleavage-anchor deltas and tail signals. The
   pilot found 136,763 deduplicated witnesses, 10.30-fold enrichment near
   external candidates, and 0.1261% retained-sidecar overhead. A production
   encoding must handle several events per molecule; the attempted one-event
   representation has a concrete losslessness counterexample.
3. **Annotation-conditional transcript equivalence classes.** Derive the set
   of compatible annotated transcripts for each UMI class using existing
   blocks, junctions, strand, and alternative placements. Equivalence classes
   are the deterministic product; transcript abundance is a separate inference
   layer. The experimental `aie query ... transcript-ecs` implementation
   exposes completeness, ambiguity, conflict, and no-compatible-transcript
   states over the retained archive quotient. Conflict and no-compatible flags
   are non-exclusive and must not be summed as a partition. Before promotion,
   run the predeclared full-read-versus-retained-representative benchmark and
   measure the fraction of singleton/ambiguous classes, especially for
   fragmented 3-prime protocols. Until that gate passes, no full-read
   equivalence, transcript abundance, isoform call, or full-transcript phasing
   claim is warranted.
4. **Paralog and mapping-ambiguity networks.** Project stored alternative
   placement patterns to locus/gene graphs and report unique, ambiguous, and
   inferred support by sample or group. A derived reverse alternative-locus to
   chunk index would accelerate these queries without adding new evidence.
5. **Sparse aligned-base/edit witnesses.** A new optional event section could
   support allele-specific expression, RNA-editing screens, and assisted donor
   demultiplexing. It requires prospective storage, ingest, bias, and privacy
   evaluation and cannot reproduce variant-aware realignment.

The implementation order selected for the current work is item 1 followed by
item 3. Item 2 remains the strongest future archive-format extension.

## 3. User-interface and workflow opportunities

### Ranked opportunities

1. **Project workspace and declarative plans.** Register archives,
   collections, annotations, sample metadata, groups, and conditions once.
   Versioned YAML/JSON plans should describe the scientific intent and compile
   into the existing typed command implementations. `plan check --explain`
   must resolve identifiers, validate assembly/annotation/scope/replicates,
   show exact predicates and defaults, estimate selected routes and I/O, and
   disclose the output schema before evidence is decoded. A fully resolved
   plan and provenance record accompanies every result and serves as its cache
   key.
2. **Scientific-intent query layer.** Accept unambiguous gene symbols, stable
   IDs, transcripts, exons, named events, and metadata predicates in addition
   to raw half-open coordinates. Always display resolved coordinates, strand,
   assembly, and annotation identity; ambiguous identifiers fail rather than
   selecting silently.
3. **Uniform typed outputs and language integration.** Standardize
   `--format text|tsv|json|arrow|mex`, `--output`, diagnostics, schema names,
   and provenance envelopes. Publish schemas and provide a Python client that
   returns Arrow, pandas, sparse matrices, and AnnData. Begin with a stable
   subprocess protocol rather than an unstable native FFI; add R after the
   result contract stabilizes.
4. **Onboarding, diagnostics, and CLI organization.** Provide `aie doctor`,
   ingest preflight, chemistry-aware recipes, packaged binaries and
   Bioconda/container distribution, a demo project, shell completions, and
   actionable errors. Move research instruments under `aie dev` or a separate
   binary so the public command tree expresses user tasks.
5. **Local interactive Explorer.** Build a localhost-only frontend for gene
   search, sample/group selection, event construction, splice/terminal views,
   routing and exclusion explanations, and export of the exact plan, command,
   or notebook fragment. It must delegate to the same plan/resolution/result
   API and contain no separate scientific logic.

### Selected delivery order and current state

1. Project/plan/check/explain and resumable execution: delivered.
2. Annotation comparison as the flagship planned workflow: delivered.
3. Scientific identifier and metadata resolution plus Python/AnnData access:
   delivered for the typed contracts described above.
4. Doctor, ingest onboarding, documentation, and public/development command
   cleanup: delivered. Packaging and release templates are built; external
   publication has not occurred.
5. Local Explorer on the shared APIs: delivered as a read-only planner and
   artifact browser.
6. Transcript-equivalence-class queries through the same plan and result
   contracts: implemented as experimental, pending the predeclared
   full-read-versus-retained benchmark.

## Evaluation

Archive changes require byte-identical scientific comparisons against the
current representation plus measured size, build time, query I/O, wall time,
and peak RSS. Capability queries require conservation properties and an
explicit distinction between exact archive-derived relations and inferred
quantities. UX changes should be tested with representative ingest, replay,
grouped splice-event, and multisample tasks, measuring first-attempt success,
time to a valid result, auxiliary files/commands, coordinate or replicate-unit
errors, and exact handoff reproducibility.
