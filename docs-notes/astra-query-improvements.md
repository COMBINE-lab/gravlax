# Astra query improvements: implementation and evaluation

Historical stage report: the later GQ implementation and release hardening
supersede statements below that the language is not implemented. See
`gq-native-execution.md` and the current GQ user guide for the implemented engine.

Branch: `astra-improvements`. Scope: a bounded implementation of the three query
directions (semantics, richer predicates, execution), plus a separate language RFC.
The textual language is not implemented. No archive/index format changed.

## What is included

- **Precise co-occurrence:** opt-in `--match-within any-placement|all-placements`.
  The entire Boolean expression is evaluated independently on each retained
  placement. Default record and explicit exact-UMI-class semantics are preserved.
  Quantifiers cannot manufacture a joint witness across multimapper alternatives.
  All is non-vacuous; terminal events are not falsely treated as placement-linked.
- **Richer geometry:** minimum single-block overlap, genomic start/end boundary
  windows, explicit junction tolerance, consecutive observed junction paths, and
  ordered subsequences. These compose in the existing Boolean expression interface.
- **Faster cooccur execution:** shared exact-junction lookup, a geometry-mask cache
  capped at 65,536 entries per chunk, and cell filtering before expensive geometry
  matching. The scalar matcher remains an independent reference for record queries
  and the default for small panels. Auto uses compiled matching at eight total
  predicates with at least seven geometry predicates.
- **Inspectable plans:** `--explain --format json` emits a serialized expression,
  resolved predicates, initial routes, cell-scope provenance, closure requirements,
  and route limitations without decoding molecule payloads.

The new geometry operations are opt-in, and do not increase archive storage.
They expose information already retained; their scientific utility is not a claim
that they make previously impossible biological inference possible.

## What is not included, and why

- Expression-based pruning of the current all-patterns result: pruning away false
  patterns would change counts and denominators. A selected-only DSL operator can
  have a different, explicitly optimized contract.
- Independent stream projection: the split-stream codec is still experimental,
  not a production archive format. The current reader still decodes complete
  selected chunks; this work does not claim fewer payload bytes for those chunks.
- General cell-level Boolean lifting or nested quantifiers in the CLI: existing
  cell aggregation is not the same operation. These are explicitly designed in
  the GQ v1 RFC, not faked by relabeling record counts.
- A new annotation-comparison or discovery engine: these already exist in
  `archivecmd/annotationcompare.rs` and `collectioncmd/search.rs`. Reuse their typed
  results and validated semantics in GQ rather than duplicate scientific code.
- A common executor for every existing command: the compiled matcher is currently
  scoped to cooccur. Existing batch/collection/replay operations are unchanged.

## Correctness and scientific checks

Full workspace/all-target/all-feature tests: **420 passed**. Strict clippy:
**passed with warnings denied**. The existing loopback-server test requires
socket permission outside the sandbox; the full run passed with that allowance.

New checks include:

- differential scalar/compiled masks across both region modes, all three placement
  modes, strands, chromosomes, compact/full geometry, tails, and cache reuse;
- counterexamples where both atoms match the record but no placement matches their
  conjunction, including mutually exclusive multimapper alternatives;
- compact-chain junction-path invariance, empty/nonempty quantifiers, unknown
  endpoint-sensitive predicates, exact threshold boundaries and coordinate overflow;
- ordered vs consecutive paths, invalid paths/thresholds, expression-depth budgets;
- end-to-end rich geometry queries under all backends and placement quantifiers;
- an archive with an intentionally invalid molecule payload: explain succeeds,
  execution fails, demonstrating that explain does not decode the payload.

`astra-query-geometry-results.json` contains 18 real-data checks, each evaluated
by scalar, compiled, and automatic backends with identical complete result data.
The candidate population is 53,931 retained unique-chain records; these are not
53,931 distinct physical molecules or deduplicated classes.

For the CD45 exon-A start boundary `[198696711,198696712)` on minus strand:
full geometry yields **76 witnessed records**, versus **72** in the compact
archive. The latter reports **21,557 unknowns**, not false absences. This high
unknown count is deliberately conservative: the current completeness test does
not further localize every omitted endpoint. Full geometry yields no unknowns.
The 12-base exon overlap query finds 9,041 witnessed records in both archives;
compact geometry leaves 18,088 additional records unresolved.

By contrast, the tested junction-only universal yields **4,339 records and zero
unknowns in both formats**. Retained and omitted members of a compact unique chain
share its junction path, so endpoint loss need not weaken junction-only logic.
This is a semantic precision improvement requiring no extra archive information.

The tested two-junction path has zero witnesses in both archives. The 1-base
tolerance test adds no witnesses over the exact junction (5,989 records each).
These are useful negative evaluations: tolerance was not enabled implicitly, and
no biological result was manufactured to demonstrate a new operator.

## Performance method

Use `astra-query-benchmark-final.json` as the authoritative benchmark, not the
earlier exploratory `astra-query-benchmark.json`. The final panel comes from the
whole archived PTPRC interval: 169 distinct observed junctions, with the most
supported 1/7/31/63 selected for panels of 2/8/32/64 total predicates including the
universe. **No duplicate predicates are used to inflate the larger panel.**

Four existing layouts are tested: indexed/fine-chunk full geometry and compact,
plus default/single-chunk full geometry and compact. Each contains 54,173 stored
records; the unique-chain population used here is 53,931 records across 10,842
observed cells. The subset selects every tenth sorted observed barcode: 1,085 cells.

Each of 32 workloads runs four backends: before-change executable, current scalar,
current compiled, and current auto. One warmup plus 11 measured repetitions per
backend are randomly interleaved with a fixed seed. This is 1,536 checked query
runs, including warmups. Use four Rayon threads and a warm filesystem cache; report
complete-process time including JSON output and peak RSS from GNU time. All result
table rows and scientific summary counts must match before timing is accepted.

The JSON records every timing/RSS sample, archive content roots, executable hashes,
host/platform, and logical result digests. It does not claim cold-cache, remote
federation, or whole-atlas performance. The archives are small and locus-specific;
the useful evidence is the controlled within-workload comparison.

### Final timing and decision

Representative full-geometry/fine-chunk timings, median complete-process ms:

| Workload | Before | Auto | Change |
| --- | ---: | ---: | ---: |
| 2 predicates, all cells | 18.00 | 18.64 | 3.5% slower |
| 8 predicates, all cells | 23.39 | 21.97 | 6.0% faster |
| 32 predicates, all cells | 41.98 | 28.47 | 32.2% faster |
| 64 predicates, all cells | 67.72 | 38.60 | 43.0% faster |
| 32 predicates, 10% cells | 36.00 | 14.27 | 60.4% faster |
| 64 predicates, 10% cells | 61.77 | 15.07 | 75.6% faster |

Across all four layouts, 32-predicate panels are 25.8–32.2% faster, and
64-predicate panels 35.8–43.0% faster. With the cell subset, those ranges become
58.9–60.4% and 73.8–75.6%. Eight-predicate unscoped panels gain 5.3–6.3%.

Two-predicate unscoped queries are **not a demonstrated speed win**: final medians
are 0.3–3.5% slower (at most 0.64 ms). They retain scalar execution; forcing the
compiled matcher is slower still. Do not describe this as accelerating every
query. The larger-panel and early-scope improvements are retained because their
benefits are substantial and repeat across all four layouts with identical rows.

Archive-size overhead is **zero bytes**. Runtime memory is not zero-cost: the
compiled path uses a bounded geometry cache, and exact peak RSS is recorded per
run rather than inferred from archive size. Scoped queries generally reduce
retained-hit memory; unscoped RSS varies with the allocator and chunk layout.
No new default ingest option or permanent access index is enabled by this work.

## Reproduction

```sh
cargo build --release -p gravlax --locked
python3 scripts/benchmark-cooccur.py \
  --before /tmp/gravlax-query-before \
  --after /tmp/gravlax-v016/target/release/aie \
  --archive /tmp/gravlax-astra-benchmark-v3/fidelity.aie \
  --archive /tmp/gravlax-astra-benchmark-v3/indexed.aie \
  --archive /tmp/gravlax-fidelity-default-for-sparse.aie \
  --archive /tmp/gravlax-astra-benchmark-v3/reader.aie \
  --junctions docs-notes/astra-query-junction-panel.json \
  --output /tmp/gravlax-query-benchmark-reproduction.json

python3 scripts/evaluate-query-geometry.py \
  --binary /tmp/gravlax-v016/target/release/aie \
  --archive /tmp/gravlax-astra-benchmark-v3/fidelity.aie \
  --archive /tmp/gravlax-astra-benchmark-v3/indexed.aie \
  --output /tmp/gravlax-query-geometry-reproduction.json
```

These commands assume the preserved local fixture archives and before executable;
they are not a standalone distribution bundle. Exact new-feature invocations are
also embedded in the geometry-results JSON. CLI syntax is documented in
`docs/src/content/docs/cli/query.md`. The proposed language, including the v1
boundary and implementation sequence, is in `gravlax-query-language-v1-rfc.md`.
