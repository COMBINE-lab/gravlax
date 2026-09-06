# GQ v1 implementation and Gravlax 0.2.0 preparation

Historical first-implementation report. Performance and executor limitations below
are superseded by `gq-execution-optimization.md`, `gq-native-execution.md` and the
release allocator/acceptance report. In particular, the original 4.6–4.9× overhead
is not the current executor's measured overhead.

Branch: `astra-improvements`. Checkout: `/tmp/gravlax-v016`.
Version metadata is synchronized to **0.2.0**; changelog status is **Unreleased**.
No commit, push, tag, package publication, or deployment was performed for this task.
Pre-existing archive/query/compression changes in this checkout were preserved.

## Review disposition and parser choice

The final specification incorporates the final assessment with one important scientific
correction: compact representatives do not generally bracket omitted ends. The extractor
orders `(start, Reverse(end))`; `[100,120)`, `[105,200)`, `[110,130)` retains the first
and last but omits the largest end. GQ uses sound start intervals, invariant junction
gaps/internal blocks, and end bounds only where provable. It returns unknown elsewhere.
Read-support bounds incorporate retained witnesses and omitted multiplicity separately.

The three domains, classical `all`, explicit transcript/alignment strand conversion,
typed annotation unions/paths, weighted-read semantics, namespaced federation units,
explicit closure allowance, and Truth-state output were retained. Catalogue counters
were not blindly reinterpreted as upper bounds on class support.

**Selected parser: winnow 1.0.4.** It was already locked in the workspace and declares
Rust 1.65, compatible with Gravlax's tested Rust 1.89 floor. The implementation combines
winnow lexical combinators with a bounded, source-spanned precedence parser. Chumsky
remains a strong alternative for multi-error recovery/editor tooling; no comparative
performance claim is made without a head-to-head measurement. pest and LALRPOP were
considered but would add a grammar-generation/lowering layer without a compelling
benefit for this bounded grammar. See the parser decision in `gq-implementation-plan.md`.

## Implemented surface

- `aie gq validate`, `explain`, and `run`, typed plans, byte-position diagnostics,
  bounded parsing/expansion, and nonrecursive local/exported functions.
- Archive/federation records, complete exact raw-UMI classes, and cells. Unique,
  diagnostic stored, grouped signature/alternative, and record quantifiers.
- Region/junction/near/endpoint/overlap predicates, fixed consecutive/subsequence
  paths, explicit strand conversion, and record-linked terminal evidence.
- Accepted read totals, exact UInt64 literals, predicate-specific support bounds,
  explicit-denominator fractions, metadata expressions, Truth tallies, summaries,
  projection, ascending sorting, and disclosed final pagination.
- Pinned project gene/transcript lookup, exon unions, transcript junction paths,
  metadata schemas, assembly assertions, source-root verification and duplicate guards.
- Separate initial routing, class/cell closure and denominator accounting; explicit
  full-scan permission and finite chunk/record/evaluation/group/result budgets.
- Versioned tables under the existing result envelope, source/plan/resource provenance,
  atomic no-clobber publication, and Python client methods.
- Bounded unique-junction enumeration and exact record/class/cell support in one scan.

The general interpreter operates directly on `LazyArchive` and decoded evidence.
It does not invoke the legacy query command separately for each predicate.

## Validation

The final 0.2.0 tree passed:

- **439 Rust tests**, with one deliberately ignored parser benchmark. Full workspace,
  all targets and features, locked dependencies, tested toolchain Rust 1.89. The existing
  explorer loopback test required permission outside the filesystem/network sandbox.
- **Strict workspace Clippy**, all targets/features, `-D warnings`.
- **91 Python tests** on Python 3.9, using the already installed vendored `tomli` for
  packaging tests that need a TOML reader.
- **32-page documentation build** with the available Node 22.19 runtime, including
  the new GQ reference and search index.
- Release-version synchronization check for `v0.2.0` and `git diff --check`.
- Real Python-client execution of the CD45 query, producing a validated
  `gravlax.gq.result.v1` bundle with five observed Truth-state rows.

New tests include exhaustive small compact read-support bounds; the omitted-end
counterexample; classical universal and Kleene laws; cross-chunk class closure;
fine/coarse chunk invariance; grouped multimapper alternatives and multiplicities;
empty populations; exact enumeration; metadata nulls and federation namespaces;
pinned annotations and exon unions; typed function options/strand rules; budgets;
UInt64 arithmetic/fractions; and location-independent logical plan digests.

## Measured execution

Every benchmark invocation checked exact Truth-count agreement with the corresponding
legacy mode. The GQ population explicitly uses a stored-unique overlap witness to match
the legacy unique-only universe; this avoids conflating different populations.

Method: one warmup, five measured repetitions, deterministic randomized interleaving,
warm filesystem cache, four Rayon threads, complete command processes including JSON
output. GQ's general per-archive interpreter is currently serial; the specialized
command can use its existing parallel execution. Results are host/workload-specific.

| Archive | Query form | GQ median (ms) | Specialized command (ms) |
| --- | --- | ---: | ---: |
| Compact indexed | Independent record marginals | 76.90 | 16.34 |
| Compact indexed | Same-observation existential | 81.40 | 17.45 |
| Compact indexed | Non-vacuous universal | 83.25 | 17.16 |
| Full fidelity | Independent record marginals | 85.74 | 18.76 |
| Full fidelity | Same-observation existential | 92.49 | 19.83 |
| Full fidelity | Non-vacuous universal | 93.63 | 20.01 |

All six comparisons agreed exactly. Both archive representations produced the same
junction-query counts: 53,936 candidate records; 9,030 positive existential results;
6,302 positive non-vacuous universal results; zero joint positives for the two
mutually exclusive junction marginals in this fixture.

The initial general executor measured roughly 98–115 ms. Profiling identified memory
allocation/copying costs. Reusing cell metadata/group keys, avoiding repeated schema
construction, boxing large compile-time feature values, borrowing records, using
small inline geometry buffers, and reusing the universe predicate reduced times by
roughly 19–22% without changing answers. It is **still about 4.6–4.9× slower** than
the specialized command on these cases. This is an expressiveness release, not a
general query-speed improvement; existing commands remain available unchanged.

A parser-only release-mode microbenchmark measured approximately **6.4 µs/document**
for the 433-byte CD45 example over 10,000 iterations. This is not a comparison against
Chumsky and excludes binding/execution; parsing is not the main cost observed above.

Reproduction:

```sh
python scripts/benchmark-gq.py --binary target/release/aie \
  --archive /tmp/gravlax-astra-benchmark-v3/indexed.aie \
  --archive /tmp/gravlax-astra-benchmark-v3/fidelity.aie \
  --repeats 5 --output /tmp/gq-repeated-benchmark.json
```

`gq-benchmark-v020.json` records per-run measurements, binary/query digests, archive
content roots, and work counts. The script refuses to replace an existing report.

## Deliberate v1 boundaries and follow-up work

- Group keys are base sample/cell metadata, not arbitrary derived evidence states;
  joint evidence states belong in `tally`. This keeps pre-filter denominators defined.
- Class closure uses an available class index; cell closure currently requires an
  explicit full scan. Record/class denominators grouped by cell metadata also need
  a disclosed full scan. Sample-only source totals do not.
- Enumeration currently uses unique junction evidence, requires full-scan permission,
  and does not enumerate alternative-only candidates or report catalogue estimates.
- Sorting is ascending output-field sorting. No editor recovery/formatter, arbitrary
  joins, recursive graphs, unbounded path regex, posterior assignment or automatic
  assembly/contig migration is included.
- The generic frontend has its own typed execution plan. It reuses native archive,
  index, annotation and output infrastructure; the existing flag-based query engine
  has not been rewritten to share every logical-IR node.
- The main remaining performance opportunities are compact scalar execution slots,
  geometry/predicate caching and safe parallel/fused reduction. These require their
  own equivalence gates; they were not enabled speculatively.

## Entry points

- Specification: `gravlax-query-language-v1-rfc.md`.
- Review disposition/parser choice: `gq-implementation-plan.md`.
- User guide: `docs/src/content/docs/cli/gq.md`.
- Runnable examples: `examples/gq/cd45.gq`, `read-support.gq`, `junction-support.gq`.
- Native implementation: `crates/aie/src/gq.rs` and `crates/aie/src/gq/`.
- Python API: `Client.gq_validate`, `gq_explain`, `gq_run`, `gq_run_to_file`.

Release publication and any commit/push remain separate actions.
