# GQ execution optimization

Implemented on `astra-improvements`, for the unreleased 0.2.0. No archive/index
format change, additional stored data, new dependency, commit, or publication.

## Outcome

The first completed interleaved benchmark found approximately **2× faster record
tallies/summaries**, **1.45× faster class summaries**, **1.33–1.35× faster distinct
counts**, and **1.69× faster read-support queries with a result limit**. All results,
denominators, unknown diagnostics, and logical expression-step counts matched the
reference executor. Peak memory for the compact CD45 example fell approximately 31%.

The final rerun additionally covers junction enumeration and confirms it is unchanged.
That rerun encountered substantially higher shared-server contention: the observed
load average afterward was 97.90 on a 128-logical-CPU host. Absolute timings increased
for both engines. Both datasets are retained below; do not compare an engine's time
from one run to its counterpart from the other, or treat either run as a latency SLA.

## What changed

1. **Physical lowering and operation fusion.** Eligible ordered filter/Truth-derive/
   tally or summarize pipelines lower to a compact kernel. Derived Truths use indexed
   slots, not per-unit named value-map entries. The kernel supports nested record,
   unique, stored, multimap and alternative quantifiers without combining their scopes.
   Native scalar leaves remain available; unsupported derived-value shapes fall back.
2. **Direct record execution and projection pruning.** Selected record ordinals and
   inline indices replace singleton unit-group trees. Unused `unit.id` strings are
   not formatted. When no metadata columns are declared and no cell fields are read,
   a sample row is reused. Required metadata validation and cell-ID bounds checks are
   not skipped. Class/cell queries still assemble complete children.
3. **Sparse aggregation.** Interned group identifiers and compact Truth-state keys
   select accumulators without building a result row each time. Count updates use
   checked integer addition; aggregate argument storage is reused. Distinct counts
   use bounded hash sets, and an already-present value remains legal at the distinct
   limit. Accumulator state moves between federation members rather than being copied
   or combined in a different order. Floating-point sums retain reference input order.
4. **Geometry microkernels.** Single-junction matching avoids an intermediate junction
   vector. Omitted-geometry proofs use inline small buffers instead of several heap
   vectors; read-total evaluation no longer allocates a record-index vector. These
   low-level improvements are shared with the reference interpreter.
5. **Stable top-k.** Adjacent `sort |> take` partitions the fully evaluated results and
   sorts only the requested prefix. Original row ordinals resolve ties, preserving the
   stable reference result. This does not truncate evidence, bypass row budgets, or
   claim a streaming-memory bound; full result materialization remains.
6. **Timing and explainability.** `--profile` reports internal post-parse timings;
   `gq explain` describes the selected physical strategy. `--engine reference`
   disables fusion, direct record iteration, projection pruning, and top-k. Dedicated
   exact junction enumeration is explicitly identified and timed under both engines.

There is no quantifier-domain early exit or predicate pushdown across universe/closure
boundaries. Even otherwise-decisive quantifiers retain the reference omitted-cause
counters and logical work charges. There is no `P | !P = true` rewrite under Kleene
logic. Read-support bounds and general non-Truth derived intermediates remain on the
interpreter, although direct iteration, projection pruning and top-k still help them.

## Where the time goes

Internal phase medians for compact record-marginal tally, first completed run:

| Phase | Reference ms | Auto ms |
|---|---:|---:|
| Preflight planning | 1.26 | 1.27 |
| Archive setup | 1.24 | 1.30 |
| Chunk read/decode/attachments | 6.19 | 6.17 |
| Universe witness and selection | 10.41 | 4.22 |
| Unit/denominator preparation | 4.19 | 0.004 |
| Evaluation and aggregation | 37.81 | 17.74 |
| Finalization | 0.008 | 0.004 |
| Entire execution | **64.97** | **31.44** |

Phase medians are computed independently and omit some teardown/dispatch time, so
they need not sum to the overall median. Resource binding/type checking, serialization
and publication are measured separately. Parsing, source reading and process startup
are outside the execution timer. Evaluation and aggregation are deliberately reported
together because fusion interleaves them; this is not a claim that aggregation is free.

The geometry-joint query still spends 38.39 of 52.92 ms in evaluation/aggregation in
that run. Read-support finalization falls from 21.48 to 13.12 ms, including sorting,
projection and disposal of unreturned materialized rows. Profiling initially showed
substantial allocator and map overhead; after optimization geometry evaluation and
the compact evaluator account for more of the remaining samples.

## Benchmarks

Host: `nomad01.umiacs.umd.edu`, AMD EPYC 9555, 128 logical CPUs; Rust 1.89.0 release
build. One warmup followed by seven measured repetitions per engine/case/archive,
randomized engine order, warm filesystem cache. `RAYON_NUM_THREADS=4` is fixed; this
round does not parallelize the GQ evaluation loop. The selected locus has 54,173
decoded records in 14 chunks in each archive. These are locus-scale fixtures, not
evidence of whole-transcriptome or cold-storage scalability.

### First completed run

This contains the full optimization set before the final enumeration timing addition.
Values are **reference → auto**, internal execution milliseconds:

| Query | Compact | Fidelity | Speedup range |
|---|---:|---:|---:|
| Record marginals | 64.97 → 31.44 | 70.44 → 35.83 | 1.97–2.07× |
| Same-observation disjunction | 70.38 → 34.97 | 76.95 → 39.69 | 1.94–2.01× |
| Nonvacuous universal | 71.18 → 35.28 | 78.16 → 40.53 | 1.93–2.02× |
| Geometry + joint tally | 100.16 → 52.92 | 107.70 → 58.92 | 1.83–1.89× |
| Truth counts + read sums | 80.05 → 40.82 | 88.99 → 46.76 | 1.90–1.96× |
| Class closure + summary | 74.03 → 50.75 | 82.85 → 56.98 | 1.45–1.46× |
| Distinct cells | 65.72 → 48.54 | 70.19 → 52.91 | 1.33–1.35× |
| Read bounds + sorted first 20 | 103.14 → 60.93 | 110.41 → 65.41 | 1.69× |

Raw repetitions, full verified results and binary/query hashes:
[first-run JSON](gq-execution-optimization-initial.json).

### Final rerun under higher shared-host load

| Query | Compact, reference → auto ms | Fidelity, reference → auto ms |
|---|---:|---:|
| Record marginals | 189.18 → 67.08 | 198.09 → 76.56 |
| Same-observation disjunction | 197.43 → 74.21 | 219.10 → 83.11 |
| Nonvacuous universal | 217.75 → 78.48 | 223.40 → 85.74 |
| Geometry + joint tally | 288.23 → 116.24 | 300.70 → 126.64 |
| Truth counts + read sums | 160.54 → 66.50 | 254.80 → 109.46 |
| Class closure + summary | 159.55 → 99.83 | 236.94 → 143.66 |
| Distinct cells | 204.72 → 145.84 | 225.88 → 121.85 |
| Read bounds + sorted first 20 | 290.65 → 135.23 | 317.45 → 150.33 |
| Exact junction enumeration | 103.78 → 104.24 | 116.92 → 118.34 |

Junction enumeration uses the same executor under both settings; the small measured
differences are not an optimization or a claimed regression. The other paths retain
their advantage, but load-dependent speedups should not be promoted as guarantees.

[Final JSON](gq-execution-optimization-benchmark.json) records all 18 comparisons.
Final binary SHA-256:
`a09a867e6a542c60ce950d23292ebcf44a4a5894b51e75bebfb72db2c9bfcd8f`.

Peak process RSS for `examples/gq/cd45.gq` on the compact archive, three runs, KiB:
auto `[34300, 34924, 34360]`; reference `[49384, 49856, 49752]`.
Medians: 34,360 versus 49,752 KiB, a **30.9% reduction**. This measures whole-process
peak memory, not just the kernel. There is zero additional archive/index storage.

## Validation and reproduction

- 441 Rust tests passed; one existing manual parser benchmark ignored.
- 91 Python tests passed; strict workspace clippy passed; 32-page documentation build
  passed; `git diff --check` passed.
- All 18 final archive/query comparisons check complete typed tables and semantic/work
  summaries on every run, including warmups. Only physical descriptions and timings
  are excluded from equality.
- Synthetic differentials cover cross-chunk class closure, cell scopes including empty
  cells, grouped multimapper alternatives, unknowns, metadata nulls, federation groups,
  wide tallies/aggregates, tied and zero/oversized top-k limits, and one-step-short
  budgets. A federation regression specifically detects floating-point reassociation:
  per-unit additions of 2.0 after a 4e16 accumulator must not become one addition of 8.0.

```sh
cargo build --release -p gravlax --bin aie --locked --offline
python scripts/benchmark-gq-execution.py \
  --binary target/release/aie \
  --archive /tmp/gravlax-astra-benchmark-v3/indexed.aie \
  --archive /tmp/gravlax-astra-benchmark-v3/fidelity.aie \
  --repeats 7 --output /tmp/gq-execution-new-run.json
```

Use a new output filename; reports are not overwritten. Each report includes full
query source, query digest, binary digest, authenticated source roots, all timings
and verified results. For a single query:

```sh
target/release/aie gq run examples/gq/cd45.gq \
  --bind pbmc=/tmp/gravlax-astra-benchmark-v3/indexed.aie --profile
```

## Remaining opportunities

The next evidence-backed targets are shared geometry evaluation across multiple
predicates, faster class/cell grouping and metadata paths, dedicated read-bound kernels,
and the still-allocation-heavy exact junction enumerator. A bounded geometry cache
would need separate measurement and must distinguish omitted geometry from literal
alignments. Parallel chunk decoding is also possible, but initial measurements put
only about 6 ms there versus 18–38 ms in evaluation/aggregation. A cost-based parallel
planner, general common-subexpression elimination, streaming top-k, and whole-genome
benchmarks are not claimed as completed by this round.
