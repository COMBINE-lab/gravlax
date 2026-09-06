# GQ native-kernel and aggregation optimization

Release follow-up: parallel decoding is now opt-in (`--parallel-decode`), following
the mixed speed and increased memory measurements below. Predicate and grouping
optimizations remain automatic. These historical measurements used parallel decode;
the release allocator assessment measures the new serial default.

## Scope and implementation

This is a second optimization round on `astra-improvements`, after physical GQ
fusion and the optional allocator experiment. The comparison baseline is the
system-allocator executable frozen immediately before this round. No archive
format, query syntax, scientific semantics, default allocator, or release version
changes are introduced here.

The automatic engine now:

1. Lowers single-interval, non-total overlap/start/end predicates and
   single-junction predicates with symmetric tolerance to the built-in shape
   matcher. This is **placement-predicate reuse**, not delegation of an entire GQ
   query to the legacy co-occurrence executor. Union/total overlap, asymmetric
   tolerance, and general paths retain the GQ evaluator.
2. Uses allocation-free full-shape validation and a shared, direct-mapped
   1,024-slot placement-result cache per fused plan/member. Full keys are compared,
   so collisions only cause recomputation. Cache storage is fixed, initialized
   lazily, and does not grow with the archive. Omitted extent-sensitive geometry
   always bypasses the cache; GQ's conservative proofs remain authoritative.
3. Reads, checks compressed payload integrity, decompresses and decodes chunks in
   windows bounded by twice the Rayon thread count. Results remain in route order.
   Cell/terminal attachment and all evaluation/aggregation remain ordered. Pending
   chunk/record budgets are checked before launching decode work.
4. Reuses one empty/sample-only aggregation key per member instead of creating a
   key and group-ID mapping for each cell. Logical grouping evaluation steps are
   still charged once for each encountered cell, as in the reference engine.
   Cell metadata validation, source/population denominators, and floating-point
   accumulation order are unchanged. Cell/metadata grouping retains the general
   path.

`--engine reference` disables these new paths along with the earlier structural
optimizations. `gq explain` describes the selected strategy. No new permission to
scan an archive is implied, and no quantifier/domain is collapsed into another.

## Rejected experiment

The initial implementation used a 65,536-entry hash table for placement-result
memoization while retaining the previous grouping machinery. It made the compact
matched queries approximately 10–18% slower and used more memory. It was replaced
with the small direct-mapped cache and allocation-free validation before the final
measurement. The initial raw report is retained separately. These experiments do
not isolate the contribution of every final component; do not attribute the whole
speedup to the shape matcher or parallel decoding alone.

## Benchmark method

The reproducible harness is `scripts/benchmark-gq-native.py`. It compares before
and after executables on ten query shapes over compact and fidelity archives.
Three shapes also run equivalent built-in co-occurrence commands: record-level
conjunction, any-placement disjunction, and nonempty all-placement disjunction.
Both archives have 14 chunks and 54,173 decoded records; the matched built-in
population contains 53,936 uniquely mapped evidence records.

Each case uses one warmup, then randomized interleaving of measured runs. Both
executables use the system allocator and four Rayon threads, pinned to the same
four allowed logical CPUs. The filesystem cache is warm. The host is heavily
loaded, so absolute milliseconds should not be compared to earlier quiet-host
reports. This is a small, two-archive benchmark, not a claim about all datasets or
storage systems.

GQ uses `post_parse_compute_ms`: binding/checking/planning, archive access,
execution and serialization, excluding source reading, parsing, startup and
publication. Built-ins use complete process wall time, including startup,
argument parsing and output. **This asymmetry favors GQ.** Similar measured
numbers are not proof of equal internal execution speed. GNU time supplies
whole-process peak RSS for both GQ builds.

Complete typed GQ results, uncertainty counts, denominators and logical work
counters are compared exactly on every repetition, excluding physical strategy
text and timings. Built-in true/false/unknown totals are compared on every matched
run. Binary hashes, archive roots, queries, raw repetitions, phase times and RSS
are recorded in the final JSON report.

## Results

The final measurement used 11 repetitions per combination: 552 invocations
including warmups, with every result comparison passing. Host load averages were
approximately 159–160. Final median milliseconds:

| Archive | Matched query | Previous GQ, post-parse | New GQ, post-parse | Built-in, whole process | GQ time reduction |
|---|---|---:|---:|---:|---:|
| Compact | Record conjunction | 98.14 | 53.23 | 54.41 | 45.8% |
| Compact | Same observation | 107.69 | 63.22 | 57.83 | 41.3% |
| Compact | Nonempty universal | 108.28 | 60.77 | 58.72 | 43.9% |
| Fidelity | Record conjunction | 104.93 | 62.09 | 58.78 | 40.8% |
| Fidelity | Same observation | 119.39 | 70.31 | 62.20 | 41.1% |
| Fidelity | Nonempty universal | 116.76 | 67.04 | 62.42 | 42.6% |

The measured GQ/built-in ratio is now 0.98–1.13, **using the asymmetric clocks
described above**. For a symmetric whole-process comparison, new GQ is still
1.23–1.37 times the built-in duration on these six cases. Thus the result closes
much of the gap, but does not establish equal internal executor speed.

Other post-parse time reductions across compact/fidelity archives:

| Query | Compact | Fidelity |
|---|---:|---:|
| Geometry tally | 50.1% | 49.7% |
| Truth/read summary | 51.5% | 46.3% |
| Class summary | 35.5% | 26.1% |
| Distinct cells by sample | 27.7% | 34.2% |

Cell-grouped summaries were essentially unchanged (+1.0% / −0.2% reduction).
Read-support fallback and junction enumeration were 0.5–2.5% slower in this run;
the earlier pilot varied around zero. There is no demonstrated speed win for
those paths. Parallel decoding adds transient memory even when the main cost is
elsewhere: matched queries added approximately 1.8–5.3 MiB peak RSS, while the
unchanged/fallback workloads added 2.6–6.9 MiB. Sample-grouped geometry, summary,
class and distinct queries instead saved approximately 3.3–8.6 MiB by eliminating
redundant per-cell grouping structures. The allocator was system in both builds.

For same-observation queries, evaluation/aggregation fell from 68.93 to 22.24 ms
(compact) and 71.22 to 26.76 ms (fidelity). Read/decode/attachment changed from
18.04 to 20.03 ms and 20.96 to 18.10 ms respectively. This strongly localizes the
overall win to evaluation/grouping; it does **not** demonstrate a consistent
decode-only speedup on this small, busy-host benchmark.

Raw artifacts:

- `gq-native-execution-benchmark.json`: final 11-repetition report with archive roots.
- `gq-native-execution-pilot.json`: revised implementation's 5-repetition pilot.
- `gq-native-execution-initial.json`: rejected large-cache experiment.

## Correctness and build validation

- Workspace tests: 445 passed, one intentionally ignored parser benchmark.
- Strict workspace/all-targets Clippy: passed.
- Opt-in mimalloc GQ tests: 24 passed, one intentionally ignored parser benchmark.
- Documentation build: all 32 pages passed. Benchmark harness syntax and
  `git diff --check` also passed.
- New tests cover native/interpreter predicate agreement, repeated placement
  evaluation, omitted-geometry cache exclusion, bounded cache storage and
  collisions, malformed/overflowing shapes, and ordered one/four-thread decoding.
- Existing differential tests retain exact logical-step counts and fail at the
  same exhausted budgets; federation tests retain non-reassociated numeric sums.

## Reproduction

Freeze the pre-change executable before rebuilding, then run:

```sh
cargo build --release --locked -p gravlax --bin aie --no-default-features
python3 scripts/benchmark-gq-native.py \
  --before /tmp/gravlax-gq-before-native \
  --after target/release/aie \
  --archive /tmp/gravlax-astra-benchmark-v3/indexed.aie \
  --archive /tmp/gravlax-astra-benchmark-v3/fidelity.aie \
  --repeats 11 --cpus 1,2,7,9 \
  --output /tmp/gq-native-new-run.json
```

Use a new output filename and CPUs allowed on the target machine. The archive and
frozen binary paths above are local benchmark artifacts, not bundled datasets.
