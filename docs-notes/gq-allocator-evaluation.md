# Allocator evaluation after GQ optimization

Historical experiment. The broader final default decision and 24-thread/scale
results are in [the 0.2.0 release assessment](gq-release-readiness.md).

## Decision

Add **opt-in mimalloc using its v2 backend**, but keep the system allocator as the
default. Both mimalloc backends improved GQ latency; neither dominated the system
allocator on peak memory. v2 offered similar speed to v3 with lower RSS on these
fixtures, making it the better optional configuration.

```sh
cargo build --release --locked -p gravlax --bin aie --features mimalloc
```

Explicit system build:

```sh
cargo build --release --locked -p gravlax --bin aie --no-default-features
```

The dependency is `mimalloc = 0.1.52` with `default-features = false` and `v2` enabled;
the lockfile selects `libmimalloc-sys = 0.1.49`, bundling native mimalloc 2.3.2 for that
backend. The comparison also tested its bundled v3.3.2 backend. The only new dependency
packages are the wrapper and native binding; the native build uses the existing `cc`
dependency. Rust 1.89 remains sufficient.

Allocator selection is in the executable, not the reusable library crates. There is
no `override` feature or `LD_PRELOAD`: native C-library allocation calls keep their own
allocator. No archive encoding, fidelity, index or query semantics change. GQ profiles
identify `system` versus `mimalloc`. No commit, push or release was performed.

## Candidate assessment

The [mimalloc Rust wrapper](https://github.com/purpleprotocol/mimalloc_rust) exposes both
backends and bundles the native build. Its local source/build integration made it a
straightforward first experiment. [Upstream mimalloc](https://github.com/microsoft/mimalloc)
documents the maintained v2/v3 lines.

[Snmalloc](https://github.com/microsoft/snmalloc) is particularly interesting for
cross-thread frees and large deallocation batches. That remains a credible candidate
for the parallel ingest/replay paths, but it was **not benchmarked** in this round.
The present result is a system/mimalloc-v2/mimalloc-v3 comparison, not evidence that
mimalloc outperforms snmalloc. GQ's current evaluation loop is serial.

## Method

- Seven GQ workloads on each of the compact and full-fidelity CD45 archives: 14 cases.
- Same already-optimized GQ engine for all allocator binaries; parsing excluded.
- One warmup and nine measured repetitions, randomized allocator order for each case.
- CPU affinity fixed to logical CPU 1; `RAYON_NUM_THREADS=4`; warm filesystem cache.
- Internal execution wall-clock timing, separately reported serialization and total
  post-parse computation, and GNU time's whole-child peak RSS.
- Complete typed tables and semantic/work summaries checked on every invocation,
  including population denominators and expression-step/unknown counters. Only timing
  fields are excluded. All 420 invocations in the three-way comparison matched.
- No allocator environment tuning, forced collection or C malloc interposition.

The host was shared and busy: `nomad01.umiacs.umd.edu`, AMD EPYC 9555, 128 logical CPUs,
load averages about 99 during the three-way comparison. Affinity does not isolate the
physical core or memory bandwidth. Compare allocators **within** each case, not absolute
times between archives or between benchmark runs. These short-lived, locus-scale
queries do not establish long-running fragmentation or whole-dataset behavior.

## Results: system versus selected mimalloc v2

All times below are median internal execution milliseconds. Memory is median peak
process RSS in MiB. “Lower time” is `1 - mimalloc_time / system_time`.

| Query | Archive | System → v2 ms | Lower time | System → v2 MiB |
|---|---|---:|---:|---:|
| Geometry tally | compact | 116.55 → 107.04 | 8.2% | 34.1 → 46.0 |
| Truth/read summary | compact | 84.52 → 73.46 | 13.1% | 34.5 → 45.0 |
| Class summary | compact | 121.11 → 104.36 | 13.8% | 33.8 → 45.0 |
| Distinct cells | compact | 102.14 → 76.71 | 24.9% | 48.8 → 59.0 |
| Per-cell groups | compact | 253.36 → 184.21 | 27.3% | 77.8 → 91.0 |
| Read support + top-k | compact | 102.99 → 62.07 | 39.7% | 79.1 → 89.0 |
| Junction enumeration | compact | 72.29 → 59.27 | 18.0% | 22.7 → 39.0 |
| Geometry tally | fidelity | 95.77 → 85.75 | 10.5% | 37.9 → 49.0 |
| Truth/read summary | fidelity | 68.53 → 58.87 | 14.1% | 37.9 → 49.0 |
| Class summary | fidelity | 97.16 → 82.81 | 14.8% | 37.5 → 49.0 |
| Distinct cells | fidelity | 154.97 → 124.40 | 19.7% | 52.2 → 63.0 |
| Per-cell groups | fidelity | 334.25 → 261.75 | 21.7% | 81.8 → 93.0 |
| Read support + top-k | fidelity | 151.46 → 88.73 | 41.4% | 82.4 → 93.0 |
| Junction enumeration | fidelity | 116.24 → 95.66 | 17.7% | 26.5 → 41.0 |

The selected backend costs roughly **10–16 MiB extra peak memory** on these GQ cases.
v3 used approximately **4–27 MiB more than v2**, with generally similar latency; the
fidelity enumeration case was slightly faster with v3. This is a useful speed option,
not a universal improvement or justification for silently replacing the default.

Source/binary hashes, individual samples, query text, source identities, result digests
and complete work summaries are in:

- [Three-way comparison](gq-allocator-v2-v3-benchmark.json).
- [Initial system/v3 experiment](gq-allocator-benchmark.json). This exploratory run
  overlapped local validation compilation as well as shared-host load; the three-way
  comparison was run after local compilation had finished.

The three-way binaries are preserved at:

- `/tmp/gravlax-gq-allocator-system`
- `/tmp/gravlax-gq-allocator-mimalloc` (v3 experiment)
- `/tmp/gravlax-gq-allocator-mimalloc-v2` (selected backend)

The system comparison executable predates only the allocator module/profile-name
addition; its query engine is the same optimized implementation. Exact binary hashes
are recorded in the reports. The ordinary `target/release/aie` is restored to a system
build at handoff, matching the unchanged default.

## Broader smoke checks

The same two allocator binaries also replayed the CD45 annotation and ingested the
complete-multiplicity source BAM. Replayed `matrix.mtx`, `features.tsv` and `barcodes.tsv`
were byte-identical. Newly ingested, access-indexed archives were byte-identical too.

Single-run observations, **not latency benchmarks**:

| Operation | System seconds / RSS KiB | v2 seconds / RSS KiB |
|---|---:|---:|
| Replay | 0.13 / 67,848 | 0.14 / 101,844 |
| Ingest | 0.62 / 41,864 | 0.61 / 80,860 |

These provide no evidence of a general ingest/replay speedup and reinforce keeping
the memory-conservative default. Scratch artifacts are in
`/tmp/gravlax-allocator-smoke.HtPFtN`; archive SHA-256 is
`d976cdcabbff8626aa3ddb609fc76460573990a689f61d5ea562fddfd87fe28d`.

## Validation and portability

- System and selected mimalloc-v2 configurations each passed **442 Rust tests**, with
  one existing manual benchmark ignored. The v3 experiment also passed the full suite.
- Added an allocator smoke test covering 4096-byte alignment, vector growth and
  shrinking, and allocation/freeing on different threads.
- Strict workspace clippy, **91 Python tests**, **37 packaging tests**, and the
  **32-page documentation build** passed. Packaging tests used the available Python
  3.13 interpreter; the host's default Python 3.9 lacks the tooling's `tomllib`.
- The portability workflow now includes system and opt-in mimalloc builds for static
  Linux musl and Windows MSVC. Those platform jobs have **not** been run locally;
  they are configured for CI. The actual allocator builds/tests here are Linux glibc.

## Reproduction

```sh
python scripts/benchmark-gq-allocators.py \
  --binary system=/tmp/gravlax-gq-allocator-system \
  --binary mimalloc-v3=/tmp/gravlax-gq-allocator-mimalloc \
  --binary mimalloc-v2=/tmp/gravlax-gq-allocator-mimalloc-v2 \
  --archive /tmp/gravlax-astra-benchmark-v3/indexed.aie \
  --archive /tmp/gravlax-astra-benchmark-v3/fidelity.aie \
  --cpu 1 --repeats 9 --output /tmp/gq-allocator-new-run.json
```

Use a new report path; the harness does not overwrite reports. For new builds, copy
the binary after each Cargo build so the subsequent configuration does not replace
the previous comparison executable. Producing the v3 experiment from the final source
requires removing the dependency's explicit `v2` selection; the final supported
`mimalloc` feature intentionally selects v2, not the wrapper's upstream v3 default.
