# Gravlax 0.2.0 release assessment

Date: 2026-09-05. This supersedes earlier branch recommendations about default
parallel decoding and allocator selection. Historical benchmark reports retain
their original measurements; they are not measurements of every final default.

## Decision and scope

Release GQ v1 and its optimized/reference executors, Python interface, examples,
documentation and diagnostics. Include the compatible reader/rANS/access
improvements, richer built-in co-occurrence predicates and compiled large-panel
matching. Include optional access indexes, finer chunks, geometry fidelity and
compression tuning. Keep structural/sparse/path compression experiments under
`aie dev`; do not introduce their experimental encodings as production codecs.

Keep serial GQ decoding and the system allocator as defaults. `--parallel-decode`
and the `mimalloc` Cargo feature remain explicit choices. Do not expand the language
with joins or general graph operations as a prerequisite for its first release.

## Final allocator reassessment

Compared separately built, source-matched system and mimalloc-v2 executables with
one warmup and nine randomized paired measurements per case. GQ measurements use
post-parse computation (preparation, execution and serialization), excluding source
parsing, process startup and publication. Native query, replay and ingest use whole
process wall time. Peak RSS covers the complete child process. All full typed
results/work summaries, native result data, replay file hashes and ingest bytes
matched before any report compaction. Large population lists are represented by
their count and digest in the reports, not omitted from equality checks.

Data: compact and fidelity CD45 archives; a real SEZ donor and the immutable
eight-donor demo federation (three selected chromosomes, 71,813 records and 24 chunks
across donors); and an explicitly synthetic 600,000-record, 60,000-cell,
three-chromosome archive with grouped multimappers. The demo is not a whole atlas.
The synthetic fixture is a scale test, not biological evidence.

Four-thread results favored mimalloc for per-cell groups and read support: roughly
35–45% lower times at about 9–32 MiB extra RSS. However, the normal 24-thread
configuration exposed a substantially larger memory cost:

| Workload, 24 threads | System → mimalloc ms | System → mimalloc MiB peak RSS |
|---|---:|---:|
| Eight-donor geometry summary | 27.7 → 25.0 | 18.7 → 73.6 |
| Eight-donor per-cell summary | 279.6 → 176.2 | 113.1 → 186.6 |
| Eight-donor read support/top-k | 78.8 → 46.0 | 91.7 → 144.6 |
| Synthetic geometry summary | 164.5 → 168.8 | 103.4 → 261.9 |
| Synthetic class summary | 292.8 → 267.7 | 137.2 → 267.9 |
| Synthetic per-cell summary | 1214.9 → 650.3 | 436.0 → 533.7 |
| Synthetic read support/top-k | 635.9 → 340.6 | 735.4 → 813.9 |
| CD45 native co-occurrence | 12.37 → 12.10 | 19.7 → 90.3 |
| CD45 replay | 74.35 → 73.44 | 66.6 → 157.8 |
| CD45 default ingest | 345.6 → 339.7 | 934.0 → 967.5 |

The allocation-heavy wins are reliable within these runs (paired bootstrap
intervals exclude parity), but memory does not increase moderately across the
workload mix. Native co-occurrence and replay gains are not statistically clear.
Turning off transparent huge pages with `MIMALLOC_ALLOW_THP=0` reduced RSS, but the
synthetic geometry case became about 22% slower (159.9 → 195.5 ms) and still used
about 90 MiB extra. Synthetic class aggregation also regressed slightly. This is
not a sufficiently consistent improvement to change the universal default.

Raw timings, paired confidence intervals, exact query text, binary hashes and
source identities:

- [Four-thread mixed operations](gq-release-allocator-benchmark.json).
- [24-thread mixed operations and synthetic scale](gq-release-allocator-scale-benchmark.json).
- [24-thread huge-page-off sensitivity test](gq-release-allocator-no-thp.json).
- Reproducer: `scripts/benchmark-release-allocator.py`; export the synthetic fixture
  with the ignored `gq::tests::export_release_scale_fixture` test and
  `GRAVLAX_GQ_SCALE_OUTPUT` set to a new output path.

These are warm-cache measurements on one shared Linux host, not cross-platform or
long-running fragmentation guarantees. Compare paired allocators within each run,
not absolute values across runs. No allocator tuning is enabled in production.

## Compatibility and release hardening

Default container/schema versions, evidence reduction, codecs and chunking stay
unchanged. A fresh default CD45 ingest produced the same 340,008-byte archive size
as 0.1.6; all 13 non-provenance compressed sections were byte-identical. The only
section change was producer-version provenance (0.1.6 → 0.2.0). Thus whole-file bytes
and content roots change, but default evidence/encoding does not. The old reader
successfully verified the new default archive, including payloads. Optional new
codecs/fidelity declarations can require the new reader and fail closed otherwise.

Release hardening includes deterministic parser mutations and oversized/deep input
tests; explicit serial/parallel/reference parity; and offline binary acceptance on
an archive and a two-member federation, including no-overwrite publication. PR CI
now exercises GQ on five native release targets, both allocators on Windows/musl,
and source, Python, documentation and container distribution checks. Extracted
Linux release archives also execute GQ, rather than only printing a version.

Local validation: 446 Rust tests pass separately with all features and with default
features (three intentionally ignored exporter/fixture tests), strict all-feature
Clippy passes, 37 packaging tests and 91 Python client tests pass, documentation
builds, and release version/package/public-source validation passes.
Platform CI remains a merge gate, not a claim of local
cross-platform execution.

[Broader acceptance evidence](gq-release-acceptance.json) records full payload
verification of all eight real capsule archives and agreement across auto serial,
auto parallel and reference execution for 20 query/data combinations, including
the synthetic scale cases. Reproduce with `scripts/validate-gq-release.py` against
the scale benchmark report. The existing row budget correctly rejects a
600,000-row support projection at a 500,000-row limit even if followed by `take 20`;
acceptance uses an explicit one-million-row limit.

## Merge and release ordering

Commit reviewed branch work, integrate the current main-branch logo without
overwriting it, and submit a PR into main. Merge only after all required checks
and reviews succeed; do not bypass protections. Prepare an annotated `v0.2.0` tag
from the validated main commit using `scripts/bump-and-publish`, then push through
the existing immutable-release workflow. PR readiness is not publication, and
platform/registry failures must remain visible rather than being waived.
