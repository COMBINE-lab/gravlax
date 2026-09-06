---
title: GQ query language
description: Compose typed evidence queries over archives, raw-UMI classes, cells, and federations.
---

GQ adds an expressive textual interface without changing the archive format. Start
with an explicit evidence unit, select a witnessed population, derive observations,
and preserve their joint true/false/unknown states.

```text
header { gq = 1, assembly = "GRCh38" }
let locus = g[1:198600000..198800000]
let RA = overlaps(g[1:198696711..198696909:-], min: 12bp)
let RO = j[1:198692373..198703297:-]
from @pbmc.records
|> within locus
|> derive {
     ra = any unique { RA },
     ro = any unique { RO },
     joint = any unique { RA & RO }
   }
|> tally {ra, ro, joint} by {sample}
```

Save this as `cd45.gq`, then:

```sh
aie gq validate cd45.gq
aie gq explain cd45.gq --bind pbmc=sample.aie
aie gq run cd45.gq --bind pbmc=sample.aie --output cd45.json
```

Unbound validation reports syntax/type checking only. Bound validation also checks
resource capabilities and route feasibility. An archive without a complete access
index may require `--allow-full-scan`. This permission is an execution option, not
query syntax. Outputs publish atomically and never overwrite existing files.

## Choose the scientific question explicitly

- `any unique {A & B}` asks for one uniquely mapped observation satisfying both.
- `any unique {A} & any unique {B}` permits separate observations of the record.
- On `.classes` or `.cells`, use `any record { ... }` to inspect complete children.
- `any stored {P}` inspects literal representatives and is diagnostic.
- `any multimap {all alternative {P}}` preserves grouped alternatives, including primary.

`all` is true on an empty domain. Use `nonempty(unique) & all unique {P}` when
nonempty evidence is required. A primary multimapper alignment is never unique.
Terminal-tail predicates are record-scoped, not alignment-scoped.

Compact archives can leave extent-sensitive questions unknown. GQ uses sound start
and splice-geometry proofs; it never assumes that the two retained ends bound every
omitted end. `where` keeps true and reports dropped unknowns; `tally` retains them.

## Counts and summaries

`reads.unique`, `reads.multimapping`, and `reads.total` count accepted observations,
not physical molecules. A multimapping signature contributes its weight once.
`support_reads(P)` returns lower/upper bounds and an exact/bounded status; `.lower`,
`.upper`, and `.status` expose individual fields.

```text
from @pbmc.records
|> within locus
|> where reads.total >= 2u64
|> summarize {records = count(), reads = sum(reads.total)} by {sample}
```

The read sum is the selected records' total weight, not the number of reads matching
a prior geometry filter. `fraction(n, denominator: d)` names its denominator and
returns null for a zero denominator. `select` projects values; `sort {field}` is
ascending; a final `take N` requires sort and discloses result truncation.

## Strand, annotations, and functions

Literal coordinates are 0-based, half-open; literal strand is alignment strand.
`tx_strand(-)` requires `library="same"` or `library="opposite"` in the header.
Annotation strand uses that same conversion. There are no implicit contig aliases.

`gene(@anno,"PTPRC")` has `.span`, union `.exons`, and `.junctions`.
`transcript(@anno,"ENST...")` also has `.junction_path`. Bind annotations through
a project with explicit assembly and annotation label; ambiguous names fail.
Use `--project PATH` to reuse project resources.

```text
fn exon_support(exon, minimum) = overlaps(exon, min: minimum)
export fn supported(exon: Region<GRCh38,alignment>, minimum: Bases)
  -> Predicate<Alignment,GRCh38,alignment> = overlaps(exon, min: minimum)
```

Functions are bounded and nonrecursive. `>>` matches consecutive junctions; `~>`
matches ordered subsequences. Unstranded paths require `path(order: genomic)`;
transcriptional paths can use `path(-)` or `path(+)` with a library rule.

## Metadata and federations

Bind a JSON federation with `--bind cohort=cohort.json`:

```json
{
  "schema_version": 1,
  "assembly": "GRCh38",
  "archives": [
    {"sample": "donor1", "path": "donor1.aie", "metadata": {"donor": "D1"}},
    {"sample": "donor2", "path": "donor2.aie", "metadata": {"donor": "D2"}}
  ]
}
```

Duplicate sample names, files, or content roots are rejected. Cells/classes are
sample-namespaced; changing federation membership never joins UMIs across samples.

`--metadata metadata.json` binds typed columns and optional cell values:

```json
{
  "columns": {
    "donor": {"type": "string"},
    "cell.type": {"type": "string", "optional": true}
  },
  "cells": {"donor1": {"AAAAAAAAAAAAAAAA": {"cell.type": "Astro"}}}
}
```

Comparisons, `in`, `not in`, arithmetic, parentheses, and Boolean composition are
typed without implicit conversions. Null comparisons yield unknown; grouping retains
a missing-value group. Grouping uses base sample/cell metadata so denominators remain
defined before evidence filters. Exact source denominators by cell metadata can require
an explicitly permitted full scan for record/class queries; explain reports that cost.

## Junction enumeration

```text
header {gq=1, assembly="GRCh38"}
from junctions(@pbmc, within: g[1:198600000..198800000])
|> support(unit: class, by: {sample})
```

The initial enumerator counts unique-evidence junction support exactly, in one bounded
full scan. Both boundaries must be in a target interval. It does not report unverified
catalogue counts or treat alternative-only mappings as unique evidence.

## Python

```python
from gravlax import Client

client = Client()
client.gq_validate("cd45.gq")
result = client.gq_run("cd45.gq", bindings={"pbmc": "sample.aie"})
table = result.table("results")
```

`gq_explain` inspects the plan; `gq_run_to_file` avoids materializing rows in Python.
GQ uses `gravlax.gq.result.v1` under the shared typed result envelope. The specialized
`aie query` commands remain available and can be faster for their fixed workflows.

## Execution planning and profiling

GQ selects its execution strategy automatically. `gq explain` reports it under
`physical_strategy`. Eligible filter/Truth-derive/tally or summarize pipelines are
fused: derived Truth values occupy compact slots, grouping uses interned identifiers,
and sparse accumulators avoid constructing a result row for every evidence unit.
Counts use checked integer updates; distinct counts use bounded hash sets. Numeric
sums retain source and unit order, including across federation members.

Compatible single-interval overlap/start/end and single-junction predicates reuse
the built-in shape matcher. A fixed 1,024-slot placement cache is shared across the
fused plan within each archive member; collisions replace entries and never change
answers. Extent-sensitive omitted geometry bypasses this cache and retains GQ's
conservative proofs. Union/total-overlap and more general patterns retain their GQ
evaluators. Empty and sample-only grouping reuse one key per member, while keeping
the reference engine's logical grouping-step charges for each encountered cell.

Chunk decoding is serial by default. `--parallel-decode` opts into bounded parallel
windows (up to twice the Rayon thread count), followed by attachment and evaluation
in original order. It can use more temporary memory and is not consistently faster
on small archives. Chunk/record budgets are checked before pending decode work.
`RAYON_NUM_THREADS` controls concurrency; floating-point aggregation stays ordered.
The Python client exposes the same option as `parallel_decode=True`.

Record queries iterate selected records directly instead of constructing singleton
unit groups. Unused unit identifiers are not formatted; queries without metadata or
cell-field references can reuse a sample row. An adjacent `sort |> take` uses stable
top-k selection. All evidence is still evaluated, and row/work budgets still apply
before pagination: top-k is not an early evidence cutoff or a streaming memory bound.

Unsupported fused shapes, including non-Truth derived intermediates such as read
support bounds, use the general interpreter. Direct record iteration and top-k can
still improve those queries. Junction enumeration retains its separate executor.
There is no predicate pushdown across quantifier or unit-closure boundaries, no
full-scan permission change, and no archive-format change.

```sh
aie gq run cd45.gq --bind pbmc=sample.aie --profile
aie gq run cd45.gq --bind pbmc=sample.aie --engine reference --profile
```

`--engine reference` disables physical fusion, projection pruning, direct record
iteration, member-constant grouping, parallel decoding, placement caching and top-k
selection; earlier shared low-level geometry/aggregation improvements
remain active. This is a diagnostic comparison path, not a different query language.
Both paths preserve logical expression-step budgets, omitted-geometry diagnostic
counts, complete class/cell closure, population denominators and three-valued results.

`--profile` emits one `gq_profile=...` JSON line to stderr, leaving the typed result on
stdout. It excludes parsing and reports resource preparation/type checking, execution,
serialization and publication separately. Regular evidence queries additionally report
preflight planning, archive setup, chunk reading/decoding, universe witnessing, unit
and denominator preparation, fused evaluation/aggregation, and final projection/sort.
Evaluation and aggregation share a loop, so their time is reported together rather
than assigned artificially to separate phases. Timings do not affect plan digests.

The repository's `scripts/benchmark-gq-execution.py` runs interleaved comparisons on
compact and fidelity archives. It checks complete typed results and semantic/work
summaries on every run and uses internal execution timing—not total subprocess time.

### Optional allocator

The executable can additionally use mimalloc for Rust allocations:

```sh
cargo build --release --locked -p gravlax --bin aie --features mimalloc
```

This selects mimalloc's maintained v2 backend through wrapper version 0.1.52. The
system allocator remains the default; use `--no-default-features` without adding
`--features mimalloc` for an explicit system build. `gq run --profile` identifies the
allocator in stderr's profile JSON. This choice does not affect archive encoding.

The GQ comparison found lower latency with mimalloc, particularly for distinct counts,
per-cell groups and read-support projections, but also higher peak memory. The v2
backend used less memory than v3 on these fixtures with similar speed. It is therefore
an opt-in tradeoff, not an unconditional recommendation for memory-constrained jobs.
The final eight-donor and synthetic-scale reassessment confirmed large gains for
allocation-heavy queries, but found much higher memory overhead at 24 threads and
little gain for some lightweight/native operations. The system default is retained.
Only the executable selects the allocator: reusable Rust libraries do not impose it
on their callers, and native C-library `malloc`/`free` calls are not overridden.

`scripts/benchmark-gq-allocators.py` compares separately built binaries using post-parse
execution timing and peak process RSS. It supports fixed CPU affinity and verifies
complete typed results and semantic/work summaries on every run.
