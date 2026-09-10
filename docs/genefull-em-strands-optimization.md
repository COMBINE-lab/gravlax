# GeneFull EM strands and optimization validation

This follow-up extends the implementation validated in `genefull-em-validation.md`.
The original manuscript results and input archives were preserved. Host:
`nomad00.umiacs.umd.edu`; Rust 1.89.0; STAR 2.7.11b. The baseline source is
`7078074` on `codex/genefull-replay`.

## Counting models and strand policies

`aie dev em --solo-strand forward|reverse|unstranded` now applies to packed
recovery EM, its eager reference, and STAR-style EM. `--gene-full --star` is
supported. Both the unique-count baseline and the ambiguous candidate sets use
the requested counting model and strand. Forward remains the default; it means
the cDNA alignment and gene have the same strand, including genes on the minus
strand. Reverse means opposite strands, and unstranded accepts both.

No archive or compiled-annotation format change is required. The biological
experiment reads the existing v1 `d2p.aie` archive directly. The GeneFull index
is derived from the annotation supplied at query time.

## Independent STARsolo comparison

`scripts/validate-em-strands.py` creates synthetic reads containing intronic,
nested-gene, antisense and genuine multiple-placement evidence. Each of the six
Gene/GeneFull × strand combinations compares the STAR-style EM contribution
against STARsolo on the same alignments, checks counting-model/strand metadata,
compares packed and eager unmasked pooled emission, and verifies that selecting
only pooled EM reproduces its complete-evaluation metrics and coverage counts.

All six combinations passed. The largest per-entry difference in the EM
increment was 0.000004 counts, within the programs' different text-output
rounding. Complete GeneFull unique-plus-EM matrices agreed within 0.000003 counts.

The Gene fixture deliberately includes equal unique-gene support for a UMI.
Existing replay keeps the lower gene id while STARsolo removes that UMI. This
produces one extra A count under forward or reverse assignment, and two under
unstranded assignment. The test reports this separately and checks the EM
increment independently; it does not label the complete Gene matrices identical.
No default Gene UMI-filtering or tie-breaking rule was changed. These synthetic
agreements do not establish universal equivalence after archive representation
or UMI-collapse differences on arbitrary input.

Rust tests additionally check strand-flipped alternative placements, missing
annotation contigs, reverse-orientation symmetry, unique-plus-ambiguous classes,
and STAR's finite convergence threshold. The GeneFull interval index is checked
against an exhaustive exon-bound overlap oracle for 6,000 query/strand cases.

## Scoring scope, mode selection and memory

`--eval-barcodes` sets the evaluation population independently of group
assignments. It neither removes other archive barcodes from prior fitting nor
enables extra inference models. If supplied with `--groups`, it overrides only
the scored population. `--modes` selects a subset of the existing nine modes,
retaining canonical execution order and the same ten updates. Paired metrics
are produced only when both participating modes are selected. Production
emission still uses pooled EM, which must be included in an explicit mode list.

`--support-memory-mib` defaults to 512 MiB and limits the retained compact
support arrays. Actual support volume triggers a switch to temporary shards;
zero forces spilling. It is not a whole-process memory cap: decoded batches,
worker buffers, spill buffers and the final packed EM state are additional
working sets. Tests force a switch partway through classes split across batches,
including deliberately underestimated input-record counts. Every packed field
matches the in-memory path, and temporary shards are removed after completion.

Metrics record the counting model, strand, evaluation-barcode source, requested
modes and support-storage scope. Their accuracy denominator remains evaluable
raw UMI classes before one-mismatch collapse. Replay matrices count collapsed
UMIs; STAR-style matrices add fractional multi-only raw-class counts to the
unique replay baseline.

## Brain benchmark protocol

`scripts/benchmark-genefull-em.py` runs three interleaved repetitions with 16
Rayon threads. The archive receives a read-through warmup. It records elapsed
time, user/system CPU time, peak process RSS, commands and binary checksums.
Intermediate binaries isolate mode selection, actual-support spilling, and the
direct block/strand index. Comparisons verify all selected metrics and coverage
counts against the nine-mode baseline, and all three replay MEX files by SHA-256.

EM scores the same 6,460 original Gene-called nuclei. Its pooled prior still
uses every archive barcode after masking. Replay emits the full raw barcode
list. All runs use the same archived alignments and v49 compiled annotation.
No nucleus-calling change is part of this benchmark.

The fixed-nucleus denominator is unchanged: 231,351 masked classes, including
74,467 truth-lost exclusions and 156,884 evaluable classes. This preserves the
manuscript's conditional pooled accuracy of 74.8999% and all its other shared
EM metrics. Existing matrix results are preserved separately from the new runs.

## Complete-command results

Medians of three interleaved repetitions; peak RSS is the peak of the complete
process, expressed in GiB (2^30 bytes). Every result equality gate passed.

| Configuration | Elapsed seconds | Peak RSS, GiB |
|---|---:|---:|
| Previous implementation, nine modes | 65.60 | 2.769 |
| Updated implementation, nine modes | 62.87 | 2.897 |
| Four modes, previous overlap/spill implementation | 32.46 | 2.741 |
| Four modes, actual-support spilling | 32.03 | 2.936 |
| Four modes, actual-support spilling and direct strand index | 32.83 | 2.861 |
| Updated implementation, pooled only | 23.24 | 2.890 |
| Previous GeneFull replay | 8.05 | 5.171 |
| Updated GeneFull replay | 7.95 | 5.184 |

Selecting the four manuscript models reduces elapsed time by 50.0% relative to
running all nine models in the previous implementation (2.00-fold throughput).
Selecting only pooled reduces elapsed time by 64.6% (2.82-fold throughput).
These comparisons explicitly change the amount of requested work; each retained
model's scientific result is identical. The same nine-mode workload is 4.2%
faster in these runs.

Actual-support spilling adds a measured peak-memory cost on this brain archive:
with four modes, the median rises from 2.741 to 2.936 GiB before the lookup change.
It enforces a candidate-volume budget when the number of candidates per archive
record is underestimated; it should not be described as a demonstrated reduction
in whole-process memory. The measured peak retained support capacity before
spilling was 529,997,216 bytes, within the 536,870,912-byte default budget.

The overlap-index change did not demonstrate a consistent whole-EM speedup:
four-mode medians were 32.03 seconds before and 32.83 seconds after. Replay
medians were 8.05 and 7.95 seconds, with overlapping individual timing ranges.
These small differences should not support a broad throughput claim.

## Isolated overlap lookup

The `gravlax-anno` example `genefull_index_benchmark` checks 200,000 deterministic
gene-local block queries per strand against the previous mixed-strand index.
It then interleaves three repetitions of 1,000,000 queries per implementation
and strand, including single-block and spliced geometries. Every candidate set
agrees; timings exclude annotation parsing, archive decoding, and EM updates.

| Strand policy | Previous lookup, seconds | Direct strand lookup, seconds | Speed ratio |
|---|---:|---:|---:|
| Forward | 0.0822 | 0.0466 | 1.77× |
| Reverse | 0.0825 | 0.0466 | 1.77× |
| Unstranded | 0.0762 | 0.0774 | 0.98× |

The index change is retained for its measured strand-specific lookup saving and
smaller per-gene index entries. It does not establish a material whole-EM speedup.

## Final checks and preserved results

The final Rust workspace suite passed 459 tests with three pre-existing ignored
tests. Clippy passed for all workspace targets with warnings denied. The 32-page
documentation site built successfully using Node 22.23.2. Changed Rust regions
were formatted, and `git diff --check` passed; unrelated historical formatting
was left intact.

A final binary rerun preserves all nine default Gene EM metric records, paired
comparisons and coverage counters exactly against the saved manuscript control.
Default Gene replay preserves all three raw MEX files byte-for-byte. The legacy
forward-strand Gene STAR-style fixture is also byte-identical to the previous
binary. The final binary passes all six EM model/strand combinations and all
24 existing STARsolo replay/representation comparisons.

A separate full-brain GeneFull control forced support spilling with a zero-byte
retained budget. All four scientific metric records and coverage counters match
the original nine-mode baseline exactly; reported peak retained support capacity
is zero. Its single observed run took 35.75 seconds and 2.789 GiB peak RSS. This
is a storage-path correctness check and an illustrative memory/time tradeoff,
not a repeated comparative benchmark.

Commands, timing logs, source/binary and input identities, raw metrics, synthetic
fixtures and MEX checksums are saved under the project run directory
`runs/genefull-em-strands-opt-20260910`. The main benchmark has 24 successful
complete-command runs. No biological ground-truth accuracy or changed nucleus
calling is inferred from these performance experiments.
