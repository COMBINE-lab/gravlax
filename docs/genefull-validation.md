# GeneFull validation

Validated 2026-09-10 with Rust 1.89.0 and STAR 2.7.11b
(`b1edc1208d91a53bf40ebae8669f71d50b994851`).

## Implementation and correctness

GeneFull indexes exon-derived gene spans separately by chromosome and strand.
A placement's candidate set is the union of spans overlapping its aligned
blocks on the accepted strand. Skipped introns alone do not create overlaps;
splice-junction concordance is not required. Existing unique-gene assignment,
best-gene selection and global UMI-class collapse are retained. No archive or
compiled-annotation format changes were made.

The indexed implementation agreed with an independent exhaustive exon-bound
and overlap calculation for 6,000 query/strand combinations. Explicit fixtures
cover half-open boundaries, introns, nested genes, disjoint isoforms, skipped
and novel introns, reused gene identifiers on distinct chromosomes/strands,
alternative placements, absent contigs, repeated class records and 1-mismatch
collapse across reducer batches.

`scripts/validate-genefull.py` runs a reproducible synthetic STARsolo experiment.
All 24 comparisons passed: Gene/GeneFull × three strand policies × streaming,
eager, compiled-annotation and direct-BAM replay. The corresponding MEX files
were byte-identical across equivalent representations, and report/MEX model
identity and incompatible flags were checked. A separate Python-client smoke
run produced a report and a readable MEX result.

The full Rust suite passed 456 tests (3 pre-existing ignored tests); Python
passed 92 tests; distribution tooling passed 43 tests. Clippy with warnings
as errors and the documentation build passed. Repository-wide `cargo fmt
--check` still flags pre-existing formatting outside the changed regions;
the base source was checked separately. New algorithm regions were formatted,
and `git diff --check` passed.

## Brain nuclei: fixed barcode set

The primary comparison uses the same 6,460 Gene-called nuclei from the
Brain_3p dataset. Raw output matrices were generated with the complete barcode
list and subsequently restricted to those nuclei. Alignments in the archive
were generated without the target annotation; reference matrices were produced
by matched annotation-aware STARsolo processing. These are end-to-end replay
comparisons, not a claim that the two alignment procedures are identical.

| Matrix | Collapsed UMIs in 6,460 nuclei | Median UMIs/nucleus | Median detected genes/nucleus |
|---|---:|---:|---:|
| STARsolo Gene | 24,280,804 | 2,900 | 1,838 |
| Replay Gene | 24,153,768 | 2,882.5 | 1,833 |
| STARsolo GeneFull | 70,815,492 | 8,353 | 3,513 |
| Replay GeneFull | 71,156,828 | 8,390.5 | 3,532 |

For Gene, the L1 difference is 220,486 UMI counts; half-L1 divided by the
STARsolo Gene total is 0.4540%. For GeneFull, L1 is 1,055,992; half-L1 divided
by the STARsolo GeneFull total is 0.7456%. Half-L1 is a matrix-distance
normalization, not a count of tracked molecules that physically moved.

Switching replay from Gene to GeneFull adds 48,631,181 counts and removes
1,628,121 across gene-by-nucleus entries, for a net gain of 47,003,060 UMIs
and a 2.94599-fold total. Nested/overlapping genes and UMI filtering mean that
some individual gene counts decrease even though the total rises.

Default Gene output matches the saved earlier matrix, features and barcodes
byte-for-byte. GeneFull replay from the full ingest BAM matches streaming
archive replay byte-for-byte for all three MEX components and all assignment
statistics. This tests the shared representation; compression can still differ
from per-read STAR geometry, and replay does not recompute genome alignments.

## Nucleus calling and downstream sensitivity

STAR 2.7.11b `soloCellFiltering` with default `EmptyDrops_CR` was applied to
each replay raw matrix separately:

| Matrix | Called nuclei | Collapsed UMIs within its own called set |
|---|---:|---:|
| STARsolo Gene | 6,460 | 24,280,804 |
| Replay Gene | 6,458 | 24,152,117 |
| STARsolo GeneFull | 6,564 | 70,923,992 |
| Replay GeneFull | 6,565 | 71,266,591 |

Replay Gene and GeneFull share 6,446 called barcodes: GeneFull adds 119 and
excludes 12 relative to Gene. Gene replay differs from the matching STARsolo
calling list by four omitted and two added barcodes. GeneFull retains all
6,564 STARsolo GeneFull calls and adds one. These calling contrasts are
separate from the fixed-nucleus comparison above.

For the fixed nuclei, a shared sensitivity analysis used log1p counts after
10,000-count library normalization, 2,000 genes ranked by pooled variance,
gene z-scores, 30 joint principal components and pooled k-means with 15
clusters (seed 2718). Adjusted Rand indices were 0.98561 for STARsolo/replay
Gene, 0.98048 for STARsolo/replay GeneFull and 0.33320 for replay Gene/GeneFull.
Thus the counting-model change affects clustering far more than replay does
under a fixed analysis. Cluster marker contrasts and marker detection counts
were saved with explicit nucleus denominators. This exploratory analysis is
not a test of cell-type ground truth or biological superiority.

## Annotation changes and the existing extension

All comparisons below use replay on the same 6,460 nuclei. Ensembl accessions
are matched after removing the version while preserving `_PAR_Y` suffixes.

| Change | Model | Before UMIs | After UMIs | Net gain / before total |
|---|---|---:|---:|---:|
| v32 → v49 | Gene | 22,642,472 | 24,153,768 | 6.6746% |
| v32 → v49 | GeneFull | 68,990,832 | 71,156,828 | 3.1395% |
| v49 → existing extended annotation | Gene | 24,153,768 | 29,415,412 | 21.7839% |
| v49 → existing extended annotation | GeneFull | 71,156,828 | 75,765,816 | 6.4772% |

The extension was reused as saved, without retuning. Its relative gain is
smaller under intron-inclusive counting. No matched v32 GeneFull STARsolo
matrix was available; the annotation-change table is an archive replay
sensitivity analysis, not an additional STARsolo fidelity claim.

## Assignment units: whole archive

| Quantity | Denominator | Gene numerator | GeneFull numerator |
|---|---:|---:|---:|
| Records with ≥1 uniquely assigned representative | 130,270,701 archive records | 30,672,329 | 81,557,405 |
| Uniquely assigned representatives | 184,680,609 representative rows | 44,325,139 | 116,961,267 |
| Classes with uniquely assigned evidence | 122,483,959 UMI classes | 29,149,820 | 78,953,348 |

Full-input collapsed UMI totals are 28,630,629 for Gene and 77,643,934 for
GeneFull. These include barcodes outside the called-nucleus set. The earlier
34% ratio divided assigned representatives by archive records; it was not a
fraction of distinct records. The corrected record fractions are 23.55% and
62.61%, respectively. Archive records can share a UMI class and are not counts
of independently identified physical molecules.

## Reproduction

Use the three scripts `validate-genefull.py`, `compare-genefull-brain.py` and
`summarize-genefull-effects.py` in `scripts/`. They require new output
directories and preserve prior results. The comparison scripts require NumPy
and SciPy; the synthetic test requires STAR 2.7.11b and the candidate `aie`.
Raw replay reports bind archive/BAM, annotation and barcode inputs to content
identities. The local experiment record retains commands, checksums, raw
outputs, per-nucleus counts, markers and calling lists.

Reference assignment code:
[STAR geneFullAlignOverlap](https://github.com/alexdobin/STAR/blob/b1edc1208d91a53bf40ebae8669f71d50b994851/source/Transcriptome_geneFullAlignOverlap.cpp)
and [exon-derived gene spans](https://github.com/alexdobin/STAR/blob/b1edc1208d91a53bf40ebae8669f71d50b994851/source/Transcriptome.cpp).
