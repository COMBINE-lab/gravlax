# GeneFull EM validation

The recovery evaluator accepts `aie dev em --gene-full`. The original manuscript
experiment below used forward-strand libraries. It builds unique and multi-gene support using the GeneFull replay
index, including ambiguity from overlapping spans at a single placement.
Masking, ten EM iterations, and the emitted pooled layer retain their existing
rules. Subsequent implementation adds `--solo-strand` to recovery and STAR-style
EM, and supports `--gene-full --star`; see `genefull-em-strands-optimization.md`
for the additional correctness checks and optimization measurements.

Tests cover intronic evidence, nested genes, mixed classes split across decoder
batches, masking without unique-evidence leakage, missing-contig alternatives,
conditional-score denominators, CLI model selection, packed/eager unmasked
agreement, model metadata, and incompatible flags. The workspace suite passed
457 tests with three pre-existing ignored tests; the updated CLI fixture also
passed separately. Clippy with warnings denied passed.

## Matched brain-nucleus experiment

The same 6,460 original called nuclei define the scored population. Both models
use 20% deterministic masking, seed 7, blend alpha 20, and ten updates. Pooled
priors use every archive barcode after masking, as in the original sample-wide
experiment. A constant group map restricts scoring without changing the base
uniform/cell/pooled/blend estimators. Every accuracy denominator is evaluable
masked UMI classes before one-mismatch collapse.

| Counting model | Masked classes | Evaluable classes | Uniform top-1 | Cell top-1 | Pooled top-1 | Blend top-1 |
|---|---:|---:|---:|---:|---:|---:|
| Gene | 72,671 | 70,278 | 50.25% | 72.66% | 88.53% | 88.10% |
| GeneFull | 231,351 | 156,884 | 43.16% | 68.93% | 74.90% | 76.13% |

GeneFull loses the unique-evidence label from the remaining candidates for
74,467 masked classes (32.19%). Conditional pooled accuracy is 74.90%; counting
these exclusions as failures gives 50.79% over all masked classes. Seeds 17
and 29 give 74.52% and 74.66% conditional pooled accuracy. The eligible truth
population changes between Gene and GeneFull, so these are not paired scores
on an identical target list or independently established biological labels.

Among GeneFull targets with maximum pooled responsibility at least 0.9,
74,650/80,412 (92.83%) have the correct top-1 label. Those targets account for
51.3% of evaluable classes. Calibration therefore depends on the population;
responsibilities should not be presented as universally calibrated probabilities.

Scoring all archive barcodes instead gives pooled accuracy 90.85% for Gene
(100,186 evaluable classes) and 76.81% for GeneFull (193,141). The original eager
Gene scores used a different masking draw and candidate tie order; its 59.4%
cell / 90.9% pooled result is historical, not the matched control above.

## Annotation and extension controls

Fresh matched STARsolo v32/v49 GeneFull matrices on these nuclei contain
68,671,764 / 70,815,492 collapsed UMIs and differ by 6,570,870 in L1 count mass.
Half-L1 divided by the v49 total is 4.6394%, compared with v49 replay deviation
0.7456%. The fresh v32 Gene result reproduces the historical fixed-nucleus
matrix exactly. GeneFull replay deviation at v32 is 0.7334%.

The saved blanket extension adds 125.86% Gene UMIs but changes GeneFull UMIs by
−0.82% on fixed nuclei. All exon intervals were preserved or expanded; the
count reduction results from assignment and UMI filtering, not span contraction.
The saved evidence-guided extension adds 4,608,988 GeneFull collapsed UMIs
(6.4772%), retaining 87.6% of its absolute Gene gain. Neither annotation was
retuned. GeneFull and Gene references use the same genome, strand, barcode,
and UMI policies within each matched comparison.
