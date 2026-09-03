---
title: Capabilities
description: What a molecular evidence index can do — annotation comparison, replay, transcript compatibility, indexed queries, discovery, differential 3′ usage, EM, and federation.
---

Every analysis on this page runs against the same `.aie` files; none requires
returning to reads or rebuilding the archive. Unless another dataset is named,
the measurements below come from four public 10x datasets spanning blood,
tumor, and brain nuclei (66.6M–383.9M reads each), compared with fresh STARsolo
runs using the same documented parameters.

## Annotation replay from fixed alignments

The index is built with no annotation, so *every* annotation is a query-time
input. One index, four annotations, no rebuild: replaying GENCODE v32, v48,
v49, and a withheld-feature panel against the same PBMC index deviates from
fresh STARsolo processing by 0.21–0.25% of UMI mass; across all four datasets
the deviation is 0.22–0.45%. That is 5–11× below the *annotation signal* the
index exists to capture — the v32→v49 release gap alone moves 2.43% of UMI
mass on blood and 4.90% on brain nuclei.

The current evidence reduction is not asserted to be minimal or universally
sufficient. Two design choices account for the observed fidelity:

- **One replay implementation.** Index-sourced and BAM-sourced replay share
  the same row abstraction. `replay-rows --from-bam` reads rows directly from
  a BAM, while the normal path reads equivalent rows from the archive.
- **Two span-extreme representatives.** Each molecule's junction chains store two
  span-extreme representatives plus a read count. One representative costs
  0.62% of UMI mass against fresh STARsolo, two cost 0.25%, and Gravlax replay
  retaining all read placements costs 0.22% — the index adopts the
  two-representative point. The remaining 0.22% includes upstream differences
  between fixed annotation-free alignments and fresh annotation-aware
  processing.

Replay is a bounded genome-ordered scan: decoded chunk batches feed a global
cell-sharded reducer and are released immediately. On two evaluated PBMC datasets this lowers
peak RSS by 3.10–3.31× against same-binary eager replay while adding only
16–24% median wall time; all outputs remain byte-identical. Against
function-matched, counts-only STARsolo, replay is 34–82× faster across the
three evaluated datasets at equal thread budgets.

Annotations can also be compiled once into a deterministic, checksummed `.aic`
artifact and supplied transparently to every existing `--gtf` consumer. The
evaluated GENCODE v49 artifact is 46.9 MB (1.41% of its 3.32 GB source GTF).
On an annotation-marked junction query in an evaluated PBMC archive, it reduces median wall time from
1.19 s to 0.09 s and maximum observed RSS from 376 MiB to 136 MiB, with
byte-identical JSON. Whole Gene replay remains byte-identical and falls from
2.94 s to 1.90 s in the
same single-run check because its annotation setup is amortized.

A second function-matched comparison isolates the storage abstraction from
SAM/BAM tags: a CRAM containing the same post-correction molecule placements
and UMI graph is 2.10–2.55× larger than `.aie`, and streaming its records back
through quantification is 6.64–13.69× slower across the four human/mouse
datasets. This is the fair context for the compact-index claim; the 9–13×
ratio to tag-preserving CRAM compares different retained information.

### RNA-velocity replay

Velocity counting is the hardest replay target: per-UMI transcript-set
intersection across reads, a tolerance-based exon/intron classifier, and UMI
correction applied only within the Gene feature's scope. `replay-rows
--velocity` reproduces the spliced/unspliced/ambiguous matrices to 0.45–6.1%
per-component deviation, with **unspliced — the component velocity analyses
lean on — always below 0.9%**. The residual concentrates in the ambiguous
component of deeply-sequenced molecules whose two span extremes stand in for
three or more reads; it scales with saturation and is purchasable by storing
additional representatives.

### Annotation optimization as a replay

Reference-optimization studies recover missing expression by extending
under-annotated 3′ ends and re-mapping every dataset. With Gravlax that entire
program is a replay: generating a 3′-extended annotation programmatically and
replaying it takes seconds. On brain nuclei — where standard counting uses
only 34% of archived molecules — counted molecules more than double, led by
canonical neuronal genes with long, under-annotated 3′ regions (NRXN1, RBFOX1,
KCNIP4, CADM2).

## Compare two annotations on fixed evidence

`aie compare-annotations` feeds one archive scan into two independent assignment
and UMI-collapse reductions. It reports signed `B - A` gene/count changes, a
complete ledger of changed UMI-class states, non-exclusive contributing causes,
and bounded molecular witnesses. Independent final collapse matters: a changed
class can alter which neighboring UMIs survive, so class transitions are not
additive contributions to the final delta.

The count deltas are exact within the archive's retained-evidence quotient under
the fixed alignment and barcode-correction policy and the requested assignment
settings. They do not claim equivalence to replaying full reads, changing barcode
correction, or performing fresh annotation-aware alignment. Source GTF/AIC bytes
are content-bound in provenance. For a direct path-based command, `--assembly`
is a caller assertion; a project plan can additionally verify assembly
compatibility between registered annotation and coordinate resources.

See [annotation comparison](/gravlax/cli/compare-annotations/) for an example,
the typed result tables, and the explanation semantics.

## Indexed molecular queries

`aie query` answers region, junction, and 3′-site questions at per-cell,
selected-cell, group, or bulk resolution without rebuilding the archive:

- **`region`** — per-cell UMI counts for molecules whose archive anchor lies in
  `chrom:start-end` (optional tracks reconstruct overlapping blocks/junctions);
- **`junction`** — per-cell molecule counts supporting an exact splice
  junction;
- **`jset`** — exact class-level inclusion-only, exclusion-only, and shared
  support for repeatable junction sets; conservative usage excludes the
  explicitly reported shared category;
- **`events`** — discover alternative-acceptor, alternative-donor, and cassette
  structures from junction coordinates and support alone, then reduce every
  event in one union decode at cell, group, or bulk resolution;
- **`splice-graph`** — construct separate strand-aware junction graphs and
  count exact UMI-class path fragments, preserving multi-junction co-support
  without claiming that a fragment is a complete transcript;
- **`junctions`** — enumerate junctions in a 0-based window from the catalogue,
  filter by index support, and optionally recover exact UMI/cell counts and
  annotation endpoint membership;
- **`apa`** — clustered 3′-end site usage in a window, strand-aware, with UMI
  and cell counts per site;
- **`apa-test`** — genome-wide site × group tests with optional
  reference-verified internal-priming filtering and BH-FDR;
- **`cohort transcript-ends`** — a recurrent cross-donor endpoint catalogue,
  exact sample/group/site matrices, reference-supported confidence tiers, and
  paired biological-sample distal-usage inference for assay-geometry audits;
- **`cohort polyasite-mixture`** — protocol-aware 3′-tag analysis using orthogonal
  PolyASite identities, reference-verified internal-priming rejection, a label-blind
  empirical fragment kernel, leave-one-donor-out deconvolution, and paired donor inference.

For fragmented 3′ protocols, `polyasite-mixture` is the preferred biological analysis:
`transcript-ends` deliberately reports aligned fragment boundaries and does not interpret them as
cleavage coordinates. On an eight-donor adult human SEZ cohort, one `polyasite-mixture` invocation
analyzes 45,705 recurrent candidates and 35.7 million expected UMIs in 22.0 s at 3.49 GiB peak
RSS, 5.18× faster than eight independent legacy `apa-test` scans. Its numerical tables are
byte-identical before and after sparse-likelihood hot-path optimization.

`query batch` accepts a strict panel of anchor-region and exact-junction
predicates, unions their chunk selections, and decodes each selected chunk
once. On a 32-query high-locality PBMC panel it reduces 48 independent
chunk decodes to 2 (95.8%) and median wall time from 0.99 s for 32 standalone
processes to 0.06 s (16.5×), with exact aggregate and top-20 per-cell results.
Peak batch RSS is 142,752 KiB (139.4 MiB).

Gene-scale index queries are typically around 0.1 s. A whole-chromosome
`junctions --with-cells` scan returning 5,089 rows takes 0.44 s median and
about 1.05 GB peak RSS; index-only locus enumeration is much smaller.

The junction catalogue, cell postings, and 3′-site table are part of the
index itself — the compression dictionaries double as query indexes.

The same strict `--cells`, `--groups`, and `--agg auto|cell|group|bulk`
interface applies to region, point-junction, junction enumeration, junction
set, event, splice-graph, batch, and transcript-equivalence-class queries. It
filters classes before aggregation rather than post-filtering displayed rows.
On an evaluated PBMC junction set, the union executor takes 0.05 s median
versus 0.09 s for two independent point-query processes, while grouped output
also takes 0.05 s and peaks at 76,576 KiB.

Across archives, `cohort splice-graph` first defines coordinate edges recurring
in a requested number of local catalogues, then performs the same exact
class-level reduction separately in every biological sample. It emits complete
sample × edge and sample × path matrices. With a strict two-condition design it
fits path-versus-rest beta-binomial contrasts at the sample level, excluding
rather than zero-filling samples below `--min-sample-umis`. A
counts-only mode supports exploratory panels without manufacturing inferential
replication.

`query events` is an event engine rather than a loop over `jset`: each junction
maps to all event sides it serves, every selected chunk is decoded once, and
each `(event, molecule class)` is globally reduced once. Candidate discovery
does not consult a GTF. An optional GTF or compiled AIC labels genes, strand,
and component annotation status without changing candidate membership or
counts. Alternative-site side names follow increasing genomic coordinate and
therefore do not imply transcript-directional 3′/5′ terminology.

## Experimental transcript compatibility

`aie query ARCHIVE transcript-ecs` derives annotation-conditional compatible
transcript sets for archived UMI classes selected by a gene or genomic window.
Alternatives for one retained record are unioned, while candidate sets across
the class's retained records and representatives are intersected globally,
including across archive chunks. Stable content-addressed EC identifiers and
cell/group/bulk counts make the result directly queryable without selecting one
transcript.

This is an exact deterministic relation over retained archive geometry, not a
claim of equivalence classes derived from every original read. In particular,
the result may differ for singleton and ambiguous classes in fragmented 3′
protocols. Counts are archived UMI-class counts—not transcript abundance,
post-replay gene-collapse counts, isoform calls, or full-transcript phasing.
`no_compatible_transcript` and `conflict` are non-exclusive flags: one class can
contain an unmatched record and also have disjoint nonempty candidate sets, so
their totals must not be summed as a partition.

See [transcript equivalence classes](/gravlax/cli/transcript-ecs/) for selectors,
scope, typed outputs, provenance, and output limits.

## Discover, then re-quantify

`aie query discover` clusters molecules unclaimed by a supplied annotation
into candidate loci, and `--emit-gtf` writes them as an annotation that feeds
straight back into `replay-rows` — the index quantifies transcription no
annotation has named yet.

Annotation history supplies an external test. Given only GENCODE v32 (the
release in production use in 2020), the conservative default recovers
**51/52 (98.1%)** sufficiently expressed v49 additions at loci free of v32
features. A complete denominator including old-annotation overlaps exposes
its intended blind spot: 164/367 (44.7%) overall.

The optional `residual-sites` mode keeps every conservative call, then adds
bounded transcript-oriented terminal-site clusters from molecules that
overlap a v32 span but are incompatible with every overlapping transcript.
Using `--residual-min-umis 75` recovers **266/367 (72.5%)**, or 81.5%
of truth UMI mass, with 26,787 candidates (2.16× the `span` result) and 99.0%
same-strand recurrence in a second PBMC dataset. Recall is 51/52 in clean loci, 52/60 for
opposite-strand overlaps, and 163/255 for same-strand overlaps. Every added
site interval is at most 1,001 bp, so broad genomic spans cannot manufacture
overlap recall. More permissive 25/50-UMI settings recover more truth and emit
more candidates.

These are later-annotation **recall** measurements, not estimates of the
biological precision of every emitted candidate. Candidates without v49
overlap remain unresolved. On a controlled withheld-gene panel the
discover→replay loop reaches 7.6% median count error — versus 26% for the
discovery clustering alone — and 79% of one PBMC sample's intergenic loci
(90% of their mass) recur independently in a second sample.

## Differential 3′-end usage between cell populations

`aie query apa --groups` emits per-site, per-population 3′ counts straight
from the index. A T-cell versus monocyte contrast on PBMC yields 15% of
testable genes shifted beyond the within-population null's 95th percentile
(5% expected), led by PTPRC (CD45), CD44, and LCP1 — canonical cell-type 3′
isoform switchers. The signal replicates across independent datasets (26 of
27 co-testable significant genes recur), and the implicated sites fall within
200 bp of annotated transcript 3′ ends. The count matrix cannot represent
this analysis at all.

## Cross-cell EM multimapper recovery

Per-cell quantifiers discard gene-informative evidence: 5–6% of countable
molecules carry only multi-gene evidence and contribute nothing to the
matrix. Because the index holds every cell's equivalence-class evidence
simultaneously — the interned paralog pattern *is* the equivalence class —
sample-wide EM needs no second pass over reads.

`aie dev em` can score recovery with a masked-evidence analysis: classes whose
unique reads identify a gene are reduced to their multimapper evidence and the
models attempt to recover that gene. In evaluated datasets, pooled cross-cell EM attains
95.7–98.3% top-1 accuracy where per-cell EM reaches 90.0–94.4% — and on
shallow brain nuclei, where per-cell priors collapse to 59.4%, pooling holds
90.9%: **sharing matters most exactly where cells are sparsest.** The
recovered layer:

- adds +6.3% of countable molecules (85% at responsibility >0.8),
  concentrated in ribosomal, mitochondrial, and histone paralog families;
- is externally corroborated (it collapses the discrepancy against a
  sequence-aware quantifier from a median 2.5× to 11%);
- is calibrated (the dominant r > 0.9 responsibility bucket is 99%
  empirically correct);
- is reproducible across datasets (per-gene recovered-mass fractions agree at
  Pearson r = 0.992).

Recovered counts are an additive, opt-in, real-valued layer (`--emit` writes
`em.mtx`); the base replay matrices are never modified. For compatibility,
`--star` implements the STARsolo `--soloMultiMappers EM` per-cell algorithm.

The default implementation streams archive batches into packed 8-byte support
records, reduces 64 deterministic cell shards, and iterates flat CSR labels.
Large inputs use temporary exact support shards, finalized one at a time. On
an evaluated 10,000-cell dataset, the hybrid evaluator took 35.86 s and
3,846,008 KiB peak RSS. `--eager` runs the older all-in-memory implementation
for direct comparison.

An optional barcode/group map adds group-only and cell/group/sample
hierarchical evaluators with negative-log-loss and multiclass-Brier output.
The `convex` evaluator normalizes cell, leave-one-cell-out group, and sample
distributions within each candidate set and mixes them with
simplex-constrained weights. The `dirichlet-proxy` evaluator uses
candidate-level posterior means; it is a screening approximation, not full
variational inference. The `depth-hybrid` evaluator transitions monotonically
from the convex to the proxy estimate as fitted unique-candidate evidence
increases. These additional modes are analysis tools only: `--emit` writes
pooled recovered counts, while `--convex-only`, `--dirichlet-only`, and
`--hybrid-only` write scores but cannot emit a recovered-count matrix. See
[`aie dev em`](/gravlax/cli/em/) for every model and option.

## Federation across indexes

`aie federate` runs one junction query across N archives and returns
per-sample, per-cell counts byte-equal to the single-index answers. Six
indexes spanning three tissues, chemistry generations from v3 to LT Chromium
X, and a 20× cell-count range — roughly 120 GB of raw reads living as 2.8 GB
of queryable indexes — answer a federated junction query in about 0.5 s, with
~94,000 cells reporting. Across PBMC samples, 89–97% of each sample's junction
support mass lies on junctions replicated in the others, so cross-sample
recurrence doubles as a confidence signal for discovery. This is the atlas
access pattern — per-sample molecule resolution under one query surface —
that coverage-domain archives cannot offer.

For catalogue-wide `cohort events`, `--min-row-informative` pushes an exact
minimum conservative denominator through the sample plan. The shallowest
sample/group scope is evaluated first, so events that cannot appear in the
requested table are discarded before deeper archives are decoded. Packed
event/class hits and a direct group reducer avoid retaining cell vectors when
only group or bulk rows are requested.

## Limitations

- **Barcode correction is fixed at ingest.** Re-correction under a different
  whitelist is the one replay the index cannot perform.
- **No read sequence or qualities.** Variant-aware and allele-specific
  analyses are out of scope, and sequence-consuming quantifiers can be
  corroborated but not replayed exactly. The format's optional edit and
  residual-sequence sections mark the extension path.
- **Chemistry coverage is still bounded.** On one public mouse 5′ v2 dataset,
  Gene replay moved 0.75% of UMI mass relative to fresh STARsolo, while the
  velocity-component differences were larger. Measurements are not yet
  available for 5′ v3 or full-length protocols.
