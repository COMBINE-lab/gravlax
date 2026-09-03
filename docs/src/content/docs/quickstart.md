---
title: Quick start
description: Align once, build the index, then replay and query forever.
---

This walkthrough takes a 10x 3′ scRNA-seq sample from raw reads to a queryable
evidence index, then replays a count matrix and runs indexed queries against
it. The workflow has three phases: **align once** (annotation-free), **ingest
once**, and **query forever**.

Before starting, check the local installation:

```sh
aie doctor
```

The commands below work with ordinary paths. For a longer-lived analysis, you
can first create a portable project and run the same operations from checked,
versioned plans:

```sh
aie project init my-analysis --name my-analysis
cd my-analysis
```

See [projects and plans](/gravlax/cli/projects/) for registering inputs and
turning a repeated workflow into a plan. The scientific operations and output
formats are identical in direct and planned runs.

For a no-download exercise first, enter `examples/demo-project` in a repository
checkout. Its tiny compile plan demonstrates checking, dry-run, execution,
identity-verified resume, and read-only Explorer browsing without requiring a
real archive.

After registering an archive and an annotation with explicit assembly/release
identity, run `aie explore`. Its read-only plan builder resolves a gene,
transcript, or exon to exact coordinates and exports the same typed plan as
YAML, JSON, a CLI command, and a Python snippet without executing it.

## 1. Align once, annotation-free

Generate the command for the library's read structure, substituting your real
paths as needed:

```sh
aie ingest recipe --chemistry 10x-3p-v3
```

Align the reads with STAR in two-pass mode against the genome with **no GTF**,
keeping secondary alignments and emitting raw barcodes and UMIs (`CR`/`UR`).
`--soloFeatures SJ` is the one solo feature defined without a gene model:

```sh
STAR \
  --runThreadN 24 \
  --genomeDir star-index-nogtf/ \
  --readFilesIn sample_R2.fastq.gz sample_R1.fastq.gz \
  --readFilesCommand zcat \
  --twopassMode Basic \
  --soloType CB_UMI_Simple \
  --soloFeatures SJ \
  --soloCBwhitelist 3M-february-2018.txt \
  --outSAMattributes NH HI AS nM CR CY UR UY \
  --outSAMtype BAM SortedByCoordinate \
  --outSAMunmapped None \
  --outFilterMultimapNmax 50 \
  --outSAMmultNmax 50 \
  --outFileNamePrefix align/
```

The essential properties of the ingest BAM:

- coordinate-sorted;
- raw barcode and UMI tags present (`CR`, `UR`, and qualities `CY`);
- secondary alignments retained (multimapping structure is evidence);
- **no annotation anywhere in the run**.

This is the only time the reads are touched. Barcode correction happens at
ingest, using an annotation-independent port of STARsolo's
`1MM_multi_Nbase_pseudocounts` corrector — so no gene model is ever needed.

## 2. Build the index once

```sh
aie ingest check align/Aligned.sortedByCoord.out.bam \
  --whitelist 3M-february-2018.txt \
  --chemistry 10x-3p-v3

aie ingest-archive align/Aligned.sortedByCoord.out.bam \
  --whitelist 3M-february-2018.txt \
  --out sample.aie
```

This scans the BAM, extracts molecule-level evidence (junction-chain span
extremes, UMI equivalence classes and 1-mismatch adjacency, interned paralog
patterns), and writes the `.aie` container with genomic chunks plus junction
and cell postings. Region and 3′-site queries derive their answers from those
chunks. Typical cost: 11–18 bits per input read — a 66.6M-read PBMC sample
becomes a 111 MB index.

After this step the BAM and FASTQ are no longer needed for the replay and
query operations below. Keep the original reads whenever future analyses may
need sequence, qualities, a different reference alignment, or different
barcode correction.

## 3. Replay a compatible annotation

Supply a compatible GTF *now* — including releases published after the index was built —
and quantify it against the archived genome alignments:

```sh
aie compile-annotation gencode.v49.annotation.gtf \
  --out gencode.v49.annotation.aic

aie replay-rows sample.aie \
  --gtf gencode.v49.annotation.aic \
  --barcodes barcodes.tsv \
  --out-dir counts-v49/
```

`--barcodes` fixes the output column order (pass a raw `barcodes.tsv` from a
reference run, or any barcode list you want columns for). The output directory
receives `matrix.mtx`, `features.tsv`, and `barcodes.tsv`.

Compilation is optional but useful when an annotation will be reused. The
checksummed `.aic` caches parsing and overlap-index construction; it does not
contain sample data and produces the same results as its source GTF.

Archive-backed Gene replay is bounded by default: decoded genomic chunks are reduced into global
cell shards and released as the scan advances. This is not a per-chunk approximation—UMI classes
that recur in different chunks are merged by the same final global sort and collapse. `--eager`
retains the full-materialization mode for profiling and direct comparison; its
output is byte-identical.

Two comparisons have different meanings. Archive-sourced replay and
`replay-rows --from-bam` are byte-identical because they consume the same
Gravlax row abstraction. A fresh annotation-aware STARsolo run may alter the
alignments themselves; across the four evaluated datasets, that end-to-end
comparison differs by 0.22–0.45% of UMI mass.

Replaying a different annotation is just another invocation — seconds, not a
pipeline:

```sh
aie replay-rows sample.aie --gtf gencode.v32.annotation.gtf \
  --barcodes barcodes.tsv --out-dir counts-v32/
```

Add `--velocity` to emit RNA-velocity (Velocyto) spliced/unspliced/ambiguous
matrices instead of Gene counts.

Resolve a biological name against an explicit reference identity before
constructing a coordinate query:

```sh
aie resolve gencode.v49.annotation.aic TP53 \
  --assembly GRCh38.p14 --annotation "GENCODE 49" --format json
```

Ambiguous symbols or identifiers return an error rather than selecting an
arbitrary first match. Compiled AIC v2 annotations preserve transcript and
source exon IDs as well as genes. Legacy AIC v1 files remain usable for replay,
but transcript/exon resolution asks you to recompile the source GTF.

To measure what changed between two annotation releases while holding the
archive, alignments, and barcode correction fixed:

```sh
aie compare-annotations sample.aie \
  --annotation-a gencode.v44.annotation.gtf \
  --annotation-b gencode.v49.annotation.aic \
  --assembly GRCh38.p14 \
  --annotation-a-label "GENCODE 44" \
  --annotation-b-label "GENCODE 49" \
  --format json --output annotation-change.json
```

This reports exact signed per-cell/gene count deltas within the retained
archive quotient, plus changed-class states, non-exclusive explanatory causes,
and bounded witnesses. It is not a counterfactual fresh realignment. See the
[annotation-comparison reference](/gravlax/cli/compare-annotations/).

## 4. Indexed queries: no annotation required

Per-cell evidence for a genomic region:

```sh
aie query sample.aie region chr6:73489308-73525587 --format json
```

Per-cell molecule counts supporting an exact splice junction:

```sh
aie query sample.aie junction chr1:155234452-155235327 --format json
```

Enumerate junctions in a window without knowing their coordinates first:

```sh
aie query sample.aie junctions chr1:155200000-155300000 \
  --min-support 5 --with-cells --format tsv
```

Compare an inclusion/exclusion junction definition across cell groups. The
same `--cells`, `--groups`, and `--agg` scope is available on region, point
junction, enumeration, and batch queries:

```sh
aie query sample.aie jset \
  --include chr11:34052636-34071725 \
  --exclude chr11:34061757-34071350 \
  --groups cell-types.tsv --format json
```

Discover all coordinate-defined alternative-donor, alternative-acceptor, and
cassette structures in a locus, then reduce them together:

```sh
aie query sample.aie events chr11:34000000-34100000 \
  --min-support 2 --min-informative 10 \
  --groups cell-types.tsv --format json
```

Run a panel of region and exact-junction predicates with shared archive and
chunk work:

```sh
aie query sample.aie batch --plan panel.tsv --top 20 \
  --format json --output panel.json
```

3′-end site usage in a window (the APA view a count matrix cannot represent):

```sh
aie query sample.aie apa chr1:198692373-198703061 --format tsv
```

Each returns in ~0.1–0.5 s; opening the index costs ~10 ms regardless of its
size.

Transcript equivalence classes can be derived for a gene or locus from
retained evidence:

```sh
aie query sample.aie transcript-ecs \
  --annotation-file gencode.v49.annotation.aic \
  --assembly GRCh38.p14 --annotation-label "GENCODE 49" \
  --feature gene:ENSG00000141510 --format json \
  --output tp53-transcript-ecs.json
```

These are compatible-transcript sets and archived UMI-class counts, not
transcript abundance, isoform calls, or full-transcript phasing. They can
differ from classes derived from every original read; see the
[transcript-equivalence reference](/gravlax/cli/transcript-ecs/).

## 5. Discover, then re-quantify

Cluster molecules unclaimed by an annotation into candidate loci, write them as
a GTF, and feed that GTF straight back into replay:

```sh
aie query sample.aie discover --gtf gencode.v49.annotation.gtf \
  --emit-gtf novel-loci.gtf
aie replay-rows sample.aie --gtf novel-loci.gtf \
  --barcodes barcodes.tsv --out-dir counts-novel/
```

The historical default treats any transcript-span overlap as claimed. For the
evaluated human 10x 3′ operating point, add bounded terminal-site candidates
from span-overlapping but transcript-incompatible evidence with:

```sh
aie query sample.aie discover --gtf gencode.v32.annotation.gtf \
  --claim-mode residual-sites --residual-min-umis 75 \
  --emit-gtf novel-loci-sensitive.gtf
```

On the evaluated human 10x 3′ data, this improves complete-denominator
later-annotation recall from 44.7% to 72.5%, emits 2.16× as many candidates as
`span`, and gives 99.0% same-strand recurrence in a second PBMC dataset. The
threshold is not asserted to transfer unchanged to 5′ or other assays.

## 6. Recover multi-gene ambiguity

Evaluate per-cell and cross-cell EM by masking evidence with known truth:

```sh
aie dev em sample.aie --gtf gencode.v49.annotation.gtf
```

Emit the pooled, additive recovered-count layer without modifying the base
matrix:

```sh
aie dev em sample.aie --gtf gencode.v49.annotation.gtf \
  --mask 0 --emit em-layer/ --barcodes barcodes.tsv
```

## 7. Across samples

One junction query federated over many indexes, with per-sample, per-cell
answers:

```sh
aie federate a.aie b.aie c.aie chr1:155234452-155235327
```

For an event table rather than one preselected junction, supply stable sample
IDs and optional per-sample group maps. The union retains only events present
in the requested number of sample catalogues; missing components are explicit
and no cross-sample p-value is invented:

```sh
aie cohort events chr1:0-248956422 \
  --sample PBMC1=pbmc1.aie --sample PBMC2=pbmc2.aie \
  --groups PBMC1=pbmc1-groups.tsv --groups PBMC2=pbmc2-groups.tsv \
  --min-samples 2 --min-informative 10 --min-row-informative 10 --format json
```

The row minimum is pushed into execution: samples with the smallest selected
scope are evaluated first, and an event failing any group or bulk denominator
is not reduced in deeper samples. This changes neither the conservative
`include_only + exclude_only` denominator nor caller sample order.

For an exact common molecular splice graph and a biological-replicate contrast,
provide one unique archive per sample in a strict four-column design:

```text
sample	condition	archive	cells
C1	control	c1.aie	.
C2	control	c2.aie	.
T1	treated	t1.aie	.
T2	treated	t2.aie	.
```

```sh
aie cohort splice-graph chr1:45550000-45571000 \
  --design experiment.tsv --contrast control:treated --format json
```

Use `--counts-only` for a single-condition or descriptive cohort. Low-depth
samples remain visible in the complete matrices but are excluded from tests;
cells and molecules are never treated as replicates.

## Next steps

- [Workflow and interfaces](/gravlax/workflow/) — projects, plans, exact
  resume, Explorer, typed results, and format availability.
- [Capabilities](/gravlax/capabilities/) — what each of these operations
  measures and how well it holds up against fresh processing.
- [CLI reference](/gravlax/cli/) — the full option set of every subcommand.
- [Python and AnnData](/gravlax/python/) — run checked plans and consume typed
  results without duplicating archive logic.
- [Distribution and integrity](/gravlax/distribution/) — release artifacts,
  verification, containers, conda templates, and publication status.
- [The .aie format](/gravlax/format/) — what is stored and the conditions
  under which it supports replay.
