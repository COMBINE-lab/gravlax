---
title: aie extend
description: Evidence-supported per-gene 3′ extension of a GTF — reference optimization as an index query.
---

Evidence-supported per-gene 3′ extension: the reference-optimization workflow of
Pool et al. (Nature Methods, 2023) as an index query. For each gene, molecules
ending downstream of the annotated 3′ end are clustered into candidate cleavage
sites; the gene's end moves to the most downstream site that is

- reachable from the annotated end without a coverage gap larger than `--evidence-gap`,
- supported by at least `--min-umis` molecules in `--min-cells` cells,
- clipped before any gene already occupying the downstream corridor, and
- (with `--genome`) not explainable as oligo(dT) internal priming on a templated A-tract.

The output GTF feeds straight back into [`replay-rows`](/gravlax/cli/replay-rows/):
discover → extend → replay, no realignment. On single-nucleus brain data this
recovers ~20% more counted UMIs in about half a minute.

## Usage

```sh
aie extend <ARCHIVE> --gtf annotation.gtf --out-gtf extended.gtf \
    --report extensions.tsv --genome genome.fa.gz
```

For automation, request the uniform report independently of the two scientific
artifacts:

```sh
aie extend sample.aie --gtf annotation.gtf --out-gtf extended.gtf \
    --report extensions.tsv \
    --report-format json --report-output extend-result.json
```

Omitting `--report-format` preserves the historical summary on standard output
and the byte format of `--report`. With `--report-format text|tsv|json`, standard
output contains only the uniform result (or is empty when `--report-output` is
used); progress and warnings go to standard error. `--report-output` is
installed atomically and never replaces an existing path.

The genome FASTA is verified against the index's stamped signature
(see [`stamp-genome`](/gravlax/cli/stamp-genome/)) before any sequence is consulted.

## Options

| Option | Default | Description |
|---|---|---|
| `--gtf <GTF>` | required | Annotation to extend |
| `--out-gtf <OUT_GTF>` | required | Extended GTF out |
| `--report <REPORT>` | — | Per-gene extension report TSV |
| `--report-format <text\|tsv\|json>` | — | Emit the versioned uniform report; no default, so legacy output remains unchanged |
| `--report-output <PATH>` | — | Atomically publish the uniform report without replacing an existing file; requires `--report-format` |
| `--genome <GENOME>` | — | Reference FASTA: drop internal-priming candidate sites |
| `--max-extend <BP>` | 10000 | Furthest a gene may extend past its annotated 3′ end |
| `--evidence-gap <BP>` | 2000 | Largest tolerated molecule-coverage gap |
| `--min-umis <N>` | 5 | Minimum molecules at the accepting site |
| `--min-cells <N>` | 3 | Minimum distinct cells at the accepting site |
| `--site-gap <BP>` | 24 | Site clustering gap |
| `--min-extension <BP>` | 50 | Extensions shorter than this are skipped |
| `--clip-any-strand` | off | Clip at the next gene on either strand (conservative) |

## Report columns

`gene_id`, `gene_name`, `chrom`, `strand`, `old_end`, `new_end`, `ext_bp`,
`site_umis`, `site_cells`, `n_sites`, `ip_dropped`, `clip` (`neighbor` or `max`).

## Uniform report contract

The outer result schema is `gravlax.extend.result.v1`. Its typed summary records
the annotation gene/transcript counts, number and total length of accepted
extensions, qualifying and internal-priming-site counts for extended genes,
neighbor-clipped extensions, GTF line count, and whether genome filtering was
enabled. Reproducibility provenance includes the archive's root and format
version from the same open archive reader, the exact input-GTF digest, the
normalized genome digest when a FASTA is used, and all extension thresholds.

The bundle contains two named tables:

- `artifacts` uses `gravlax.extend.artifacts.v1`. It is a set keyed by
  `(artifact_kind, path)` and records the extended GTF plus the optional legacy
  per-gene TSV, with byte count and an explicit record unit.
- `extensions` uses `gravlax.extend.genes.v1`. It is a set keyed by
  `gene_index`; physical row order is intentionally unspecified. Each row is a
  typed per-gene extension with identity, locus/strand, old and new 3′ boundary,
  extension length, accepting-site support, qualifying-site count,
  internal-priming drops, and corridor clipping reason.

The GTF and optional legacy TSV remain the primary artifacts. They are written
before the uniform report and are not one cross-file transaction with it. A
late report-publication race can therefore leave complete primary artifacts and
no uniform report. The report itself is staged beside its destination and
installed atomically; producer failure leaves no partial report, and an
existing destination is never overwritten. The uniform report path must differ
from both `--out-gtf` and `--report`.
