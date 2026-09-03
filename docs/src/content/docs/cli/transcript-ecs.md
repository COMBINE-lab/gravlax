---
title: Transcript equivalence-class queries
description: Derive annotation-conditional compatible-transcript sets from retained archive evidence.
---

`aie query ARCHIVE transcript-ecs` derives the set of annotated transcripts
compatible with each archived UMI class in a selected gene or genomic window.
It uses the archive's retained alignment blocks, junctions, strand, and
alternative placements rather than reconstructing reads.

This is an exact query over the evidence retained in the archive, not a claim
that the result equals equivalence classes derived from every original read.
The distinction is especially important for singleton and ambiguous classes in
fragmented 3′ protocols.

```sh
aie query sample.aie transcript-ecs \
  --annotation-file gencode.v49.annotation.aic \
  --assembly GRCh38.p14 \
  --annotation-label GENCODE-v49 \
  --feature gene:ENSG00000141510 \
  --format json -o tp53-transcript-ecs.json
```

Exactly one selector is required. `--feature` must resolve unambiguously to a
gene through the bound annotation identifier resolver. `--locus` accepts a
0-based, half-open `contig:start-end` window and selects every transcript whose
span overlaps that window. Missing archive or annotation contigs fail instead
of being reported as a biological zero.

For a direct path-based invocation, `--assembly` is a caller assertion retained
in provenance. The annotation digest binds the parsed bytes but does not prove
their assembly. A [project plan](/gravlax/cli/projects/) can additionally verify
declared compatibility between registered annotation and archive resources.

## Compatibility semantics

Alternative placements for one retained record are unioned: a transcript is a
candidate when any retained alternative is concordant with it. Candidate sets
from separate retained records and representatives of the same UMI class are
then intersected globally, including when those records occur in different
archive chunks. Strand compatibility follows `--solo-strand`.

Selection and compatibility are deliberately separate. A UMI class belongs to
the query when at least one retained alignment block overlaps the selected
gene/window. If such a class has no compatible selected transcript, it remains
in the result as `no_compatible_transcript` with a null `ec_id`; it is not
silently discarded or interned as an ordinary empty equivalence class. A class
whose nonempty record-level candidate sets have an empty intersection is also
flagged as a conflict.

The diagnostic flags are not one mutually exclusive outcome partition. In
particular, `no_compatible_transcript` and `conflict` can both be true when one
retained record has no compatible transcript and other nonempty record sets are
disjoint. Ambiguity can overlap incompleteness as well. Do not sum these flag
counts. Only assigned/unassigned and complete/incomplete are documented
complementary pairs over the scoped UMI classes.

The catalog's `ec_id` is content-addressed from the exact annotation digest and
the sorted transcript IDs. It is stable across output ordering and cell
filtering. GTF and AIC v2 annotations are supported. AIC v1 files lack stable
transcript identifiers and fail with a request to recompile the source
annotation.

These results describe compatibility within the archive's retained-evidence
quotient. Counts are counts of archived UMI classes. They are not post-replay
gene-collapse counts, transcript-abundance estimates, isoform calls, or
evidence of full-transcript phasing, and the sets are not asserted to equal
full-read-derived equivalence classes. Completeness and retained-representative
flags make these limits explicit per class.

## Cell scopes and aggregation

The command uses the shared query scope:

- `--cells cells.txt` selects a headerless list of archive barcodes;
- `--groups groups.tsv` selects a headerless `barcode<TAB>group` mapping;
- `--agg auto|cell|group|bulk` controls count-row aggregation.

`auto` emits group rows with `--groups` and cell rows otherwise. Scope and
aggregation are applied only after per-cell equivalence classes have been
derived. UMI classes are never merged or deduplicated across cells.

When a scope file is used, result scope and provenance retain its path, a
BLAKE3 digest of the exact bytes parsed, and a second canonical digest of the
resolved sorted barcodes and barcode-to-group assignments. The all-cells scope
is likewise bound to its resolved archive-barcode population. Thus two
same-sized but different cell lists or group mappings do not share query
provenance.

The current implementation scans the full archive to preserve UMI-class
intersection across chunks. The result records `archive_access:
full_archive_scan` and `chunk_pruning_applied: false`; the feature/locus and
cell scope filter the derived classes, not the archive read itself.

## Output

Text is the default and provides a concise summary plus catalog and count
tables. JSON is one `gravlax.result-envelope.v1` result with result schema
`gravlax.query.transcript-ecs.v1`. Its `data` contains the exact selection,
selected transcript IDs, scope, scientific semantics, summary, and versioned
typed tables:

| JSON member / TSV selector | Table schema | Contents |
|---|---|---|
| `catalog` / `catalog` | `gravlax.query.transcript-ecs.catalog.v1` | Stable EC IDs, transcript/gene sets, ambiguity, class and cell counts |
| `counts` / `counts` | `gravlax.query.transcript-ecs.counts.v1` | Cell/group/bulk archived UMI-class counts and outcome/completeness totals |
| `membership` / `membership` | `gravlax.query.transcript-ecs.membership.v1` | One diagnostic row per scoped UMI class |

Membership is omitted unless `--emit-membership` is supplied. TSV emits one
typed table and therefore requires `--table`; `--table` is rejected for text
and JSON, and the membership table additionally requires
`--emit-membership`.

```sh
# A locus-selected catalog for a pipeline
aie query sample.aie transcript-ecs \
  --annotation-file gencode.v49.annotation.gtf \
  --assembly GRCh38.p14 --annotation-label GENCODE-v49 \
  --locus chr17:7668402-7687550 \
  --format tsv --table catalog > catalog.tsv

# Per-class diagnostic membership for selected cells
aie query sample.aie transcript-ecs \
  --annotation-file gencode.v49.annotation.aic \
  --assembly GRCh38.p14 --annotation-label GENCODE-v49 \
  --feature TP53 --cells cells.txt --emit-membership \
  --format tsv --table membership > membership.tsv
```

JSON provenance includes the assembly and annotation labels, observed
annotation BLAKE3 digest, exact query parameters, and the authenticated archive
root read by the same archive reader used for derivation. Rooted archives use
the canonical `aie-directory-root-v2:<root>` identity shared with collections,
projects, and federation. An optional
`--annotation-digest blake3:<64 lowercase hex>` verifies the expected
annotation before analysis. Legacy AIE v1 input has no rooted commitment and
is reported with an explicit nonportable path locator and provenance warning.

`--max-ecs` (default 100,000) and `--max-memberships` (default 1,000,000) are
hard limits: excess output fails without truncation. The membership limit is
applied only when membership was requested. The counts table also has a fixed
1,000,000-row safety limit, recorded as `max_count_rows` in provenance.
`-o/--output` fully renders and atomically installs a new file; existing paths
are never replaced, and a validation or cap failure leaves no partial result.

## Options

| Option | Default | Description |
|---|---|---|
| `--annotation-file <GTF/AIC>` | required | Exact annotation resource used for transcript compatibility |
| `--assembly <ASSEMBLY>` | required | Caller-asserted reference assembly identity; project plans can additionally verify compatibility |
| `--annotation-label <LABEL>` | required | Annotation release or immutable label |
| `--annotation-digest <DIGEST>` | — | Expected `blake3:<64 lowercase hex>` content digest |
| `--feature <GENE>` | one selector required | Gene stable ID or symbol; explicit `gene:` qualification is accepted |
| `--locus <LOCUS>` | one selector required | 0-based, half-open transcript-selection and evidence-overlap window |
| `--solo-strand <STRAND>` | `forward` | `forward`, `reverse`, or `unstranded` |
| `--cells <FILE>` | — | Headerless archive-barcode scope |
| `--groups <FILE>` | — | Headerless `barcode<TAB>group` scope |
| `--agg <LEVEL>` | `auto` | `auto`, `cell`, `group`, or `bulk` |
| `--emit-membership` | off | Include the optional per-class membership table |
| `--max-ecs <N>` | `100000` | Hard catalog-row limit; zero is rejected |
| `--max-memberships <N>` | `1000000` | Hard requested-membership-row limit; zero is rejected |
| `--format <FORMAT>` | `text` | `text`, `json`, or `tsv` |
| `--table <TABLE>` | — | Required for TSV: `catalog`, `counts`, or `membership` |
| `-o, --output <FILE>` | stdout | Atomically create a new result file without clobbering |
