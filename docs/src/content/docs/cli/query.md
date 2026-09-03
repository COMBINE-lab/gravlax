---
title: aie query
description: Indexed region, junction, APA, discovery, and transcript-compatibility queries against an .aie archive.
---

Indexed queries against an `.aie` archive. `region`, `junction`, `junctions`,
and `apa` are annotation-free unless an optional GTF overlay is requested;
`discover` takes an annotation only to define what counts as *claimed*;
`transcript-ecs` uses a content-bound annotation to define its selected
transcript universe. Most indexed commands decode only the chunks they touch.
Transcript equivalence instead discloses a full archive scan so all retained
records of a UMI class are intersected globally. Gene-scale indexed queries
return in roughly 0.1 s; whole-chromosome and full-archive queries are naturally
larger scans.

## Usage

```sh
aie query <ARCHIVE> <COMMAND> [OPTIONS] ...
```

## Shared cell and group scopes

`region`, `junction`, `junctions`, `jset`, `events`, `splice-graph`, `batch`,
and `transcript-ecs` share one strict scope contract:

- `--cells cells.txt` selects a headerless list containing one archive barcode
  per line;
- `--groups groups.tsv` selects a headerless `barcode<TAB>group` mapping;
- `--agg auto|cell|group|bulk` controls the output reduction. `auto` means
  group rows with `--groups` and cell rows otherwise.

The two scope files are mutually exclusive. Empty scopes, malformed or
duplicate rows, unknown archive barcodes, and `--agg group` without
`--groups` fail before evidence is reported. Group order follows first
appearance in the mapping. For older scope-aware commands, omitting all three
options preserves the original unscoped output schemas and bytes; scoped
`junction`, `junctions`, and `batch` results use their v2 schemas.
`transcript-ecs` always uses its own typed envelope and applies scope only after
per-cell compatibility classes have been derived.

## Uniform output for scientific queries

Every query in this page has an opt-in uniform output contract. This includes
`batch`, `region`, `junction`, `junctions`, `jset`, `events`, `splice-graph`,
`apa`, `apa-test`, and `discover`. Omitting `--format` preserves each command's
historical default, `--tsv`, and `--json` output, including its bytes and
command-specific behavior. Select the uniform contract with
`--format text|tsv|json`; do not combine it with a legacy `--tsv` or `--json`
flag. Diagnostics and timing stay on stderr, so stdout contains only the
selected result representation. `transcript-ecs` already uses a typed result
contract and retains its established format and table-selection flags.

```sh
aie query sample.aie junction chr1:155234452-155235327 \
  --groups cell-types.tsv --agg group --format json \
  --output junction-groups.json
```

The JSON form uses `gravlax.result-envelope.v1`. Text and TSV carry the same
typed summary, provenance, warnings, table schemas, and selection metadata.
Multi-table results write their named tables sequentially, so they do not need
to materialize a second copy of all rows merely to change presentation.

| Command | Result schema | Named tables |
|---|---|---|
| `batch` | `gravlax.query.batch.result.v1` | `queries`, `counts` |
| `region` | `gravlax.query.region.result.v1` | `counts` |
| `junction` | `gravlax.query.junction.result.v1` | `counts` |
| `junctions` | `gravlax.query.junctions.result.v1` | `junctions`, optional `counts` |
| `jset` | `gravlax.query.jset.result.v1` | `junctions`, `counts` |
| `events` | `gravlax.query.events.result.v1` | `events`, `components`, `counts` |
| `splice-graph` | `gravlax.query.splice-graph.result.v1` | `nodes`, `edges`, `paths`, optional `group_counts` |
| `apa` | `gravlax.query.apa.result.v1` | `sites`, optional `group_counts` and `group_test` |
| `apa-test` | `gravlax.query.apa-test.result.v1` | `genes` |
| `discover` | `gravlax.query.discover.result.v1` | `candidates` |

Each table declares whether rows are a set, multiset, or sequence, together
with a key and ordering only when those claims are scientifically meaningful.
The normalized `region` and `junction` `counts` tables use schemas
`gravlax.query.region.counts.v1` and `gravlax.query.junction.counts.v1` with
these fields:

| Field | Meaning |
|---|---|
| `aggregation` | `cell`, `group`, or `bulk` |
| `entity` | Visible barcode, group name, or `bulk` |
| `umis` | Scope-filtered UMI-class count |
| `cells` | Nonzero contributing cells for group/bulk rows; null for cell rows |
| `selected_cells` | Cells in the group/bulk scope; null for cell rows |

The table declares set semantics and the key `(aggregation, entity)`. It does
not promise physical row order. For cell output, `--top N` nevertheless needs
a reproducible rule for deciding *which* rows survive truncation: UMI count
descending, then visible barcode ascending. This comparator defines the
selected subset, not a general ordering requirement. Every table reports exact
`available_rows`, `emitted_rows`, and `truncated` values. Under the uniform
contract, `--top 0` consistently means all rows. The legacy unscoped `region`
path retains its historical exception in which `--top 0` emits no cell rows.

Scientific totals live in the typed result summary rather than provenance.
The junction summary names catalogue-wide quantities
`archive_supporting_children` and `archive_posting_chunks`, while `umis` and
`cells` reflect the selected scope. Provenance records the archive access mode,
aggregation, selection policy, and cell scope. A supplied scope file is bound
both by its content digest and by a canonical digest of its resolved archive
mapping. Rooted v2 archives contribute their content identity; legacy v1
archives emit an explicit warning that their path is not a portable content
identity. Uniform batch plans, cohort graph designs, annotations, and APA group
mappings are parsed from the same captured file snapshot whose canonical
`blake3:<hex>` digest is recorded; provenance never reopens their pathname to
describe potentially different bytes.

`--output` (or `-o`) requires `--format`. It stages and installs the complete
result atomically without replacing an existing path, including a dangling
symlink; its current `Flush` durability is not an fsync/crash-durability
promise. The destination parent must already be a directory. Uniform `region`
output rejects `--plot` and `--export-prefix`; uniform `apa` rejects `--plot`;
and uniform `discover` rejects `--emit-gtf`. Those additional files would not
be part of the same transaction, so run a side-artifact command separately.
For uniform APA results, group files use the same strict malformed-row,
duplicate-barcode, and unknown-barcode checks as other scoped queries; legacy
APA parsing remains unchanged for compatibility.

## `batch` — share work across a query panel

Run many anchor-region and exact-junction predicates while opening the archive
once and decoding the union of selected chunks once. The strict input is a
three-column TSV:

```text
id	kind	locus
promoter-a	region	chr1:1000000-1020000
splice-a	junction	chr1:1012345-1016789
```

```sh
aie query sample.aie batch --plan panel.tsv --top 20 > panel.json
```

`kind` is exactly `region` or `junction`; identifiers must be unique and may
not contain whitespace. Coordinates are 0-based and half-open. The plan is
limited to 100,000 predicates. Unknown contigs, malformed rows, duplicate
identifiers, and unsupported kinds fail before any chunk is decoded.

The unscoped JSON result uses schema `gravlax.query.batch.v1`; scoped results
use `gravlax.query.batch.v2`. Both preserve plan order and include the
independent and unique chunk-decode counts. A junction absent
from the archive catalogue is represented with `present: false` and zero
counts rather than aborting the panel. Aggregate counts and returned per-cell
rows have the same class-deduplication semantics as their standalone commands;
`--top 0` returns every cell.

The first release deliberately batches only predicates with mature indexed
standalone semantics. Batched junction enumeration, APA, discovery, and track
rendering remain separate commands.

Use `--format text|tsv|json` for the normalized `queries` and `counts`
tables, and `-o/--output` for atomic no-clobber publication.

## `region` — per-cell evidence in a genomic window

Molecules whose archive anchor lies in `chrom:start-end`; per-cell UMI counts to
stdout. Optional plot/export tracks reconstruct aligned-block overlap and junctions.

```sh
aie query sample.aie region chr6:73489308-73525587 --top 20
```

| Option | Default | Description |
|---|---|---|
| `--top <N>` | `20` | Print the top N cells; with `--format`, `0` means all |
| `--plot <SVG/PNG>` | — | Render strand-split coverage and junction arcs |
| `--export-prefix <PREFIX>` | — | Write plus/minus bedGraph and junction BED12 tracks |
| `--gtf <GTF>` | — | Add a gene underlay to `--plot` |
| `--tsv` / `--json` | off | Emit scoped machine-readable output |
| `--format <FORMAT>` | — | Opt into uniform `text`, `tsv`, or `json` output |
| `-o, --output <PATH>` | stdout | Atomically publish uniform output without replacing a path |

## `junction` — per-cell support for an exact splice junction

Per-cell molecule counts supporting the junction `chrom:donor-acceptor`
(0-based, exact).

```sh
aie query sample.aie junction chr1:155234452-155235327 --top 20
```

| Option | Default | Description |
|---|---|---|
| `--top <N>` | `20` | Print only the top N cells; `0` means all cells |
| `--tsv` | off | Emit a header and machine-readable per-cell rows |
| `--json` | off | Emit one JSON object including all selected cells |
| `--format <FORMAT>` | — | Opt into uniform `text`, `tsv`, or `json` output |
| `-o, --output <PATH>` | stdout | Atomically publish uniform output without replacing a path |

## `jset` — conservative inclusion/exclusion usage

Count UMI classes supporting an inclusion junction set, an exclusion junction
set, or both, with the same cell/group scope used by the point-query paths:

```sh
aie query sample.aie jset \
  --include chr11:34052636-34071725 \
  --exclude chr11:34061757-34071350 \
  --groups cell-types.tsv --json
```

Repeat `--include` and `--exclude` for multi-junction definitions. Each UMI
class is placed in exactly one of `include_only`, `exclude_only`, or `both`,
even when it occurs in several selected chunks or supports several junctions
on one side. The reported usage is
`include_only / (include_only + exclude_only)`; `both` is shown explicitly but
excluded from the denominator, and a zero denominator is JSON `null` / TSV
`NA`. Duplicate loci, loci present on both sides, and unknown chromosomes are
errors. A junction absent from the archive is retained in metadata with
`present: false`. The executor unions posting lists and decodes every selected
chunk once.

| Option | Default | Description |
|---|---|---|
| `--include <LOCUS>` | required | Inclusion junction; repeat to define a set |
| `--exclude <LOCUS>` | required | Exclusion junction; repeat to define a set |
| `--top <N>` | `20` | Cell rows to emit; `0` means all |
| `--tsv` / `--json` | off | Emit machine-readable count rows or a complete object |
| `--format <FORMAT>` | — | Uniform `text`, `tsv`, or `json` tables |
| `-o, --output <PATH>` | stdout | Atomically publish uniform output without replacing a path |

## `events` — discover and reduce splice-event sets

Discover coordinate-defined event candidates from the archive junction
catalogue, then count all candidates with one union decode:

```sh
aie query sample.aie events chr1:45550000-45571000 \
  --event-type cassette --min-support 2 --min-informative 10 \
  --groups cell-types.tsv --json
```

Supported types are `alt-acceptor`, `alt-donor`, and `cassette`; repeat
`--event-type` to select several or omit it for all three. Alternative-site
events pair junctions sharing one endpoint. Cassette events require two
observed flanks and an observed skipping junction. The lower genomic
alternative is deterministically named the inclusion side; this is a stable
coordinate convention, not a transcript-directional claim.

Every component must meet `--min-support`. `--min-informative` is applied
after scope selection and exact reduction. `--max-events` is a hard safety
limit and fails rather than truncating. `--gtf` accepts a GTF or compiled AIC
and adds gene/strand/annotation labels without altering discovery or counts.
The JSON schema is `gravlax.query.events.v1`.

| Option | Default | Description |
|---|---|---|
| `--event-type <TYPE>` | all | Repeatable event-type selector |
| `--min-support <N>` | `2` | Minimum catalogue support for every component |
| `--min-informative <N>` | `1` | Minimum scoped include-only + exclude-only classes |
| `--max-events <N>` | `100000` | Hard candidate limit |
| `--gtf <GTF/AIC>` | — | Count-neutral event labels |
| `--top <N>` | `20` | Per-event cell rows; `0` means all |
| `--tsv` / `--json` | off | Long table or versioned object |
| `--format <FORMAT>` | — | Uniform event, component, and count tables |
| `-o, --output <PATH>` | stdout | Atomically publish uniform output without replacing a path |

## `splice-graph` — exact molecular path fragments

Build a strand-aware junction graph within a locus and retain the exact set of
selected junctions co-supported by each archive UMI class:

```sh
aie query sample.aie splice-graph chr1:45550000-45571000 \
  --groups cell-types.tsv --min-path-umis 2 --json
```

Edges are archived splice junctions. A path is a *molecular path fragment*:
positive evidence that one UMI class supports all listed junctions on one
strand. It is deliberately not called a transcript or isoform, and a
single-edge fragment does not assert that the molecule lacked other exons.
Forward and reverse evidence are separate directed graphs. Multimappers use
the archived anchor placement without attempting to resolve alternatives.

The executor unions catalogue posting lists and decodes each selected archive
chunk once. Repeated representatives and chunks cannot count one UMI class
twice in a path. Edge UMI/cell counts are then derived from the retained path
fragments, so they conserve the path counts containing that edge. With
`--groups`, every edge and path includes exact per-group UMI and cell counts.
The versioned JSON schema is `gravlax.query.splice-graph.v1` and records these
lower-bound semantics explicitly; it emits no population p-value.

| Option | Default | Description |
|---|---|---|
| `--min-support <N>` | `1` | Minimum strand-combined catalogue support for a junction |
| `--min-path-umis <N>` | `1` | Minimum scoped UMI classes for an exact path fragment |
| `--max-paths <N>` | `100000` | Hard candidate-path limit; never truncates |
| `--json` | off | Emit the complete versioned graph object |
| `--format <FORMAT>` | — | Uniform node, edge, path, and group tables |
| `-o, --output <PATH>` | stdout | Atomically publish uniform output without replacing a path |

## `transcript-ecs` — transcript compatibility

Derive deterministic sets of annotated transcripts compatible with archived
UMI classes selected by one gene or one 0-based, half-open window:

```sh
aie query sample.aie transcript-ecs \
  --annotation-file gencode.v49.aic \
  --assembly GRCh38.p14 --annotation-label "GENCODE 49" \
  --feature gene:ENSG00000141510 --groups cell-types.tsv \
  --format json -o tp53-transcript-ecs.json
```

The command unions alternative placements within a retained record and
intersects candidate sets across a class's retained records globally. Results
are exact for this retained archive quotient, but can differ from equivalence
classes built from every original read. It
does not estimate transcript abundance, call an isoform, or phase a complete
transcript. `no_compatible_transcript` and `conflict` are non-exclusive flags
and therefore are not a count partition. A direct `--assembly` value is a caller
assertion; project-plan compatibility checks can additionally verify registered
annotation and coordinate resources.

See the dedicated [transcript-equivalence-class reference](/gravlax/cli/transcript-ecs/)
for selector semantics, typed catalog/count/membership tables, output caps, and
provenance.

## `junctions` — enumerate junctions in a window

Enumerate catalogue junctions whose two endpoints lie in a 0-based,
half-open interval. This does not require knowing an exact junction first.
The index-only path reads catalogue/posting metadata; `--with-cells` decodes
the union of selected posting chunks once and adds exact, class-deduplicated
UMI and cell counts.

```sh
aie query sample.aie junctions chr11:35138870-35232402 \
  --min-support 20 --with-cells --tsv
```

| Option | Default | Description |
|---|---|---|
| `--either` | off | Include a junction when either endpoint is in the window; default requires both |
| `--min-support <N>` | `1` | Minimum index supporting-child count; this is not a read or UMI count |
| `--with-cells` | off | Add exact class-deduplicated UMI, cell, and per-cell counts |
| `--min-cells <N>` | `0` | Minimum exact cell count; nonzero implies `--with-cells` |
| `--gtf <GTF/AIC>` | — | Mark exact, donor, and acceptor annotation membership |
| `--tsv` | off | Emit machine-readable rows |
| `--json` | off | Emit one JSON object with rows and optional cell counts |
| `--format <FORMAT>` | — | Uniform junction and optional count tables |
| `-o, --output <PATH>` | stdout | Atomically publish uniform output without replacing a path |

Format v1 does not store a strand bit in the junction catalogue, so GTF
membership flags are explicitly either-strand. Coordinate and cell-count
answers remain exact.

## `apa` — 3′-end site usage in a window

Clustered molecule 3′-most coordinates (strand-aware) with UMI and cell
counts per site. The count matrix cannot represent this at all.

```sh
aie query sample.aie apa chr1:198692373-198703061 --tsv
```

| Option | Default | Description |
|---|---|---|
| `--site-gap <BP>` | `24` | Site clustering gap in bp |
| `--strand <+/->` | both | Restrict to one strand |
| `--tsv` | off | Emit TSV rows instead of a summary |
| `--groups <TSV>` | — | Two-column TSV (barcode, group): emit per-site per-group UMI counts — differential 3′ usage between cell populations, straight from the archive |
| `--genome <FASTA>` | — | Verify the stamped reference and flag internal-priming sites |
| `--drop-ip` | off | With `--genome`, omit flagged sites instead of only marking them |
| `--permute <N>` | `0` | Cell-label permutations for the site × group G-test |
| `--seed <N>` | `1` | Permutation seed |
| `--plot <SVG/PNG>` | — | Render a 3′-site lollipop plot |
| `--format <FORMAT>` | — | Uniform site, group-count, and test tables |
| `-o, --output <PATH>` | stdout | Atomically publish uniform output without replacing a path |

The `--groups` mode is the differential-APA analysis: label cells by
population (e.g. T cells vs monocytes from your clustering) and compare
per-site usage between groups.

## `discover` — unannotated transcription

Cluster molecules unclaimed by `--gtf` into candidate loci.

```sh
aie query sample.aie discover --gtf gencode.v49.gtf \
  --emit-gtf novel-loci.gtf
```

| Option | Default | Description |
|---|---|---|
| `--gtf <GTF/AIC>` | required | The GTF or compiled `.aic` annotation defining *claimed* molecules |
| `--merge-gap <BP>` | `1000` | Molecule clustering gap |
| `--min-umis <N>` | `10` | Minimum class-deduplicated UMIs per candidate |
| `--claim-mode <MODE>` | `span` | Claiming rule: `span`, `strand-span`, `compatible`, or `residual-sites` |
| `--residual-min-umis <N>` | `10` | Residual-channel support in `residual-sites`; the span channel still uses `--min-umis` |
| `--solo-strand <STRAND>` | `forward` | `forward`, `reverse`, or `unstranded`; ignored by historical `span` mode |
| `--tsv` | off | Emit TSV rows (chrom, start, end, strand, umis, cells) instead of a summary |
| `--emit-gtf <PATH>` | — | Also write candidates as a GTF (single-exon genes) — this file feeds straight back into `replay-rows`, closing the discover→replay loop |
| `--format <FORMAT>` | — | Uniform `text`, `tsv`, or `json` candidate table |
| `-o, --output <PATH>` | stdout | Atomically publish uniform output without replacing a path |

The modes answer different questions:

- `span` is the historical compatibility default. Any transcript-span overlap
  claims a molecule. It is deliberately conservative and independent of
  library-strand configuration.
- `strand-span` claims only same-library-strand overlaps.
- `compatible` claims only evidence exon/junction-concordant with a transcript;
  it is useful as a diagnostic but can join residual evidence into broad
  transitive components.
- `residual-sites` retains every `span` candidate and adds a second channel
  from span-overlapping evidence incompatible with every overlapping
  transcript. That channel clusters transcript-oriented terminal bases in
  non-transitive `--merge-gap` windows; its intervals are bounded site
  extents, never full splice spans.

On evaluated human 10x 3′ datasets, `--residual-min-umis 75` raises
complete-denominator GENCODE v32→v49 recall from 164/367 (44.7%) to 266/367
(72.5%), with 81.5% UMI-weighted recall, 2.16× as many candidates as `span`,
and 99.0% same-strand recurrence in a second PBMC dataset. The default is not
changed: 75 is a measured human 3′ operating point, not a universal
threshold for 5′ or other protocols. Later-annotation overlap measures recall;
it is not the biological precision of all emitted candidates.

```sh
aie query sample.aie discover --gtf gencode.v32.gtf \
  --claim-mode residual-sites --residual-min-umis 75 \
  --emit-gtf novel-loci.gtf
```

The emitted GTF is the input to the **discover → replay** loop: replaying it
quantifies the candidate loci with the full assignment and collapse machinery,
which measures substantially more accurately than the discovery clustering
alone (see [Capabilities](/gravlax/capabilities/#discover-then-re-quantify)).

## Internal-priming filtering and testing (`apa`, `apa-test`)

Supplying `--genome` to `apa` activates the internal-priming filter: sites whose
downstream sequence is A-rich (≥12 A in 20 nt, or an 8-A run within 140 nt, in
transcript orientation) are flagged — these replicate across datasets because they
are sequence-templated, so only the genome exposes them. The FASTA is verified
against the index's stamped signature first. A uniform `apa` or `apa-test`
result with `--genome` therefore requires a stamped archive and fails before
consulting sequence when that identity is absent; stamp the archive with
`aie stamp-genome`. The legacy output paths retain their historical
warning-and-continue behavior for unstamped archives. With `--groups`, a site × group
G-test is reported (`--permute N` adds a label-permutation p-value); `--drop-ip`
excludes flagged sites.

`apa-test` runs the differential analysis genome-wide: an annotation-free
site × group table per gene, a multinomial G-test, and Benjamini–Hochberg FDR
across genes (about 11k genes in ~14 s on a 1k-cell dataset).

```sh
aie query sample.aie apa chr1:198629899-198759346 \
    --groups populations.tsv --genome genome.fa.gz --permute 1000
aie query sample.aie apa-test --gtf gencode.gtf --groups populations.tsv \
    --genome genome.fa.gz > apa-test.tsv
```

For a typed `gravlax.query.apa-test.result.v1` result with a `genes` table,
replace shell redirection with `--format tsv -o apa-test.tsv`; JSON and text
use the same scientific summary and table semantics. The group mapping is
strict under the uniform contract and is bound by both source and resolved
mapping digests.

### `apa-test` options

| Option | Default | Description |
|---|---|---|
| `--gtf <GTF/AIC>` | required | Supplies gene spans and identifiers; it does not supply site positions |
| `--groups <TSV>` | required | Two-column `barcode<TAB>group` mapping |
| `--genome <FASTA>` | — | Verify the stamped reference and filter internal-priming sites |
| `--site-gap <BP>` | `24` | Maximum distance for clustering adjacent terminal coordinates |
| `--min-site-umis <N>` | `5` | Grouped UMIs required for a site to enter a gene's table |
| `--min-gene-umis <N>` | `20` | Grouped UMIs required across sites before testing a gene |
| `--tail-extend <BP>` | `2000` | Include molecules ending this far past the annotated 3′ end |
| `--permute <N>` | `0` | Cell-label permutations per gene; `0` uses only the chi-square approximation |
| `--seed <N>` | `1` | Permutation seed |
| `--format <FORMAT>` | — | Uniform `text`, `tsv`, or `json` result |
| `-o, --output <PATH>` | stdout | Atomically publish uniform output without replacing an existing file; requires `--format` |

## Protocol-aware replicated 3′-end cohorts (`cohort polyasite-mixture`)

In fragmented 3′-tag libraries, an aligned fragment boundary is generally upstream of the
RNA cleavage site and must not itself be called a polyadenylation site. The production cohort
command therefore takes site identities from a same-assembly PolyASite catalogue, rejects
genomic internal-priming candidates, and learns the assay's transcript-upstream fragment-distance
kernel directly from unambiguous molecules. Each donor is deconvolved with a kernel learned from
the other donors, so neither that donor nor its group labels train its observation model.

```sh
aie cohort polyasite-mixture --design donors.tsv --gtf gencode.aic \
    --genome genome.fa.gz --polyasite atlas.bed.gz \
    --group-contrast astro_nsc:mature_neuron --out-dir polyasite-mixture
```

The strict design header is `sample<TAB>condition<TAB>archive<TAB>groups`; each archive is one
biological sample and each groups file is `barcode<TAB>group`. The reducer scans every archive
once, assigns each UMI class at most once per uniquely claimed same-strand terminal region, fits
the candidate mixture by EM, constructs a recurrent donor-level catalogue, and emits complete
sample × group × site expected counts. Inference uses paired donors, an exact sign-flip
calibration, BH correction, effect-size/concordance thresholds, and leave-one-donor-out sign
stability. `--shuffle-seed` supplies a within-donor negative control. On the evaluated eight-donor
human 10x 3′ cohort, 45,705 recurrent sites and 35.7 million assigned expected UMIs are analyzed
in 22.0 s at 3.49 GiB peak RSS—5.18× faster than eight independent legacy `apa-test` scans.

`fragment-kernel.tsv` exposes the learned observation model; `sites.tsv`, `genes.tsv`, and
`summary.json` expose all fitted counts, thresholds, held-out predictive checks, and runtime
semantics. Output overflow and mixed reference identities are hard errors. The model estimates
catalogued terminal-site usage; it does not phase an entire transcript isoform.

The directory remains the primary scientific artifact. Add
`--report-format text|tsv|json` to emit the typed
`gravlax.cohort.polyasite-mixture.result.v1` bundle after the directory is
complete; add `--report-output result.json` to install that report atomically
without replacing an existing file. With no report output, the bundle is the
only stdout content and progress remains on stderr. The bundle normalizes the
wide artifact matrices into `samples`, `sites`, `site_counts`, `genes`,
`gene_usages`, `fragment_kernel`, `heldout_kernel`, and `artifacts` tables. Rows
have explicit set or sequence semantics and keys; physical row order is not
scientifically meaningful except for the increasing fragment-kernel bins. No
table is capped or silently truncated.

Uniform provenance binds the exact design and group-file snapshots, a digest
of the annotation snapshot that was parsed, the normalized PolyASite catalogue,
and every archive's authenticated v2 root (or a full-file digest for a legacy
v1 archive). It also records an `aie-genome-blake3-v1` identity computed from
the same normalized FASTA traversal used for internal-priming checks; the
archives' stamped genome signature is recorded separately after sequence
verification. These identities are inputs; fitted counts and diagnostics stay
in the typed summary and tables.

## Raw endpoint audit (`cohort transcript-ends`)

`aie cohort transcript-ends` constructs one common
cross-sample endpoint-site catalogue and preserves exact sample × group × site
counts before testing. It is useful as an assay-geometry and evidence audit, but broad fragment
endpoint clusters from fragmented 3′ libraries are not cleavage-site calls; use
`polyasite-mixture` for biological polyadenylation-site inference.

```sh
aie cohort transcript-ends --design donors.tsv --gtf gencode.aic \
    --genome genome.fa.gz --polyasite atlas.bed.gz \
    --group-contrast astro_nsc:mature_neuron --out-dir terminal-atlas
```

The reducer scans each archive once, assigns an endpoint only when exactly one
same-strand gene window claims it, and deduplicates each archive UMI class once
per gene. Sites recur by donor, not by cell. Internal-priming sequence, a
strand-correct PAS-motif search, and same-strand PolyASite distance define
separate evidence fields; an aligned endpoint is never silently renamed a
cleavage site. `sites.tsv` contains a zero-explicit packed count matrix and
`genes.tsv` contains donor-paired distal-usage tests, exact sign-flip
calibration, BH q-values, effect-size/concordance guards, and every donor's raw
usage pair. `summary.json` records planning, accuracy, inference, and runtime
semantics. `--shuffle-seed` preserves group sizes while shuffling labels within
each donor for a negative-control run. The output directory must not exist;
site overflow is a hard error rather than truncation.

`--report-format text|tsv|json` provides the corresponding typed
`gravlax.cohort.transcript-ends.result.v1` report, and `--report-output` makes
the report an atomic, no-clobber sidecar. It includes `samples`, `sites`,
`site_counts`, `mixture_sites`, `mixture_site_counts`, `genes`, `gene_usages`,
`fragment_kernel`, and `artifacts`. The site-count tables are lossless normalized
views of the wide TSV matrices: endpoint counts remain unsigned integers and
mixture estimates remain floating point. Gene rows carry an explicit
`endpoint`, `polyasite_only`, or `polyasite_mixture` analysis label. Omitting
the report flags preserves the established directory files and legacy pretty
JSON stdout.

The directory and report have deliberately separate commit boundaries. The
directory is written incrementally and is not transactional; `summary.json` is
its last completion marker. Only after that marker is written does Gravlax
stream the report. A `--report-output` sidecar is installed atomically, but it
is not the directory's completion marker and a late report-publication failure
does not roll back an already complete scientific directory. Destination and
format errors are checked before archive reduction, and a site-cap failure
publishes neither the directory nor the report.

## Plots and IGV export

`region` and `apa` accept `--plot <out.svg|out.png>` (PNG renders in-tool, no
browser needed; a font is bundled so headless hosts work):

- `region --plot` — a sashimi-style portrait: per-strand molecule coverage with
  junction arcs weighted by support, and a gene underlay with `--gtf`.
- `apa --plot` — 3′-site lollipops on a genomic axis (log-scale UMIs); flagged
  internal-priming sites are drawn hollow, and `--groups` adds a per-site
  usage-share band.
- `aie dev em --mask 0.2 --plot` — a reliability diagram of the masked recovery run
  (per-mode empirical accuracy by responsibility decile).

`region --export-prefix <p>` writes IGV-ready tracks instead of (or alongside)
a picture: `<p>.plus.bedgraph` and `<p>.minus.bedgraph` (per-strand molecule
coverage) and `<p>.junctions.bed` (BED12, rendered as a junction track).

```sh
aie query sample.aie region chr1:198690000-198760000 \
    --plot ptprc.png --gtf gencode.gtf --export-prefix ptprc
aie query sample.aie apa chr1:198629899-198759346 \
    --groups populations.tsv --genome genome.fa.gz --plot apa.png
```
