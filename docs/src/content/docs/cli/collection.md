---
title: aie collection
description: Build, inspect, and query a content-addressed routing index over independent .aie archives.
---

An `.aicollection` is an authenticated routing sidecar over independent `.aie`
archives. It stores source identities and interval, junction, and optional
local-shape routes; molecule evidence remains in the source archives and exact
counts are recomputed there.

## Build and inspect

```sh
aie collection build \
  --sample donor-a=/data/donor-a.aie \
  --sample donor-b=/data/donor-b.aie \
  --shape-routes \
  --out atlas.aicollection \
  --format json --output atlas.build.json

aie collection inspect atlas.aicollection \
  --verify-routes --format json --output atlas.inspect.json
```

Sample IDs are sorted before encoding. Sources must have identical chromosome
dictionaries and, unless `--allow-unstamped` is used, the same stamped genome
identity. Duplicate IDs, resolved paths, inodes, and encoded archive content
are rejected. `--shape-routes` derives source-root-bound exact intron-span
routes without reading molecule chunks.

An incremental build adds an immutable layer:

```sh
aie collection build --base atlas.aicollection \
  --sample donor-c=/data/donor-c.aie \
  --shape-routes --out atlas-plus-c.aicollection
```

The child records its parent's canonical path and authenticated root. Queries
verify every layer, reject changed parents and cycles, and remap layer-local
archive ordinals into one sample order.

`inspect` authenticates all collection payloads and checks every source's
recorded filesystem and content identity. `--verify-routes` reconstructs stored local-shape
routes from their root-bound source dictionaries. `--verify-content` also
verifies complete source content.

### Build and inspection options

| Command | Argument or option | Default | Description |
|---|---|---|---|
| `build` | `--sample <ID=ARCHIVE>` | repeatable | Add a named source archive; a build needs a sample or `--base` |
| `build` | `--source-digest <ID=BLAKE3>` | — | Require the named source's native v1 file digest or v2 directory root to match |
| `build` | `--base <COLLECTION>` | — | Extend an existing sidecar without rescanning its source indexes |
| `build` | `--out <PATH>` | required | New `.aicollection` destination |
| `build` | `--allow-unstamped` | off | Allow sources without a genome digest; chromosome dictionaries must still match |
| `build` | `--shape-routes` | off | Add exact source-bound intron-span routes from shape dictionaries |
| `build` | `--json` | off | Emit the command-specific JSON build summary |
| `inspect` | `<COLLECTION>` | required | Collection to inspect |
| `inspect` | `--verify-routes` | off | Decode and reconstruct every stored shape-route block |
| `inspect` | `--verify-content` | off | Re-hash every source archive and reconstruct routed shapes |

## Exact queries

```sh
aie collection junction atlas.aicollection chr2:1234567-1250000 \
  --top 0 --format json

aie collection region atlas.aicollection chr2:1200000-1300000 \
  --format tsv

aie collection jset atlas.aicollection \
  --include chr2:1234567-1250000 \
  --exclude chr2:1234567-1260000 \
  --format text
```

- `junction` returns exact per-sample point-junction UMI and cell counts.
- `region` returns exact per-sample anchor-window molecule, UMI, and cell
  counts.
- `jset` classifies UMI classes as `include_only`, `exclude_only`, or `both`;
  usage is `include_only / (include_only + exclude_only)`.

## Atlas-wide reverse search

`find-events` searches the collection without requiring a locus. For splice
entities it scans the authenticated collection catalogue, applies
catalogue-safe recurrence bounds, and then opens only the routed source chunks
needed to recompute exact, class-deduplicated support. The reported unit is an
archive UMI class: one corrected cell barcode paired with one retained raw UMI
value. This search does not apply the archive's annotation-dependent 1-mismatch
edges, so these values are not final gene-level corrected-UMI counts. For
terminal tails it reads each archive's sparse `index.tail` first and decodes only the
tail-bearing chunks; it does not scan ordinary molecule chunks that cannot
contain a retained terminal event.

```sh
aie collection find-events atlas.aicollection \
  --kind junction --kind cassette --kind terminal-tail \
  --design donors.tsv \
  --groups cell-groups.tsv \
  --require-group neuron --require-group astrocyte \
  --min-support 2 --min-samples 4 --min-donors 4 \
  --min-umi-classes 20 --min-side-umi-classes 2 \
  --min-group-umi-classes 5 \
  --terminal-cluster-bp 25 \
  --annotation gencode.v49.aic \
  --assembly GRCh38.p14 --annotation-label GENCODE-v49 \
  --novel-only \
  --format json --output recurrent-events.json
```

`donors.tsv` is strict and has this header:

```text
sample	donor
sample-a	donor-1
sample-b	donor-2
```

Without `--design`, each collection sample is treated as a distinct donor and
that choice is recorded in provenance. `cell-groups.tsv` is also strict:

```text
sample	barcode	group
sample-a	AAACCCAAGAAACACT	neuron
sample-a	AAACCCAAGAAACCAT	astrocyte
```

When `--groups` is present, listed cells define the scope; unlisted cells do
not contribute. `--require-group NAME` requires the pooled exact count for
that group to reach `--min-group-umi-classes`. `--min-samples`, `--min-donors`,
and `--min-umi-classes` are checked again against exact source-derived counts after
routing. Alternative-splicing entities also require both alternatives to
reach `--min-side-umi-classes`, so a one-sided catalogue pattern is not
reported as an observed event. Every component junction must independently
reach that threshold as well; in particular, a cassette event cannot survive
when only one inclusion flank has exact support. Cassette `include_only` is a
population union of classes witnessing either inclusion flank without the
skip junction. The marginal component checks may be satisfied by different
classes, so this command does not claim that one class spans the full inclusion
path; use `query cooccur` for that molecule-level question. These are
recurrence and group-presence predicates, not a
differential or paired-donor significance test. The sparse `counts` table
retains sample, donor, and group labels for a downstream contrast.

Splice-event support is deliberately restricted to unique-read chain
representatives. The collection junction catalogue is a conservative routing
superset that also contains BAM-designated primary multimapper anchors, but a
candidate with no unique-chain support is removed by exact filtering.
Multimapper primary placements and their pattern-dictionary alternatives are
not counted or searched. This avoids making a result depend on an arbitrary
primary alignment and does not claim complete all-placement event discovery.
The summary and provenance record
`evidence_placement_policy=unique_chain_representatives_only`,
`multimapper_placements_included=false`, and
`multimapper_alternatives_available_to_search=false`.

The supported entities are `junction`, `alt-acceptor`, `alt-donor`, `cassette`,
and `terminal-tail`; omitting `--kind` selects all five. In `forward` and
`reverse` modes, splice entities are split by exact evidence strand before
filtering. Consequently, donor/acceptor names follow transcript orientation on
both strands rather than genomic left/right names. In `unstranded` mode,
`junction` and `cassette` evidence from both alignment orientations is
coalesced before any UMI, sample, donor, or group threshold is applied. The
single entity reports strand `.`; its forward/reverse alignment-UMI columns
retain the informative-class orientation totals for auditing. As with
`exact_umi_classes`, both-side-only classes are excluded from those totals.

Alternative donor/acceptor events use neutral coordinate sides rather than a
potentially misleading exon-inclusion label. `side_a` is the lexicographically
lower `(donor, acceptor)` junction and `side_b` is the higher junction,
independent of transcript strand. Their counts are
`side_a_only_umi_classes` and `side_b_only_umi_classes`; they are not PSI or a
biological inclusion/exclusion assignment. Cassette events retain biological
`include` flank and `exclude` skip roles.

Terminal-tail entities are deterministic, strand-aware single-linkage
clusters: consecutive exact cleavage anchors are joined when their gap is no
larger than `--terminal-cluster-bp` (25 by default). `0` reports one entity per
exact anchor. The entity row gives the half-open cluster bounds and its summit,
chosen by UMI support with the lower anchor breaking ties. Clustering never
hides its exact coordinate components: `components` and `terminal_anchors`
retain every exact anchor, while `terminal_counts` retains its nonzero
sample/group UMI-class counts and `terminal_anchors` reports signal maxima and
per-anchor annotation compatibility. With `--novel-only`, compatible exact
anchors are discarded before clustering; clusters and every recurrence/support
threshold are then recomputed from incompatible anchors only. Thus nearby
annotated support cannot cause a novel anchor to disappear or prop up its
counts. The archive
remains authoritative for every retained per-molecule signal. `--min-support`
is a splice-catalogue pruning bound and does not filter terminal tails.

Terminal-tail evidence is optional per archive. `capabilities` reports both an
aggregate and one explicit availability/denominator row per source archive.
Recurrence and support for a terminal entity use only archives and donors that
declare the typed terminal-tail capability; an unsupported archive is
unavailable evidence, not a measured zero. A capable archive with no event is
a measured zero. Tail
extraction has the archive's fixed `forward-cdna-terminal-softclip-v1`
strand/edge semantics, so requesting `terminal-tail` with any
`--solo-strand` other than `forward` is rejected. Before loading `index.tail`
or any molecule chunk, `--max-terminal-events` checks the sum of the rooted
capability declarations and fails rather than truncating an oversized search.
The combined splice-candidate and terminal-cluster budget is consumed before
each terminal entity is built, so separated tail anchors cannot allocate past
`--max-candidates`.

`--annotation` accepts an uncompressed GTF or compiled AIC and requires a
caller-declared `--assembly` and immutable or
descriptive `--annotation-label`; `--annotation-digest` can additionally bind
the exact file bytes. The result records the observed annotation identity, the
collection's stamped genome algorithm/digest (when present), and compatibility
status `caller_declared_unverified`. Exact contig-name overlap is required but
does not prove reference-sequence identity. A zero-overlap annotation fails,
and candidates and tail routes on a collection contig absent from the
annotation are omitted rather than confidently called novel or treated as
zero. The summary and provenance report exact omitted splice-candidate,
terminal-route, and declared terminal-event counts plus the sorted contig
names. Names such as `1` and `chr1` are never silently equated.

Each classified entity reports the number of compatible transcripts and
overlapping gene IDs. Annotation-gap flags are deliberately nonexclusive:

- `missing_junction`: both splice endpoints are known on the requested strand,
  but their exact junction is absent from every eligible transcript;
- `boundary`: one or both splice endpoints are absent on the requested strand,
  or a terminal cluster has no exact same-strand annotated 3-prime boundary;
- `strand`: an exact opposite-strand junction exists, or only opposite-strand
  annotation overlaps the evidence;
- `overlap`: the exact component junctions of a multi-junction side are known
  in a same-strand locus, but no single transcript contains the required path,
  or alternative sides are annotated only in disjoint genes.

`gap_primary_class` uses the deterministic precedence `strand`, `boundary`,
`overlap`, then `missing_junction`. This makes the four primary classes
structurally distinct while retaining nonexclusive flags where, for example,
opposite-strand context and a boundary discrepancy co-occur. For a terminal
cluster, compatibility means that at least one retained exact anchor equals an annotated 3-prime transcript
boundary on the same strand. `--novel-only` retains only entities whose required
evidence sides cannot all be explained within a common annotated gene by
eligible transcripts; different transcripts of that gene may explain the two
alternatives. `--solo-strand`
specifies the STARsolo alignment/transcript strand relationship used for splice
decisions and donor/acceptor naming. Unstranded splice evidence cannot orient
alternative donor versus alternative acceptor, so `unstranded` accepts only
explicit `junction` and `cassette` searches. Reverse and unstranded searches
must explicitly omit `terminal-tail`; the no-`--kind` default includes it. For
an unstranded entity, transcripts on either strand are eligible and the
annotation classifier does not claim an opposite-strand gap.

The uniform result always contains `capabilities`, `entities`, `components`,
`counts`, `terminal_anchors`, and `terminal_counts` tables. The two terminal
tables are empty when tails were not requested or no capable evidence passed
the predicates. `counts` contains nonzero sample/group rows; an absent row is
an exact logical zero only within the capability denominator declared for that
entity kind. Catalogue support is always an upper bound used for planning, not
a reported molecular count. The
`catalogue_min_component_route_upper_bound`,
`catalogue_sample_route_upper_bound`, and
`catalogue_donor_route_upper_bound` columns are all planning bounds formed from
component routes; `exact_samples` and `exact_donors` are the authoritative
observed recurrence. `entities.exact_umi_classes` is pooled support;
the two `*_alignment_umi_classes` columns audit informative-class orientations
and may overlap when one informative class has evidence in both orientations;
both-side-only classes are excluded. In `counts`, junctions
use `support_umi_classes`; alternative-site events use the neutral side-A/B
columns; cassette events use `include_only_umi_classes` and
`exclude_only_umi_classes`. Both use `both_umi_classes`, which is excluded from
`informative_umi_classes`. `terminal_counts.umi_classes` provides the same
archive-class unit per exact anchor. `entities.exact_cells` and `counts.cells`
count cells contributing at least one informative class, so a class supporting
both alternatives does not inflate those cell counts. Every `components` row
reports its marginal `exact_umi_classes` count.

`--min-support` uses catalogue support only as a safe pruning bound. An upper
bound never substitutes for a source-derived count. `--verify-content`
performs the stronger source check before execution. `find-events` emits every
retained sparse count row; it has no `--top` truncation or implicit
differential inference. `--max-candidates-considered` bounds every attempted
splice-event definition, including definitions later rejected by recurrence;
`--max-candidates` bounds retained routed definitions/clusters. Both fail rather
than truncate. `--max-routed-entries` bounds the sum of candidate-to-archive
target associations and routed chunk postings before the exact plan
materializes. `--max-exact-match-attempts` then bounds target-list expansion
while decoded molecule junctions are matched. Per-chunk `(entity, UMI class)`
hits are combined incrementally rather than retained as one row per match. The
annotation interval index restricts classification to local transcripts, and
`--max-annotation-comparisons` bounds the indexed transcript-comparison work.
Observed route, match-attempt, and annotation-comparison counts are recorded in
summary and provenance.

### Query options

| Command | Argument or option | Default | Description |
|---|---|---|---|
| `junction` | `<COLLECTION> <LOCUS>` | required | Collection and exact `chrom:donor-acceptor` junction |
| `region` | `<COLLECTION> <LOCUS>` | required | Collection and 0-based, half-open `chrom:start-end` window |
| `jset` | `<COLLECTION>` | required | Collection to query |
| `jset` | `--include <JUNCTION>...` | required | One or more inclusion-side exact junctions |
| `jset` | `--exclude <JUNCTION>...` | required | One or more exclusion-side exact junctions |
| `junction`, `jset` | `--min-support <N>` | `0` | Skip decoding when the collection-wide support upper bound is below this value |
| `find-events` | `<COLLECTION>` | required | Collection to search genome-wide |
| `find-events` | `--kind <KIND>` | all kinds | Repeat to select `junction`, `alt-acceptor`, `alt-donor`, `cassette`, or `terminal-tail` |
| `find-events` | `--design <TSV>` | sample = donor | Bind collection samples to biological donors |
| `find-events` | `--groups <TSV>` | all cells as `bulk` | Define the exact cell scope and cell groups |
| `find-events` | `--require-group <NAME>` | none | Require exact support in each named group |
| `find-events` | `--min-support <N>` | `2` | Minimum catalogue upper bound for every component |
| `find-events` | `--min-samples <N>` | `1` | Minimum archives with exact scoped support |
| `find-events` | `--min-donors <N>` | `1` | Minimum biological donors with exact scoped support |
| `find-events` | `--min-umi-classes <N>` | `1` | Minimum pooled exact archive UMI classes; `--min-umis` is an alias |
| `find-events` | `--min-side-umi-classes <N>` | `1` | Minimum exact archive UMI classes on each alternative-splicing side; `--min-side-umis` is an alias |
| `find-events` | `--min-group-umi-classes <N>` | `1` | Required exact archive UMI classes in every `--require-group`; `--min-group-umis` is an alias |
| `find-events` | `--terminal-cluster-bp <N>` | `25` | Maximum consecutive-anchor gap in a strand-aware terminal-tail cluster; `0` preserves one entity per anchor |
| `find-events` | `--max-terminal-events <N>` | `10000000` | Fail before tail-index or molecule decoding if rooted capability declarations exceed this event count |
| `find-events` | `--annotation <GTF_OR_AIC>` | — | Classify gaps from an uncompressed GTF or compiled AIC; requires assembly and label |
| `find-events` | `--assembly <NAME>` | — | Caller-declared assembly for the supplied annotation |
| `find-events` | `--annotation-label <LABEL>` | — | Release or immutable descriptive label for the supplied annotation |
| `find-events` | `--annotation-digest <DIGEST>` | — | Require an exact `blake3:<64 lowercase hex>` annotation digest |
| `find-events` | `--novel-only` | off | Keep only annotation-incompatible entities |
| `find-events` | `--solo-strand <MODE>` | `forward` | Alignment/transcript strand relationship; `unstranded` coalesces junction/cassette orientations into strand `.` before filtering |
| `find-events` | `--max-candidates <N>` | `100000` | Fail rather than truncate above this derived-candidate limit |
| `find-events` | `--max-candidates-considered <N>` | `1000000` | Fail when attempted splice definitions exceed this bound, even when recurrence filters reject them |
| `find-events` | `--max-routed-entries <N>` | `10000000` | Fail before exact-plan materialization when candidate/archive associations plus chunk postings exceed this bound |
| `find-events` | `--max-exact-match-attempts <N>` | `25000000` | Fail before expanding a hot candidate-target list when exact molecule matching would exceed this bound |
| `find-events` | `--max-annotation-comparisons <N>` | `10000000` | Fail when indexed annotation transcript-comparison work would exceed this bound |
| `junction`, `region`, `jset` | `--top <N>` | `5` | Per-sample cell rows to return; under `--format`, `0` returns all |
| `junction`, `region`, `jset` | `--explain` | off | Include per-sample routing decisions and compressed-byte estimates |
| all queries | `--verify-content` | off | Re-hash every source archive before querying |
| `junction`, `region`, `jset` | `--json` | off | Emit the command-specific JSON representation |

## Uniform output

All collection subcommands accept opt-in `--format text|tsv|json` and
`-o, --output PATH`. Omitting `--format` preserves the historical text and
`--json` presentations exactly. The legacy `--json` flag and uniform
`--format` are mutually exclusive.

Uniform JSON uses `gravlax.result-envelope.v1`; text and TSV carry the same
summary, provenance, schemas, and selection metadata. Results contain these
named tables:

| Command | Tables |
|---|---|
| `build` | `archives`, `source_io`, `source_sections` |
| `inspect` | `layers`, `chromosomes`, `archives`, `shape_route_blocks` |
| `junction` | `samples`, `cells` |
| `region` | `samples`, `cells` |
| `jset` | `requests`, `samples`, `cells` |
| `find-events` | `capabilities`, `entities`, `components`, `counts`, `terminal_anchors`, `terminal_counts` |

Sample and cell tables declare set semantics and keys; their physical order is
not contractual. Layer, chromosome, and request tables are explicit sequences.
For bounded cell output, the summary records the per-sample ranking comparator
and the cell table records exact aggregate `available_rows`, `emitted_rows`,
and `truncated` metadata.

Scientific totals live in `data.summary`. Provenance records source content
identities, every collection-layer root, invocation parameters, and the access
strategy. Rows are borrowed and streamed directly, with no output-only global
sort or result-sized second copy.

`--output` requires `--format`. Its parent must already exist; an occupied
destination is rejected before expensive work. The complete result is staged
beside the destination and installed atomically without overwriting an
existing file. Machine-readable stdout contains only the selected result;
diagnostics use stderr.

Collection files are derived routing indexes rather than replacements for
their source archives. Keep every source path reachable for inspection and
queries.
