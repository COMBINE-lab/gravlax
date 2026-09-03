---
title: aie cohort
description: Exact cross-archive event, graph, and replicated transcript-end analyses.
---

The `events` and `splice-graph` modes support the opt-in uniform result
contract with `--format text|tsv|json` and atomic no-clobber
`-o/--output`. Omitting `--format` preserves every historical default,
`--tsv`, and `--json` byte contract. Diagnostics remain on stderr. Uniform
output cannot be combined with `events --sparse-dir`, because the sparse
directory is a separate multi-file artifact rather than part of the result
transaction.

## `events` — descriptive event tables

`aie cohort events` builds the union of coordinate-defined splice-event
catalogues from two or more archives, retains events present in at least
`--min-samples` catalogues, and executes the same exact class-level reducer in
each archive.

```sh
aie cohort events chr1:0-248956422 \
  --sample donor-a=donor-a.aie --sample donor-b=donor-b.aie \
  --groups donor-a=donor-a-groups.tsv --groups donor-b=donor-b-groups.tsv \
  --min-support 2 --min-samples 2 --min-informative 10 \
  --min-row-informative 10 --json
```

Sample and group arguments use `ID=PATH`; IDs must be unique and every group
ID must name a supplied sample. Samples without a group map emit bulk rows.
Stamped archives must have the same genome signature, and stamped and
unstamped identities cannot be mixed.

The output schema `gravlax.cohort.events.v1` preserves caller sample order and
deterministic event order. A component absent from one archive is represented
by `present: false`; evidence is not imputed. V1 intentionally reports
descriptive counts and usage only. Samples, donors, chemistries, and depths are
not treated as interchangeable replicates for a convenience p-value.

`--min-row-informative N` is an exact output predicate, not a statistical
contrast. With grouped samples, every group row in every sample must contain at
least `N` conservative informative classes (`include_only + exclude_only`);
with ungrouped samples, the same requirement applies to every bulk row. An
event retained by the cross-sample catalogue recurrence may pass on observed
evidence even when `present: false` in one local catalogue; no evidence is
imputed. Gravlax evaluates the smallest scope first and removes denominator
failures before decoding deeper samples, but restores caller sample order in
the output. The default `0` disables this pushdown and preserves the original
output byte-for-byte.

When pushdown is active, `planning.candidate_events` is the recurrent catalogue
size and `planning.post_row_candidate_events` is the number surviving the row
threshold. The usual summed `--min-informative` filter is then applied to those
survivors.

| Option | Default | Description |
|---|---|---|
| `--sample <ID=ARCHIVE>` | required, ≥2 | Named archive; repeat |
| `--groups <ID=TSV>` | — | Optional barcode/group map; repeat |
| `--event-type <TYPE>` | all | Repeatable event selector |
| `--min-support <N>` | `2` | Per-component support within each catalogue |
| `--min-samples <N>` | `2` | Catalogues in which all components are present |
| `--min-informative <N>` | `1` | Minimum informative classes summed across samples |
| `--min-row-informative <N>` | `0` | Minimum `include_only + exclude_only` in every emitted sample/group row; `0` disables pushdown |
| `--max-events <N>` | `100000` | Hard union-catalogue limit |
| `--gtf <GTF/AIC>` | — | Count-neutral gene/strand labels |
| `--tsv` / `--json` | off | Long sample/group table or versioned object |
| `--format <FORMAT>` | — | Uniform `samples`, `events`, `components`, and `counts` tables |
| `-o, --output <PATH>` | stdout | Atomically publish uniform output without replacing a path |

The uniform result schema is `gravlax.cohort.events.result.v1`. Samples retain
caller order; event, component, and count tables declare set keys without a
gratuitous physical ordering promise. A grouped sample emits both its exact
total row and every group row. For each count category, the total is exactly
the sum of those group rows. Ungrouped samples emit one bulk row. Archive
content identities are bound to sample IDs, preserving legitimate replicate
multiplicity even when two archive byte streams are identical.

## `splice-graph` — sample-level path contrasts

`aie cohort splice-graph` constructs a common strand-aware junction graph,
reduces each archive independently, emits complete zero-explicit sample × path
and sample × edge matrices, and can test path usage across biological
replicates:

```sh
aie cohort splice-graph chr1:45550000-45571000 \
  --design experiment.tsv --contrast control:treated --json
```

The design is a strict tab-separated file. Its header and four columns are
exactly:

```text
sample	condition	archive	cells
C1	control	archives/c1.aie	cells/c1.txt
C2	control	archives/c2.aie	.
T1	treated	archives/t1.aie	cells/t1.txt
T2	treated	archives/t2.aie	.
```

Paths are resolved relative to the design file. `.` selects every archived
cell. Sample and condition IDs may contain ASCII letters, digits, `.`, `_`, or
`-`. Sample IDs and resolved archive paths must be unique, so the same archive
cannot masquerade as two replicates. The two design conditions must exactly
match `--contrast A:B`; effects are fitted usage in B minus fitted usage in A.
Use `--counts-only` instead of `--contrast` to emit matrices without inference.

The common coordinate-edge catalogue retains an edge when at least
`--min-edge-samples` local catalogues meet `--min-support`. Gravlax then reduces
all observed evidence against those edges in every sample. Falling below the
local discovery threshold does not force a zero: the reducer still recovers
observed evidence. The resulting exact matrix explicitly records genuine
zero-count rows, while low-depth samples remain reported and are marked
ineligible for inference.

For each strand and exact junction-set path fragment, the model compares that
path with all other common-graph path fragments on the same strand. A sample
enters only when that strand has at least `--min-sample-umis` fragments. The
test is a beta-binomial likelihood-ratio contrast with condition-specific means
and a shared alternative concentration; its unit is one biological sample,
never a cell or molecule. Paths must pass the replicate, path-UMI, supporting-
sample, and nonzero-comparator thresholds. P values use a one-degree-of-freedom
asymptotic likelihood-ratio tail and are Benjamini–Hochberg adjusted
across tested paths in the locus. These are path-fragment usage tests, not
complete-transcript tests or a substitute for an experimental design with
biological replication.

The output schema is `gravlax.cohort.splice-graph.v1`. It records the complete
design, reference digest, graph, thresholds, per-sample scope and planning
counts, every path/edge row, fitted model parameters, all tested paths, and the
reason each untested path was skipped.

| Option | Default | Description |
|---|---|---|
| `--design <TSV>` | required | Strict sample/condition/archive/cells design |
| `--contrast <A:B>` | required for inference | Ordered two-condition contrast; effect is B minus A |
| `--counts-only` | off | Emit exact matrices without inference; conflicts with `--contrast` |
| `--min-support <N>` | `1` | Local catalogue support for edge recurrence |
| `--min-edge-samples <N>` | `2` | Catalogues required for a common coordinate edge |
| `--min-sample-umis <N>` | `10` | Same-strand graph fragments required for sample eligibility |
| `--min-replicates <N>` | `2` | Eligible biological samples required per condition |
| `--min-path-umis <N>` | `5` | Path UMIs required across eligible samples |
| `--min-path-samples <N>` | `2` | Eligible nonzero samples required for a path |
| `--max-paths <N>` | `100000` | Hard union-path limit; never truncates |
| `--json` | off | Emit the full versioned result object |
| `--format <FORMAT>` | — | Uniform sample, graph, matrix, test, and skipped-test tables |
| `-o, --output <PATH>` | stdout | Atomically publish uniform output without replacing a path |

The uniform schema is `gravlax.cohort.splice-graph.result.v1`. Named tables are
`samples`, `nodes`, `edges`, `paths`, `path_counts`, `edge_counts`, and, for a
contrast, `tests` and `skipped_tests`. Matrix tables are complete and
zero-explicit; node/edge/path IDs are sequential topology identifiers, while
the scientific set tables use explicit keys. Hard candidate caps remain
errors—none is represented as successful truncation. The design file and any
cell scopes are content-bound, and each already-open archive contributes its
rooted identity to provenance.

## Replicated 3′-end analyses

`transcript-ends` and `polyasite-mixture` produce multi-file scientific
directories rather than replacing those artifacts with stdout. Their optional
uniform interface is therefore an operation report:

| Option | Default | Description |
|---|---|---|
| `--design <TSV>` | required | Strict `sample<TAB>condition<TAB>archive<TAB>groups` design; paths are resolved relative to the design file |
| `--gtf <GTF/AIC>` | required | Gene and terminal-region annotation; observed sites still come from archive evidence or PolyASite |
| `--genome <FASTA>` | required | Stamped reference used to identify or reject internal-priming sequence |
| `--polyasite <BED/BED.GZ>` | required | Same-assembly PolyASite catalogue |
| `--group-contrast <A:B>` | required | Paired within-donor group contrast; effect is B minus A |
| `--out-dir <PATH>` | required | New output directory; an existing path is rejected |
| `--site-gap <BP>` | `24` | Maximum distance for clustering endpoint sites or merging adjacent catalogue coordinates |
| `--tail-extend <BP>` | `2000` | Extension past annotated 3′ regions in transcript direction |
| `--min-site-umis <N>` | `10` | Donor support required for an endpoint or fitted site |
| `--min-site-samples <N>` | `3` | Donors required for a site to enter the recurrent catalogue |
| `--motif-min-samples <N>` | `4` | `transcript-ends` only: donors required for a motif-only high-confidence site |
| `--min-group-gene-umis <N>` | `20` | Per donor, minimum gene-site UMIs required in each contrasted group |
| `--min-samples <N>` | `6` | Eligible paired biological samples required for a gene test |
| `--min-distal-umis <N>` | `20` | Most-distal-site UMIs required across eligible sample/group rows |
| `--max-sites <N>` | `1000000` | Hard site-catalogue limit; exceeding it returns an error rather than truncating |
| `--shuffle-seed <N>` | — | Shuffle the two group labels within donors reproducibly for a negative-control analysis |

| Option | Default | Description |
|---|---|---|
| `--report-format <text|tsv|json>` | — | Stream a typed summary plus normalized scientific and artifact tables |
| `--report-output <PATH>` | stdout | Atomically install the report without replacing an existing file; requires `--report-format` |

Without these flags, both commands preserve their existing directory files and
legacy pretty-JSON stdout. A report is emitted only after the directory's
`summary.json` completion marker. The directory itself is incremental rather
than transactional; the optional report file is a separate atomic sidecar and
cannot roll the directory back. The [query and evidence
guide](/gravlax/cli/query/) gives the scientific interpretation, named tables,
and provenance rules for both modes.
