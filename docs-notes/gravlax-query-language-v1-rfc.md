# Gravlax Query Language v1 specification

Target: Gravlax 0.2.0. Implemented and validated; the release has not been published.
This document supersedes the original proposal and both review rounds.

## Motivation and scope

GQ combines a pipeline over evidence units with explicit quantifiers over archived
observations. The native executor is justified by scientific semantics and direct
archive access, not by a claim that these questions cannot be expressed relationally.
The seqproc influence is named domain objects, bounded composition, version negotiation,
and a validated execution graph; Boolean OR is commutative, not ordered fallback.

## Evidence model

An archive contains serialized **records** associated with corrected cells and exact
raw-UMI **classes**. A unique-evidence record groups chains at a locus; multimapping
signatures in one record may have different direct anchors. A class is not a proven
physical molecule: raw-UMI collisions remain possible. IDs include a sample namespace.

A **chain** retains junction-invariant unique alignment observations, represented by
one or two geometries plus an accepted-observation weight. A **representative** is a
stored geometry, not a weight-one replacement for the whole chain. In compact archives,
the extrema minimize/maximize `(start, Reverse(end))`. They do **not** independently
bracket all ends. Full-fidelity singleton geometries retain exact per-geometry weights.

A **multimapping signature** has one weight and a grouped set of alignment alternatives,
including its primary. Primary alignments are not unique alignments. A **cell** contains
all its records. No implicit coordinate, unit, or strand-frame conversion is permitted.

## Language by example

```text
header { gq = 1, assembly = "GRCh38" }
let locus = g[1:198600000..198800000]
let RA = overlaps(g[1:198696711..198696909:-], min: 12bp)
let RO = j[1:198692373..198703297:-]

from @pbmc.records
|> within locus
|> derive {
     ra = any unique { RA },
     ro = any unique { RO },
     joint = any unique { RA & RO }
   }
|> tally {ra, ro, joint} by {sample}
```

`any unique { A } & any unique { B }` permits different observations of a record
to witness A and B. `any unique { A & B }` requires one observation to witness both.
On `.classes` or `.cells`, write `any record { ... }` to descend into records.
`within` selects a positive witnessed population; it never crops the children of a
selected class or cell. The universe uses retained unique representatives and retained
multimapper alternatives. Unwitnessed omitted geometries do not silently enlarge it.

## Types, quantifiers, and uncertainty

Scalar types include Truth, String, Number, Count (UInt64), Bases, null, homogeneous
lists, and read-support Bounds. Coordinate values are Region/RegionSet, Pattern,
PatternSet, GeneFeature and TranscriptFeature under the declared assembly.

Use `12bp` for base quantities and `12u64` for exact Count literals (for example
`reads.total >= 2u64`). Unsuffixed numbers are Number values. `fraction(n,
denominator: d)` explicitly converts matching numeric inputs to a floating ratio;
a zero or null denominator gives null, never zero.

| Expression | Domain / meaning |
| --- | --- |
| `any unique { P }`, `all unique { P }` | Accepted uniquely mapped observations, interpreted through the quotient |
| `any stored { P }`, `all stored { P }` | Literal retained representatives; diagnostic |
| `any multimap { all alternative { P } }` | Existential signature with universal retained alternatives, including primary |
| `any record { P }`, `all record { P }` | Complete records of the current class/cell |
| `nonempty(D)` | Explicit nonempty-domain test |

`all D { P } == !any D { !P }`, including empty domains. Existing CLI non-vacuous
`all-placements` corresponds to `nonempty(D) & all D { P }`, not bare `all`.
`placement` is reserved and produces a migration diagnostic, not an alias.

Truth has values `true`, `false`, and `unknown`. Negating unknown yields unknown;
false AND unknown is false; true OR unknown is true. `P | !P` need not be true.
`where` retains only true and reports dropped unknowns. `is_unknown` and `is_null`
are distinct operations. None of these answers establishes biological absence.

### Sound proofs under compact retention

Junction and fixed-path predicates are chain-invariant under the archive's chain
grouping model. Retained witnesses establish an existential true or universal false.
At multiplicity at most two all observations are represented. A singleton geometry
has exact weight. Otherwise, evaluate omitted observations conservatively:

- Starts lie between retained starts. Containment/disjointness proves endpoint tests.
- Ends are bracketed by the retained ends only when all starts agree. Otherwise
  use structural lower bounds and the coordinate limit, not the two stored ends.
- Invariant introns cannot gain aligned overlap. First-block start bounds and fixed
  internal blocks provide safe overlap envelopes; unknown terminal extents remain unknown.
- Compose atom proofs in three-valued logic, without assuming arbitrary expressions
  are monotone. This is sound but not a complete decision procedure for all tautologies.

Counterexample to the review's proposed end bracketing: `[100,120)`, `[105,200)`,
`[110,130)` retains the first and last. The omitted end 200 is outside both.

`reads.unique`, `reads.multimapping`, and `reads.total` count accepted observations
of the current record/class/cell. Multimapping weight is counted once per signature,
not once per alternative. `sum(reads.total)` after a filter sums selected units'
totals, not just matching reads.

`support_reads(P)` is record-scoped and returns `{lower, upper, status}`. Unique
bounds include known retained witnesses plus conservative omitted weight; a signature
contributes its weight to the lower bound when all alternatives satisfy P and to
the upper bound when at least one might. `status` is `exact` iff bounds coincide.
Projections `.lower`, `.upper`, `.status` are available.

## Geometry, paths, and strand frames

Coordinates are 0-based, half-open; `end_in` tests the exclusive end boundary itself.
Literal strands denote **alignment** strand; omitted strand is unconstrained.
Quoted chromosome names support punctuation. There is no automatic `chr` aliasing.

```text
g[chr1:100..200:+]
j[chr1:125..225:+]
overlaps(G, min: 12bp, blocks: one)
overlaps(G, min: 12bp, blocks: total)
start_in(G)
end_in(G)
near(J, left: 2bp, right: 0bp)
terminal(G)
```

Overlap uses unions: `one` is the maximum target-union overlap of one aligned block;
`total` is aligned-block-union overlap with the target union. Duplicated annotation
exons cannot multiply overlap. `terminal` requires record-level tail capability;
placing it inside a unique/stored/alternative quantifier is a type error.

`header {gq=1, assembly="GRCh38", library="opposite"}` declares the alignment-to-
transcript strand relationship. `tx_strand(-)` (also `tx_strand("-")`) and annotation
features require this rule. Use `same` for matching frames. The rule is recorded in
provenance, not inferred from a filename, assay label, or strand on a literal.

`>>` means consecutive junctions; `~>` permits intervening junctions. Bare composition
requires agreeing literal strands and follows that strand's orientation. Explicit
`path(order: genomic) { ... }` uses increasing coordinates. `path(-) { ... }` uses
descending transcript order and requires a library rule. A conflicting written order
is an error. Lowering to an increasing genomic matcher is explicit, not syntax repair.
Mixed chromosomes/strands and mixed `>>`/`~>` in one path are rejected in v1.

## Definitions and pinned annotations

Local nonrecursive functions permit inference:

```text
fn supported(exon, minimum) = overlaps(exon, min: minimum)
```

Exported coordinate/predicate signatures declare unit, assembly, and strand requirements:

```text
export fn supported(exon: Region<GRCh38,alignment>, minimum: Bases)
  -> Predicate<Alignment,GRCh38,alignment> = overlaps(exon, min: minimum)
```

No recursion, mutable state, loops, arbitrary filesystem calls, joins or general UDFs.
Definitions expand under a finite node budget; arguments are values, not source fragments.

Project annotation resources must have assembly and an immutable annotation label;
their exact consumed bytes are digested. `gene(@gencode,"PTPRC")` resolves one gene,
with `.span`, `.exons`, `.junctions`. `transcript(@gencode,"ENST...")` also offers
`.junction_path`. Ambiguous or missing labels and unpinned resources fail. A gene's
junction union is not an invented transcript. Annotation strands are transcript-frame
and convert through the declared library rule.

## Pipeline and results

`within` is required as the first stage (`within all` is explicit). `derive` names
typed observations. `tally` accepts Truth expressions and emits observed combinations
only. `summarize` supports `count`, `count_distinct`, `sum`, and `count_true`,
`count_false`, `count_unknown`. Counts are exact UInt64, not floating-point reads.

Grouping uses base sample/cell metadata, retaining a missing-value group. This v1
restriction makes pre-filter denominators well-defined. Metadata columns are declared
as string, number, truth (Boolean inputs), or homogeneous list types. Explicit null
values are nullable; absent columns require an optional declaration. No coercions or
regular expressions. Comparisons, `in`, `not in`, arithmetic, `!`, `&`, `|`, and
parentheses are supported. Numeric division by zero yields null; count arithmetic
is checked and fails on overflow or division by zero.

`select` projects values. `sort {field, ...}` is ascending, with a deterministic full-row
tie-break. `take N` requires sort, occurs after result computation, and discloses
truncation. It does not reduce the scientific aggregation population or bypass budgets.

GQ uses the existing result envelope and a `results` table with schema
`gravlax.gq.table.v1`. Tally has grouping fields, ordered `<name>_state` string columns,
and UInt64 `count`. Unknown Truth is the string `"unknown"`, never bare null.
Metadata includes unit, ordered state fields, source/population totals per group,
dropped unknown counts, execution paths, and pagination. Read bounds retain their status.

## Compiler and execution contract

```sh
aie gq validate examples/gq/cd45.gq
aie gq explain examples/gq/cd45.gq --bind pbmc=sample.aie
aie gq run examples/gq/cd45.gq --bind pbmc=sample.aie --output cd45.json
```

Unbound validate checks syntax/types, explicitly reporting that resources were not
checked. Bound validate and explain check capabilities/routes without decoding records.
Run uses native archive access, not per-predicate subprocesses. No query rewrites
that discard false tally states, intersect independent class witnesses, or crop
children across unit boundaries are enabled.

`--project` reuses project aliases. `--bind name=archive.aie` explicitly associates the
header's assembly assertion with that archive; this is not proof of its genome identity.
A binding may instead point to a version-1 JSON federation with `assembly` and
`archives: [{sample, path, metadata}]`. Relative paths resolve beside that manifest.
Duplicate sample names, paths, or content roots are rejected. Sample namespaces are
independent; a federation is not a way to split one UMI class across different samples.

`--metadata` supplies `{columns: {"cell.type": {type:"string",optional:true}},
cells: {sample: {barcode: {"cell.type":"Astro"}}}}`. Sample-wide metadata can live
on federation members. No implicit donor inference occurs.

Initial positive routing and unit closure are separate. Class indexes close selected
classes; cell closure currently needs a full scan. Missing complete routes require
`--allow-full-scan`. Exact source-unit denominators grouped by cell metadata also
require a disclosed full scan for record/class queries; sample-only denominators use
archive totals. No hidden fallback. Chunk, record, expression-step, terminal-event,
group and output-row budgets fail instead of truncating scientific counts.

Source text, normalized logical plan, bound metadata/annotations, archive roots,
engine version, coordinate frame, strand rule and execution policy are recorded.
Output files publish atomically without overwriting an existing result.

## Junction enumeration

```text
from junctions(@pbmc, within: g[1:198600000..198800000])
|> support(unit: class, by: {sample})
```

v1 enumerates unique-evidence junctions with both boundaries in one target interval,
then counts distinct supporting records/classes/cells under that same domain. It
uses one bounded full scan, not a scan per candidate. MM-only candidate enumeration
is outside this initial contract. Counts are exact retained-evidence support; no
catalogue statistic is silently relabeled an upper bound. Candidate/group budgets
apply; explain discloses the required full-scan allowance.

## Normative decisions and boundaries

1. Three domains, no implicit primary/unique policy and no `placement` alias.
2. Classical universal quantification; explicit nonempty requirements.
3. Sound start/junction proofs; never assume independently bracketed ends.
4. Literal alignment strand; explicit transcript conversion and path orientation.
5. Accepted-observation weights and predicate support bounds are distinct.
6. Annotation unions and transcript paths are typed, pinned and unambiguous.
7. Grouped Truth tables retain unknowns and disclose denominators.
8. Full scans require execution authority; budget exhaustion is an error.
9. Native execution, no scientifically unsafe optimizer rewrites.
10. Bounded unique-junction enumeration is included; general graph/event discovery,
    arbitrary joins, statistical models, sequence/quality queries, regex paths,
    posterior assignment and automatic assembly migration are not v1 features.

The original proposal is preserved separately for review history. Parser selection
and acceptance work are recorded in `gq-implementation-plan.md`.
