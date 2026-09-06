# Gravlax Query Language: original proposal (superseded)

Status: design proposal, not an implemented textual language. Working short name
**GQ**, extension `.gq`. This proposal accompanies the cooccur execution work on
`astra-improvements`; approval of the syntax should precede parser implementation.

## Recommendation

Use two small, composable algebras:

1. A relational pipeline for choosing evidence units, selecting metadata,
   deriving observations, filtering, and aggregating.
2. A typed evidence-pattern algebra inside explicit quantifiers for geometry,
   ordered paths, and retained-placement alternatives.

Do not make users write JSON execution plans, and do not turn biological evidence
into arbitrary SQL joins. The language should make the scientific distinction
between *same placement*, *same record*, *same raw-UMI class*, and *same cell*
as easy to express as a region overlap.

The useful seqproc lesson is concise, named domain objects compiled into a
validated execution graph. Its EFGDL 2 source documentation also demonstrates
version negotiation, bounded layout algebra, explicit ambiguity policies, and
`validate`/`explain`/`run`. Borrow those engineering principles, not its literal
FASTQ grammar. In particular, EFGDL ordered fallback is appropriate for read
recognition; GQ Boolean disjunction must be commutative. PRQL is a useful precedent
for the outer pipeline, but GQ should compile to Gravlax's evidence engine, not SQL.

Sources inspected: [official seqproc repository](https://github.com/COMBINE-lab/seqproc),
especially `website/src/content/docs/efgdl/{version-2,layout-algebra,annotations-and-ambiguity}.md`,
and the [PRQL language book](https://prql-lang.org/book/).
The supplied bioRxiv v2 URL was rate-limited; the offered local PDF path was not
visible from the active terminal on nomad01. The paper-specific comparison remains
pending; this proposal does not claim to have read that PDF.

## What a query looks like

Illustrative CD45 query, using the coordinates in our existing case study:

```text
header { gq = 1 }

let locus = g[1:198692346..198703372:-]
let RA = overlaps(g[1:198696711..198696909:-], min: 12bp)
let RO = j[1:198692373..198703297:-]

from @pbmc.records
|> within locus
|> derive {
     ra = any placement { RA },
     ro = any placement { RO },
     joint = any placement { RA & RO }
   }
|> tally {ra, ro, joint} by {sample, cell.type}
```

This returns the observed truth-state combinations and their record counts,
including false and unknown states, per sample and cell type. It does not silently
discard the denominator. `@pbmc` is a bound archive or federation resource, not an
uncontrolled filename or a URL to execute. A cell metadata resource supplies
`cell.type`; an unavailable column is a validation error.

To ask whether the same *cell* has both forms, change the unit and quantify records:

```text
from @pbmc.cells
|> within locus
|> derive {
     ra = any record { any placement { RA } },
     ro = any record { any placement { RO } }
   }
|> tally {ra, ro} by {sample, cell.type}
```

These examples express observations, not a claim that a physical molecule contains
both forms, nor a classification of cells from an uncalibrated absence of coverage.

The four expressions below are deliberately different:

```text
any placement { A & B }                         # one placement witnesses both
any placement { A } & any placement { B }       # possibly different placements of this record
any record { any placement { A & B } }          # one record has a joint placement witness
any record { any placement { A } }
  & any record { any placement { B } }          # possibly different records of this cell/class
```

There is no implicit promotion from a placement predicate to a record/cell
predicate. Omitting the quantifier is a useful type error, not shorthand that
changes scientific meaning.

## Patterns and functions

| Form | Meaning |
| --- | --- |
| `g[chr1:100..200:+]` | Typed genomic region, 0-based half-open |
| `j[chr1:125..225:+]` | Exact intron-boundary pattern, smaller coordinate first |
| `overlaps(G, min: 12bp)` | At least 12 bases in one aligned block of one placement |
| `start_in(G)`, `end_in(G)` | Genomic left start / exclusive right end boundary |
| `near(J, left: 2bp, right: 0bp)` | Explicit per-boundary tolerance; no implicit merging |
| `J1 >> J2` | Consecutive observed junctions on one placement |
| `J1 ~> J2` | Ordered subsequence; intervening observed junctions allowed |
| `terminal(G)` | Retained terminal-tail cleavage anchor in G; record-level |
| `any placement { P }` | At least one retained placement satisfies P |
| `all placement { P }` | Nonempty retained placement set, every member satisfies P |
| `any record { P }` | At least one record of the current cell/class satisfies P |

Junction literals have a pattern type; matching that pattern inside a placement
predicate is its canonical interpretation. Region literals are values and require
an operation such as `overlaps`; this avoids conflating anchor selection with
aligned-block overlap. Functions may accept typed region or pattern arguments.

`>>` and `~>` use **transcriptional order** when a strand is specified. On minus
strand, a multi-junction path is therefore written in descending genomic order;
each individual `j[...]` still uses smaller-first coordinates. Unstranded paths
use increasing genomic order. Mixed-strand or cross-contig paths fail validation
in v1. The current CLI `path:`/`subpath:` interface uses increasing genomic order
on both strands: the future GQ compiler must make this lowering explicit.

Precedence is path composition, `!`, `&`, `|` (tightest to loosest);
quantifier braces establish scope. Pattern composition takes junction patterns,
not Boolean expressions; invalid combinations fail type checking. Pattern
alternation uses `|` and is unordered set union, never first-match fallback.
Path composition is left-associative. Parentheses are always available.

Named, nonrecursive parameterized definitions provide reusable templates:

```text
fn inclusion(exon, left, right) =
  overlaps(exon, min: 12bp) & (left >> right)
```

Types can be inferred here; the arguments resolve to Region, JunctionPattern,
JunctionPattern and the result to PlacementPredicate. Explicit type annotations
are allowed in exported definitions. This is not an arbitrary computation language:
no recursion, mutable variables, loops, filesystem calls, or general-purpose UDFs.

An observed path is not a reconstructed full transcript. `>>` never creates an
edge from junctions witnessed on separate reads. v1 has no unbounded regex `*` or
recursive graph traversal. A finite path matcher/automaton is sufficient for its
fixed, ordered patterns without exponential expansion into all alternatives.

## Units, universes, and unknowns

Sources expose `.records`, `.classes`, and `.cells`. A class means a corrected
cell plus an exact raw UMI value, not error-corrected physical-molecule identity.
Do not call this source `.molecules` or silently collapse one-mismatch UMI edges.
Every unit ID includes its archive/sample namespace.

`within G` defines a positive **population**, not a destructive crop of its
children: records with witnessed aligned-block overlap of G; classes/cells with
at least one such record. `within junction J` is a more selective positive universe.
The unit's children remain logically complete for subsequent quantifiers. An
explicit region constraint inside a quantifier can narrow those children.
This preserves the current cross-chunk UMI-class closure contract.

For cell queries, this is a covered-cell denominator, not all annotated cells.
`within all` requests the complete metadata-selected population and requires an
explicit full-scan allowance where necessary. Results report both population size
and source cell-scope size. Units whose universe membership is not witnessed are
not candidates; explain/results disclose the witnessed-universe contract and
whether its membership is complete. It must not be called the biological universe.

Evidence predicates return `Truth = true | false | unknown`:

- true: positive retained evidence, or a valid finite-scope logical proof;
- false: refuted within the declared retained-evidence scope;
- unknown: the archive representation does not resolve the proposition.

Use three-valued logic: `!unknown = unknown`, `false & unknown = false`,
`true | unknown = true`. `all` is non-vacuous: it is false on an empty domain.
Omitted compact geometry can leave a geometry-sensitive unwitnessed existential
or unrefuted universal unknown; exact junction paths remain invariant within a
compact chain. No operator turns lack of retained evidence into biological
absence. An optional tail observable is not proof of a complete transcript end.

`where P` keeps true rows and reports dropped unknown counts. `where is_unknown(P)`
inspects unresolved observations. `tally {P, Q}` preserves all joint truth states,
including unknown, and is the preferred contrast-building operation.
`count()` counts the current unit; reads and weighted UMIs are not synonyms.
Named aggregation supports `count()`, `count_distinct(cell)`, `count_true(P)`,
`count_false(P)`, and `count_unknown(P)`. Fractions must name a denominator, with
zero denominators producing null, not zero.

Missing capability is a compile/plan error by default. Representation-limited
unknown is distinct from an absent entire terminal stream. v1 should not offer
an implicit "missing means false" mode. Explicit optional-capability expressions
can be added after their result/provenance contract is settled.

## Protocol-aware composition

This record query is meaningful if terminal-tail evidence is present:

```text
from @sez.records
|> within gene_span
|> where any placement { splice_pattern } & terminal(distal_window)
|> summarize {records = count(), cells = count_distinct(cell)}
     by {sample, donor, cell.type}
```

It asserts record-level co-occurrence. This must fail validation:

```text
any placement { splice_pattern & terminal(distal_window) }
```

Reason: retained terminal events are linked to a serialized record, not to a
specific placement. More expressive syntax must not manufacture missing linkage.
Likewise, `all placement` is not "all possible alignments" or posterior probability.

## Federation and reproducibility

Exactly the same query runs over one archive or a federation. The resolver binds
resource aliases, assembly, archive roots, metadata schema/digests, and pinned
annotation resources before planning. Coordinates from incompatible assemblies
cannot be combined. Contig alias conversion must be declared and validated.

Filter sample/cell metadata before evidence matching; execute predicates near each
archive; merge only typed results. Namespace all barcodes/classes by sample and
reject accidental duplicate archive membership. Donor is explicit metadata, not
inferred from a filename or silently equated with sample.

Named annotations can resolve features, e.g. `gene(@gencode, "PTPRC")`; ambiguous
gene labels and unpinned resources are errors. Existing annotation-comparison,
transcript-compatibility, and discovery outputs should be bindable as typed tables
in v1. Do not reimplement their scientific inference inside the textual parser.
Direct general-purpose joins to those tables are beyond v1; first support
key-validated metadata/feature lookup.

The result records language version, source text digest, normalized logical plan
digest, engine version, content-bound resources, coordinates, placement policy,
retention capabilities, unit identity, and execution limits. Optimizer hints affect
the physical plan only; they cannot change tolerances, ambiguity, or unit semantics.

## Compiler and execution contract

`aie gq validate query.gq`, `aie gq explain query.gq`, `aie gq run query.gq` are
proposed entry points, with resource bindings resolved from the existing project
system or explicit `--bind` arguments. Reuse the existing uniform result tables
and plan runner rather than creating a second workflow framework.

The compiler stages are: parse with source spans; bind names/resources; type and
capability check; normalize/hash common expressions; build a logical evidence
plan; choose safe routes/projections; execute; emit typed results and provenance.
The common plan is the API shared by textual queries, existing CLI flags, and
future programmatic clients. Do not make shell argument generation the long-term IR.

Types include Region<Assembly>, JunctionPattern<Assembly>, PathPattern<Assembly>,
Predicate<Unit>, Truth, and Table<Row>. Predicate effects include required streams,
geometry completeness, placement policy, and closure requirements. A geometry
filter cannot be moved across a quantifier or unit conversion without a proof
that the meaning is unchanged.

Explain must show resolved coordinates, scientific unit, quantifier nesting,
required capabilities, initial routes, data-dependent closure, full-scan reasons,
unknown semantics, projection, cache bounds, and which optimizations were used.
It must distinguish estimates from actual I/O, never label an initial route count
as the final UMI-class cost. The current cooccur explain implements a first subset.

High-value optimizations include shared exact-junction lookup, common subexpression
elimination, once-per-geometry path matching, early metadata filters, cached shape
decoding, and fused aggregation. For a selected-only query, additional safe route
pruning becomes possible. For `tally` over the full universe, throwing away false
patterns is incorrect. Intersecting chunk postings for a class/cell conjunction
is also incorrect when its witnesses can live in different chunks.

Finite compilation budgets apply to AST depth/nodes, expanded definitions, path
states, and predicate count. Runtime budgets apply to chunks, decoded records,
geometry cache, terminal events, groups, witnesses, and result rows. Exceeding a
budget is an error, not hidden truncation. Pagination is separate from scientific
aggregation; `take` is allowed only after a deterministic order and exposes truncation.

## The v1 boundary

Ship:

- versioned documents, named bindings, typed parameters, bounded nonrecursive
  predicate functions, and source-located diagnostics;
- archive/federation records, exact-UMI classes, and cells with explicit closure;
- region/junction/endpoint predicates, minimum block overlap, independent left/right
  junction tolerance, fixed consecutive/subsequence paths, record-level terminals;
- explicit any/all placement and any record quantifiers with three-valued logic;
- `within`, metadata `where`, `derive`, evidence `where`, `tally`, grouped
  `summarize`, projection, deterministic sort and disclosed result pagination;
- pinned feature resolution, typed existing result resources, validate/explain/run,
  bounded execution and complete provenance.

Do not ship in v1: recursive graph queries, arbitrary joins/UDFs, unbounded regex,
free-form statistical modelling, implicit probabilistic assignment, automatic
annotation/assembly migration, or new mandatory archive indexes. Querying sequences,
qualities, or placement-linked tails remains impossible without the required data.

Build in four increments after syntax review:

1. Parser, formatter, typed AST, source-span diagnostics, golden example corpus.
2. Shared logical IR and local executor; differential equivalence to current CLI;
   adversarial tests for quantifier scopes, unknowns, paths, and retention modes.
3. Cell/class closure and federated execution with namespace-safe typed merging;
   single-archive vs partitioned-federation equivalence tests and memory bounds.
4. Parameterized templates, grouped summaries, feature binding, stable explain and
   provenance; benchmark real workflows before enabling new physical rewrites.

The acceptance criterion is not merely parsing these examples: the compiler must
reject the tempting but invalid same-placement terminal query, preserve unknowns,
and produce identical logical results across scalar/compiled execution, chunk
layouts, and federation partitions of the same nonduplicated samples.

## Grammar sketch and diagnostics

The parser should use an ordinary expression parser plus a small pipeline grammar,
not a grammar that enumerates every CLI subcommand. A deliberately abbreviated
grammar (literals, argument lists and metadata expressions omitted) is:

```ebnf
document    = header, { binding | function }, pipeline ;
binding     = "let", name, "=", expression ;
function    = "fn", name, "(", parameters, ")", "=", expression ;
pipeline    = "from", source, { "|>", stage } ;
stage       = "within", universe
            | "where", expression
            | "derive", named_fields
            | "tally", fields, [ "by", fields ]
            | "summarize", named_fields, [ "by", fields ]
            | "select", fields | "sort", order_fields | "take", integer ;
expression  = disjunction ;
disjunction = conjunction, { "|", conjunction } ;
conjunction = unary, { "&", unary } ;
unary       = "!", unary | path ;
path        = primary, { (">>" | "~>"), primary } ;
primary     = name | literal | call | "(", expression, ")"
            | quantifier, "{", expression, "}" ;
quantifier  = ("any" | "all"), "placement" | "any", "record" ;
```

Newlines and trailing commas are allowed at natural document/field boundaries;
`#` starts a comment. Lexical rules for genomic literals are separate from numeric
expressions, avoiding the common `chr1:100-200` subtraction ambiguity. Typed
parameters bind values, not source-code fragments, so substitution cannot inject
operators or change the evidence unit.

Examples of required diagnostics:

- `terminal(...)` inside `any placement`: "terminal requires Record; found Placement;
  move it outside this quantifier to assert record-level co-occurrence".
- A minus-strand path written in the wrong order: show both offending junctions,
  the declared transcriptional order, and the normalized genomic path.
- `@cohort.classes` without a viable bounded closure plan: show the full-scan or
  index requirement, not a misleading initial-route-only estimate.
- `where RA & RO` at record level when RA/RO are placement predicates: offer the
  two **different** explicit quantifier forms; do not automatically choose one.
