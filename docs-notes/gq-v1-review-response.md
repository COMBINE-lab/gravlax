# Response to the GQ v1 design review

This document evaluates the [GQ v1 review](/nfshomes/nomad/.codex/attachments/5e9ae10f-5a90-4623-8dc6-a290c23a662d/gqv1review.md)
against the [proposed language RFC](/tmp/gravlax-v016/docs-notes/gravlax-query-language-v1-rfc.md)
and the current archive/query implementation.

**Status:** design assessment and proposed disposition, not implementation
authorization. The recommendations below have not yet been incorporated into the
RFC or code. All new language syntax is illustrative, not implemented syntax.

## Overall assessment

The review is largely sound. Most recommendations should be accepted, but **the
quantifier-domain proposal and read-weight proposal need refinement before
adoption**. The original position on non-vacuous `all` should also change.

The two-algebra architecture, explicit evidence-unit boundaries, population-not-crop
semantics, three-valued truth, record-typed terminal evidence, and provenance and
execution-budget requirements remain the foundation of the design.

## Evaluation of the eight recommendations

### 1. Split the domains currently called “placement”: accept, but revise the proposal

The review identifies a genuine specification gap. However, three corrections
matter:

- **Primary is not uniquely mapped.** A multimapping read still has a primary
  alignment. `placements = primary` would not establish the conservative domain
  the review intends.
- Representatives are not guaranteed to come from the same physical molecule.
  They belong to an archive evidence group; raw-UMI collisions remain possible.
- The existing CLI already defaults to unique-chain evidence. Including
  multimapper alternatives requires an explicit choice. This is a gap in the
  proposed language specification, not an unsafe CLI default. See the
  [current placement policies](/tmp/gravlax-v016/crates/aie/src/querycmd/cooccur.rs:73).

There is a deeper distinction the review’s naming does not resolve:

> Quantifying over the **stored representatives** is different from quantifying
> over the **accepted alignments those representatives summarize**.

If we literally ask whether any stored representative satisfies a fully evaluable
predicate, inspecting them all gives a definite answer. Missing interior geometry
causes uncertainty only when the intended domain includes the observations omitted
by the quotient.

The language should therefore distinguish three domains explicitly:

| Domain | Meaning |
| --- | --- |
| Unique-alignment evidence | Accepted uniquely mapped alignment observations, evaluated through the archive’s possibly incomplete representation |
| Stored representatives | The geometries physically retained; principally useful for diagnostics |
| Multimapper alternatives | Competing retained alignment hypotheses, grouped by multimapping signature |

Illustrative syntax:

```text
any unique { P }

any multimap {
  all alternative { P }
}
```

The nested form is important: it distinguishes “one multimapping group has only
alternatives satisfying P” from “every alternative of every group satisfies P.”

The alternative set should include the primary placement. “Other alternatives”
would be a separate operation, not the default interpretation.

**Decision:** add the evidence-model section and explicit domains; reject a global
header that silently widens “representatives” to include alternatives.

### 2. Define region strand: accept, with an additional distinction

Strand on a region must have an operational meaning, not be decorative annotation.

The specification must also distinguish **alignment strand from transcript
strand**. Those cannot be conflated when protocol conventions require
opposite-strand interpretation.

Specify:

- Whether a literal constrains alignment or transcript strand.
- That omitted strand imposes no strand constraint.
- That protocol-aware conversion is explicit and recorded in provenance.
- That annotation-derived strand is not silently treated as alignment strand.

The review is right about the danger; simply saying “placement strand must match”
does not fully solve the annotation/chemistry interaction.

### 3. Declare path orientation explicitly: accept, without requiring redundant declarations

The explicit form is useful:

```text
path(-) { J1 >> J2 }
```

An explicitly genomic-order form should also be supported:

```text
path(order: genomic) { J1 >> J2 }
```

Retain `J1 >> J2` when its orientation is unambiguously established by its types.
Requiring an additional declaration in that case adds ceremony without additional
safety.

Remove the original silent convention switch: transcriptional order when stranded,
genomic order otherwise. An unresolved orientation should require clarification.

The compiler must validate the supplied order, not silently reorder the pattern.

### 4. Expose read multiplicity: accept, but distinguish totals from matching support

This is an important omission. Counting representatives is not counting reads.

However, an unqualified `reads` field hides several choices. Prefer:

```text
reads.unique
reads.multimapping
reads.total
```

Unique counts sum chain multiplicities. Multimapping counts sum signature weights
**once per signature**, not once per alternative. These correspond to accepted
input observations, not all sequenced reads. The
[ingest implementation](/tmp/gravlax-v016/crates/aie/src/rows.rs:3831) retains these
multiplicities.

The crucial limitation is:

> Total read multiplicity does not determine how many reads satisfy a geometry
> predicate.

Suppose a compact chain represents ten reads. One retained representative overlaps
a window; the other does not. The total is ten, but the matching count may not be
recoverable.

Similarly:

```text
where any unique { P }
|> summarize { reads = sum(reads.unique) }
```

counts **all unique-read weight in the selected units**, not just reads satisfying
P—and, under population-not-crop semantics, potentially includes reads outside the
universe window.

Predicate-specific read support must be provably exact, return bounds, or be
rejected as unavailable. It must never assign an entire chain’s weight to a
matching representative.

### 5. Return typed annotation features: accept

Make this more specific than generic `Feature<Assembly>`:

```text
let ptprc = gene(@gencode, "PTPRC")  # GeneFeature<Assembly>

ptprc.span
ptprc.exons
ptprc.junctions
```

Their types should distinguish a region, a region set, and a junction-pattern set.

Two additional rules are necessary:

- Overlapping annotation intervals must not cause duplicate base counting.
- The union of a gene’s junctions is **not** a set of valid annotated transcript
  paths. Transcript membership/order must remain available separately when that
  is the question.

This keeps annotation resolution explicit without inventing transcript
compatibility from a union.

### 6. Distinguish single-block and total overlap: accept

Accept the proposed option, with its mathematical meaning explicit:

- Single-block: maximum overlap between one aligned block and the target region
  union.
- Total: overlap between the union of aligned blocks and the target region union.

Using unions prevents duplicate counting when annotation exons overlap across
transcripts.

The parameter adds little syntactic complexity, but its implementation and
completeness analysis are not literally free. Compact archives may leave either
result unresolved.

### 7. Justify non-vacuous `all`: change the original design

The review’s objection is persuasive. For a composable language, favor **classical
`all`**, true on a known-empty domain, preserving:

```text
all D { P } == !any D { !P }
```

When nonempty evidence is required, state it:

```text
nonempty(D) & all D { P }
```

Domains cannot universally be assumed nonempty: a record may lack multimapper
alternatives, and an explicitly selected cell population may contain cells without
relevant records.

This is a proposed **DSL** decision. The existing CLI’s
`--match-within all-placements` remains non-vacuous; the compiler must not silently
equate the two on empty domains.

### 8. Reorganize as a specification: accept

The proposed structure is better. Operational notes belong in a work log. Keep
short motivation-level acknowledgments of seqproc and pipeline-language influences,
plus normal references.

One incidental claim in the review is too strong: these semantics are not
inherently impossible to represent relationally or in SQL. The justification for
a native engine is explicit scientific semantics and efficient execution over this
representation—not logical inexpressibility elsewhere.

## Answers to the five open questions

### 1. Enumeration

Include **bounded catalogue-based junction enumeration in v1**, rather than defer
all discovery. It corresponds to an existing capability and makes the language
immediately useful beyond confirmation.

General event-pattern generation and combinatorial cassette discovery can follow
in v1.1.

Importantly, enumerated patterns cannot simply feed ordinary `tally`: that would
count pattern rows, not supporting biological evidence units. An explicit
support-evaluation operation is needed, conceptually:

```text
from junctions(@cohort, within: G)
|> support(unit: class, by: {sample, cell.type})
```

The enumerator produces typed candidates; `support` evaluates them against the
archives with explicit unit semantics and budgets. Catalogue statistics must not
masquerade as exact class support.

### 2. Result schema

Define `tally` as:

> Grouped counts of the joint truth states of named expressions.

For `tally {ra, ro} by {sample, cell.type}`, illustrative rows would be:

| sample | cell.type | ra_state | ro_state | evidence_units |
| --- | --- | --- | --- | ---: |
| S1 | T cell | true | false | 120 |
| S1 | T cell | true | unknown | 8 |
| S1 | T cell | false | true | 42 |

Recommended contract:

- Explicit `"true"`, `"false"`, `"unknown"` state values.
- `UInt64` counts.
- A row key containing the grouping fields and all state fields.
- Observed combinations only; no automatic Cartesian expansion of zero-count
  combinations.
- Separate group/population summaries identifying denominators and source scope.
- Evidence-unit identity recorded in the result metadata.

Use the existing result envelope and named-table machinery with a new versioned GQ
table schema. The existing
[table contract](/tmp/gravlax-v016/crates/gravlax-output/src/table.rs:10) supports
string-valued states and integer counts.

This requires no second wire-format reader. Schema-aware clients may need to
recognize a new table identifier. Avoid using bare JSON null for evidence unknown,
so it remains distinguishable from missing metadata.

### 3. Explicit allowance for class closure

**Yes, whenever execution requires unrestricted full-scan fallback—not only for
`within all`.**

Make this an execution policy, supplied through command options or a bound
execution manifest, rather than changing the query’s mathematical meaning:

- Indexed closure proceeds within declared budgets.
- Full-scan fallback requires explicit allowance.
- Data-dependent expansion beyond a budget fails.
- The engine never substitutes incomplete closure or approximate answers.

Explain must disclose initial routing separately from subsequent closure.
Permission to scan does not override resource limits.

### 4. Metadata expressions

This belongs in v1’s actual grammar and type specification. Include:

- String, numeric, Boolean, null, and homogeneous list literals.
- `==`, `!=`, `<`, `<=`, `>`, `>=`.
- `in` and `not in`.
- Boolean `!`, `&`, `|`.
- Parentheses and typed numeric arithmetic.
- `is_null` and explicit handling of nullable metadata.

For example:

```text
where cell.type in ["T cell", "B cell"] & qc.umis >= 500
```

Require compatible operand types; do not silently convert `"500"` into a count.

A missing column is a schema error unless explicitly declared optional. A missing
value in a nullable column is different: comparisons can yield unknown, which
filtering must account for. Grouping must retain a missing-value group rather than
silently dropping it.

Leave arbitrary regular expressions out of v1.

### 5. Exported function signatures

**Require explicit parameter and return types on exported functions.** Allow
inference for file-local definitions.

Exported signatures should expose assembly/strand requirements and whether they
produce a placement-, record-, or higher-unit predicate. Those constraints are
part of the reusable scientific contract, not merely implementation details.

## Smaller edits and acceptance tests

Accept:

- An expected assembly declaration in the header, checked against bound resources.
  A matching label alone is not proof of reference compatibility.
- `all record { P }` in v1, using the same classical universal semantics.
- Explicitly defining `near` as widening the matching set—not merging junction
  identities.
- Diagnostics for missing strand when a library function requires it.

**Keep `tally`**, define it immediately, and restrict its expression fields to
Truth in v1. Categorical metadata already fits naturally in `by {…}`.

The proposed tests need one correction:

```text
any unique { A } & any unique { B }
```

corresponds to ordinary record-level co-occurrence, whereas:

```text
any unique { A & B }
```

corresponds to the placement-local existential mode. Requiring both to agree with
default `cooccur` would enforce the very semantic conflation the language is
intended to prevent. Comparisons must also use corresponding predicate and
placement policies.

Accept the optimizer property tests, with two qualifications:

- Compare scientific results under stable logical-unit mappings, not physical
  chunk IDs or provenance bytes.
- Include unknown-valued counterexamples. For example, rewriting `P | !P` to
  `true` is invalid under the proposed three-valued logic.

## Bottom line

Accept the review’s direction, but make the next revision more precise about
quantifier domains, strand frames, and what multiplicity can actually measure.
Those are substantive correctness issues; the remaining recommendations mostly
improve clarity and completeness.
