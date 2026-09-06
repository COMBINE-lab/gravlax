# GQ v1 implementation plan for Gravlax 0.2.0

Status: implemented and validated for an **unpublished 0.2.0** release on
`astra-improvements`. See `gq-v020-implementation-results.md` for measured outcomes,
the implemented boundary, and remaining limitations. This is not a release announcement.

## Review disposition

Adopt the final review's three evidence domains, classical universal quantification,
explicit strand frames and path order, typed metadata, multiplicity bounds, typed
annotation features, namespaced federation units, bounded closure, and typed results.
Keep the distinction between a record's marginal conjunction and a conjunction on
one observation. No header option changes unique evidence into multimapper primaries.

Correct the proposed bracketing rule before implementing it. In `rows.rs`, compact
representatives minimize and maximize `(start, Reverse(end))`, not independent start
and end coordinates. For reads `[100,120)`, `[105,200)`, `[110,130)`, the representatives
are the first and last; the omitted end is outside both retained ends. End and
overlap decisions cannot use the proposed enclosing/contained geometry proof.
Start intervals and invariant junctions can prove answers. Otherwise use conservative
three-valued abstract interpretation, combined with actual retained witnesses.
Never infer arbitrary Boolean predicate monotonicity from atom monotonicity.
Read-support bounds also include retained witnesses, not only all-or-nothing chains.
Singleton weighted geometry is exact; two different representatives do not reveal
their individual multiplicities when the chain weight exceeds two.

Describe records as serialized evidence groups rather than universally one locus:
multimapping signatures can have different direct anchors within the same record.
Catalogue count units must be verified individually; they are not automatically
upper bounds on distinct class support. Exact support is a separate reduction.

## Parser decision

Use **winnow 1.0.4**, the version already locked by this workspace, with a spanned
token layer and bounded precedence parser. Its package declares Rust 1.65, below
our Rust 1.89 requirement. Keep parsing independent of resources and execution.

| Candidate | Fit and tradeoff |
| --- | --- |
| winnow | Flexible combinators for lexical boundaries and committed failures; explicit control over recursion/work limits. Existing dependency. Selected. |
| Chumsky | Particularly attractive for rich recovery, multiple syntax diagnostics, and an eventual editor. Recovery is not a scientific execution policy: recovered documents must never execute. Strong alternative, not rejected on performance grounds. |
| pest | Readable separate PEG grammar, but generated parse trees still need typed lowering and scientific diagnostics. |
| LALRPOP | Established LR generator; additional grammar/build layer not needed for this bounded expression/pipeline grammar. |

This is an architectural evaluation, **not a measured comparative speed claim**.
Benchmark the selected parser separately from archive execution; dependency choice
must not be conflated with evidence-query performance.

Primary references: [winnow 1.0.4](https://docs.rs/winnow/1.0.4/winnow/),
[Chumsky](https://docs.rs/chumsky/latest/chumsky/),
[pest book](https://pest.rs/book/), [LALRPOP book](https://lalrpop.github.io/lalrpop/).

## Ordered increments and acceptance

1. Spanned parser, bounded AST, type checker, clean specification, golden examples.
   Reject obsolete `placement`, invalid unit scope, ambiguous path orientation,
   recursion, and malformed/oversized documents before reading evidence payloads.
2. Native local executor: unique/stored/signature/alternative domains, geometry,
   classical `all`, conservative omitted-geometry interpretation, read bounds,
   terminal attachment, typed metadata and grouped reductions. Differential tests
   must map each expression to its actual existing CLI semantics.
3. Class/cell closure, archive federation, pinned bindings and typed result envelope.
   Fail closed on missing capabilities, assembly mismatch, duplicate sources and
   budget exhaustion. Full-scan permission is an execution flag, not query text.
4. Functions, annotation features, expression ergonomics, Python support, docs and
   workflow benchmarks; junction enumeration and exact support are the final increment.
5. Full workspace tests and strict lint, release-version synchronization and 0.2.0
   notes only when acceptance criteria are satisfied. Publishing is a separate action.

Do not enable unvalidated optimizer rewrites. In particular, `P | !P` need not be
true, and child cropping, class-posting intersections, and elimination of false
tally combinations are not generally semantics-preserving optimizations.
