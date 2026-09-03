# Scientific-intent and result-output contract

This note describes the interfaces shipped in Gravlax 0.1.0.

The identifier resolver is used by direct commands, plans, Explorer, and the
Python client. The shared result envelope is available across supported direct CLI
query, federation, cohort, collection, archive-lifecycle, compilation, export,
and extension surfaces, and across the corresponding subset represented by
source-plan v1 step kinds. Existing command-specific `--json` and `--tsv` modes
remain unchanged; the shared contract is selected
explicitly with `--format`, or with `--report-format` on commands that expose a
separate operation report. This distinction is command-specific rather than a
rule for every artifact producer. Native Arrow IPC is not yet a general CLI
output, and no R client is provided. This note describes the contract and the
commands that use a different output format.

## Identifier resolution

`anno::intent::IntentResolver` turns biological identifiers into genomic loci without allowing the
reference context to remain implicit. Construction requires an `AnnotationIdentity` containing a
non-empty assembly and annotation release; an optional content digest can bind the label to a
specific artifact. Digests use the canonical `blake3:<64 lowercase hex characters>` form and
`IntentResolver::from_path` verifies the digest before exposing any resolutions;
`from_annotation` refuses an unverified digest on an in-memory model. Every `ResolvedFeature`
repeats this identity.

Accepted input is either unprefixed or explicitly typed:

```text
TP53
gene:ENSG00000141510
transcript:ENST00000269305
exon:ENSE00003753508
```

Resolution rules are deterministic:

1. Exact stable identifiers and exact, case-sensitive gene symbols are considered first.
2. An unversioned stable identifier can match the version carried by the selected annotation
   (`ENST00000269305` can match `ENST00000269305.9`). A supplied version is exact: `.8` never
   silently resolves to `.9`.
3. Multiple symbol matches, multiple stable-ID versions, or cross-kind matches are errors with
   structured candidates. Nothing is chosen by first occurrence.
4. Coordinates are always zero-based, half-open and carry contig and strand.
5. Exons reused by multiple transcripts resolve once with every known parent and locus.

The GTF compiler retains transcript IDs and each source `exon_id` with that source record's
exact interval. Assignment exons remain separately merged where they overlap or touch. Compiled
annotation format v2 serializes these dictionaries and reads them back losslessly. The reader
still accepts v1: gene IDs and symbols remain resolvable, while transcript/exon queries return a
structured `IdentifierMetadataUnavailable` error. Partial source metadata is also reported as
unavailable for unknown identifiers rather than incorrectly claiming they do not exist. It does
not invent stable IDs for legacy files. Recompile a source GTF to create a v2 artifact.

## Result envelope

`gravlax-output` defines the outer JSON contract `gravlax.result-envelope.v1`:

```json
{
  "$schema": "gravlax.result-envelope.v1",
  "result_schema": "gravlax.query.region.result.v1",
  "producer": {"name": "aie", "version": "0.1.0"},
  "provenance": {
    "archives": ["aie-directory-root-v2:<64-hex-root>"],
    "assembly": "GRCh38.p14",
    "annotation": "GENCODE 49",
    "annotation_digest": "blake3:...",
    "parameters": {
      "aggregation": "cell",
      "archive_access": "range-index-selected archive chunks",
      "archive_path": "sample.aie",
      "archive_version": 2,
      "cell_scope": {"source": "all", "aggregation": "cell", "selected_cells": 1},
      "selection_policy": {
        "requested_top": 20,
        "top_zero_means_all": true,
        "comparator": "umis descending, entity ascending (barcode)"
      }
    }
  },
  "warnings": [],
  "data": {
    "summary": {"molecules": 20, "umis": 17, "cells": 1},
    "tables": [{
      "name": "counts",
      "schema": {
        "id": "gravlax.query.region.counts.v1",
        "fields": [
          {"name": "aggregation", "data_type": "string", "nullable": false},
          {"name": "entity", "data_type": "string", "nullable": false},
          {"name": "umis", "data_type": "uint64", "nullable": false},
          {"name": "cells", "data_type": "uint64", "nullable": true},
          {"name": "selected_cells", "data_type": "uint64", "nullable": true}
        ],
        "semantics": {
          "row_semantics": "set",
          "key": ["aggregation", "entity"]
        }
      },
      "selection": {"available_rows": 1, "emitted_rows": 1, "truncated": false},
      "rows": [["cell", "AAAC...-1", 17, null, null]]
    }]
  }
}
```

There is deliberately no automatic timestamp: identical inputs and parameters can produce
byte-identical results when the command also promises canonical presentation order; the envelope
does not by itself promise byte-identical reruns. Command-specific schemas must be versioned. Field types are `string`,
`int64`, `uint64`, `float64`, `boolean`, and `json`; rows are checked for width, nullability, finite
floats, and exact logical types before emission.

### Command coverage

| Surface | Output |
| --- | --- |
| `aie resolve` | Single typed table in `gravlax.result-envelope.v1`; text/TSV/JSON and no-clobber `--output`. |
| `aie compare-annotations` | One typed JSON bundle with four tables, or one selected typed TSV table; atomic no-clobber file output. |
| `aie query … transcript-ecs` | One typed JSON bundle with catalog/counts and optional membership, or one selected typed TSV table; hard non-truncating row caps and atomic no-clobber file output. |
| Other `aie query` surfaces | `region`, `junction`, `batch`, `junctions`, `jset`, `events`, `splice-graph`, `apa`, `apa-test`, and `discover` provide opt-in streaming text/TSV/JSON bundles with typed summaries and tables. |
| `aie federate` | Opt-in archive/count bundle with source identities, routing information, and exact selection metadata. |
| `aie cohort` | `events`, `splice-graph`, `transcript-ends`, and `polyasite-mixture` provide typed reports; directory-producing endpoint analyses keep `summary.json` as their primary completion marker. |
| `aie collection build`, `inspect`, `region`, `junction`, `jset` | Opt-in streaming text/TSV/JSON bundles with typed summaries, explicit table semantics, exact cell-selection metadata, content/root provenance, and atomic no-clobber `--output`; omitting `--format` preserves each command's default presentation. |
| Archive lifecycle | `ingest-archive`, `replay-rows`, `stamp-genome`, and `seal-archive` provide opt-in `--report-format`; `inspect-archive` provides `--format`. Reports bind archive/artifact identities without replacing the primary artifact. |
| Artifact producers | `compile-annotation`, `export-molecule-bam`, and `extend` provide typed opt-in operation reports. Their report file is published independently of the primary AIC, BAM, GTF, or auxiliary files. |
| Plans and Explorer | Source plans opt into uniform output/report blocks; resolved-plan v6 discloses schemas and publication behavior, and Explorer exports those same plans without executing them. |
| `aie doctor`, `aie ingest check`, projects, plans, and completions | Stable purpose-specific machine contracts where applicable, but not result envelopes. |
| `aie dev` | Command-specific evaluation output; these subcommands do not emit a shared scientific-result envelope. |
| Python | Strict project, resolved-plan v3-v6, named-bundle, envelope, and MEX models; explicit `result_raw()` for command-specific JSON and streaming `run_to_file()` for large output. |

### Performance measurements

Benchmarks alternated the typed `--format` or `--report-format` interface with
the corresponding command-specific output while checking scientific
equivalence, wall time, peak RSS, and the unchanged bytes of default stdout or
primary artifacts. Selected commands were also checked against the executable
from before typed output was introduced. The measurements below use repeated,
alternating runs on the same inputs; they describe those workloads and systems,
not every possible filesystem or future schema.

| Evaluated surface | Typed / command-specific median time | Typed / command-specific peak RSS |
| --- | ---: | ---: |
| Region and junction queries on the evaluated archive | 0.919–1.009 across scan-, scope-, JSON-, TSV-, and file-output cases | at most 1.010× (raw maximum 1.009335×) |
| Batch, junction enumeration and sets, events, splice graph, and federation on the evaluated archive | 0.937–1.003 | 0.901–1.013× |
| Collection junction with local shape routes | 0.978–0.996 | increase of less than 4 MiB |
| Archive inspection | 0.993 | 1.010× |
| Gene replay with an operation report | 1.029 | 1.086× |
| Archive extension with a 3.3 GB GTF | 0.982 | 0.999× |
| Eight-donor PolyASite mixture | 1.028 | 1.004× |
| Eight-donor transcript-end analysis | 1.020 | 0.991× |

Every comparison produced byte-identical scientific output. On these workloads,
typed streaming output did not impose a material performance penalty. In
particular, no command adds a global row sort merely for the shared contract.

Direct annotation-dependent commands record `--assembly` as a caller
assertion. A checked project plan can additionally verify it against registered
coordinate-resource metadata; matching labels or annotation bytes alone do not
prove assembly compatibility.

### Row semantics and ordering

A uniform result contract does **not** inherently require a globally
deterministic row order. A complete count, delta, transition, or catalog table
can represent a relation, set, or multiset while its physical rows stream in
producer or parallel traversal order. Writers must not add an implicit global
sort merely to make a table typed.

Each migrated table should therefore document separately:

- whether rows are a set, multiset, or sequence;
- the key or uniqueness fields, when one exists;
- whether presentation order is unspecified or explicitly defined; and
- any selection/ranking comparator used for `--top`, pagination, truncation,
  or bounded witnesses.

Determinism remains necessary for the *selected subset* of a ranked or capped
result, including tie-breaking; for content-addressing a mathematical set; for
bounded audit witnesses that promise reproducibility; and for tables whose row
ordinals are semantic. Matrix feature/barcode arrays must remain aligned with
matrix indices. None of these requirements implies sorting every complete
table.

Integrity and semantic identity are also distinct. A stored artifact digest
authenticates the exact bytes, including their physical row order, and is what
plan resume rechecks. An optional future logical-result digest for an unordered
table would require canonical row encodings with multiplicity-preserving
canonicalization; a commutative XOR-style shortcut would be insufficient.

## Format behavior

| Format | Contract |
| --- | --- |
| text | Human-readable schema/identity prelude and pipe-separated typed rows. |
| TSV | Comment metadata, one header row, then rows. Null is `\N`; backslash, tab, CR and LF use C-style escapes. |
| JSON | The envelope above, streamed row by row without materializing the result. |
| MEX | Non-overwriting feature-by-barcode directory with `matrix.mtx`, `features.tsv`, `barcodes.tsv`, and envelope `metadata.json` written last as its completion marker. Indices in Matrix Market are one-based. |
| Arrow | Dependency-neutral, bounded `ColumnarBatch` stream delivered to an `ArrowBatchSink`. Primitive types map directly; logical JSON maps to canonical UTF-8 with `gravlax.logical_type=json`. |

The core library does **not** yet ship an Arrow IPC encoder. This is an intentional dependency
boundary rather than a placeholder data model: batches validate the same schema and nullability as
the text writers, and a small adapter can map each column one-to-one to `arrow-rs`, PyArrow, Polars,
or an in-process consumer. A future CLI integration that promises `--format arrow` must provide
such a sink; the generic row writer returns an explanatory unsupported-format error instead of
silently emitting another representation.

MEX is separate from ordinary table writing because it needs dimensions, feature/barcode labels,
unique coordinates, and sparse nonzero semantics. Existing destinations are always rejected.

## Publication and failure semantics

`--output` on envelope-producing commands renders a complete sibling
temporary file and installs it atomically at a new destination. Validation,
digest, cap, or rendering failure leaves no advertised result, and an existing
destination is rejected. MEX is a new directory whose metadata envelope is
written last as its completion marker; readers reject a directory without that
marker or with inconsistent dimensions, labels, coordinates, or paths.

Diagnostics belong on standard error. Machine-readable stdout or `--output`
bytes contain only the selected result representation. Timings and other
environment-dependent diagnostics must not silently enter a deterministic
scientific result.

## Streaming and memory behavior

The interface is intended to add schema and provenance without changing the
scientific computation or forcing result materialization. Flat results should
use fallible streaming row producers, borrowed or directly encoded values, one
buffered writer, and bounded Arrow record batches. Multi-table bundles require
a bounded bundle writer rather than collecting every table in memory. Matrix
outputs require transactional streaming MEX support.

Output follows natural traversal order unless a command's selection semantics
requires ranking. A gratuitous global sort would change an
`O(n)` bounded-memory output path into `O(n log n)` time and usually `O(n)`
additional memory. Likewise, mechanically representing every row as an owned
`Vec<ScalarValue>` is a convenience for small tables, not the large-result
implementation target.

Cross-format scientific equivalence, peak memory independent of result
cardinality for streaming tables, and measured end-to-end and serialization
throughput characterize an implementation of this interface. Atomic
publication durability is a separate policy: a same-filesystem atomic install
is cheap, while mandatory `fsync` can be observable on local or network storage
and is not part of the format contract.
