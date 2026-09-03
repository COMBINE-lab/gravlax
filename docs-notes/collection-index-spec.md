# Federated collection index

The `.aicollection` v4 format is a deterministic, authenticated routing layer over independent
Gravlax archives. It contains archive identities, genomic interval routes, junction routes, and
optionally local shape routes. It does **not** contain molecule records, annotations, group labels,
or query results. The source `.aie` files remain the authoritative molecular evidence: deleting
and rebuilding a collection cannot change an exact answer.

## Building a collection

```sh
# Content-addressed root layer.
aie collection build \
  --sample donor-a=/data/donor-a.aie \
  --sample donor-b=/data/donor-b.aie \
  --shape-routes \
  --out atlas.aicollection \
  --format json --output atlas.build.json

# Add archives without rewriting the root layer.
aie collection build \
  --base atlas.aicollection \
  --sample donor-c=/data/donor-c.aie \
  --shape-routes \
  --out atlas-plus-c.aicollection
```

Sample IDs are sorted before encoding, so the order of `--sample` arguments does not affect a root
layer. Each input must have the same chromosome dictionary and, by default, the same stamped
reference-genome digest. `--allow-unstamped` permits chromosome-only compatibility checks when a
genome digest is unavailable. Duplicate IDs, paths, inodes, and encoded archive contents are
rejected so that the same evidence cannot be counted twice.

The builder reads `meta`, `chroms`, `rans.tables`, `index.chunks`, `index.junctions`, and
`index.jpost`. It does not read any `cN` molecule chunk. With `--shape-routes`, it also reads the
`shapes` dictionary and derives the exact routes described below; this option requires rooted `.aie`
v2 sources. A build without `--shape-routes` remains fully functional and uses the source shape
dictionary during point and junction-set execution.

## Source identity

A rooted `.aie` v2 source is identified by its authenticated directory root. The directory already
commits the exact compressed digest of every source section, so collection construction and normal
queries can establish content identity without scanning molecule payloads. The manifest also stores
the scheme-independent encoded-section identity, permitting duplicate detection across a legacy v1
container and its byte-preserving v2 seal.

Seekable `.aie` v1 sources remain supported for collections without shape routes. Because v1 has no
authenticated directory, deriving its native and encoded-section identities requires a complete
source scan. `aie seal-archive` can copy a v1 archive into the rooted v2 container without changing
its compressed evidence payloads.

For each archive, the manifest records its canonical path, byte length, nanosecond modification and
change times, device, inode, archive format, native content identity, encoded-section identity,
chromosome identity, and optional genome identity. Queries open the recorded source and check this
filesystem identity. For v2 they also authenticate and compare the directory root using the same
open file description that supplies subsequent section reads. A normal v1 query uses the filesystem
guard; `--verify-content` performs its complete digest scan. For v2, that option verifies every
source payload. These checks fail closed if a source is missing, replaced, truncated, or
incompatible.

## Immutable layers

A root collection contains one immutable layer. An incremental build writes only the new archives
and their routes, plus a canonical path and authenticated root digest for its parent. It does not
decode, copy, or recompress the parent's junction or shape-route sections.

Queries follow at most 32 parents, verify every parent digest, reject cycles or changed parents,
and form one sample order by ID. Routes use layer-local archive ordinals and are remapped while the
chain is opened; inserting a sample whose ID sorts before an older sample therefore does not
invalidate an existing layer. A compact root and an equivalent incremental chain produce the same
sample rows and totals, although their bytes need not be identical.

## Container and sections

Every v4 layer begins with `GRVLXCOL`, a little-endian version and section count, and a 32-byte
BLAKE3 root. The root covers the fixed header and complete variable-length directory. Each
directory entry contains a section name, raw and compressed lengths, and a digest of the raw
payload. Payloads are independent zstd frames laid out in directory order. Opening a layer verifies
the complete directory before using its offsets; reading a selected section performs bounded
decompression and checks its authenticated payload digest.

Sections are canonical and ordered:

- `manifest` records the parent binding, reference and chromosome identities, archive identities,
  interval chunk metadata, junction cardinalities, and optional shape-route bindings;
- `j.<chromosome-id>.<donor-bin>` stores junction rows for a 16-Mb donor-coordinate bin;
- `s.<layer-local-archive>.<first-span>` stores one optional local shape-route block.

A junction row is sorted by `(chromosome, donor, acceptor)` and contains a layer-local sample
presence bitmap, a sum of supporting-child counts, and delta-coded source chunk postings for each
present archive. The support sum is an upper bound used only for pruning. Cell-, group-, and
UMI-resolved counts are always recomputed from source chunks.

## Local shape routes

Junction geometry in a source archive is represented by a shape identifier and a placement
coordinate. Loading the complete shape dictionary for every donor is avoidable because the exact
intron length is known from a point-junction or junction-set request.

For each archive, `--shape-routes` derives a sorted table keyed by exact intron span. A row contains
the unique candidate pairs

```text
(shape_id, donor_offset)
```

where `donor_offset` is the end of one shape block and
`acceptor_offset = donor_offset + intron_span` is the start of the next. The table preserves every
adjacent-block intron tuple needed for an exact junction predicate; it is not a replacement for the
full shape dictionary and is not filtered by observed counts. Rows are partitioned into blocks of
at most 256 distinct spans and delta-coded in canonical order.

Each archive's route binding commits:

- the layer-local archive ordinal;
- the source `.aie` v2 directory root;
- the source directory's committed digest for the compressed `shapes` payload;
- the source shape count; and
- the exact names and span bounds of all route blocks.

For a requested junction, the collection first selects archives and source chunks from the `j.*`
route. If an archive has a shape route, the query reads only the block containing the requested
span and tests both genomic boundaries against its candidate `(shape_id, donor_offset)` pairs.
Chain representatives and multimapper placements are handled by the same exact predicate. The
route supplies no molecule, UMI, or cell count: source chunk records, UMI classes, and cell-of-class
blocks remain authoritative. Collections without a route, and older layers within a mixed chain,
fall back to the complete source shape dictionary with identical query semantics.

Shape routes accelerate `collection junction` and `collection jset`. `collection region` is
unchanged because it routes directly through anchor intervals and source chunk metadata.

## Query surfaces and accounting

```sh
aie collection junction atlas.aicollection chr2:1234567-1250000 --explain --json
aie collection region atlas.aicollection chr2:1200000-1300000 --json
aie collection jset atlas.aicollection \
  --include chr2:1234567-1250000 \
  --exclude chr2:1234567-1260000 \
  --json
```

The default text and command-specific `--json` presentations remain byte-compatible when
`--format` is omitted. All five collection surfaces (`build`, `inspect`, `junction`, `region`, and
`jset`) also accept the opt-in uniform contract:

```sh
aie collection junction atlas.aicollection chr2:1234567-1250000 \
  --top 0 --format json --output junction.json
aie collection region atlas.aicollection chr2:1200000-1300000 --format tsv
aie collection jset atlas.aicollection \
  --include chr2:1234567-1250000 \
  --exclude chr2:1234567-1260000 \
  --format text
```

`--format text|tsv|json` emits the same versioned result envelope and named typed tables in every
encoding. JSON can carry the whole bundle directly; text and TSV use explicit table sections.
`--output PATH` requires `--format`, writes through a sibling staging file, and atomically installs
the completed result without replacing an existing path. An occupied or invalid destination is
rejected before collection or source scanning begins. Without `--output`, machine output is the
only content written to stdout.

Uniform collection results use these result/table families:

| Command | Result schema | Named tables |
| --- | --- | --- |
| `build` | `gravlax.collection.build.result.v1` | `archives`, `source_io`, `source_sections` |
| `inspect` | `gravlax.collection.inspect.result.v1` | `layers`, `chromosomes`, `archives`, `shape_route_blocks` |
| `junction` | `gravlax.collection.junction.result.v1` | `samples`, `cells` |
| `region` | `gravlax.collection.region.result.v1` | `samples`, `cells` |
| `jset` | `gravlax.collection.jset.result.v1` | `requests`, `samples`, `cells` |

Every table declares whether its rows are a set or sequence, its logical key, and any guaranteed
ordering. Sample and cell tables are logical sets: their current emission order is not part of the
contract and the writer does not sort merely to make bytes deterministic. Request, layer, and
chromosome tables are sequences with explicit integer order fields. A `top` limit applies per
sample; `top=0` means all cells when `--format` is used. The cell table's exact aggregate availability,
emitted-row count, and truncation state are recorded as selection metadata, while the per-sample
selection comparator is recorded in the typed summary.

Scientific query totals and category semantics live under `data.summary`; source identities,
collection-layer roots, invocation parameters, and access strategy live under `provenance`.
Archive identities come from the same opened manifests/readers used to plan and execute the query.
The serializers stream borrowed rows directly and never introduce result-wide sorting or a second
row materialization.

- `collection junction` returns exact per-sample point-junction UMI totals, supporting-cell counts,
  and the requested number of top-cell rows.
- `collection region` returns exact per-sample anchor-window molecule, UMI, and cell totals while
  retaining one decoded source chunk at a time.
- `collection jset` returns exact `include_only`, `exclude_only`, and `both` UMI-class categories.
  Requested junctions sharing a sidecar block or source chunk are decoded together.

`junction --min-support N` may skip all source work when the global upper bound is below `N`.
`jset --min-support N` requires every requested component to meet the bound. An absent sample is
reported explicitly as zero; an upper bound never substitutes for source-derived counts.

Versioned JSON output separates the principal I/O terms:

- `source_archive_identity_bytes_read`: source bytes read to establish or check identity;
- `source_archive_execution_bytes_read`: source bytes read to compute the answer;
- `source_archive_bytes_read`: the sum of source identity and execution reads;
- `collection_sidecar_bytes_read`: authenticated directories and selected sidecar payloads;
- `shape_route_sidecar_payload_bytes_read`: the shape-route subset of sidecar reads; and
- `total_logical_bytes_read`: source and sidecar bytes combined.

The plan also reports layers, opened and pruned archives, decoded chunks, route blocks, routed and
fallback archives, and timing for collection loading, identity checking, route planning, source
execution, and the complete command. Build JSON separately reports total source I/O and the
payload/source bytes attributable to deriving shape routes.

## Inspection and integrity

```sh
# Verify every sidecar payload and ordinary source identity guard, using command-specific JSON.
aie collection inspect atlas.aicollection

# Emit the uniform inspection bundle atomically.
aie collection inspect atlas.aicollection --format json --output atlas.inspect.json

# Also reconstruct every shape route from its bound source dictionary.
aie collection inspect atlas.aicollection --verify-routes

# Verify complete rooted-v2 source content as well as route reconstruction.
aie collection inspect atlas.aicollection --verify-content
```

Ordinary inspection authenticates every layer directory, checksum-verifies every sidecar section,
and validates source identities. `--verify-routes` additionally reads each bound source `shapes`
section and requires structural equality between the stored route blocks and a fresh deterministic
derivation. `--verify-content` includes that reconstruction and verifies every rooted-v2 source
payload; a legacy-v1 source instead receives a complete full-file and encoded-section digest scan.

Readers reject unknown versions, noncanonical section names or ordering, undeclared or missing
route sections, malformed varints, unsorted or duplicate coordinates and postings, invalid
presence bitmaps, out-of-range shape/chunk/archive identifiers, arithmetic overflow, truncation,
trailing bytes, unsafe allocation or compression ratios, checksum failures, and source/root/shape
binding mismatches.

## Scope

Collection files use local filesystem paths for sources and parents. Relocatable object-store
locators, chain compaction, and direct group-map aggregation are not part of the v4 format. Group
maps and annotations remain query-time inputs; the `cohort` commands provide group- and
sample-scoped biological reductions. Readers accept older v2/v3 collection layers, but a v2 layer
cannot be used as an incremental v4 base because it lacks the cross-container encoded identity.
