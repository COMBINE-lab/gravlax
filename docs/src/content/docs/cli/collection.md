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

`--min-support` uses catalogue support only as a safe pruning bound. An upper
bound never substitutes for a source-derived count. `--explain` adds each
sample's routing decision. `--verify-content` performs the stronger source
check before execution. `--top 0` means all cell rows under the uniform
contract.

### Query options

| Command | Argument or option | Default | Description |
|---|---|---|---|
| `junction` | `<COLLECTION> <LOCUS>` | required | Collection and exact `chrom:donor-acceptor` junction |
| `region` | `<COLLECTION> <LOCUS>` | required | Collection and 0-based, half-open `chrom:start-end` window |
| `jset` | `<COLLECTION>` | required | Collection to query |
| `jset` | `--include <JUNCTION>...` | required | One or more inclusion-side exact junctions |
| `jset` | `--exclude <JUNCTION>...` | required | One or more exclusion-side exact junctions |
| `junction`, `jset` | `--min-support <N>` | `0` | Skip decoding when the collection-wide support upper bound is below this value |
| all queries | `--top <N>` | `5` | Per-sample cell rows to return; under `--format`, `0` returns all |
| all queries | `--explain` | off | Include per-sample routing decisions and compressed-byte estimates |
| all queries | `--verify-content` | off | Re-hash every source archive before querying |
| all queries | `--json` | off | Emit the command-specific JSON representation |

## Uniform output

All five collection subcommands accept opt-in `--format text|tsv|json` and
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
