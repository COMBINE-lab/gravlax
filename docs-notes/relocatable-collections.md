# Relocatable collections

## Contract

The expected archive commitment is authoritative; a pathname is a locator. Rooted
archives copied or moved to new paths, filesystems or machines reuse the original
collection and shape routes. No archive or collection format version changes, no
index regeneration, and no rewritten commitments are involved.

`--locations` supplies a bounded, versioned JSON list of scheme-qualified identity
and local-path pairs. Relative paths are relative to the manifest, so a directory
bundle can move intact. Unmapped identities retain embedded hints. Overrides never
silently fall back after a missing file or mismatching commitment. The resolver also
handles parent collection layers. Expected identities come from authenticated
collection bytes, not from the location manifest. Extra entries may describe other
members of a shared bundle; duplicate keys and unknown schemes/fields fail closed.

Historical size, inode, device and timestamps remain readable in existing collection
v2/v3/v4 manifests for wire compatibility, but no longer determine source acceptance
or cross-layer duplicate identity. Current-source build duplicate checks remain;
content-based duplicate prevention remains across layers. Inspection, actual query
opens and shape-route reconstruction share one content guard. Same-operation checks
use fresh metadata snapshots and can represent pre-epoch timestamps, rather than
requiring the serialized build-time stat tuple to match.

Opening a rooted archive authenticates its existing header/directory and compares
its root and encoded-section identity. Consumed sections retain ordinary payload
verification. Unread payloads are not represented as fully audited; `--verify-content`
remains explicit. Legacy unrooted v1 sources require full-file hashing on each
preflight and actual source open; they cannot provide the zero-additional-scan path.

All collection CLI commands accept the option (build uses it for the existing base
chain). Python event-search methods expose `locations=...`. Collection region,
junction and junction-set plan steps accept a registered metadata resource named by
`locations`, which is bound in plan provenance. Uniform collection results separately
report committed identities, resolved source paths and the location-file digest.
No implicit directory search, remote resolver, persistent cache or recompressed
evidence-equivalence claim is introduced.

## Validation

The executable acceptance test constructs a three-archive, two-layer synthetic
collection, copies the whole bundle, makes the original directory unavailable and
then runs inspection, full-content verification, shape-route reconstruction, region,
junction, junction-set and event search through the location manifest. It checks
complete science/work outputs, archive and sidecar I/O counts, and all original
archive/collection file hashes. It also verifies identical replacement at an
unchanged path and explicit wrong-root rejection without falling back.

The same test passed between the actual `/tmp` and `/dev/shm` filesystems with the
release-mode executable. See [cross-filesystem evidence](relocatable-collection-cross-filesystem.json).
Unit tests additionally cover historical inode/path collisions across layers,
arbitrary saved stat fields, pre-epoch current timestamps, legacy full-file identity,
same-open-file pathname replacement, unchanged directory commitments with corrupted
consumed payloads, malformed/oversized location manifests and wrong parent bindings.

The existing portability matrix now invokes relocation acceptance alongside GQ on
five release targets and both allocator variants on Windows/musl. Extracted Linux
release archives inherit the same check. Local validation includes the complete
workspace tests, strict all-feature Clippy, Python tests, distribution tests, plan
resolution and a documentation build. Hosted platform results are reported by PR CI,
not inferred from the local Linux run.

## Paired real-data benchmark

The immutable eight-donor SEZ demo capsule was copied into two locations on the same
filesystem. Four donors form the base layer and four the extension; both carry shape
routes. Original-path and relocated-path queries use the same release executable.
One warmup and nine paired repetitions use randomized order and four Rayon threads.
The benchmark verifies the entire typed result and work summary on every invocation,
excluding timing fields and locator-specific provenance. This is selected-locus demo
data, not a whole-atlas benchmark. Timing is full-process wall time, including startup,
argument parsing and output publication, on a noisy shared Linux host.

| Query | Original → relocated median ms | Archive bytes read, both | Chunks decoded, both |
|---|---:|---:|---:|
| chr1 region | 73.6 → 73.8 | 420,135 | 8 |
| chr12 region | 62.2 → 61.7 | 372,313 | 8 |
| Observed junction | 56.3 → 51.8 | 372,313 | 8 |
| Observed junction set | 67.6 → 64.8 | 400,113 | 16 |

All paired timing confidence intervals include parity and are wide; these numbers
do not establish a speedup or a tight bound on latency overhead. The stronger result
is exact equality of archive/sidecar reads and decoding work. Relocation adds only
location-manifest parsing and lookup, not archive processing. Location-file metadata
I/O is included in elapsed time but not in the existing archive/sidecar I/O counters.

Raw timings, paired intervals, source checksums, binary identity, selected junctions
and verified work summaries: [benchmark report](relocatable-collection-benchmark.json).

Reproduce using:

```sh
python3 packaging/collection_relocation_smoke.py --binary target/release/aie
python3 packaging/collection_relocation_smoke.py --binary target/release/aie \
  --relocation-dir /dev/shm --report cross-filesystem.json
python3 scripts/benchmark-collection-locations.py --binary target/release/aie \
  --demo-dir /path/to/verified/demo-data-v1 --repeats 9 --output benchmark.json
```

Both report destinations must be new files. No release version is bumped by this
change; it is proposed separately against main on `relocatable-federation`.
