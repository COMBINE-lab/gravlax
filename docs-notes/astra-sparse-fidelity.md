# Compact base plus sparse fidelity corrections

## Outcome: exact, but larger on the CD45 fixture

This experiment preserves the compact archive byte-for-byte and adds a rooted
correction capsule. It does not beat the existing full-fidelity representation,
even after entropy-coding the correction columns. It is not enabled for ingest
or production queries.

These are **actual file sizes**, including the separate capsule's manifest,
framing, directory, roots, references and coded multiplicities—not projected
payload-only or entropy-bound sizes.

| Indexed, 4,096-record target | Total bytes | Additional bytes over compact | Premium |
|---|---:|---:|---:|
| Compact archive alone | 384,486 | — | — |
| Existing full-fidelity archive | 466,431 | 81,945 | 21.31% |
| Base + sparse corrections, absolute shape IDs | 484,434 | 99,948 | 26.00% |
| Base + sparse corrections, shape-ID deltas | 481,713 | 97,227 | 25.29% |
| Base + sparse corrections, deltas + adaptive local rANS | 480,042 | 95,556 | **24.85%** |

The best sparse candidate is 13,611 bytes larger than the existing full-fidelity
archive. It also fails to improve on the earlier 16.20% matched-encoding premium
of the structural position-delta prototype (which remains a projection, not a
production codec).

Default, single-chunk encoding confirms that this is not just fine-chunk overhead:

| Default chunking, no access index | Total bytes | Fidelity premium |
|---|---:|---:|
| Compact | 340,008 | — |
| Existing full fidelity | 419,312 | 23.32% |
| Best compact + sparse capsule | 425,891 | 25.26% |

All candidates use zstd level 19 and the same accepted source evidence. This is
one narrow-locus dataset under two layouts, not two independent datasets.

## Representation tested

The compact base supplies molecule/class/cell identities, multimapper records,
shape dictionaries, junction-chain representatives and total chain read counts.
The capsule binds both source roots in its manifest and records its column codec.

For each compact chain:

1. Deduplicate identical retained representatives.
2. Store only geometries missing from that retained set. Positions are deltas;
   shape IDs are absolute or differences from the preceding shape ID.
3. Default each nonfinal geometry's count to one and infer the final count from
   the existing total chain count.
4. Store explicit count corrections wherever that default differs from the full
   evidence. The encoder checks total-count conservation and the decoder rejects
   an impossible or zero remainder.

Only chains needing extra geometry or explicit count corrections get a capsule
record. Their addresses are delta-coded flattened chain ordinals within each
original chunk. Seven columns separate routes, extra-geometry counts, count-patch
counts, positions, shapes, count-patch indexes and replacement counts. The entropy
variant chooses varints or local rANS by comparing each candidate column's final
zstd size, including its table; the completed section is then actually compressed.
This selection is a tested heuristic, not an exhaustive globally optimal coder.

Reconstruction emits the exact full-fidelity singleton chains in their original
canonical order. Molecule boundaries, classes, cell identity and multimapper
records stay unchanged. Neither read sequences nor qualities are added.

## Why “sparse” was not cheap here

Corrections touch **23,304 chains** and carry **52,974 omitted geometries**, across
all 14 fine chunks. That is approximately 32% of the 72,600 compact chains. The
omitted geometry is substantial; its locations and multiplicities still have to
be encoded. The immutable-base approach also cannot recover any bytes already
spent on the compact representation. These observations explain why an exception
layer is not automatically a smaller representation; they are not a proof that
every possible sparse codec must lose.

## Routing and integrity

When both inputs have access indexes, the experiment authenticates and semantically
verifies their postings against decoded records. It constructs the exact additional
postings required to turn base routes into full-fidelity routes. Nonempty additions
are stored as `index.access.supplement` in the capsule, and their space is charged.
They must be unioned with the base routes; a supplement is not a standalone index.

This fixture needs **no additional postings**: its base and fidelity chunk-route
sets are identical. Thus the measured size penalty is not an unpaid or enlarged
geometry index. Without an access index the manifest explicitly requires full
scans for fidelity geometry queries. A production query adapter is not implemented.

Capsules are experimental rooted containers, not independently usable `.aie`
archives. Their generic container integrity is checked, but `aie doctor` and
ordinary queries do not understand this experiment's evidence schema.

## Validation

- Both original archives are fully payload-verified.
- Paired dictionaries, molecule identities, chunk identities, multimapper records,
  retained geometry inclusion, and per-chain total multiplicities are checked.
- Every reconstructed `MolRec` is compared for exact equality with the full-fidelity
  input, including geometry ordering and weights.
- The written capsules are reopened, their roots/payloads verified, and every
  chunk reconstructed and compared again from the stored bytes.
- Tests cover missing middle geometry, count corrections, inferred counts, duplicate
  extremes, sparse record routes, same-class separate loci, truncated streams,
  entropy envelopes, and exact supplementary-index unions.
- **406 tests pass** across 23 workspace test targets, all features/targets enabled.
  Strict clippy with `-D warnings` passes.

This demonstrates reconstruction fidelity, not independent biological accuracy.
The experiment does not establish production query latency or native decoding cost.

## Reproduce and locate artifacts

```sh
cargo build --release --locked
target/release/aie dev sparse-bench \
  /tmp/gravlax-astra-benchmark-v3/indexed.aie \
  /tmp/gravlax-astra-benchmark-v3/fidelity.aie \
  --out /tmp/gravlax-sparse-new-run --level 19
```

The output directory must not exist. The tool creates three experimental capsules
and `report.json`; it never rewrites either source archive. The measured files are
under `/tmp/gravlax-sparse-fidelity-v2` and `/tmp/gravlax-sparse-default-v1`. These
are temporary local artifacts, not a published bundle. Checked-in reports are
`astra-sparse-indexed-results.json` and `astra-sparse-default-results.json`.

The default-layout full-fidelity input was generated from the same complete
multiplicity BAM/whitelist with `ingest-archive --geometry-fidelity`, omitting
the access-index and fine-chunk flags. It is at
`/tmp/gravlax-fidelity-default-for-sparse.aie`.

Measured binary SHA256:
`e62132dc1fbbc020d70aa1a0f5003076674ebffb645cdf9139e399ad3e589f99`.
