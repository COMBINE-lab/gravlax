---
title: aie stamp-genome
description: Stamp the reference genome's BLAKE3 signature into an index.
---

Stamps (or re-stamps) a reference-genome signature into an index's metadata so
sequence-consulting analyses — the internal-priming filter in
[`query apa`](/gravlax/cli/query/), [`aie extend`](/gravlax/cli/extend/) — can verify
that they are using the same bound reference. A FASTA is refused if it lacks any
chromosome named by the archive. If a different signature is already present,
the command reports that it is replacing it.

The signature is per-contig BLAKE3 over the uppercased bases, so it is invariant
to line wrapping, case, and gzip framing. It adds ~8 KB. All evidence streams are
copied compressed, byte-for-byte — stamping cannot perturb the data, and replayed
matrices remain byte-identical. For a logical
`gravlax.molecular-evidence.v2` archive, `meta.genome_reference_binding` records
the exact FASTA identity, normalized signature, the `stamp-genome` action, and a
caller-declared relationship. The original `alignment.provenance` section is
never rewritten: a later stamp cannot prove which reference produced the BAM.
Older archives update only `meta.genome_sig` and retain legacy/unattributed
binding semantics.

Indexes built with `ingest-archive --genome` are stamped from the start;
`stamp-genome` retrofits existing indexes.

## Usage

```sh
aie stamp-genome sample.aie --genome GRCh38.primary_assembly.genome.fa.gz
# or write to a new file instead of replacing in place:
aie stamp-genome sample.aie --genome genome.fa --out stamped.aie
```

| Option | Description |
|---|---|
| `--genome <GENOME>` | Reference FASTA (plain or gzipped) to bind for sequence-consulting queries; the command verifies its bytes and archive-contig coverage, while its relationship to the original alignment is caller-declared |
| `--out <OUT>` | Write here instead of replacing the input in place |
| `--report-format <FORMAT>` | Opt in to a versioned `text`, `tsv`, or `json` operation report |
| `--report-output <PATH>` | Atomically publish the report without replacement; requires `--report-format` |

With uniform reporting enabled, `gravlax.archive.stamp-genome-report.v1`
records the exact source and output archive identities, the raw FASTA content
identity, the normalized per-contig genome signature, and per-section byte
accounting. It also distinguishes a completed rewrite from the no-op case in
which the same genome signature was already present. In that no-op case,
supplying `--out` still publishes an exact byte-for-byte copy at the requested
new path; without `--out`, the source remains untouched. The report is separate
from the archive and is preflighted before genome hashing or archive rewriting;
omitting the report flags preserves the legacy stdout and archive behavior.
