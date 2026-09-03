---
title: aie export-molecule-bam
description: Export the exact post-correction molecule abstraction as a sequence-free BAM.
---

`export-molecule-bam` writes the archived molecule abstraction to a conventional
BAM container without inventing nucleotide UMI strings that the archive no
longer stores:

```sh
aie export-molecule-bam sample.aie \
  --fai GRCh38.fa.fai \
  --out sample.molecules.bam
```

The FASTA index supplies standards-compliant reference lengths. Placements use
ordinary BAM alignment fields; Gravlax-local tags retain the cell, global UMI
class, representative group, alternative, weight, and anchor state. Separate
unmapped records preserve one-mismatch UMI-class edges. Generic BAM tools do
not interpret this molecule model automatically, so this is an interchange and
capability-matched storage form rather than an ordinary read-level alignment.

| Argument or option | Description |
|---|---|
| `<ARCHIVE>` | Input `.aie` archive |
| `--fai <PATH>` | FASTA index supplying reference names and lengths for BAM `@SQ` records |
| `--out <PATH>` | New sequence-free molecule BAM; an existing file is not replaced |
| `--report-format <FORMAT>` | Emit a `text`, `tsv`, or `json` operation report in addition to the BAM |
| `--report-output <PATH>` | Write that report to a new file instead of standard output; requires `--report-format` |

Mapped records use these local tags:

| Tag | Meaning |
|---|---|
| `CB` / `XC` | Corrected barcode and dense cell ID |
| `XI` | Opaque global UMI-class ID; no nucleotide UMI is invented |
| `XM` | Dense molecule ID |
| `XW` | Signature read weight |
| `XK` | Record kind: `C` for a chain representative, `M` for a multimapper alternative |
| `XG` / `XA` | Group ID and representative/alternative index within that group |
| `XP` | Multimapper anchor flag |
| `NH` | One for a chain, or the number of multimapper alternatives |

After mapped placements, unmapped `XK:E` records store each one-mismatch UMI
edge using its smaller `XI` endpoint and larger `XJ` endpoint. Import these
records with `replay-rows --from-molecule-bam`; ordinary read-level tools do not
reconstruct the molecule graph.

## Uniform export report

The historical status line remains the default. Add
`--report-format text|tsv|json` for a versioned
`gravlax.molecule-bam.export.result.v1` report:

```sh
aie export-molecule-bam sample.aie \
  --fai GRCh38.fa.fai \
  --out sample.molecules.bam \
  --report-format json \
  --report-output sample.molecules.export.json
```

The typed summary records molecule, UMI-edge, BAM-record, archive-byte, and
output-byte counts. Its one-row `artifacts` table identifies the BAM and its
exact size and BLAKE3 identity. The report binds the archive identity obtained
from the same open reader used for export—the authenticated directory root for
v2 or a complete full-file digest for legacy v1—and binds the exact FASTA-index
bytes parsed to build the BAM header. Paths remain invocation locators rather
than substitutes for those content identities.

`--report-output` is optional. Without it, the selected representation is the
only content written to standard output and operational diagnostics go to
standard error. A report file is staged beside its destination and atomically
installed without replacing an existing path. The BAM remains the primary
artifact; report publication does not wrap the BAM itself in a cross-file
transaction.
