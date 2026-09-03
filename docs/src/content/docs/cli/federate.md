---
title: aie federate
description: One junction query across N .aie indexes — the atlas access pattern.
---

One junction query across N archives — the atlas access pattern at molecule
resolution. Returns per-sample, per-cell counts byte-equal to running the
single-index query against each archive separately.

## Usage

```sh
aie federate [OPTIONS] <ARCHIVES>... <LOCUS>
```

```sh
aie federate pbmc1k.aie pbmc5k.aie pbmc10k.aie chr1:155234452-155235327
```

## Arguments

| Argument | Description |
|---|---|
| `<ARCHIVES>...` | Two or more `.aie` archives |
| `<LOCUS>` | Junction, written `chrom:donor-acceptor` (0-based, exact) |

## Options

| Option | Default | Description |
|---|---|---|
| `--top <N>` | `5` | Top cells reported per archive; under `--format`, `0` means all |
| `--format <FORMAT>` | — | Opt into uniform `text`, `tsv`, or `json` output |
| `-o, --output <PATH>` | stdout | Atomically publish uniform output without replacing a path |

## Uniform result contract

Omitting `--format` preserves the historical human-readable output bytes and
its historical `--top 0` behavior. With `--format`, the result schema is
`gravlax.federate.junction.result.v1`. Its `archives` table is a sequence in
caller archive order; its `counts` table is a sequence ordered by archive and
rank. The selected top-N subset within each archive uses UMI count descending,
then the visible barcode ascending. This tie-break selects reproducibly without
claiming that barcode order has scientific meaning.

The typed summary carries the cross-archive UMI and cell totals. Each archive
row also distinguishes `present`, `junction_absent`, and `chromosome_absent`
and reports exact available/emitted/truncated row counts. Rooted archive
identities are bound to their input positions, so repeated or byte-identical
archive inputs remain explicit rather than being collapsed in provenance.
Diagnostics remain on stderr. `--output` requires `--format`, checks its parent
before querying, and atomically installs the complete result without replacing
an existing file.

## Notes

- Archives are independent files: no shared catalogue, no merge step, no
  reprocessing. Any set of indexes — different tissues, chemistries, and cell
  counts — federates directly.
- Each archive is opened lazily and only the chunks holding the junction are
  decoded, so a federated query over many samples completes in seconds.
