---
title: aie resolve
description: Resolve biological identifiers against an explicit assembly and annotation.
---

`aie resolve` turns gene symbols and stable gene, transcript, or exon IDs into
zero-based, half-open genomic loci. The assembly and annotation release are
required: an identifier is never interpreted against an unnamed reference.

```sh
aie resolve gencode.v49.annotation.aic TP53 transcript:ENST00000269305 \
  --assembly GRCh38.p14 \
  --annotation "GENCODE 49" \
  --format json \
  --output resolved.json
```

The annotation input can be a source GTF or a compiled `.aic`. Prefix an
identifier with `gene:`, `transcript:`, or `exon:` when its kind is known;
unprefixed input searches all supported kinds. Matching is case-sensitive.
Unversioned stable IDs may match the version present in the selected release,
but a supplied version must match exactly.

Duplicate symbols, multiple stable-ID versions, cross-kind matches, missing
identifiers, or unavailable legacy `.aic` metadata return an error for the
entire request. No partial output file is created.

## Arguments and options

| Argument or option | Default | Description |
|---|---|---|
| `<ANNOTATION_FILE>` | required | Source GTF or compiled `.aic` used for resolution |
| `<IDENTIFIER>...` | required | One or more symbols or stable IDs; optional `gene:`, `transcript:`, or `exon:` prefix constrains the kind |
| `--assembly <ASSEMBLY>` | required | Exact reference assembly identity recorded with the result |
| `--annotation <RELEASE>` | required | Exact annotation release or immutable label |
| `--annotation-digest <DIGEST>` | — | Require the observed annotation identity to equal `blake3:<64 lowercase hex>` |
| `--format <FORMAT>` | `text` | `text`, `tsv`, or typed-envelope `json` |
| `-o, --output <PATH>` | stdout | Atomically create a new result file; existing files are not replaced |

Current AIC v2 files retain exact gene, transcript, and source exon identifier
dictionaries. AIC v1 remains readable and can resolve the gene IDs and symbols
it contains, but it cannot answer transcript/exon requests. In that case the
command reports identifier metadata as unavailable and asks you to recompile
the source GTF; it does not guess that an identifier was absent from the
original annotation.

## Output

`--format text|tsv|json` selects the representation and `--output` writes a
new file instead of standard output. Existing files are never overwritten.
All three formats carry assembly, annotation, annotation-content digest, and
the versioned result schema. JSON is a typed
`gravlax.result-envelope.v1` result with schema
`gravlax.annotation.resolve.v1`, ready for the Python client:

```python
from gravlax import Client

resolved = Client().resolve(
    "gencode.v49.annotation.aic",
    ["TP53", "transcript:ENST00000269305"],
    assembly="GRCh38.p14",
    annotation="GENCODE 49",
)

print(resolved.table.records())
print(resolved.provenance.annotation_digest)
```

Each row includes the requested value, resolved kind and stable ID, match
basis, parent gene/transcript IDs, contig, strand, and exact locus. Features
with multiple loci produce one row per locus.
