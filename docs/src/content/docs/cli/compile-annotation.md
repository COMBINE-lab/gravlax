---
title: aie compile-annotation
description: Compile a GTF once into a deterministic, checksummed annotation artifact for fast reuse.
---

Parse and index an uncompressed GTF once, then reuse the compiled `.aic`
artifact anywhere the CLI accepts `--gtf`.

```sh
aie compile-annotation gencode.v49.annotation.gtf \
  --out gencode.v49.annotation.aic

aie replay-rows sample.aie \
  --gtf gencode.v49.annotation.aic \
  --barcodes barcodes.tsv \
  --out-dir counts-v49/
```

The input GTF remains the scientific source of record. The `.aic` is a
derived cache containing the gene and contig dictionaries, transcript exon
models, and the augmented transcript-overlap index replay would otherwise
rebuild for every process. It is independent of any `.aie` archive, so one
compiled annotation can be shared by every sample aligned to a compatible
reference coordinate system.

Current AIC v2 artifacts also retain every stable transcript ID and each source
`exon_id` with its exact source interval. Those dictionaries support
ambiguity-detecting gene, transcript, and exon resolution and transcript-
compatibility queries. They are deliberately separate from the overlapping or
touching exon intervals merged for assignment, so adding identifier metadata
does not change replay semantics.

## Usage

```sh
aie compile-annotation <INPUT.gtf> --out <OUTPUT.aic>
```

The command refuses to overwrite a destination and installs a completed file
with an atomic, no-replace link. Recompiling the same GTF produces identical
bytes.

Add `--report-format text|tsv|json` to replace the historical status line with
a versioned `gravlax.annotation.compile.result.v1` report. The report contains
a typed scientific summary and a one-row `artifacts` table. Provenance binds
the exact GTF bytes consumed, and the artifact row carries the AIC v2 payload
identity already committed by its checksummed header; elapsed time is a diagnostic
on standard error rather than part of the reproducible result.
`--report-output FILE` publishes that report atomically without replacing an
existing file:

```sh
aie compile-annotation genes.gtf \
  --out genes.aic \
  --report-format json \
  --report-output genes.compile.json
```

Omitting `--report-format` preserves the legacy status output.

## Validation and compatibility

Compiled annotations have their own eight-byte magic, explicit format
version, declared little-endian payload length, and BLAKE3 payload checksum.
The reader rejects future versions, truncation, trailing data, invalid UTF-8,
invalid dictionary references or strands, malformed exon lists, hostile
allocation counts, and an overlap index inconsistent with the serialized
transcripts.

An `.aic` is accepted transparently by existing `--gtf` arguments. It does
not change assignment, replay, velocity, query-annotation, or discovery
semantics. On the evaluated GENCODE v49 artifact, GTF and `.aic` inputs
produced byte-identical Gene replay matrices, discovery output, and annotated
junction JSON.

The format version is deliberately separate from the `.aie` archive version.
Regenerate the `.aic` with the current Gravlax binary if a future release
changes the compiled-annotation format.

Legacy AIC v1 files remain readable. They support replay and the gene metadata
they actually contain, but lack the exact transcript/exon identifier
dictionaries. Transcript or exon resolution therefore returns a structured
metadata-unavailable error and asks for recompilation rather than guessing an
identifier or claiming it is absent from the source annotation.
