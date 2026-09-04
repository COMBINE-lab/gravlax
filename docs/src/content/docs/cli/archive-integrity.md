---
title: Archive identity and sealing
description: Inspect authenticated archive identities and migrate legacy v1 containers to rooted v2.
---

`inspect-archive` reports the identity of an `.aie` archive. `seal-archive`
copies a legacy seekable v1 container into the authenticated v2 container
without decompressing, recompressing, or changing its encoded section
payloads.

## Inspect an archive

```sh
aie inspect-archive sample.aie
aie inspect-archive sample.aie --verify-content
```

Normal v2 inspection authenticates the directory/root and reports the
scheme-independent encoded-section identity. Payload digests are checked when
sections are selected. `--verify-content` additionally reads and verifies
every compressed payload. A v1 archive has no committed directory, so even
ordinary identity inspection requires one complete file scan.

Inspection also reports the archive's logical evidence capabilities. For a
current ingest it includes `gravlax.molecular-evidence.v2`, the full parsed
alignment-provenance manifest, and terminal-tail availability plus its event,
molecule, and routed-chunk counts. It also reports the current genome-reference
binding, including whether it was established at ingest or by `stamp-genome`.
It verifies an embedded junction catalogue's length, BLAKE3 digest, and parsed
row count against the manifest. Older archives state that alignment provenance
and terminal tails are **unavailable**, and label any genome signature as
legacy/unattributed;
inspection never guesses “one pass” or reports a missing tail section as zero.
An unknown logical schema, a declared section that is absent, an undeclared
capability section, or a side-section identity mismatch is an error.

The human output and `--json` object remain available; both include these
capability fields.
For the shared result contract, use:

```sh
aie inspect-archive sample.aie --format json
aie inspect-archive sample.aie --verify-content \
  --format tsv --output sample.identity.tsv
```

`--format` accepts `text`, `tsv`, or `json`. `--output` requires `--format`,
is checked before the archive is opened, and atomically installs a complete
file without replacing an existing path. The
`gravlax.archive.inspect-report.v1` summary contains the native and
encoded-section identities, format and file sizes, verification mode, and
bytes read, together with the same logical-schema, alignment-provenance,
terminal-tail, and genome-reference-binding fields. The streamed `sections`
table carries exact raw and compressed byte counts. Machine-readable stdout
contains no progress text.

## Seal a legacy archive

```sh
aie seal-archive legacy.aie --out rooted.aie
aie seal-archive legacy.aie --out rooted.aie \
  --report-format json --report-output seal-report.json
```

The destination archive is always installed without replacement. Sealing
validates every legacy frame, copies its exact compressed bytes, builds and
authenticates the v2 directory, verifies the completed output, and checks
that its `aie-encoded-sections-v1` identity equals the source identity before
installation.

Omitting report flags preserves the established human output and `--json`
record. `--report-format text|tsv|json` selects the uniform
`gravlax.archive.seal-report.v1` result; `--report-output` publishes it
atomically without replacement. The report path is preflighted before the
archive conversion starts. Scientific and byte-accounting totals are result
data, while source/output paths and invocation choices remain provenance.

The archive and operation report have separate commit points. A report
publication error can occur after the already-verified archive has been
successfully installed; retrying must use a new report path and must not
re-run sealing against the occupied archive destination.
