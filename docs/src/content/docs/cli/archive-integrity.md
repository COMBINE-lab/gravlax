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

The established human output and `--json` object remain available unchanged.
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
bytes read. The streamed `sections` table carries exact raw and compressed
byte counts. Machine-readable stdout contains no progress text.

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
