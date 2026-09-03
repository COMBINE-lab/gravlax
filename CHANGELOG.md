# Changelog

This file records user-visible changes in each Gravlax release.

## [0.1.1] - 2026-09-03

This is the first complete Gravlax distribution. Version 0.1.0 established the
Rust packages on crates.io, but its immutable GitHub release contains only the
distribution manifest: native archives, installers, the vendored source
archive, and Python packages were not published for that version. Install
version 0.1.1 when using a packaged release.

### Distribution changes

- Publish native archives and installers for 64-bit GNU and musl Linux, Intel
  and Apple Silicon macOS, and 64-bit Windows, together with checksums, build
  attestations, an SPDX dependency inventory, and a vendored source archive.
- Publish the `gravlax-client` wheel and source distribution to PyPI from the
  exact files attached to the immutable GitHub release.
- Pin portable release code generation explicitly: `x86-64` for 64-bit Linux
  and Windows, `penryn` for Intel macOS, and `apple-m1` for Apple Silicon.
  Intel macOS artifacts target macOS 10.12 or newer; Apple Silicon artifacts
  target macOS 11.0 or newer.
- Require successful validation and artifact assembly before creating the
  GitHub release, and keep Python 3.10 packaging validation compatible with
  the client's supported Python versions.

There are no archive-format, result-format, or command-interface changes from
version 0.1.0.

## [0.1.0] - 2026-09-03

This is the first public release of Gravlax. It introduces compact,
molecule-resolved evidence archives for single-cell RNA-seq and lets an
annotation be changed or compared without realigning the source reads.

### Highlights

- Build authenticated `.aie` v2 archives from annotation-free alignments;
  validate, inspect, seal, stamp, extend, replay, and combine those archives.
- Replay Gene, GeneFull, and Velocyto-style matrices against a GTF or compiled
  AIC annotation while retaining fixed alignment and barcode-correction
  decisions.
- Query regions, junctions, junction sets, splice events, splice graphs, 3′
  endpoints, and transcript-compatibility classes at cell, group, bulk, or
  multi-sample scope where supported.
- Compare two annotations on the same retained evidence and report signed
  count changes, class transitions, non-exclusive causes, and bounded
  witnesses.
- Define reusable projects and analysis plans with content-bound inputs,
  explicit biological intent, safe resume, and a loopback-only Explorer for
  building and inspecting plans.
- Request typed text, TSV, or JSON results from supported commands without
  changing their established default output; operation reports and directory
  bundles use the same schema and provenance model.
- Control the `aie` executable and validate typed results from Python with the
  dependency-light `gravlax-client` source package.
- Install the `aie` executable from crates.io. Native archives, installers,
  the vendored source archive, and the PyPI package first become available in
  version 0.1.1.

### Compatibility notes

- New archives use `.aie` v2. Seekable v1 archives remain readable and can be
  sealed into authenticated v2 containers without recompressing section
  payloads.
- Transcript-equivalence results describe compatibility with the retained
  archive representatives; they are not transcript-abundance estimates or
  full-isoform phasing.
- The Python distribution does not embed the Rust executable. Install `aie`
  separately and keep its version aligned with `gravlax-client`.
- Gravlax 0.1 is an initial public interface. Result and archive formats are
  versioned so readers can reject incompatible future changes explicitly.

[0.1.1]: https://github.com/COMBINE-lab/gravlax/releases/tag/v0.1.1
[0.1.0]: https://github.com/COMBINE-lab/gravlax/releases/tag/v0.1.0
