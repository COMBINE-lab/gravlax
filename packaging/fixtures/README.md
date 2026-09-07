# GQ smoke fixture

The three `.aie.b64` files are synthetic, authenticated archives with four records, three
cells and three UMI classes each. They differ in chunking (100, 200 and 1,000 bases),
giving distinct roots. They contain compact unique geometry and a grouped
multimapper signature, not biological or personal data.

Regenerate to a new temporary file using the ignored
`gq::tests::export_release_smoke_fixture` test with `GRAVLAX_GQ_SMOKE_OUTPUT` set,
then base64-encode it. The generator is in `crates/aie/src/gq/tests.rs`.
Set `GRAVLAX_GQ_SMOKE_CHUNK_BP=1000` for the coarse fixture.
Set `GRAVLAX_GQ_SMOKE_CHUNK_BP=200` for the middle fixture.

`packaging/gq_smoke.py` decodes the fixture into a temporary directory, runs actual
queries through the supplied executable, verifies exact counts on a two-member
federation, and checks no-clobber result publication. It needs only Python's
standard library and never downloads data.

It also runs `packaging/collection_relocation_smoke.py`, which copies a three-archive,
two-layer collection, makes the original paths unavailable, and verifies unchanged
results, archive-I/O counters and file hashes using relative content-root locations.
The standalone relocation script accepts `--relocation-dir` to select a different
filesystem and `--report` to record its evidence.
