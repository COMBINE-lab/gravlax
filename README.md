<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="assets/logo-dark.svg">
    <img src="assets/logo-light.svg" alt="Gravlax" width="380">
  </picture>
</p>

**Align once and query forever — a compact molecular-evidence index for annotation replay in
single-cell RNA-seq.**

Gravlax builds a compact, molecule-resolved evidence index (`.aie`) from raw 10x reads and a
genome — **no annotation touches the build**. Supplying a compatible GTF at query time replays
Gene, GeneFull, and Velocyto-style quantification from fixed genome alignments and fixed barcode
correction. Archive-sourced and BAM-sourced Gravlax replay are byte-identical; compared with a
fresh annotation-aware STARsolo run, the archive's two-representative encoding differs by
0.22–0.45% of UMI mass across the four evaluated datasets. It occupies 11–18 bits per input read
and opens in ~10 ms. Against a function-matched CRAM containing the same post-correction molecule
placements and UMI graph, it is 2.10–2.55× smaller and 6.64–13.69× faster to quantify across the
four evaluated datasets. Tag-preserving BAM/CRAM retain additional read-level information and are
reported separately.

Because the evidence — not the interpretation — is what is stored, reanalysis becomes a query:

- **Fast, bounded annotation replay** — quantify compatible past or future GTFs in seconds,
  without realignment or materializing the archive's molecule table.
- **Indexed queries** — cell/group-scoped region, exact-junction, junction-set, automatic
  splice-event, and junction-enumeration queries, plus 3′-site queries, from milliseconds to
  chromosome-scale scans.
- **Compiled annotations and query panels** — compile a GTF once into a guarded `.aic`, then
  share archive opens and chunk decodes across batched region/junction predicates.
- **Paired annotation comparison** — replay two bound annotations independently and report exact
  signed count changes, class transitions, non-exclusive causes, and bounded witnesses within the
  retained archive quotient and fixed alignment/barcode policy.
- **Experimental transcript compatibility** — derive deterministic annotation-conditional
  transcript sets for archived UMI classes, with explicit ambiguity, conflict, and completeness;
  these are not abundance estimates or full-transcript phasing.
- **Discover → replay** — find unannotated loci from the index alone, then re-quantify them. A
  bounded residual-site mode raises complete-denominator later-annotation recall from 44.7% to
  72.5% at 2.16× the candidate volume.
- **EM multimapper recovery** — packed, cell-sharded cross-cell EM over the stored paralog
  evidence, with exact disk-backed support shards for large archives (3.67 GiB peak RSS and
  35.86 s for the evaluated 10k-cell masked-evidence analysis).
- **Federation** — discover one coordinate-defined splice-event catalogue and reduce it exactly
  across named samples and groups, or build a common molecular splice graph with complete
  per-sample matrices and replicate-aware path-usage contrasts.
- **Content-addressed federation** — deterministic `.aicollection` layers bind each source by its
  authenticated archive root. Optional local shape routes accelerate exact point-junction and
  junction-set queries without copying molecule evidence; region queries use the collection's
  interval routes. Source archives remain authoritative and immutable.

## Build

```sh
cargo build --release
# binary: target/release/aie
```

Requires Rust 1.98+. No system dependencies beyond a C toolchain (for zstd).

Check the installation and create an optional portable workspace before a new
analysis:

```sh
aie doctor
aie project init my-analysis --name my-analysis
cd my-analysis
```

Projects register stable input names in `aie-project.yaml`. Versioned YAML/JSON
plans can then be validated with `aie plan check --explain`, resolved into an
exact content-addressed snapshot, and run through the existing typed commands.
`aie explore` provides a loopback-only, read-only scientific plan builder: it
resolves gene, transcript, and exon identifiers against registered annotation
identity, explains the evidence route and exclusions, and exports synchronized
plan YAML/JSON, CLI, and Python forms. It also browses exact stored plans and
results. Direct path-based commands remain fully supported.

The [workflow and interfaces
guide](docs/src/content/docs/workflow.md) is the consolidated entry point for
projects, typed step-to-step dataflow, identity-checked resume, biological
identifier resolution, Explorer, Python/AnnData, and the exact boundary of the
shared result contract. The small
[`examples/demo-project`](examples/demo-project/README.md) exercises the
project, plan, resume, and Explorer paths without downloading a dataset.

## Quick start

```sh
# 1. Align annotation-free (STAR two-pass, no GTF, secondaries kept) → BAM
aie ingest recipe --chemistry 10x-3p-v3
# 2. Build the index once:
aie ingest check align.bam --whitelist 3M-february-2018.txt --chemistry 10x-3p-v3
aie ingest-archive align.bam --whitelist 3M-february-2018.txt --out sample.aie
aie inspect-archive sample.aie --format json                               # typed content identity

# 3. Query forever:
aie compile-annotation gencode.v49.gtf --out gencode.v49.aic
aie resolve gencode.v49.aic TP53 --assembly GRCh38.p14 \
  --annotation "GENCODE 49" --format json                                # typed identifier result
aie replay-rows sample.aie --gtf gencode.v49.aic \
  --barcodes barcodes.tsv --out-dir counts/                                # replay matrices
aie compare-annotations sample.aie --annotation-a gencode.v44.gtf \
  --annotation-b gencode.v49.aic --assembly GRCh38.p14 \
  --annotation-a-label "GENCODE 44" --annotation-b-label "GENCODE 49"      # exact archive counterfactual
aie query sample.aie transcript-ecs --annotation-file gencode.v49.aic \
  --assembly GRCh38.p14 --annotation-label "GENCODE 49" \
  --feature gene:ENSG00000141510 --format json                             # experimental compatibility sets
aie query sample.aie region chr1:1000000-2000000 --format json             # per-cell evidence
aie query sample.aie junction chr2:1234567-1250000 --format json           # junction support
aie query sample.aie junctions chr2:1200000-1300000 --with-cells --format tsv
aie query sample.aie jset --include chr2:1234567-1250000 \
  --exclude chr2:1234567-1260000 --groups cell-types.tsv --format json     # group splice usage
aie query sample.aie events chr2:1200000-1300000 \
  --groups cell-types.tsv --min-informative 10 --format json               # discover + reduce events
aie query sample.aie splice-graph chr2:1200000-1300000 \
  --groups cell-types.tsv --min-path-umis 2 --format json                  # molecular path fragments
aie query sample.aie batch --plan panel.tsv --top 20 \
  --format json --output panel.json                                        # query panel
aie federate a.aie b.aie c.aie chr2:1234567-1250000 --format json         # across samples
aie cohort events chr2:1200000-1300000 \
  --sample A=a.aie --sample B=b.aie --groups A=a-groups.tsv \
  --min-row-informative 10 --format json                                   # exact cohort event table
aie cohort events chr2:1200000-1300000 \
  --sample A=a.aie --sample B=b.aie --min-row-informative 10 \
  --sparse-dir cohort-events/                                              # zero-reconstructible tables
aie cohort splice-graph chr2:1200000-1300000 \
  --design experiment.tsv --contrast control:treated --format json        # sample-level path test
aie collection build --sample A=a.aie --sample B=b.aie \
  --shape-routes --out atlas.aicollection                                 # exact local-shape routes
aie collection inspect atlas.aicollection --verify-routes --format json   # reconstruct routes
aie collection junction atlas.aicollection chr2:1234567-1250000 \
  --explain --format json --output junction.json                           # uniform bundle
aie collection region atlas.aicollection chr2:1200000-1300000 --format tsv
aie collection jset atlas.aicollection --include chr2:1234567-1250000 \
  --exclude chr2:1234567-1260000 --format json
aie dev em sample.aie --gtf gencode.v49.gtf                               # multimapper experiment
```

See the [annotation-comparison](docs/src/content/docs/cli/compare-annotations.md)
and [transcript-equivalence-class](docs/src/content/docs/cli/transcript-ecs.md)
references for their typed tables and scientific limits. The comparison is
exact only for the archive's retained evidence quotient under the fixed
alignment and barcode-correction policy. Transcript equivalence classes are
derived from the representatives retained in the archive, not from every
source read, and do not estimate abundance or phase complete isoforms. On direct
path-based commands, `--assembly` is a caller assertion recorded in provenance;
a project plan can additionally verify it against registered coordinate-resource
compatibility.

At the direct CLI boundary, supported query, federation, cohort, and collection
commands can emit the shared `gravlax.result-envelope.v1` contract.
Result-streaming commands use `--format text|tsv|json` with atomic
no-clobber `--output`; commands that advertise a separate operation report use
the parallel `--report-format` and `--report-output` interface. This includes
`ingest-archive`, `replay-rows`, `stamp-genome`, `seal-archive`,
`compile-annotation`, `export-molecule-bam`, and `extend`; `inspect-archive` and
artifact-producing `collection build` instead use `--format`, so the
distinction is explicitly per command. Omitting these opt-in flags preserves
each command's default output and primary-artifact behavior. Project plans
expose the contract for the source-plan v1 step kinds they support, not every
direct command. The `aie dev` subcommands and operational checks use their own
machine-readable formats. The Python client distinguishes command-specific and
shared result formats explicitly; native Arrow IPC is not yet a general CLI
output and no R client is currently provided. See [Python and
AnnData](docs/src/content/docs/python.md) and the [result-contract
note](docs-notes/scientific-intent-and-output-contract.md).

Ordinary archive-backed Gene replay streams bounded genomic chunk batches by default. The
diagnostic `--eager` flag retains the full-materialization reference path for result comparison
and profiling; both modes perform the same global UMI-class aggregation and emit identical bytes.

Discovery uses conservative span claiming by default.
For human 10x 3′ data, the evaluated higher-recall mode preserves those calls and adds bounded
terminal-site clusters from span-overlapping evidence that is incompatible with every overlapping
transcript:

```sh
aie query sample.aie discover --gtf gencode.v32.gtf \
  --claim-mode residual-sites --residual-min-umis 75 \
  --emit-gtf novel-loci.gtf
```

The threshold is an evaluated operating point, not a universal default for other chemistries.

For an explicit post-correction interchange form and a function-matched BAM/CRAM container
comparison, export the molecule abstraction without inventing nucleotide UMI values:

```sh
aie export-molecule-bam sample.aie --fai GRCh38.fa.fai --out sample.molecules.bam
aie replay-rows sample.molecules.bam --from-molecule-bam --gtf gencode.v49.gtf \
  --barcodes barcodes.tsv --out-dir replay-from-molecule-bam
```

The custom-tag contract and its interoperability limits are documented in
[`docs-notes/molecule-bam.md`](docs-notes/molecule-bam.md).

New archives use the authenticated `.aie` v2 container. Its root commits the directory and the
BLAKE3 digest of every compressed section, so ordinary operations authenticate the directory at
open and verify only the payloads they select. A complete audit is available with
`aie inspect-archive sample.aie --verify-content`. Legacy seekable v1 archives remain readable and
can be converted without recompressing their section payloads:

```sh
aie seal-archive legacy.aie --out rooted.aie --json
```

A collection build reads archive metadata and indexes, but never molecule chunks. For rooted v2
sources, identity comes directly from the authenticated directory; `--shape-routes` additionally
reads the shape dictionary to derive exact, source-bound route blocks. Incremental builds write a
new immutable layer rather than rewriting their parent. Collection support totals are pruning upper
bounds, never substitutes for cell/group counts. Group maps and annotations stay query-time inputs;
use the cohort commands for direct group-scoped reductions. Query JSON separates source identity,
source execution, collection-sidecar, shape-route, and total logical bytes. The format, integrity
guards, and deliberate limits are documented in
the [`aie collection` reference](docs/src/content/docs/cli/collection.md) and
[`docs-notes/collection-index-spec.md`](docs-notes/collection-index-spec.md).

Run `aie <command> --help` for the full option set of each subcommand.

## Repository layout

| Crate | Responsibility |
|---|---|
| `crates/evidence-io` | `.aie` container: chunked streams, static rANS + zstd coding, lazy open |
| `crates/ingest` | annotation-free BAM → molecule evidence (UMI classes + edges, paralog patterns) |
| `crates/anno` | GTF parsing and annotation compilation (exon models, junction sets) |
| `crates/replay` | Reserved library boundary (current replay implementation is in `crates/aie`) |
| `crates/eval` | Reserved library boundary (current evaluation commands are in `crates/aie`) |
| `crates/aie` | the `aie` CLI |

`docs-notes/format-spec.md` documents the on-disk `.aie` format; the optional derived collection
format is documented separately in `docs-notes/collection-index-spec.md`.

## Citation

Manuscript in preparation. This repository is the reference implementation.

## License

BSD 3-Clause. See [LICENSE](LICENSE).
