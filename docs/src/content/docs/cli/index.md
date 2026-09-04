---
title: CLI overview
description: The aie command-line interface at a glance.
---

Gravlax ships one binary, `aie`. Run `aie --help` for the command list and
`aie <command> --help` for the full option set of any subcommand. `-h` is the
short spelling of `--help`; `aie --version` (or `-V`) prints the installed
version.

## User-facing commands

| Command | Purpose |
|---|---|
| [`project`](/gravlax/cli/projects/) | Create a portable workspace and register named input resources |
| [`plan`](/gravlax/cli/projects/) | Check, explain, snapshot, and run versioned YAML/JSON analyses |
| [`doctor`](/gravlax/cli/doctor/) | Diagnose the installation, project, and selected archives or annotations |
| [`explore`](/gravlax/cli/explore/) | Browse exact project artifacts and export checked scientific plans from a loopback-only, read-only UI |
| [`resolve`](/gravlax/cli/resolve/) | Resolve gene, transcript, and exon identifiers against an explicit reference identity |
| [`ingest`](/gravlax/cli/ingest/) | Generate chemistry-specific STAR recipes and preflight BAM/whitelist inputs |
| [`ingest-archive`](/gravlax/cli/ingest-archive/) | Build the `.aie` index from a compatible tagged, coordinate-sorted BAM |
| [`compile-annotation`](/gravlax/cli/compile-annotation/) | Compile a GTF once into a checksummed, reusable `.aic` artifact |
| [`export-molecule-bam`](/gravlax/cli/export-molecule-bam/) | Export the exact post-correction molecule abstraction for interchange or a function-matched BAM/CRAM baseline |
| [`replay-rows`](/gravlax/cli/replay-rows/) | Quantify a compatible GTF from an index using Gene or Velocyto semantics |
| [`compare-annotations`](/gravlax/cli/compare-annotations/) | Compare two bound annotations on one fixed archive and explain count changes |
| [`query`](/gravlax/cli/query/) | Indexed region/junction queries, Boolean same-record predicates, APA, discovery, and transcript compatibility |
| [`query … transcript-ecs`](/gravlax/cli/transcript-ecs/) | Derive annotation-conditional transcript compatibility sets for archived UMI classes |
| [`federate`](/gravlax/cli/federate/) | One junction query across N indexes |
| [`cohort`](/gravlax/cli/cohort/) | Coordinate-defined splice events across named indexes and groups |
| [`collection`](/gravlax/cli/collection/) | Build a content-addressed federation and reverse-search events across samples, donors, and cell groups |
| [`extend`](/gravlax/cli/extend/) | Propose evidence-supported per-gene 3′ annotation extensions |
| [`stamp-genome`](/gravlax/cli/stamp-genome/) | Bind the reference-genome signature used by sequence-consulting analyses |
| [`seal-archive`](/gravlax/cli/archive-integrity/) | Copy a legacy v1 archive into an authenticated v2 container without recompressing its sections |
| [`inspect-archive`](/gravlax/cli/archive-integrity/) | Report archive identities and optionally verify every compressed payload |
| [`completions`](/gravlax/cli/completions/) | Generate Bash, Zsh, or Fish completions from the installed command graph |

## Conventions

- **Loci** are written `chrom:start-end` using 0-based, half-open intervals.
  Exact junction loci use `chrom:donor-acceptor`.
- **Whitelist** files are the 10x barcode whitelists (one 16 bp barcode per
  line, e.g. `3M-february-2018.txt`).
- **Barcode lists** (`--barcodes`) define output column order for emitted
  matrices; pass the `barcodes.tsv` you want columns aligned to.
- Matrix outputs are Matrix Market (`matrix.mtx` + `features.tsv` +
  `barcodes.tsv`), the same layout STARsolo emits.
- Scientific result commands select typed text/TSV/JSON with `--format` and
  optionally publish through `--output`. Commands whose primary product is an
  archive, matrix, BAM, annotation, or revised GTF use the parallel
  `--report-format`/`--report-output` interface. Omitting these flags preserves
  historical presentations.
- The default thread budget is 24; heavy stages (ingest compression, replay
  decode) parallelize automatically.

The [workflow and interfaces guide](/gravlax/workflow/) explains when to use
direct commands, a checked project plan, Explorer, or Python, and explains the
uniform-output contract alongside byte-compatible historical formats.
