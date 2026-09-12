---
title: Installation
description: Install a Gravlax release or build the aie command from source.
---

Gravlax is implemented in Rust and ships as a single binary, `aie`.

## Install a release

The recommended installation is the Bioconda package:

```sh
conda install -c bioconda gravlax
aie --version
```

The package is named `gravlax`; the command it installs is `aie`. Builds are
available for `linux-64`, `linux-aarch64`, `osx-64`, and `osx-arm64`.

Prebuilt Linux, macOS, and Windows archives and installers are attached to
every [GitHub release](https://github.com/COMBINE-lab/gravlax/releases) and are
described on the [releases and distribution page](/gravlax/distribution/):

```sh
curl --proto '=https' --tlsv1.2 -LsSf \
  https://github.com/COMBINE-lab/gravlax/releases/latest/download/gravlax-installer.sh | sh
```

The same version is published to crates.io, so Cargo can also install it:

```sh
cargo install gravlax
```

The Python client is a separate distribution, `gravlax-client` on PyPI; it
drives an `aie` executable rather than embedding one. See
[Python and AnnData](/gravlax/python/).

## Source-build requirements

- **Rust 1.89 or newer** — install via [rustup](https://rustup.rs) if needed.
- A C toolchain (for the vendored zstd). No other system dependencies.

## Build from source

```sh
git clone https://github.com/COMBINE-lab/gravlax.git
cd gravlax
cargo build --release
```

The executable is `target/release/aie`. Install it under the Cargo binary
prefix (normally `~/.cargo/bin`):

```sh
cargo install --path crates/aie
```

Make sure Cargo's binary directory is on `PATH`, then verify both the binary
and the local workspace:

```sh
aie --version
aie doctor
```

`aie doctor` reports setup checks independently and tells you how to correct
anything it finds. Missing STAR or samtools is only a warning: neither is
needed on a machine that only replays and queries existing archives. See the
[`aie doctor` reference](/gravlax/cli/doctor/) for archive validation and JSON
output.

## Shell completions

Generate completions from the installed binary so they always match its
command set:

```sh
# Bash, for this session
source <(aie completions bash)

# Fish, for future sessions
mkdir -p ~/.config/fish/completions
aie completions fish > ~/.config/fish/completions/aie.fish
```

Persistent Bash and Zsh setup is covered on the
[shell completions page](/gravlax/cli/completions/).

## First project

Projects are optional—the direct commands continue to accept ordinary
paths—but they provide named inputs, checked plans, and exact resolved-plan
snapshots:

```sh
aie project init my-analysis --name my-analysis
cd my-analysis
aie project show
```

Continue with the [workflow and interfaces guide](/gravlax/workflow/), the
detailed [projects and plans reference](/gravlax/cli/projects/), or the
[direct-command quick start](/gravlax/quickstart/).

## Companion tools

Gravlax consumes an annotation-free alignment, so you will also want:

- **[STAR](https://github.com/alexdobin/STAR)** — to produce the one-time,
  annotation-free BAM the index is built from (see the
  [quick start](/gravlax/quickstart/)).
- **samtools** — convenient for sorting and inspecting the ingest BAM.

Neither is needed after the index is built: every replay and query runs from
the `.aie` file alone.

## Repository layout

| Crate | Responsibility |
|---|---|
| `crates/evidence-io` | `.aie` container: chunked streams, static rANS + zstd coding, lazy open |
| `crates/ingest` | annotation-free BAM → molecule evidence (UMI classes + edges, paralog patterns) |
| `crates/anno` | GTF parsing and annotation compilation (exon models, junction sets) |
| `crates/replay` | Reserved library boundary; current replay implementation is in `crates/aie` |
| `crates/eval` | Reserved library boundary; current evaluation commands are in `crates/aie` |
| `crates/aie` | the `aie` CLI |
