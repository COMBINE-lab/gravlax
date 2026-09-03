---
title: Releases and distribution
description: Install native releases or the Python client, verify downloads, and prepare a Gravlax release.
---

Gravlax releases use one version for the Rust crates, the `aie` executable,
the `gravlax-client` Python distribution, and the documentation. A release is
available only after it appears on the relevant service; before the first
release is published, use the [source installation](/gravlax/installation/).

## Install the `aie` command

After the release is available on crates.io, install the Rust package named
`gravlax`:

```sh
cargo install gravlax
aie --version
```

The package name is `gravlax`; the installed command remains `aie`.

GitHub Releases also provides cargo-dist installers. On Linux or macOS:

```sh
curl --proto '=https' --tlsv1.2 -LsSf \
  https://github.com/COMBINE-lab/gravlax/releases/latest/download/gravlax-installer.sh | sh
```

On Windows PowerShell:

```powershell
powershell -ExecutionPolicy ByPass -c "irm https://github.com/COMBINE-lab/gravlax/releases/latest/download/gravlax-installer.ps1 | iex"
```

The installers select from these native targets:

| Target | Platform |
|---|---|
| `x86_64-unknown-linux-musl` | 64-bit Linux, statically linked C runtime; preferred portable Linux build |
| `x86_64-unknown-linux-gnu` | 64-bit Linux compatible with the Ubuntu 22.04 build environment |
| `x86_64-apple-darwin` | Intel macOS |
| `aarch64-apple-darwin` | Apple Silicon macOS |
| `x86_64-pc-windows-msvc` | 64-bit Windows |

Each native archive contains `aie` (or `aie.exe`) and its release metadata.
Generate Bash, Zsh, or Fish completions from the installed binary with
`aie completions`. Platforms not listed here can build from source with Rust
1.98 or newer.

## Install the Python client

The Python distribution is named `gravlax-client` and its import name is
`gravlax`:

```sh
python -m pip install gravlax-client
python -c "import gravlax; print(gravlax.__version__)"
```

The Python package controls and reads results from a compatible `aie`
executable; it does not embed the Rust program. Install the executable
separately and confirm that `aie --version` is available on `PATH`.

## Verify a downloaded release

Native cargo-dist archives have SHA-256 checksum files on the GitHub release
page. The release also includes a checksum manifest for the Python packages,
vendored source archive, and SPDX software bill of materials. Verify the
downloaded file before installing it. For example, on Linux:

```sh
sha256sum --check <downloaded-checksum-file>.sha256
```

GitHub build attestations can additionally bind an artifact to this repository
and its release workflow:

```sh
gh attestation verify <downloaded-artifact> --repo COMBINE-lab/gravlax
```

Checksums detect changed bytes. Attestations identify the GitHub Actions run
that produced them; you should still confirm the repository and version you
intended to install.

## Containers

The repository's `Dockerfile` builds a minimal, non-root runtime image with the
`aie` command:

```sh
docker build --file Dockerfile --tag gravlax:local .
docker run --rm gravlax:local --version
```

No container-registry location is documented until an image is published
there.

## Bioconda

After a release recipe has been accepted by Bioconda, install it with:

```sh
conda install -c bioconda gravlax
```

The Bioconda package installs `aie` and its shell completions. It does not
include `gravlax-client`, which remains a separate PyPI package. Until the
recipe is visible on Bioconda, use a native installer, crates.io, or a source
build.

## Maintainer release procedure

The repository uses cargo-dist for native archives and installers. The
`scripts/bump-and-publish` command keeps the Rust, Python, documentation, and
local Conda versions synchronized and creates the annotated release tag. The
release manager needs Rust 1.98, Python 3.11 or newer, Node.js 22.19 or newer,
GitHub CLI, and the normal Cargo and npm tooling.

1. Add dated release notes to `CHANGELOG.md`, commit all release changes, and
   work from a clean `main` branch.
2. Preview the release without changing files or tags:

   ```sh
   ./scripts/bump-and-publish 0.1.0 --dry-run --check-history
   ```

3. Update coordinated versions when needed, run the package and test suite,
   and create the local annotated tag:

   ```sh
   ./scripts/bump-and-publish 0.1.0 --prepare
   ```

   For the initial `v0.1.0` tag, where the repository already declares
   version `0.1.0`, add `--allow-current`.
4. For the first crates.io publication only, publish the dependency-ordered
   crate set with a crates.io token, then configure crates.io Trusted
   Publishing for later tagged releases:

   ```sh
   CARGO_REGISTRY_TOKEN=<token> \
     ./scripts/bump-and-publish 0.1.0 --publish-crates \
       --confirm-publish v0.1.0
   ```

5. Revoke the initial crates.io token. Configure each crate's Trusted Publisher
   for top-level workflow `release.yml` and environment `crates-io`. Configure
   the pending PyPI Trusted Publisher for `gravlax-client`, top-level workflow
   `publish-python.yml`, and environment `pypi`.
6. Enable immutable GitHub Releases; add an active `refs/tags/v*` ruleset that
   prevents tag updates and deletion; and configure the protected GitHub
   environments `release`, `crates-io`, and `pypi` before pushing the first
   tag. Restrict `release` and `crates-io` to selected tags matching `v*`, and
   restrict `pypi` to the selected branch `main`; do not allow every ref.
   Limit creation or bypass of protected release tags to the designated release
   maintainers.
7. Push `main` and its prepared tag atomically:

   ```sh
   ./scripts/bump-and-publish 0.1.0 --push
   ```

   Do not push release tags directly; the helper checks the repository,
   protected release settings, ancestry, and annotated tag before sending the
   branch and tag together.

8. After the complete tag workflow succeeds, dispatch the protected Python
   publication and approve its `pypi` environment deployment:

   ```sh
   ./scripts/bump-and-publish 0.1.0 --dispatch-python
   ```

The tag starts the generated cargo-dist workflow. It builds the five native
targets and installers, builds the Python wheel/source distribution, vendored
source archive, checksums, and SBOM, publishes the Rust crates, and then creates
one immutable GitHub release containing the complete artifact set and changelog.
The separate `--dispatch-python` step verifies that immutable release and
publishes its exact wheel and source distribution through PyPI Trusted
Publishing. Safe retries skip an already-complete matching version or upload
only a missing release file after an interrupted partial publication; unexpected
filenames or bytes are rejected.

Publication is permanent. The helper refuses mismatched versions, an
uncommitted worktree, a pre-existing tag, missing release notes, or a push from
a branch other than `main`. `--check-only` validates the already-declared
version without creating a tag.

After the GitHub release is public and its source archive verifies, render and
submit the Bioconda recipe by following `packaging/bioconda/README.md` in the
tagged source tree.
