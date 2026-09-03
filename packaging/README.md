# Distribution maintenance

Gravlax has one version across its Rust workspace, Python package,
documentation package, and conda recipe. An annotated `vMAJOR.MINOR.PATCH` tag
starts the cargo-dist 0.32.0 workflow. cargo-dist builds and attests the native
archives and installers; a reusable job adds the Python distributions, the
offline source archive, and an SPDX SBOM. Rust publication finishes before the
complete GitHub release is announced. Python publication is a separate,
protected manual step that consumes that release.

The release manager needs Python 3.11 or newer (the scripts use `tomllib`),
Node.js 22.19 or newer, the repository-pinned Rust 1.98.0 toolchain,
cargo-dist 0.32.0, GitHub CLI, and the normal Cargo and npm tooling.

The release workflow is generated from `[workspace.metadata.dist]` in
`Cargo.toml`. Do not edit `.github/workflows/release.yml` by hand. After changing
cargo-dist configuration, install exactly version 0.32.0 and run:

```sh
dist generate --mode=ci
dist generate --mode=ci --check
dist plan --tag=v0.1.0
```

The generated native assets are named `gravlax-TARGET.tar.gz` (or `.zip` on
Windows), with per-file SHA-256 checksums and shell and PowerShell installers.
Linux releases include both an Ubuntu 22.04 GNU build and a portable
`x86_64-unknown-linux-musl` build; users on older enterprise or HPC systems
should select the musl archive.
The separately built, vendored source asset deliberately retains the stable
`gravlax-VERSION-source.tar.gz` name consumed by the Bioconda recipe. The custom
source, Python, checksum, and SBOM files join cargo-dist's artifact set before
the GitHub release is created, so the release becomes immutable with its full
asset set already attached.

## One-time service setup

Enable release immutability in the repository's GitHub settings before creating
v0.1.0; the setting applies only to releases created afterward.

Create an active repository ruleset for `refs/tags/v*` that prevents tag
updates and deletion. Permit only the designated release maintainers to create
or bypass protected release tags. This keeps the tag bound to the commit that
cargo-dist builds while the release workflow is running; immutable releases
protect the resulting release only after it has been created.

Create three protected GitHub environments and add required reviewers if the
repository plan supports them:

- `release` protects the source, Python, checksum, and SBOM artifact build and
  should allow only selected tags matching `v*`;
- `crates-io` protects the Rust Trusted Publishers and should allow release
  only from selected tags matching `v*`;
- `pypi` controls the Python Trusted Publisher and must allow only the selected
  branch `main`, because the manual workflow is dispatched from `main` after
  checking out and verifying the immutable tag.

Do not leave these environments open to every branch or tag. Registry Trusted
Publishing identifies a workflow and environment, so the environment's own
deployment restrictions are part of the publication boundary.

On PyPI, create a pending Trusted Publisher for project `gravlax-client`, owner
`COMBINE-lab`, repository `gravlax`, top-level workflow `publish-python.yml`, and
environment `pypi`. PyPI does not authorize reusable workflows as Trusted
Publishers, so this workflow is a separate manual dispatch that consumes the
immutable GitHub release. PyPI can create a new project through a pending
publisher, so no PyPI API token is needed for the first release.

crates.io cannot configure a Trusted Publisher until a crate has been published
once. For v0.1.0 only, create a narrowly scoped `publish-new` token and publish
the five crates locally in dependency order using the command below. Revoke the
token immediately. Then configure a Trusted Publisher for each of
`gravlax-evidence-io`, `gravlax-output`, `gravlax-anno`, `gravlax-ingest`, and
`gravlax`, all bound to repository `COMBINE-lab/gravlax`, workflow
`release.yml`, and environment `crates-io`. The crates.io OIDC check identifies
the top-level caller through GitHub's `workflow_ref` claim even though the
publishing steps live in reusable `publish-crates.yml`. Future releases use
temporary OIDC credentials and need no stored crates.io token.

Crate-name availability is first-come, first-served. Recheck all five names
immediately before the bootstrap publish.

## First release

Work from a clean, up-to-date `main`. The preparation command validates version
consistency, the publishable dependency graph, package contents, release notes,
the public source policy, tests, and documentation. It creates only a local
annotated tag:

```sh
./scripts/bump-and-publish 0.1.0 --dry-run --check-history
./scripts/bump-and-publish 0.1.0 --prepare --allow-current
```

The history check must pass before the repository is made public. It examines
the current tracked tree and every commit in `HEAD` ancestry for the former
private-project marker. Rewriting published Git history is a separate,
reviewable operation; the release script will not do it implicitly.

Before changing repository visibility, fetch every GitHub branch and tag and
audit every commit returned by `git rev-list --remotes=origin --tags`, not just
the release checkout. The current history requires rewriting the ancestry of
`main` and rewriting or deleting the remote `uniform-io` branch. Review all
other GitHub refs after the rewrite as well. The release script intentionally
does not scan unrelated local refs and never automates a force-push.

Bootstrap crates.io with the temporary token, then remove the token from the
shell and revoke it on crates.io:

```sh
CARGO_REGISTRY_TOKEN=... \
  ./scripts/bump-and-publish 0.1.0 \
  --publish-crates --confirm-publish v0.1.0
unset CARGO_REGISTRY_TOKEN
```

After adding the five crates.io Trusted Publishers, atomically push `main` and
the tag:

```sh
./scripts/bump-and-publish 0.1.0 --push
```

Do not push a release tag directly. The helper verifies the canonical remote,
clean history, protected release settings, branch ancestry, and annotated tag
before it sends `main` and the tag together. CI independently rejects a
lightweight tag or a tag whose commit differs from the workflow commit.

The tag workflow publishes the immutable GitHub release. Its crates.io job sees
the already-published v0.1.0 crates and performs no upload. After the entire tag
workflow succeeds, dispatch the top-level Python workflow and approve its
`pypi` environment deployment:

```sh
./scripts/bump-and-publish 0.1.0 --dispatch-python
```

It accepts only a stable, non-draft, non-prerelease, immutable release whose tag,
commit, and package versions agree, then verifies and publishes the exact wheel
and sdist attached to that release. A retry skips a complete matching PyPI
version; after an interrupted partial upload, it publishes only the missing
matching file. Any unexpected filename or checksum fails closed. Do not manually
upload a second build of any artifact.

## Later releases

Add a dated `## [VERSION] - YYYY-MM-DD` section to `CHANGELOG.md` and commit it.
Then run:

```sh
./scripts/bump-and-publish 0.2.0 --dry-run --check-history
./scripts/bump-and-publish 0.2.0 --prepare
./scripts/bump-and-publish 0.2.0 --push
./scripts/bump-and-publish 0.2.0 --dispatch-python
```

`--prepare` updates the Rust manifests and lockfile, Python package and runtime,
documentation manifest and lockfile, and local conda versions together, runs
the release checks, commits the bump, and creates an
annotated local tag. `--push` accepts only a clean `main` containing
`origin/main`, and requires `origin` to be the canonical GitHub repository over
HTTPS or SSH; it then pushes the branch and tag atomically. crates.io and PyPI
run in separate protected jobs and workflows, so a credential for one registry
is never available to the other.

## Local conda package

The local recipe builds the current checkout on Linux or macOS:

```sh
conda build -c conda-forge packaging/conda/local
```

`GRAVLAX_VERSION` may be set when validating a future version. The recipe
defaults to the workspace's current `0.1.0` version; the distribution tests
fail when those versions diverge.

## Bioconda submission template

`bioconda/meta.yaml.in` is deliberately not a claim that Gravlax is published
in Bioconda. It has conspicuous tokens and cannot be submitted as-is. The
release workflow vendors the locked Cargo dependencies and checks that Cargo
can resolve them offline before it creates the source asset. After the GitHub
release has produced that immutable asset, render a recipe whose checksum is
computed from those exact bytes:

```sh
python packaging/bioconda/render_recipe.py \
  --version 0.1.0 \
  --source-archive dist/gravlax-0.1.0-source.tar.gz \
  --output /tmp/bioconda-recipes/recipes/gravlax/meta.yaml
conda render -c conda-forge -c bioconda /tmp/bioconda-recipes/recipes/gravlax
bioconda-utils lint --packages gravlax /tmp/bioconda-recipes
```

Review Bioconda's current contributor policy before opening a submission. A
successful local render or lint is not publication.
