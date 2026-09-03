# Publishing Gravlax through Bioconda

The Bioconda package installs the `aie` command-line program and shell
completions. The separately versioned `gravlax-client` Python package is
distributed through PyPI rather than bundled into this Conda package.

The files in this directory prepare a recipe; they do not publish anything.
Submit only a rendered `meta.yaml`, never `meta.yaml.in`.

## Prerequisites

Wait until the public, non-draft GitHub release exists. For version `0.1.0`, it
must contain these two files:

- `gravlax-0.1.0-source.tar.gz`
- `SHA256SUMS`

The source archive contains the tagged repository, `Cargo.lock`, and vendored
Cargo dependencies. The rendered recipe pins the SHA-256 digest of those exact
bytes, and its build does not access the network.

Download and verify the published files in a new directory:

```sh
mkdir -p /tmp/gravlax-0.1.0-release
gh release download v0.1.0 \
  --repo COMBINE-lab/gravlax \
  --pattern 'gravlax-0.1.0-source.tar.gz' \
  --pattern SHA256SUMS \
  --dir /tmp/gravlax-0.1.0-release
cd /tmp/gravlax-0.1.0-release
sha256sum --check --ignore-missing SHA256SUMS
```

The checksum command must report
`gravlax-0.1.0-source.tar.gz: OK` before continuing.

## Create the submission

Fork `bioconda/bioconda-recipes`, then work from its current `master` branch:

```sh
git clone git@github.com:YOUR-GITHUB-ACCOUNT/bioconda-recipes.git
cd bioconda-recipes
git remote add upstream https://github.com/bioconda/bioconda-recipes.git
git fetch upstream
git switch -c add-gravlax-0.1.0 upstream/master
```

From that checkout, render the recipe using the script in the tagged Gravlax
source tree:

```sh
python /path/to/gravlax/packaging/bioconda/render_recipe.py \
  --version 0.1.0 \
  --source-archive /tmp/gravlax-0.1.0-release/gravlax-0.1.0-source.tar.gz \
  --output recipes/gravlax/meta.yaml
```

The renderer refuses to overwrite an existing recipe. It creates
`recipes/gravlax/meta.yaml` with the release asset's SHA-256 digest.

Inspect both files and confirm that the source URL downloads successfully
without GitHub authentication.

## Test the recipe

Bioconda recommends testing from the root of the `bioconda-recipes` checkout
with its released tooling:

```sh
mamba create -n bioconda-review -c conda-forge -c bioconda \
  --strict-channel-priority bioconda-utils
conda activate bioconda-review
bioconda-utils lint recipes/ config.yml --packages gravlax
bioconda-utils build recipes/ config.yml \
  --packages gravlax \
  --docker \
  --mulled-build-and-test \
  --force
```

The build compiles from the vendored dependency tree, installs `aie`, strips
the executable, generates Bash, Zsh, and Fish completions, and records
third-party licenses. The tests run `aie --version`, `aie --help`, `aie doctor
--help`, completion generation, and checks for all three installed completion
files. The mulled test repeats the commands in a minimal runtime container.

The recipe requests the native `rust >=1.98` package and documents a narrow
`should_use_compilers` lint exception. At the time of the 0.1.0 release,
Conda-forge's Rust 1.98 package is available on its main channel, while the
corresponding compiler activation package on that channel is still Rust 1.97.
Remove the exception and switch to `{{ compiler('rust') }}` once the activation
package reaches Rust 1.98 or newer.

The recipe also documents a `missing_run_exports` exception. Gravlax installs a
standalone executable rather than a shared library or other ABI consumed by
downstream builds, so propagating a version pin into downstream environments
would add constraints without protecting binary compatibility.

When lint and build both succeed, commit the recipe, push the branch to your
fork, and open a pull request against `bioconda/bioconda-recipes:master`. A new
upstream version starts at build number `0`; recipe-only corrections for the
same upstream version increment the build number.

Current Bioconda references:

- <https://bioconda.github.io/contributor/workflow.html>
- <https://bioconda.github.io/contributor/guidelines.html>
- <https://bioconda.github.io/contributor/building-locally.html>
- <https://bioconda.github.io/contributor/linting.html>
