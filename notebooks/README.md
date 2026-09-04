# Gravlax preprint demonstrations

These notebooks are executable Google Colab entry points for three claims:

1. reinterpret one evidence archive with two annotations;
2. discover recurrent molecule-supported events across donors; and
3. query one junction across a collection, then inspect same-record evidence
   predicates in one source archive.

They intentionally contain no guessed data URLs and no saved scientific output. Before running a
notebook, publish a demo capsule, instantiate `demo-manifest.template.json`, and set the first cell's
`MANIFEST_URL` and `MANIFEST_SHA256`. Every downloaded byte (including the CLI and Python wheel) is
checked against SHA-256 before execution. Source archives are downloaded individually and each
notebook builds its path-bound `.aicollection` locally; a downloaded collection could not retain
valid source paths in Colab. Missing locators, non-HTTPS URLs, malformed digests,
and identity mismatches stop the notebook.
The installed `aie` executable and Python package must both report the exact version declared by
the manifest before any scientific command runs.

Each notebook ends with a small deterministic SVG visualization constructed directly from its
typed live result tables. The plotting code uses the Python standard library plus the IPython
display API already provided by the Colab notebook runtime; it installs no mutable plotting
dependency and carries no cached values.

The event-discovery notebook consumes the frozen `find-events` tables
`capabilities`, `entities`, `components`, `counts`, `terminal_anchors`, and
`terminal_counts`. Its event support is counted from uniquely mapped chains and
exact raw-UMI-value classes; one-mismatch edges and multimapper alternatives are
not folded into those counts. The co-occurrence notebook displays observed and completeness masks as well as
the three-valued expression result; `unknown` is never silently converted to absence.
When event discovery supplies an annotation, the story also supplies its
assembly and release label (and may pin its observed BLAKE3 digest). Those
values are recorded caller declarations; exact contig names are checked, but a
matching assembly is not inferred from names alone.

The v1 manifest and data assets must live under the exact
`https://github.com/COMBINE-lab/gravlax/releases/download/<fixed-tag>/...` URL form. The archive
roots reported by Gravlax remain the scientific content identities; transport SHA-256 protects
downloads before the archive reader opens them.

Do not publish the notebooks as one-click demonstrations until all required fields in the template
manifest are populated with real assets and its digest has been recorded in the notebook links.

## Prepare the immutable demo capsule

The capsule is built in two phases. The first phase derives the small public data assets from
hash-pinned inputs. The second phase binds those assets to official Gravlax release artifacts and
immutable HTTPS URLs. Both commands refuse to replace an existing output directory.

The intended first capsule is a small locus-restricted resource made from all eight independent
public SEZ donors. Its validated 0-based half-open windows are LRRC7
`chr1:70122500-70142500`, FNBP1 `chr9:129907500-129926000`, and the terminal-tail/junction
example `chr12:18714000-18854000`. The authoritative candidate build retained 163,659 selected
alignment records across the eight donors; the exact count is recorded by every build. The builder
remains generic over any
non-overlapping set of windows, so these coordinates belong in the hash-pinned build specification
rather than in notebook code.

Each source declaration in the build specification has this form:

```json
{
  "path": "path-below-source-root",
  "bytes": 123,
  "sha256": "64-lowercase-hex-digits",
  "provenance": {
    "accession": "public accession or release",
    "source_url": "https://stable-source-page"
  }
}
```

`sha256` and `bytes` identify the exact local object consumed by the build. For a derived BAM,
whitelist, group map, junction catalogue, decompressed annotation, or log, `source_url` identifies
the stable public sample/study or upstream release rather than claiming that the URL serves those
exact derived bytes; the remaining provenance fields must describe that derivation or role.

The top-level build specification uses schema `gravlax.demo-capsule-build.v1` and supplies:

- a public `capsule_id` and one-line `summary` (for this capsule: “A curated
  locus-restricted demonstration from eight independent public SEZ donors.”);
- the exact future release `software.version` and full 40-character source revision;
- the assembly, reference FASTA declaration, and curated 0-based half-open `windows`;
- a path-free alignment declaration with STAR version, normalized public command description,
  `junction_discovery` set to `per-library-two-pass`, index identity, and chemistry;
- one donor entry per sample, with public sample/donor IDs, archive resource and filename,
  accessions, and declarations for the tagged BAM, whitelist, headerless barcode/group map,
  exact pass-1 `_STARpass1/SJ.out.tab`, and `Log.out` (the resolved defaults and
  effective command);
- `before` and `after` GTF declarations, each with its output AIC resource, output AIC filename,
  subset-GTF filename, and public annotation label;
- the three story objects defined by `demo-manifest.schema.json`; and
- a notice for every third-party source with its stable source URL, terms, terms URL, and
  citation.

For the generated resources, the event story must name `collection_groups` and `design`; the
drilldown story must name `drilldown_groups`. Both collection stories map every staged sample to
its archive resource. `archive_sample` chooses which donor's headerless group map is emitted for
the same-molecule drilldown.

The candidate's broad cell-group labels were frozen without using FNBP1 counts. Its maps include
only cells already present in the prior FNBP1-targeted archive: A 867/4,275 fully labeled cells, B
594/1,892, C 644/2,888, D 3,199/7,259, E 392/1,159, F 771/2,966, G 3,998/14,213, and H
407/1,974. They cover the FNBP1-support numerators used by the event query but are not full
cell-group denominators. The chr12 mechanics drilldown is evaluated only within donor A's
867-cell labeled subset. These statements are required in `group_map_scope`, preserved in
`BUILD-RECORD.json`, and rendered into the generated README and release notes.

The validated event-discovery headline deliberately searches against the older v32 annotation,
not v49. Configure its story with `annotation` pointing to the v32 AIC, `kinds` set to
`["cassette"]`, `novel_only` set to `true`, and thresholds `min_donors: 4`, `min_samples: 4`,
`min_umi_classes: 8`, `min_side_umi_classes: 2`, and `min_support: 2`. With the eight-donor
capsule, also require the `neuronal` and `astrocyte` groups with
`min_group_umi_classes: 2`. The leading event is FNBP1
`cassette:chr9:-:129908999-129915965,129919242-129923843,129908999-129923843`: 44 informative
UMI classes in 44 cells across six donors and six groups, with no compatible v32 transcript and
ten compatible v49 transcripts. The annotation-reinterpretation story should likewise use v32 as
`annotation_before` and v49 as `annotation_after`. In donor C, lock the LRRC7 gene
`ENSG00000033122` to a minimum signed v49-minus-v32 UMI delta of 2,000 (the validated result is
+2,003, from 750 to 2,753) without remapping. This window lies within the extension of the
ENST00000651989 terminal exon between v32 and the v49 canonical/MANE Select transcript.

The chr12 same-molecule story is a query-mechanics positive control, not a biological discovery:
its raw support is low-complexity oligo-like. Preserve that interpretation in the required
`story_note`, require a true pattern containing `tail`, `exon`, and `splice`, and lock at least 20
selected molecule records. The donor-A group map contains the 867 scoped cells used by this
drilldown. Set `emit_membership` to `true`: this frozen story supplies the bounded
`max_memberships` option, whose CLI contract requires membership emission.
Use `chr12:18715181-18853140` (without a strand suffix) for the federated `junction` field; retain
`:+` in the named co-occurrence predicates, where strand is part of each predicate locus.

Use a candidate binary to tune the specification, but treat those outputs as disposable. Run the
final data build with the `aie` executable extracted from the official cargo-dist archive; the
finalizer requires its exact executable SHA-256, not merely a matching version string:

```sh
python packaging/build_demo_capsule.py \
  --spec /path/to/demo-capsule-build.json \
  --source-root /path/to/hash-pinned-inputs \
  --aie /path/to/aie \
  --samtools /path/to/samtools \
  --output-dir /path/to/gravlax-demo-data-build
```

Before filtering, the builder checks the declared STAR version and `--twopassMode Basic` against
each exact `Log.out`, and requires the resolved `sjdbGTFfile` and `sjdbFileChrStartEnd` values to
both be `-`. This verifies that no GTF or explicit junction list was supplied to these alignment
runs. It treats the normalized public command and the STAR index manifest as caller declarations;
the builder records but does not consume or independently verify the original index components. The
sanitized BAM retains `@HD`/`@SQ` plus a truthful public scope comment and contains no manufactured
STAR `@PG`, so the archive cannot mislabel a rewritten command as original header evidence. The
builder filters complete alignment records to the selected windows, coordinate-sorts the retained
records as required by archive ingest,
and invokes ingest from a private temporary working directory using only public relative
basenames. It records the exact source identities but publishes no BAM, read sequence, base
quality, or read name. Molecule/read-derived records and alternative placements wholly outside the
windows are absent, so multimapper multiplicity and completeness are local to the selected windows
and must not be interpreted genome-wide. Each `.aie` does embed the donor's complete genome-wide
aggregate STAR pass-1 junction catalogue as exact root-bound alignment provenance. That catalogue
contains junction coordinates and aggregate support columns, but no sequences, cell barcodes, or
UMI identities. Do not publish the source BAMs, reference, or STAR logs as separate capsule assets.

Do not finalize against version 0.1.4: those official artifacts predate the archive and query
features exercised here. After the next version has been released through the normal cargo-dist
and PyPI workflows, download—not rebuild—the official cargo-dist archive for the deployment
platform and the exact `gravlax-client` wheel, rebuild the data with the extracted official
executable, and then bind the staged data to those bytes:

```sh
python packaging/finalize_demo_capsule.py \
  --build-dir /path/to/gravlax-demo-data-build \
  --output-dir /path/to/gravlax-demo-data-v1 \
  --data-base-url https://github.com/COMBINE-lab/gravlax/releases/download/demo-data-v1 \
  --aie-asset /path/to/official-gravlax-release.tar.gz \
  --aie-url https://github.com/COMBINE-lab/gravlax/releases/download/vX.Y.Z/official-gravlax-release.tar.gz \
  --aie-member target-triple/aie \
  --python-wheel /path/to/gravlax_client-X.Y.Z-py3-none-any.whl \
  --python-wheel-url https://files.pythonhosted.org/packages/.../gravlax_client-X.Y.Z-py3-none-any.whl \
  --source-repository /path/to/gravlax-checkout
```

The v1 finalizer accepts data and `aie` assets only at fixed-tag
`github.com/COMBINE-lab/gravlax/releases/download/...` URLs and the wheel only at its
`files.pythonhosted.org/packages/...` content URL. Each URL basename must equal the local asset
filename, and the executable and wheel must both report the build-record version. Finalization
also requires the exact `vX.Y.Z` software tag in the supplied source repository and verifies that
the local tag peels to `BUILD-RECORD.json`'s full source revision. Before finalizing, independently
confirm the tag against the canonical remote, then fetch it:

```sh
test "$(git -C /path/to/gravlax-checkout remote get-url origin)" = \
  https://github.com/COMBINE-lab/gravlax.git
gravlax_remote_tag="$(git -C /path/to/gravlax-checkout ls-remote origin \
  refs/tags/vX.Y.Z 'refs/tags/vX.Y.Z^{}')"
gravlax_remote_revision="$(printf '%s\n' "$gravlax_remote_tag" | \
  awk '$2 == "refs/tags/vX.Y.Z^{}" { print $1 }')"
if test -z "$gravlax_remote_revision"; then
  gravlax_remote_revision="$(printf '%s\n' "$gravlax_remote_tag" | \
    awk '$2 == "refs/tags/vX.Y.Z" { print $1 }')"
fi
test "$gravlax_remote_revision" = FULL_SOURCE_REVISION
git -C /path/to/gravlax-checkout fetch origin tag vX.Y.Z
```

This is a required publication gate: for an annotated tag it compares the peeled `^{}` row, and
for a lightweight tag it compares the direct row. The finalizer then repeats the peel check
against the fetched local tag and records that local tag/revision binding. It deliberately makes
no network request, so its check alone is not proof of the canonical remote; retain the successful
remote-gate log with the release record.
The finalizer and standalone verifier are the authoritative v1 URL-policy checks; the JSON Schema
and notebook loader separately enforce HTTPS transport and hash-pinned bytes. Finalization copies
data rather than hard-linking it, writes `demo-manifest.json`, `README.md`,
`RELEASE-NOTES.md`, and `FINALIZATION-RECORD.json`, and writes a `SHA256SUMS` covering every other
file in the flat capsule. It never publishes a path-bound `.aicollection`.

Verify the finalized directory before upload:

```sh
python packaging/verify_demo_capsule.py /path/to/gravlax-demo-data-v1 \
  --aie /path/to/extracted/aie \
  --aie-asset /path/to/official-gravlax-release.tar.gz \
  --python-wheel /path/to/gravlax_client-X.Y.Z-py3-none-any.whl
```

This re-hashes the flat capsule, verifies every archive root and root-bound provenance manifest,
checks both software identities, rebuilds collections locally, and asserts the manifest's frozen
scientific invariants for all three stories.

## Deploy the capsule

Before creating either the next software release or `demo-data-v1`, enable repository-level
**Release immutability** under **Settings > Releases > Enable release immutability**. It applies
only to releases created after it is enabled. Confirm the setting before creating the draft; the
equivalent administrative API calls are:

```sh
gh api --method PUT repos/COMBINE-lab/gravlax/immutable-releases
gh api repos/COMBINE-lab/gravlax/immutable-releases
```

When a draft is published, GitHub release immutability locks both its tag and attached assets and
generates a release attestation. A tag ruleset remains useful defense in depth: add one for
`refs/tags/demo-data-*` that permits the release maintainers to create a tag but prevents later
updates or deletion. See GitHub's
[enablement procedure](https://docs.github.com/en/code-security/how-tos/secure-your-supply-chain/establish-provenance-and-integrity/prevent-release-changes)
and [immutable-release model](https://docs.github.com/en/code-security/concepts/supply-chain-security/immutable-releases).

Create an annotated
non-SemVer `demo-data-v1` tag at the exact feature-release commit recorded by the capsule, verify
both the tag object and its peeled commit locally, and push that tag explicitly:

```sh
git fetch origin
git tag -a demo-data-v1 FULL_SOURCE_REVISION \
  -m "Gravlax immutable demo data v1"
test "$(git cat-file -t demo-data-v1)" = tag
test "$(git rev-parse demo-data-v1^{})" = "$(git rev-parse FULL_SOURCE_REVISION^{commit})"
git push origin refs/tags/demo-data-v1
```

If the hosting settings cannot protect a newly created tag immediately, keep the release draft
private to maintainers during this short creation window, install the rule before publication,
and re-check that the remote tag has not moved. Create the GitHub release only from the verified
existing tag, leave it non-latest, and upload every file from the finalized flat directory without
`--clobber`:

```sh
gh release create demo-data-v1 \
  --repo COMBINE-lab/gravlax \
  --verify-tag \
  --draft \
  --latest=false \
  --title "Gravlax immutable demo data v1" \
  --notes-file /path/to/gravlax-demo-data-v1/RELEASE-NOTES.md

gh release upload demo-data-v1 /path/to/gravlax-demo-data-v1/* \
  --repo COMBINE-lab/gravlax
```

Download the draft data assets into a fresh empty directory. The official `aie` bundle and Python
wheel are not duplicated in the data release: separately download them into another fresh
directory from the exact `software.aie.url` and `software.python_wheel.url` recorded in the
downloaded manifest, verify their manifest SHA-256 values, and extract `aie` from the recorded
member path. Then run the verifier against only these cleanly downloaded bytes. Compare the
`digest` reported for every GitHub demo-release asset with the matching row in the downloaded
`SHA256SUMS`; compare `SHA256SUMS` itself with the local finalized copy. If anything differs,
delete the unpublished draft and use a new capsule tag rather than replacing bytes under a
published URL.

```sh
mkdir -p /path/to/fresh-demo-download /path/to/fresh-software-download
gh release download demo-data-v1 \
  --repo COMBINE-lab/gravlax \
  --dir /path/to/fresh-demo-download
gh api repos/COMBINE-lab/gravlax/releases/tags/demo-data-v1 \
  --jq '.assets[] | [.name, .digest] | @tsv'
curl --fail --location MANIFEST_SOFTWARE_AIE_URL \
  --output /path/to/fresh-software-download/official-gravlax-release.tar.gz
curl --fail --location MANIFEST_PYTHON_WHEEL_URL \
  --output /path/to/fresh-software-download/gravlax_client-X.Y.Z-py3-none-any.whl
python packaging/verify_demo_capsule.py /path/to/fresh-demo-download \
  --aie /path/to/fresh-software-download/extracted/member/path/aie \
  --aie-asset /path/to/fresh-software-download/official-gravlax-release.tar.gz \
  --python-wheel /path/to/fresh-software-download/gravlax_client-X.Y.Z-py3-none-any.whl
```

Once this clean-download verification succeeds, publish the draft without making it the latest
software release:

```sh
gh release edit demo-data-v1 \
  --repo COMBINE-lab/gravlax \
  --draft=false \
  --latest=false
```

Finally, record the SHA-256 of the published `demo-manifest.json`, place its immutable release URL
and digest in the first cell of each notebook, run all three notebooks from a clean Colab runtime,
and commit those locator-only notebook changes. A DOI-backed repository snapshot can mirror the
same byte-identical flat capsule later; changing any byte requires a new capsule version and new
manifest digest.
