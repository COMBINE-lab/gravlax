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

The manifest should itself live at an immutable URL, such as a versioned GitHub release asset or a
versioned data-repository record. The archive roots reported by Gravlax remain the scientific
content identities; transport SHA-256 protects downloads before the archive reader opens them.

Do not publish the notebooks as one-click demonstrations until all required fields in the template
manifest are populated with real assets and its digest has been recorded in the notebook links.
