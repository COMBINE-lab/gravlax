---
title: aie explore
description: Resolve biological intent and build inspectable plans in a local, read-only UI.
---

The Gravlax Explorer is a local scientific plan builder and artifact browser:

```sh
cd experiments/pbmc
aie explore
```

Open the address printed by the command, normally
`http://127.0.0.1:8787/`. Select a different project or port with:

```sh
aie explore --project ../other-project --port 8788
```

Explorer first shows every named project resource and whether it resolves. An
annotation resource displays its registered assembly and release; an archive
can display its declared coordinate assembly. These labels are scientific
compatibility claims, while content identities are computed from the files.

The plan builder supports six common evidence questions:

- a full-archive comparison of two annotation releases;
- regional molecule evidence for a gene, transcript, or exon;
- coordinate-defined splice-event usage in a resolved feature;
- support for one exact splice junction;
- inclusion/exclusion usage for two junction sets;
- strand-aware terminal-boundary (APA) evidence in a resolved feature.

Gene symbols and stable gene, transcript, and exon identifiers are resolved by
the same `IntentResolver` path used by `aie plan check`. Resolution is exact and
case-sensitive. Ambiguous identifiers fail with a request for a typed prefix
such as `gene:`, `transcript:`, or `exon:`. A successful result shows the exact
0-based half-open coordinates, strand, assembly, annotation release, and BLAKE3
identity of the annotation snapshot.

Choose a named archive and, where supported, a named cells or groups resource.
The backend then constructs one typed plan and compiles it with the normal plan
resolver. The UI displays the evidence route, scope and exclusions, output
schema, resource identities, exact known selected-input sizes, and any available
I/O bounds from that resolution. Four synchronized exports are available:

- versioned plan YAML;
- the identical plan as JSON;
- the resolved `aie` command with shell-safe arguments;
- a Python/notebook snippet using `gravlax.Client` with the same argument list.

Explorer's result-producing plans explicitly select uniform JSON and a
project-contained `results/explorer-<view>.json` destination. The YAML and JSON
therefore preserve `uniform_output`, while the resolved CLI and Python exports
contain equivalent `--format json --output <final-project-path>` arguments that
can be run directly. The detailed resolved step retains the private staging
path used transactionally by `aie plan run`, but that path is never copied into
a standalone command or Python export. Annotation comparison already uses its
native typed output contract and receives the same explicit result destination.
Explorer never creates that file; it only previews the checked plan.

## Preview an annotation comparison

Choose **Compare annotation releases**, then select the before (A) and after
(B) annotation resources. Both annotations must have registered assembly and
release labels, and the normal plan resolver verifies their content digests,
shared assembly, and compatibility with the selected archive. The preview
shows those A/B labels, exact digests, and verified or unverified compatibility
notes before presenting any export.

The gene-key, solo-feature strand, witness limits, and intentional identical
A/A-control option map directly to the source plan v1 `compare-annotations`
step. The exported JSON result contract includes the outer comparison schema
and the count-delta, class-transition, contributing-cause, and witness table
schemas.

Annotation comparison is deliberately a full-archive operation: Explorer does
not offer locus, cell, group, or output filters for this view. Final signed
B-minus-A counts come from two complete independent assignment and UMI-collapse
reductions. Class transitions, non-exclusive causes, and bounded molecule
witnesses explain changed states; they are not additive or uniquely
attributable components of the nonlinear final-count delta. Explorer only
previews this work and exports its exact plan—it never starts the scan.

Explorer never runs any of these exports and never saves them into the project.
Copy the plan text into a file only after reviewing it, then use `aie plan check`
or `aie plan run` explicitly. For example:

```sh
aie plan check plans/explorer-region.yaml --project . --explain
aie plan run plans/explorer-region.yaml --project .
```

Existing source plans, content-addressed resolved plans, and results remain
available in the artifact browser. Selecting one previews text locally;
**Download exact file** returns the stored bytes rather than a reconstruction.

## Project metadata required for name resolution

Register an annotation with explicit scientific identity:

```sh
aie project add genes references/gencode.v49.annotation.gtf \
  --kind annotation \
  --assembly GRCh38.p14 \
  --annotation-label "GENCODE 49"
```

Coordinate-bearing archives should also carry a declared assembly when known.
Plan checking rejects a declared mismatch. If an older archive has no declared
assembly, Explorer reports compatibility as unverified rather than inferring an
assembly from chromosome names alone.

Terminal-boundary plans use the resolved feature strand by default. They expose
aligned molecule boundaries as evidence locations, not asserted cleavage sites.
Without a registered genome input, the exported plan does not filter possible
internal priming. Cell-list selection is supported by region, event, junction,
and junction-set views; the current APA query supports all cells or a groups
resource.

## Security and scope

Explorer is intentionally narrow:

- it binds to `127.0.0.1` only, with no option to widen the listening address;
- it accepts only loopback `Host` headers;
- it exposes only known project, plan, resolved-plan, and result locations;
- it rejects path traversal and symlinks that escape the project;
- it accepts only read-only `GET` and `HEAD` requests;
- it bounds request headers, query size, field count, and individual values;
- it resolves identifiers and compiles previews in memory, without temporary
  plan files;
- it cannot execute a plan or alter any artifact.

The browser renders all project and identifier text through safe text nodes; it
does not interpret resource names, identifiers, errors, or plan text as HTML.

The UI is embedded in the `aie` binary. It does not require Node.js, a package
manager, a hosted service, or an internet connection. Press `Ctrl-C` in the
terminal to stop it.

When Gravlax runs on a remote machine, keep Explorer bound there and use an SSH
loopback tunnel rather than exposing a server port:

```sh
ssh -L 8787:127.0.0.1:8787 research-host
# on research-host:
aie explore --port 8787
```

Then open `http://127.0.0.1:8787/` on the local computer.
