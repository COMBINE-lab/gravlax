---
title: aie doctor
description: Check an installation, project, and selected Gravlax artifacts.
---

`aie doctor` is the quickest way to find setup problems before starting a long
analysis. It reports each check separately and gives an actionable next step
for warnings and failures.

```sh
aie doctor
```

The default diagnosis checks the running executable, available compute,
project discovery, workspace write access, and the optional alignment tools on
`PATH`. When a project is found, every internal and explicitly external named
resource is resolved and checked for availability. STAR and samtools are only
required on machines that perform their respective alignment or BAM-inspection
steps; their absence does not prevent
archive replay or queries.

The workspace test creates and immediately removes one empty, uniquely named
probe file. It does not change a project manifest, plan, archive, or result.

## Validate inputs

Pass one or more paths to add format-specific checks:

```sh
aie doctor sample.aie annotation.aic
```

For `.aie` files, the normal check authenticates the archive directory/root
and reads the small metadata and chromosome sections through the production
archive reader. Payloads remain lazy, just as they do during ordinary queries.
Use a full check when copying or auditing an archive:

```sh
aie doctor sample.aie --verify-content
```

This reads and verifies every compressed payload. A compiled `.aic`
annotation is fully decoded and structurally validated.

## Projects and automation

Project discovery normally walks upward from the current directory. Select a
different project explicitly with either its directory or manifest:

```sh
aie doctor --project experiments/pbmc
aie doctor --project experiments/pbmc/aie-project.yaml
```

For scripts and support reports, request the stable
`gravlax.doctor.v1` JSON document:

```sh
aie doctor sample.aie --json
```

Failures return a non-zero exit status. Warnings are advisory by default; add
`--strict` when a warning should also fail a CI or deployment check.
