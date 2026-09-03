# Gravlax demo project

This tiny project exercises the project/plan workflow without downloading a
dataset. Its one-step plan compiles a two-transcript demonstration annotation
into `results/demo.aic`.

It is intentionally a workflow demonstration, not an archive or scientific
benchmark. The repository's
[`workflow and interfaces` guide](../../docs/src/content/docs/workflow.md)
explains how the same plan, resume, Explorer, and Python contracts scale to
registered archives and annotations.

From this directory:

```sh
aie doctor --project .
aie project show
aie plan check plans/compile-demo.yaml --explain
aie plan run plans/compile-demo.yaml --dry-run
aie plan run plans/compile-demo.yaml
aie explore
```

The real run first writes an exact resolved-plan snapshot beneath
`.aie/resolved-plans/`, then creates `results/demo.aic`. The Explorer is
read-only and shows the source plan, snapshot, and result as separate exact
artifacts.

Repeat the completed run with `aie plan run plans/compile-demo.yaml --resume`;
the step is skipped only after its exact completion and output digest are
verified. For an intentional recomputation, remove only `results/demo.aic`
and run with `--resume` again. Keep `.aie/`: it contains the immutable resolved
snapshot and completion provenance.
