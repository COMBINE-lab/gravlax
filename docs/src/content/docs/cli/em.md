---
title: aie dev em
description: Cross-cell EM multimapper recovery from the stored paralog evidence.
---

EM multimapper recovery over the archive's paralog-pattern evidence. Two modes
share one command:

- **Evaluation** (default): a masked-evidence protocol scores recovery
  accuracy on your own data — classes whose unique reads pin their gene are
  masked down to their multimapper evidence and must be recovered — comparing
  pooled cross-cell EM, per-cell EM, and a uniform baseline.
- **Emission** (`--mask 0 --emit`): run EM on the real (unmasked) evidence and
  write the additive fractional recovered-counts layer.

## Usage

```sh
aie dev em [OPTIONS] --gtf <GTF> <ARCHIVE>
```

```sh
# score recovery accuracy (masked evaluation):
aie dev em sample.aie --gtf gencode.v49.gtf

# emit the recovered-counts layer:
aie dev em sample.aie --gtf gencode.v49.gtf \
  --mask 0 --emit em-layer/ --barcodes barcodes.tsv
```

## Options

| Option | Default | Description |
|---|---|---|
| `--gtf <GTF>` | required | Annotation defining the gene candidates |
| `--mask <FRAC>` | `0.2` | Fraction of mixed classes to mask for the labeled evaluation; `0` switches to emission |
| `--seed <SEED>` | `7` | Masking RNG seed |
| `--alpha <ALPHA>` | `20` | Blend-mode global-prior weight |
| `--groups <TSV>` | — | Two-column barcode/group map; adds group and hierarchical evaluation modes and fixes the scored cell set |
| `--group-alpha <ALPHA>` | `20` | Hierarchical group-prior pseudo-count mass |
| `--global-alpha <ALPHA>` | `5` | Hierarchical whole-sample pseudo-count mass |
| `--convex-cell-weight <W>` | `0.10` | Candidate-normalized convex weight on the target-cell distribution |
| `--convex-group-weight <W>` | `0.45` | Candidate-normalized convex weight on the leave-one-cell-out group distribution; the sample weight is `1 - cell - group` |
| `--convex-group-prior <K>` | `20` | Candidate-level sample pseudo-count mass shrinking the group distribution toward the sample distribution |
| `--convex-only` | off | With `--groups`, evaluate only the convex mode; useful for parameter grids and incompatible with `--emit` |
| `--dirichlet-cell-prior <K>` | `16` | Candidate-level group-posterior mass borrowed by the target cell in the posterior-mean Dirichlet proxy |
| `--dirichlet-group-prior <K>` | `20` | Candidate-level sample-posterior mass borrowed by the leave-one-cell-out group in the posterior-mean Dirichlet proxy |
| `--dirichlet-only` | off | With `--groups`, evaluate only the posterior-mean Dirichlet proxy; incompatible with `--convex-only` and `--emit` |
| `--hybrid-depth-scale <D>` | `8` | Half-transition depth for the monotone fixed-convex/Dirichlet-proxy evaluator |
| `--hybrid-depth-power <P>` | `8` | Positive Hill power controlling the transition sharpness |
| `--hybrid-only` | off | With `--groups`, evaluate only the monotone depth hybrid; incompatible with the other diagnostic-only modes and `--emit` |
| `--collapse-groups` | off | Assign all archive cells to one group while scoring only `--groups` barcodes (the pooled-invariance control) |
| `--metrics-json <JSON>` | — | Write top-1, truth probability, negative log loss, multiclass Brier, calibration counts, and fixed evidence-depth strata |
| `--candidate-genes-out <TXT>` | — | Write genes occurring in masked evaluation candidate sets |
| `--candidate-genes-only` | off | Stop after writing `--candidate-genes-out` |
| `--emit <DIR>` | — | With `--mask 0`: write `em.mtx` (real-valued, additive) into this directory; requires `--barcodes` |
| `--barcodes <BARCODES>` | — | Barcode list defining emitted column order |
| `--star` | off | Use the STARsolo-compatible `--soloMultiMappers EM` design (per-cell, intersection candidate sets, STAR's init/zeroing/convergence) and emit `UniqueAndMult-EM.mtx` into `--emit` |
| `--eager` | off | Use the historical full-materialization implementation as a semantic/performance reference |
| `--plot <SVG/PNG>` | — | With a masked run, write a per-mode reliability diagram |

## Sharing models

The masked evaluation reports four base models over identical candidate sets:

- **uniform** assigns equal responsibility to every candidate;
- **cell** estimates expression only from the target cell;
- **pooled** estimates one sample-wide expression vector;
- **blend** adds `--alpha` pseudo-counts distributed according to the
  sample-wide expression vector to each cell's local evidence.

With `--groups`, it also reports a group-only model, an additive hierarchy,
and a candidate-normalized convex model. For a target
candidate set `C`, the convex model forms cell and sample probability vectors
by normalizing their current gene abundances within `C`. It subtracts the
target cell from its group's candidate counts and shrinks that leave-one-cell-
out vector toward the sample vector with `--convex-group-prior`. The final
responsibility is

```text
q(g | c,C) = w_cell p_cell(g | C)
           + w_group p_group,-c(g | C)
           + (1 - w_cell - w_group) p_sample(g | C).
```

Weights must be finite, non-negative, and sum to at most one. Empty cell or
leave-one-cell-out group components fall back to the sample distribution;
cells absent from the group map transfer the group weight to the sample. This
makes pooling strength comparable between common and rare candidate sets and
prevents whole-transcriptome pseudo-count dilution. `--convex-only` skips the
eight other modes during grids. The emitted production layer still uses
`pooled`; supplying groups or experimental parameters does not change emission.

The `dirichlet-proxy` mode is an approximation to a full
hierarchical model. It first obtains the same leave-one-cell-out group
posterior and then forms the cell posterior mean

```text
p_group(g | C) = (n_group,-cell(g) + k_group p_sample(g | C))
                 / (N_group,-cell + k_group)
p_cell(g | C)  = (n_cell(g) + k_cell p_group(g | C))
                 / (N_cell + k_cell).
```

Consequently, cell and group weights change automatically with candidate-set
evidence depth. It provides evidence-dependent pooling without latent
concentrations or variational inference. It is not
the full hierarchical Dirichlet model. Machine-readable metrics stratify every
mode by the mode-independent fitted unique mass over the target candidates:
`0-1`, `(1,4]`, `(4,16]`, and `16+`.

The `depth-hybrid` mode combines the convex and proxy predictions using
that same pre-mode evidence depth `d`:

```text
s(d) = d^P / (d^P + D^P)
q_hybrid = (1 - s(d)) q_convex + s(d) q_proxy.
```

`D` is `--hybrid-depth-scale` and `P` is `--hybrid-depth-power`. The mixing weight is zero at no evidence,
one-half at `d=D`, monotone, and approaches one without a discontinuous
threshold. Its candidate-independent weight preserves normalization. Both
constituents evolve from the hybrid mode's shared EM state; the mixing weight
is fixed before mode fitting and therefore cannot feed back through the model.
The defaults are `D=8, P=8`. This mode produces evaluation scores only;
`--emit` always uses the pooled model.

## What pooling buys

Pooled EM replaces the per-cell rate estimate with the sample-wide one —
possible in a single pass only because the index holds every cell's
equivalence-class evidence at once. In masked evaluation it attains
95.7–98.3% top-1 accuracy versus 90.0–94.4% for the per-cell design, with the
largest gains on the sparsest cells (see
[Capabilities](/gravlax/capabilities/#cross-cell-em-multimapper-recovery)).

The default path does not materialize the archive or allocate a candidate
vector per target. It streams into 8-byte support records, builds 64
deterministic cell shards, and stores candidates as flat `u32` labels with
`u64` CSR offsets. Large inputs spill exact batch-compacted support records to
temporary shard streams and finalize one shard at a time; small inputs retain
the in-memory path. On an evaluated 10,000-cell dataset this reduced peak RSS from
7,653,508 to 3,846,008 KiB at 35.86 s, with byte-identical metrics and stdout.
Temporary shards are removed on success or error.

## The emitted layer

`em.mtx` is an **additive, opt-in layer** of real-valued recovered counts: the
base replay matrices are never modified. Responsibilities are calibrated, so
they can be consumed as probabilities; thresholding at responsibility > 0.8
keeps the layer's high-confidence core. Use `--star` to request the
STARsolo-compatible update scheme; cross-tool byte identity has not yet been
established.

The built-in masked analysis is useful but does not prove that pooled or
group-aware sharing is unbiased: it defines truth using mixed classes whose
unique evidence is then hidden. Group-aware results should therefore be read
together with size-preserving label shuffles and the one-group comparison.
Across evaluated PBMC datasets, the benefit of group-aware estimators varied
with evidence depth and much of the apparent improvement could be reproduced
after shuffling group labels. The `convex`, `dirichlet-only`, and `hybrid-only`
modes are consequently diagnostic model comparisons, not count-emission
choices. `--emit` remains pooled.
