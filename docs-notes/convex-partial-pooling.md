# Candidate-normalized convex partial pooling

Status: implemented as an experimental masked-evaluation mode; production emission remains
pooled.

For target class `k` in cell `c`, let `C_k` be its candidate genes and let the current EM
abundances be `pi_cg`, `Pi_hg`, and `Pi_g`. The sample component is

```text
p0(g | C_k) = (Pi_g + epsilon) / sum[j in C_k](Pi_j + epsilon).
```

The cell component is `pi_cg / sum(C_k)` when that denominator is nonzero and otherwise `p0`.
To prevent a target cell from re-entering through its group prior, define

```text
n_hg,-c = max(Pi_hg - pi_cg, 0)
ph(g | C_k) = (n_hg,-c + kappa p0(g | C_k))
              / (sum[j in C_k] n_hj,-c + kappa).
```

An empty leave-one-cell-out group falls back to `p0`. Responsibilities are the convex mixture

```text
q(g | c,C_k) = lambda_c pc(g | C_k)
             + lambda_h ph(g | C_k)
             + (1 - lambda_c - lambda_h) p0(g | C_k).
```

The CLI rejects weights outside the simplex rather than silently renormalizing them. A cell absent
from the group map transfers `lambda_h` to the sample component. At `(lambda_c,lambda_h)=(0,0)`,
the mode reproduces pooled candidate probabilities exactly. `kappa` is measured in effective
candidate-set counts, so it retains the same interpretation for common and rare genes.

The packed implementation retains the existing 64 deterministic shards and flat CSR candidate
labels. It makes one pass over each target to obtain the three candidate denominators and one pass
to update responsibilities; it stores no per-target probability vector. `--convex-only` avoids
recomputing reference modes during a grid.

## Validation boundary

The model is not a production recommendation until its weights and `kappa` are selected on
development data and confirmed on datasets not used in the additive scale audit. The comparison
baseline is the best candidate-normalized no-group arm (`lambda_h=0`), not untuned pooled EM.
Group labels must retain the existing candidate-gene exclusion, real-versus-size-preserving-
shuffle control, and stable-partition audit. Proper scores are primary; top-1 is a guardrail.

## Reserved full hierarchical model

If fixed convex weights leave repeatable group-specific calibration errors while real labels
still beat shuffled labels, the next model is a hierarchical Dirichlet construction:

```text
theta_0 ~ Dirichlet(eta)
theta_h ~ Dirichlet(kappa_0 theta_0)
theta_c ~ Dirichlet(kappa_h theta_h)
P(class k | theta_c) = sum[g in C_k] theta_cg.
```

Concentrations would be learned by empirical Bayes or variational inference, with packed
equivalence classes as sufficient optimization units. Do not implement it merely because it is
more formal. Escalate only if the convex model demonstrates reproducible biological group signal
but fails because one global mixture cannot adapt to group depth, candidate class, or uncertainty.
Stop instead if untouched datasets show less than a meaningful increment over the tuned no-group
arm or if shuffled labels explain the gain.

## Posterior-mean screening proxy

Before implementing the full model, `dirichlet-proxy` tests its distinctive adaptive-pooling
mechanism with two conjugate posterior means:

```text
ph(g | C_k) = (n_hg,-c + kappa_0 p0(g | C_k)) / (N_h,-c + kappa_0)
pc(g | C_k) = (n_cg + kappa_h ph(g | C_k)) / (N_c + kappa_h).
```

Unlike the fixed convex model, the effective cell weight is `N_c/(N_c+kappa_h)`, and the group
weight similarly depends on leave-one-cell-out group depth. The packed evaluator reuses the same
two target passes and flat counts. It reports proper scores in four strata defined before mode
fitting by unique candidate mass: `0-1`, `(1,4]`, `(4,16]`, and `16+`. This is a falsification
gate: failure means that learned latent concentrations have no demonstrated mechanism to rescue;
success only licenses a full variational or empirical-Bayes implementation and untouched testing.

### D0 screening result

The locked 25-configuration D0 screen selected `kappa_h=64` and `kappa_0=80` over mask seeds 7,
17, and 29. Against the previously selected fixed convex model, the posterior-mean proxy lowered
mean negative log loss from 0.058604 to 0.058107 (0.848%), lowered multiclass Brier from 0.027329
to 0.026780, and increased top-1 by 0.053 percentage point. Real labels beat size-preserving
shuffled labels on every seed; the shuffled control explained 19.9% of the proxy's gain over its
matched fixed-convex control. These aggregate results are encouraging evidence that both
adaptive weighting and real group structure contain signal.

The predeclared depth gate nevertheless failed. Relative to fixed convex, proxy loss was 0.625%,
1.035%, and 0.532% worse in the `0-1`, `(1,4]`, and `(4,16]` strata, respectively, while it was
3.405% better in `16+`. The result is therefore **marginal**, not a license to implement the full
hierarchical Dirichlet model: one concentration pair improved the dominant high-evidence
population by relying more on cell evidence, but its stronger shrinkage was inferior for shallow
targets. Production remains pooled. A future escalation needs an independently validated model
that can retain the fixed convex behavior at low evidence while adapting at high evidence; model
formality alone is not sufficient.

## Monotone depth-gated hybrid

The next diagnostic retains both already selected D0 constituents and changes only their mixing
rule. Let `d` be the mode-independent initial fitted unique mass over a target's candidate set and
let `D>0` be a transition scale:

```text
s(d;D) = d^2 / (d^2 + D^2)
q_hybrid = (1-s) q_fixed-convex + s q_Dirichlet-proxy.
```

The fixed component uses `(lambda_c,lambda_h,kappa)=(0.2,0.6,80)` and the proxy uses
`(kappa_h,kappa_0)=(64,80)`. The squared Hill gate is monotone, equals 0.5 at `d=D`, remains close
to fixed convex for shallow targets, and smoothly approaches the proxy for deep targets. The gate
is candidate-independent within a target, so the result remains normalized. Both constituent
probabilities are recomputed from the hybrid mode's shared current EM state; only `d` is frozen
before fitting. `--hybrid-only` avoids running the eight reference modes during a scale grid.

### D0 hybrid result

The locked ten-scale grid selected the interior value `D=64`. Mean negative log loss was
0.057808, 0.513% below the selected proxy and 1.357% below fixed convex; multiclass Brier also
improved to 0.026720, while top-1 was only 0.006 percentage point below the better constituent.
The hybrid was within 0.0041% of fixed convex in every shallow stratum and improved the proxy's
`16+` loss by 0.240%. It recovered 120.2% of the post-hoc hard-stratum oracle gain because smooth
per-target mixing and a shared iterative state slightly outperformed choosing one fitted
constituent per coarse depth bin.

Seven of eight locked criteria passed. Real-group hybrid beat shuffled-group hybrid on every
seed, but the hybrid's increment over the proxy was 0.513% with real groups and 0.684% with
shuffled groups. The resulting predeclared shuffle-explained ratio was 1.343, above its 0.50
limit, so the formal verdict is **marginal**, not PASS. A descriptive factorial view separates
the two effects: real groups still improved hybrid loss by 0.654%, while the depth gate improved
loss under both real and shuffled labels. That makes the failed ratio a poor measure of the
depth mechanism, but it remains binding because it was locked in advance. Any D4 confirmation
therefore needs a separately frozen factorial hypothesis rather than silently deleting this
control after seeing D0. Production remains pooled.
