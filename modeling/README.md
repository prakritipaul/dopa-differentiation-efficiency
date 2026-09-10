# D11 -> D52 modeling: CV/modeling plan

> **Project-level findings summary with literature context: [`FINDINGS.md`](../FINDINGS.md)**


Predicting D52 differentiation efficiency from D11 features (138 cell
lines, from `metadata_eda/cohort/qualifying_cell_line_pool_min10_per_timepoint.csv`
and `metadata_eda/cohort/d52_diff_efficiency_label.csv`). Designed to generalize
to a D30 -> D52 baseline model later (see Architecture).

## Layout

```
modeling/
  *.py           the package: folds, features, models, harness,
                 run_* entry points, feature_importance
  tests/         pytest suite (30 fast + 1 slow integration test)
  fold_data/     per-fold feature tables and fold assignments
                 (`_full` = baseline basis, `_qualonly` = restricted-fit variant)
  results/       CV grids, nested per-fold selections, feature-importance tables
  docs/          PCA_interpretation.md, feature_importance_plan.md,
                 pool_correction_investigation.md
  plots/         diagnostic figures
  archive/       superseded first-pass artifacts, kept for provenance only --
                 nothing reads these
```

Writers create their own subdirectory, so a fresh clone or a redirected
`out_dir` works without manual `mkdir`.

## First-pass prototype scope (decided while implementing)

Timed a single fold's PCA fit empirically: `compute_gene_stats` 11.2s +
`extract_hvg_matrix` 7.6s + `run_pca` 29.1s (k=10) / 23.5s (k=5) = 48-60s
per fold. The two streaming passes don't depend on `k` at all, and the
PCA step itself only saves ~12% going from k=10 to k=5 (randomized SVD
cost scales mostly with matrix size -- 250K cells x 2000 HVGs -- not much
with target rank in this range). This changed two decisions:

- **`n_pcs` stays at 10, not reduced to 5** -- reducing it barely helps
  runtime (~12% off the PCA step only) so isn't worth losing PCs 6-10 for
  a first pass.
- **Repeats = 5, not 10 or 3** -- 3 is too thin to see a real distribution
  of results (the stated goal); 10 is the better final number eventually,
  but for a first-pass prototype the marginal precision isn't worth the
  extra ~50 min, especially once we noticed LOCO/LODO (below) is the
  actual dominant cost, not repeats.
- **This first pass runs `plain` + `donor_grouped` only; LOCO and LODO are
  deferred**, not dropped. Total folds = `2 x 5 x n_repeats` (plain +
  donor_grouped) `+ 138` (loco) `+ 20` (lodo) -- LOCO alone (138 folds) is
  larger than even the full 10-repeat plain+donor_grouped schemes
  combined (100 folds), so it dominates total cost far more than the
  repeats count does. Validating the harness end-to-end on the cheaper
  schemes first, before paying for LOCO/LODO's ~2.4 hr, is the more
  sensible order for troubleshooting.
- Estimated cost for this first pass: 50 folds (`5 x 5` plain + `5 x 5`
  donor_grouped), k=10, ~48-50s/fold once the process is warm ≈ **~40-45
  min**, run as a single background job (avoids paying Python/import
  startup cost per fold).

## Full results (258 folds: plain 50, donor_grouped 50, LOCO 138, LODO 20)

Headline (donor_grouped, nested CV, **no pool-correction** -- see
"Pool-correction dropped" below):

The pre-registration named a specific model per task (**ridge** /
**logistic_l2**), so those are the headline numbers. The L1 variants
scored better, but quoting them as "the" headline would be exactly the
cherry-picking the pre-registration exists to prevent -- they are
reported below as secondary.

| task | model | key metrics |
|---|---|---|
| regression (**headline**) | ridge | MAE=0.138, RMSE=0.176, **R2=0.653 ± 0.021** |
| regression (secondary) | lasso | MAE=0.135, RMSE=0.172, R2=0.668 ± 0.015 |
| classification (**headline**) | logistic_l2 | **ROC-AUC=0.943 ± 0.008**, PR-AUC=0.965, balanced_acc=0.906 |
| classification (secondary) | logistic_l1 | ROC-AUC=0.951 ± 0.008, PR-AUC=0.970, balanced_acc=0.913 |

D11 features are clearly predictive of D52 differentiation efficiency,
even under the strict donor-grouped (never-seen-donor) test. Going from
the 5-repeat first pass to the full 10-repeat + LOCO/LODO run barely
moved these numbers (secondary lasso R2 0.667->0.668, secondary L1 AUC
0.949->0.951) -- the estimates were already stable.

**All four schemes, nested CV:**

| scheme | regression R2 (lasso / ridge) | classification ROC-AUC (L1 / L2) |
|---|---|---|
| LOCO (leave-one-line-out) | 0.691 / 0.683 | 0.931 / 0.939 |
| plain 5-fold | 0.680 / 0.677 | 0.947 / 0.945 |
| LODO (leave-one-donor-out) | 0.679 / 0.665 | 0.953 / 0.958 |
| donor_grouped 5-fold | 0.668 / 0.653 | 0.951 / 0.943 |

**Donor-leakage effect is real but small (~0.012 R2), and consistent
across two independent comparisons**: plain vs. donor_grouped
(0.680->0.668) and LOCO vs. LODO (0.691->0.679) give the same gap. The
scheme ordering is mechanistically sensible -- LOCO is most optimistic
(137/138 lines in training AND siblings allowed), donor_grouped strictest
(siblings excluded AND ~28 lines held out), LODO in between (prevents
donor leakage but still trains on ~131-137 lines). **Classification is
essentially flat across all four schemes** (0.93-0.96) -- no meaningful
donor-leakage penalty at all.

Caveat on LOCO: its metrics pool all 138 single-line fold predictions (a
one-line fold can't support per-fold metrics), so it has no variance
estimate and isn't perfectly comparable to the repeat-averaged schemes.

**Caveat on the regression R2 (checked, see `docs/pool_correction_investigation.md`
"Follow-up"): mostly reflects correctly separating success from failure,
not fine-grained precision.** The outcome is bimodal (42 failures at
0-0.185, 96 successes at 0.221-~0.92, real gap at the threshold).
Splitting predictions by true clump: within-success R2 is modest but real
(~0.17-0.22); within-failure R2 is strongly negative (~-14, i.e. far worse
than a trivial baseline -- absolute error ~0.135 is nearly the whole
failure clump's true range). Trust this model for success/failure
classification; don't trust the regression output as a precise efficiency
estimate, especially for predicted failures.

Visual confirmation in `plots/plot_regression_diagnostics.png` /
`plots/plot_classification_diagnostics.png` (predicted-vs-true and residual
scatter for regression; ROC curve, confusion matrix, and predicted-
probability histogram for classification, all donor_grouped/nested/
repeat 0). Two things visible there that the numbers alone didn't show:
- Classification's predicted probabilities are almost perfectly bimodal
  (failures near 0, successes near 1, almost nothing in between) --
  directly mirrors the gap in the true label and is *why* AUC~0.94 is
  achievable even though the regression signal is much weaker.
- Within the success clump, regression residuals show a **systematic
  bias, not just noise**: over-predicts around true~0.3-0.4, under-
  predicts around true~0.8-0.9 -- the model compresses predictions
  toward the middle of the success range (shrinkage/regression-to-the-
  mean) rather than tracking the full spread.

## Which features matter (feature importance)

`feature_importance.py` -> `results/feature_importance_table_{task}_{model}.csv`.
Design and caveats in `docs/feature_importance_plan.md`. The single
configuration these tables are fitted at comes from the **flat** grid, not
nested -- see "Flat vs nested: which is used for what" for why. Results below are on
the full 258-fold data; they were essentially unchanged from the 5-repeat
first pass (selection frequencies firmed up slightly, LOCO deltas moved
by <0.005), so they look stable.

> **These tables are the BASELINE basis only.** For the current numbers on
> both bases, including the qualifying-only variant, see
> **"Latest results"** at the bottom of this file -- that section is
> authoritative. The narrative and interpretation here still hold; only
> read the specific figures from the bottom table.

### What the "one-SE rule" is

A standard way to pick a hyperparameter config that avoids chasing noise.
Instead of taking the raw best-scoring config (the argmax across all 55
`(k, regularization)` combos), you:

1. Find the best mean cross-validated score.
2. Take every config whose score is within **one standard error** of that
   best score -- i.e. everything statistically indistinguishable from the
   winner given how noisy the estimate is.
3. Among those, pick the **simplest** one.

Rationale: with n=138 and 55 candidates, the literal argmax is often a
noise peak -- some config got lucky on these particular folds. Anything
within 1 SE of it is, on the evidence, just as good, so preferring the
simplest of them gives a more robust and more interpretable model at no
real cost in performance. "Simplest" here is defined in advance (not
after seeing results) as: fewest PCs first, then strongest regularization
(larger `alpha` for ridge/lasso; smaller `C` for logistic, since `C` is
inverse regularization strength).

**The one-SE rule is applied ONLY in flat CV, never inside nested CV.**
`select_winning_config_one_se` lives in `feature_importance.py` and is called
only from there, on the flat grid. Nested CV's inner loop
(`harness.run_nested_cv`) takes the **raw argmax** of the mean inner score:

```python
mean_inner = float(np.mean(inner_scores))
if mean_inner > best_score:                 # plain argmax, no tolerance band
    best_score, best_combo = mean_inner, (k, param_combo)
```

This is not a correctness problem -- nested CV stays unbiased either way,
because the selection happens inside the outer fold and never sees the
outer-test lines. But it is worth knowing for two reasons:

1. **The inner loop is the noisier of the two selections**, not the cleaner
   one. It ranks the same 55 candidates using only 3 inner folds of ~36
   validation lines each, versus flat's 10 repeats over all 138. If argmax on
   55 candidates is a noise peak anywhere, it is there.
2. **It partly explains the hyperparameter instability documented below** --
   the modal nested config winning only ~9-10 of 50 outer folds is exactly
   what raw argmax over 55 near-tied candidates on small inner folds
   produces. A one-SE (or any tolerance) rule in the inner loop would
   concentrate those selections considerably.

Applying it there would change reported performance, so it has not been done
on a whim; it is a live option rather than an oversight.

Implementation note: the tolerance band is SD_across_repeats /
sqrt(n_repeats), with n_repeats read from the fold-assignment data. (An
earlier version used the raw SD, which is ~3x too wide at 10 repeats --
see "Corrections from the audit". Fixed; selections happened not to
change.) Effect here: classification
went from the raw argmax's k=4 to k=2 -- a materially simpler model for
statistically indistinguishable performance.

**Regression (lasso, one-SE-selected k=5):**

| Feature | Univariate ρ | Coef | Sel. freq | LOCO Δ | Perm Δ | SHAP | Technical covariate |
|---|---|---|---|---|---|---|---|
| phat_FPP | 0.19 | -0.045 | 0.88 | -0.0015 | 0.035 (grp) | 0.033 | High (η²=0.76) |
| **phat_NB** | **-0.73** | +0.067 | 1.00 | +0.0017 | 0.035 (grp) | 0.052 | **Low (η²=0.04)** |
| phat_P_FPP | 0.34 | reference | -- | -0.0007 (grp) | 0.035 (grp) | -- | High (0.65) |
| PC1 | 0.56 | 0.004 | 0.80 | -0.0024 | 0.005 | 0.003 | Moderate (0.49) |
| **PC2** | 0.67 | **+0.312** | 1.00 | **+0.0186** | **0.176** | 0.247 | **High (η²=0.76)** |
| PC3 | -0.22 | -0.251 | 1.00 | +0.0141 | 0.146 | 0.189 | Moderate (0.50) |
| PC4 | 0.53 | 0.044 | 1.00 | -0.0009 | 0.015 | 0.035 | Moderate (0.45) |
| PC5 | -0.02 | -0.052 | 1.00 | +0.0016 | 0.009 | 0.040 | High (0.79) |

**Classification (logistic_l1, one-SE-selected k=2** -- notably simpler
than the raw argmax's k=4, the one-SE rule working as intended):

| Feature | Univariate ρ | Coef | Sel. freq | LOCO Δ | Perm Δ | SHAP | Technical covariate |
|---|---|---|---|---|---|---|---|
| phat_FPP | 0.20 | 0.0 (regularized out) | 0.00 | 0.0 | 0.277 (grp) | 0.0 | High (0.76) |
| **phat_NB** | **-0.73** | **-1.535** | 0.90 | **+0.051** | **0.277 (grp)** | **1.190** | **Low (η²=0.04)** |
| phat_P_FPP | 0.41 | reference | -- | **+0.064** (grp) | 0.277 (grp) | -- | High (0.65) |
| PC1 | 0.58 | 0.0 (regularized out) | 0.10 | 0.0 | 0.000 | 0.0 | Moderate (0.49) |
| PC2 | 0.59 | +0.648 | 0.82 | +0.004 | 0.080 | 0.514 | High (0.76) |

("grp" = grouped: permutation always shuffles the 3 proportions as one
block, since permuting one alone implies an impossible third coordinate.)

**LOCO vs. permutation disagree for the regression proportions, and the
disagreement is the finding**: permutation Δ=0.035 (the fitted model does
rely on them) but LOCO Δ≈0 (refitting without them costs nothing). That's
the correlated-feature signature -- the proportions carry real signal,
but it's redundantly available in the PCs, so a refit compensates. For
classification both measures agree they're essential (perm 0.277, LOCO
0.064, both largest in the table), meaning there the proportions carry
something the PCs don't replicate. Running only one of the two measures
would have told a misleading story either way.

Three findings:
1. **`phat_NB` is the standout trustworthy feature** -- strongest
   univariate signal of anything (ρ=-0.73), near-always selected, largest
   individual LOCO contribution for classification, and by far the
   cleanest technically (η²=0.04, the only "Low" feature in either
   table). The one result to lean on biologically.
2. **PC2 is the strongest PC for regression but is heavily
   pool-confounded** (η²=0.76) -- largest coefficient, permutation and
   LOCO delta for regression, but its predictive power may be
   substantially technical rather than biological. Treat PC2-based claims
   cautiously. (Note: an earlier version of this README claimed PC2 also
   dominated SHAP; that was a units bug -- see "Corrections from the
   audit" below. Corrected, `phat_NB` has the larger SHAP value for
   classification, 1.19 vs PC2's 0.51.)
3. **Proportions matter greatly for classification, almost not at all
   for regression** -- dropping both costs +0.064 AUC (largest single
   effect in either table) but costs regression ~nothing (-0.0007).

Caveat: `phat_NB`'s univariate correlation is negative (-0.73) but its
*regression* coefficient is positive (+0.067) -- a suppression effect
from multivariate adjustment (PC2/PC3 absorb the shared signal). Don't
read that sign in isolation; the univariate direction and the
classification coefficient (-1.53) agree that more NB at D11 -> worse D52
outcome.

### Why the tables stop at PC5 (regression) / PC2 (classification)

`k` is a tuned hyperparameter, so the tables only list features actually
*in* the selected model. The flat-CV grid tested k=0..10 exhaustively;
larger k didn't score better, so the one-SE rule took the simpler model.

**Known design limitation** (raised by Codex): we only ever test PC
*prefixes* (PC1..k), never arbitrary subsets. PCs are ordered by variance
in the *predictors*, not by outcome relevance -- so a predictive PC8 sitting
behind noisy PC6/PC7 could be missed, since reaching it requires accepting
k=8 and dragging the noise in with it.

Checked this directly -- univariate association (no model involved) for
every PC, including the excluded ones:

| PC | ρ (regression) | ρ (classification) | pool η² | In model? |
|---|---|---|---|---|
| PC1 | 0.563 | 0.584 | 0.490 | both |
| PC2 | **0.670** | **0.593** | 0.764 | both |
| PC3 | -0.221 | -0.370 | 0.498 | regression only |
| PC4 | 0.525 | 0.477 | 0.448 | regression only |
| PC5 | -0.018 | -0.015 | 0.788 | regression only |
| PC6 | -0.101 | -0.089 | **0.909** | no |
| PC7 | 0.275 | 0.364 | 0.531 | no |
| PC8 | 0.260 | 0.197 | **0.963** | no |
| PC9 | 0.061 | -0.147 | **0.946** | no |
| PC10 | -0.203 | -0.157 | 0.461 | no |

**Conclusion: little is being missed.** No excluded PC approaches the
included ones' univariate strength (best excluded is PC7 at ρ≈0.28-0.36
vs. PC2's 0.67). More tellingly, **PC6/PC8/PC9 are nearly pure batch
signal** (η² = 0.91 / 0.96 / 0.95 -- PC8's means 96% of its variance is
explained by which pool a line came from). Excluding them is a feature,
not a loss.

One nuance: **PC5 is included despite ~zero univariate correlation**
(-0.018). Since the one-SE rule takes the simplest config within
tolerance and still chose k=5 over k=4, PC5 must contribute
multivariately despite being marginally useless -- classic suppressor
behavior (correlating with noise in the other PCs so the model can cancel
it out). Its own LOCO (+0.0016) and permutation (0.0085) deltas are
small, consistent with a minor supporting role rather than a driver.

## Fitting population: the qualifying-only PCA variant

The D11 HVG/PCA basis was originally fit on every cell except the fold's
held-out lines and pool11. That INCLUDES ~25,400 cells from 39 cell lines
that are not among the 138 and never appear in any train or test set, plus
~1,100 cells from non-qualifying pools of qualifying lines. Not leakage
(the held-out mask applies by `cell_line`, so every cell of a held-out
line is dropped regardless of pool), but not the study population either.

Both variants are now built and kept side by side:

| | fit cells (one fold) | outputs |
|---|---|---|
| baseline | 243,066 | `*_full.csv`, unsuffixed tables |
| qualifying-only | 217,211 | `*_qualonly.csv` |

Fold assignments are byte-identical between the two runs (verified: 35,604
rows each), so the comparison is like-for-like and any difference is
attributable to the fitting population alone.

**Result — pre-registered headline (donor_grouped, nested):**

| | baseline | qualifying-only |
|---|---|---|
| regression, ridge, R2 | 0.6531 +/- 0.0211 | **0.6672 +/- 0.0124** |
| classification, logistic_l2, ROC-AUC | 0.9470 +/- 0.0070 | 0.9483 +/- 0.0157 |

Regression improves modestly (+0.014 R2, ~0.67 SD) and — more notably —
its across-repeat SD nearly halves. Classification is a wash: AUC up a
hair, balanced accuracy down a hair, SD doubled. So the extra 25,400
foreign cells were adding noise to the regression basis rather than
helping it.

### Do NOT compare the two importance tables row by row

The components REORDER between bases, so "PC2" does not mean the same
thing in each. Correlating per-line coordinates across the 138 lines:

| | closest match | \|r\| | |
|---|---|---|---|
| baseline PC1 | qualonly PC1 | 0.995 | same component |
| baseline PC2 | qualonly **PC3** | 0.985 | same component |
| qualonly PC2 | baseline PC2 | 0.564 | **a different component** |

Read positionally, the qualonly table looks like the fit restriction
cleaned the pool signal out of the strongest PC (eta2 by pool
0.764 -> 0.114). It did not. The pool-associated axis MOVED from position
2 to position 3 and became slightly *more* pool-associated
(eta2 0.764 -> 0.832). What rose into position 2 is a genuinely different
component that happens to be largely pool-independent (eta2 0.114) and
strongly outcome-associated (Spearman rho -0.745 vs binary `success`,
-0.809 vs continuous `diff_efficiency`).

Explained variance is near-identical across bases
(0.0539/0.0265/0.0259/0.0203 -> 0.0539/0.0268/0.0254/0.0203), i.e. PC2 and
PC3 are near-tied, which is exactly the regime where ordering is unstable.
**Align components by correlation before interpreting anything
biologically.**

## Corrections from the audit

After the analysis was complete, an independent Codex correctness audit
(implementation-vs-intent, not methodology) checked every stated design
intention across `folds.py`, `features.py`, `harness.py`,
`feature_importance.py`, and `run_experiment.py`. Most passed. Three real
defects were found and fixed; all three are recorded here rather than
quietly patched, since two of them changed reported numbers.

1. **SHAP mixed unit systems (changed a conclusion).** Coefficients come
   from a model fit on standardized features (units: per 1 SD), but they
   were multiplied by *raw* `(x - mean)` deviations. Every contribution
   was therefore scaled by that feature's raw SD, systematically
   inflating wide-scale features (PCs, raw SD ~1-9) and deflating narrow
   ones (proportions, raw SD ~0.05) -- exactly the comparison the column
   exists to make. Fixed to use standardized deviations. **This reversed
   the classification SHAP ranking**: `phat_NB` went 0.068 -> 1.190 and
   PC2 went 2.438 -> 0.514, so `phat_NB` now leads rather than trailing
   PC2 by ~36x. It strengthens the `phat_NB` finding and weakens the PC2
   one; the earlier README claim that PC2 dominated SHAP was wrong.
2. **The "one-SE rule" used one SD, not one SE.** ~3x too wide at 10
   repeats, so it selected simpler models than a real one-SE rule would.
   Fixed to SD/sqrt(n_repeats). Verified the tolerance genuinely narrowed
   (lasso: 5 configs within tolerance -> 1), but the *selections* happened
   not to change (k=5 regression, k=2 classification both still qualify),
   so no downstream numbers moved.
3. **`select_headline` didn't filter by model, and the reporting
   compounded it.** The pre-registration named ridge / logistic_l2, but
   the function returned all models and the better-scoring non-registered
   models (lasso R2=0.668, logistic_l1 AUC=0.951) were being quoted as
   "the" headline -- the exact cherry-picking the pre-registration exists
   to prevent. Fixed to filter on the pre-registered model; headline
   numbers are now ridge R2=0.653 and logistic_l2 AUC=0.943, with the L1
   variants reported as secondary.

Feature-importance tables are now produced for **all four** models
(`results/feature_importance_table_{task}_{model}.csv`) rather than just the L1
pair, for the same reason. The L1 tables remain the more informative ones
(only L1 produces sparsity, so `regularized_to_zero` /
`selection_frequency` are meaningful), but both are reported.

## Second correctness review (variant plumbing + test pollution)

A second independent Codex review, run specifically on the qualifying-only
work, found two more defects. Recorded here because one of them means a
previous entry in this README overstated what had been fixed.

All of these share one root cause: **adding a `--suffix` variant required
threading "which variant am I?" through several layers, and it was
threaded through some but not others.** Every resulting failure produced
plausible, non-crashing, wrong output.

1. **Mixed-basis importance tables.** `--suffix` redirected the
   fold-features input and the output filenames, but
   `load_full_fit_features()` always read 005's baseline global PCA. So a
   single `_qualonly` table had fold-level columns (LOCO, permutation,
   selection frequency) from one basis and full-fit columns (coefficient,
   SHAP, univariate, technical covariate) from another. The tell was that
   the variant's `univariate_spearman` column was byte-identical to
   baseline. Fixed by parameterizing `005` (`--restrict-to-qualifying`,
   `--suffix`) and threading `pca_csv` through every full-fit consumer.
2. **That fix was incomplete, and was described here as complete.**
   `compute_shap_like` was still calling `load_full_fit_features()` with
   no argument, so `shap_mean_abs` stayed on the baseline basis. The edit
   had been applied with a string replacement whose pattern did not match,
   and string replacement fails silently. Re-applied with an assertion
   that the edit landed plus a check that no no-arg call survives.
3. **The test suite was writing into the repo.** `run_all` saves
   `results/nested_selections_{task}.csv` to the package directory, and
   `test_run_experiment.py`'s fixture patched its inputs but not that
   path. Every `pytest` run overwrote the committed
   `nested_selections_regression.csv` with 32 rows of synthetic fixture
   data (`k=1, alpha=1.0`, and `repeat=1` for LOCO/LODO, which is
   structurally impossible) — and that corrupted file had been committed.
   `run_all` now takes an injectable `out_dir`, the fixture redirects it,
   and a regression test asserts `run_all` never touches committed files.
4. **`_infer_n_repeats` read a hardcoded table.** It always loaded
   `fold_data/fold_features_D11_full.csv` (falling back silently to
   `fold_features_D11.csv`), not the table the run was using, so the
   one-SE tolerance could come from an unrelated experiment. Harmless in
   fact here — both runs share fold assignments, so `n_repeats` was 10
   either way — but nothing enforced that. Now required and passed
   explicitly; a missing file is a loud error.
5. **`logistic_l1` was not reproducible.** `saga` is a stochastic solver
   and no `random_state` was set, so coefficients drifted run to run
   (measured max|diff| ~5.7e-4), producing spurious 4th-decimal diffs on
   re-runs. `models.SEED` now pins both logistic factories.

A guard now prevents defect 1 from recurring silently: with `--suffix`
set and no explicit basis, the matching suffixed basis is REQUIRED, and
its absence is a hard error naming the command that generates it.

## D30 readiness — NOT ready (read before starting D30 -> D52)

The earlier claim in "Architecture" that feature extraction is simply
"timepoint-parameterized" is too optimistic. What is genuinely reusable:
`folds.py` (fold construction is timepoint-independent — same 138 lines,
same D52 labels), `models.py`, most of `harness.py`, and
`run_experiment.run_all` given a correct feature table.

What must change first, ranked by risk. The loud failures are fine; the
SILENT ones are the danger, and they are marked.

1. **Cell-type schema differs and is hardcoded.** D11 has 3 types
   (`FPP, NB, P_FPP`); **D30 has 7** (`DA, Epen1, FPP, P_FPP, Sert,
   U_Neur1, U_Neur2`). `harness.ALL_PROPORTION_COLS` /
   `MODEL_PROPORTION_COLS` hardcode the D11 three. Fails loudly on a
   missing column, but the compositional reference category must be
   chosen and documented for 7 types, not patched ad hoc.
2. **SILENT: extraction is hardcoded to D11.** `run_feature_extraction.py`
   passes `"D11"` to `load_cell_metadata` and
   `compute_pca_features_for_fold` and names outputs
   `fold_features_D11{suffix}.csv`. Changing only the output filename
   would produce a plausible file named D30 containing D11 features.
3. **SILENT: resume has no provenance.** The incremental-write resume keys
   only on `(scheme, repeat, fold)`. It does not record timepoint,
   fitting population, or PCA settings, so reusing a destination after
   changing any of those silently keeps stale rows. Write a provenance
   header/sidecar and refuse to resume across a settings change.
4. **SILENT: full-fit inputs are D11-specific.** `005` is hardwired to
   `DAY11_FILE` and `d11_*` output names; `feature_importance.py`
   hardcodes `D11_PCA_CSV`, `D11_PROPORTIONS_CSV` and derives the variant
   basis as `d11_pca_coords_per_line{suffix}.csv`. Passing D30 fold
   features *unsuffixed* would combine D30 fold-level values with D11
   global PCA and D11 proportions. (A suffixed call would fail loudly on
   the missing `d11_...{suffix}` file.) Make these explicit parameters.
5. **`run_variant.py` input is hardcoded** to `fold_features_D11{suffix}`.
6. **Decide whether excluding pool11 is still justified at D30.** It is
   currently inherited from a D11 sequencing-depth observation without
   re-examination.

Note the qualifying `(cell_line, pool)` combos ARE the same 159 across
138 lines at both timepoints (qualification already requires >=10 cells at
D11 and D30 and D52), but the per-combo cell COUNTS are not: zero of 159
combos has identical D11 and D30 counts (D11 range 25-14,640; D30 range
13-13,112; correlation 0.748). Sampling precision therefore differs
substantially, which is a further reason D11 proportions must never be
reused as D30 features.

## What PC1 is: a D11 proliferation axis

> Full write-up, loadings tables and literature sources:
> **`modeling/docs/PCA_interpretation.md`**. Summary below.

`005_d11_pca_features.py` now also writes
`metadata_eda/pca/d11_pca_gene_loadings{suffix}.csv` -- all 10 PCs x 2000 HVGs,
long format, sorted by |loading| within each PC.

PC1 (5.4% of HVG variance, Spearman rho = **+0.563** vs D52 efficiency, so
higher PC1 = higher efficiency) is a **cell-cycle / proliferation axis**:

| positive loadings (high in high-efficiency lines) | negative loadings |
|---|---|
| HMGB2, PTTG1, NUSAP1, UBE2C, CENPF, CKS2, TOP2A, PLK1, CCNB1/CCNB2, KPNA2, CDC20, BIRC5, CKS1B, AURKA/AURKB, CCNA2, TPX2, CDK1, KIF2C, SMC4, CDKN3, MKI67 | RPL12, RPL10, RPS3, RPS12, RPS28, RPL37A (ribosomal proteins); EIF3E, EIF3L, EIF4A2 (translation initiation); GAPDH, COX7C, UQCRB, TOMM7, ETFB (housekeeping/OXPHOS); SNHG8, MIAT, EPB41L4A-AS1, APOE, CCND2, SLC2A1 |

The positive side is a textbook G2/M signature. So lines whose D11
cultures are **more proliferative differentiate better by D52**. This is
consistent with the line-level correlations: PC1 tracks `phat_P_FPP`
(proliferating FPP) at rho = **+0.798** and is negatively correlated with
`phat_NB` (rho = -0.552), the neuroblast fraction that independently
predicts *worse* outcome.

Top 25 genes carry 29.6% of PC1's total weight (loadings are unit-norm);
median |loading| is 0.0012, so PC1 is dominated by a compact, coherent
gene set rather than being diffuse.

**Replicates across bases**: PC1 loadings correlate r = **+0.9995**
between the baseline and qualifying-only bases, with 25/25 top-gene
overlap. Unlike PC2/PC3, PC1 is stable and safe to compare across the two.

Caveats:
- **Loading sign is arbitrary in general.** It is interpretable here only
  because PC1 was oriented against the outcome (rho > 0). Do not carry the
  sign convention to other PCs without re-checking.
- Loadings are weights on **standardized** expression (per 1 SD of that
  gene, after the +/-`SCALE_CLIP` clip), so they are comparable across
  genes of different absolute expression -- but they are not fold-changes.
- **PC1 is moderately pool-associated (eta2 = 0.49).** A
  proliferation-vs-ribosomal/housekeeping contrast is also the classic
  shape of a library-size / cell-quality technical axis, so some of this
  variance is likely technical. The biological reading is plausible and
  matches the `phat_P_FPP` correlation, but PC1 is not as technically
  clean as `phat_NB` (eta2 = 0.04).

## Pool-correction dropped

Investigated (with an independent Codex review) why pool-correction cut
regression R2 by ~0.10-0.12 but barely touched classification AUC. No
single settled explanation emerged -- see `docs/pool_correction_investigation.md`
for the concise summary (my reasoning, Codex's oracle-test refutation of
my leading theory, its alternative explanation, and two real
implementation issues found along the way). Decision: drop pool-correction
from the active pipeline rather than resolve this now -- it added
significant complexity for an unclear net benefit at this prototype
stage. `features.py`/`harness.py` still support it if revisited.

## Setup
- 138 lines, one row each. Regression (`diff_efficiency`) + classification (`>=0.2`), same features/folds for both.
- Features: 3 D11 proportions (always in) + up to 10 PCs, `k=0..10` tuned by truncating one PCA fit (no refit per k).
- Models, full parallel arms: Ridge + Lasso (regression); L2 + L1 logistic (classification).

### Hyperparameters: 55 configurations per model

Two are tuned, jointly: `k` (11 values) x regularisation strength (5 values).

| hyperparameter | values | n |
|---|---|---|
| `k` -- number of PCs | 0, 1, 2, ... 10 | **11** |
| regularisation strength | see per-model grid below | **5** |

`k=0` means proportions only. The 2 modelled proportions are always in; `k`
controls how many PCs join them. PCs are hierarchical, so `k` truncates
columns from one PCA fit -- no refit per `k`.

| model | task | parameter | values |
|---|---|---|---|
| ridge | regression | `alpha` | 0.01, 0.1, 1, 10, 100 |
| lasso | regression | `alpha` | 0.001, 0.01, 0.1, 1, 10 |
| logistic_l2 | classification | `C` | 0.01, 0.1, 1, 10, 100 |
| logistic_l1 | classification | `C` | 0.01, 0.1, 1, 10, 100 |

**`C` is sklearn's inverse regularisation strength** -- the logistic
equivalent of `alpha`, running the opposite direction (`C` ~ 1/`alpha`):

- `alpha` **up** -> more regularisation -> coefficients shrink
- `C` **up** -> *less* regularisation -> coefficients grow freely

That inversion is why "simplest" in the one-SE rule means *larger alpha,
smaller C*: both mean more heavily regularised.

Lasso's `alpha` grid sits one decade lower (0.001-10) because L1 shrinks
coefficients to exactly zero, so it needs smaller values before it stops
zeroing out everything.

## CV — both variants run, neither "primary"
- **Plain split**: ordinary split over the 138 rows (grouping by cell_line is a no-op here — each row is already one whole line).
- **Donor-grouped split**: lines sharing a donor forced into the same fold (136/138 lines have a donor-sibling; one donor has 18 lines). See "What donor_grouped actually does" below.
- **LOCO + LODO**: full leave-one-line-out and leave-one-donor-out runs, as comparability checks vs. the paper's (non-donor-grouped) LOOCV.
- **Flat CV + nested CV**: both run (flat = grid search, report best mean CV score; nested = outer loop for honest estimate, inner loop tunes `k`/regularization). Inner tuning metric: PR-AUC (classification), MAE (regression). They are used for **different jobs** -- see below.
- **Repeated K-fold**: report mean ± SD **per repeat** — never pool predictions across repeats (double-counts each line); pooling *within* one repeat is fine.
- **Weighting**: equal weight per line (not per donor).



### Flat vs nested: which is used for what

Both are computed for every scheme, but only one is ever quoted as
performance. Per task, `results/results_{task}.csv` holds:

| tuning | rows | what a row is | used for |
|---|---|---|---|
| flat | 440 | one `(scheme, model, k, param)` combination | the robustness grid, **and picking the config for the importance tables** |
| nested | 8 | one `(scheme, model)` | **every reported performance number** |

(440 = 4 schemes x 2 models x 11 `k` x 5 regularisation values; 8 = 4 schemes
x 2 models.)

**Every performance figure in this README is nested.** R2 0.653, AUC 0.947,
the all-schemes table -- all nested. Flat scores are optimistically biased,
because the same predictions are used both to choose the winning
configuration and to report its score. They are never quoted as performance.

**But the feature-importance tables take their configuration from the FLAT
grid** (`feature_importance.py`, `select_winning_config_one_se`, which filters
`tuning == "flat"`). That looks inconsistent and is deliberate:

> Nested CV does not produce *a* configuration. It produces 50 of them -- one
> per outer fold -- and they disagree substantially (the modal choice wins
> only ~9-10 of 50; see "Hyperparameter selection is unstable"). There is
> nothing to read a coefficient off.

So to fit one set of coefficients and SHAP values you need one fixed model,
and the one-SE rule scans the flat grid to choose it defensibly. Using flat
here is safe because **the flat score is never reported** -- it is used only
to *rank* configurations, not to claim performance.

The consequence to keep in mind: **the importance tables and the headline
numbers describe slightly different models.** The "Configurations selected"
table at the bottom lists both side by side (modal nested vs. one-SE flat)
precisely so the gap is visible rather than implicit.


### How nested CV actually runs (one repeat)

```
138 cell lines
│
├─ Outer fold 0: hold out 28 lines
│  ├─ Use the remaining 110 for inner CV
│  │  ├─ Inner fold 0: ~74 train / ~36 validate -> score every candidate
│  │  ├─ Inner fold 1: ~74 train / ~36 validate -> score every candidate
│  │  └─ Inner fold 2: ~74 train / ~36 validate -> score every candidate
│  ├─ Average the 3 inner scores per candidate; pick the best (raw argmax --
  │     the one-SE rule is NOT applied here, see "What the one-SE rule is")
│  ├─ Discard the three inner models
│  ├─ Refit ONE model on all 110 outer-training lines
│  └─ Predict the 28 outer-test lines ONCE
│
├─ Outer fold 1: 28 held-out lines predicted once
├─ Outer fold 2: 27
├─ Outer fold 3: 27
└─ Outer fold 4: 28
                 ───
Final out-of-fold predictions: 28+28+27+27+28 = 138, each line exactly once
```

The 28 outer-test lines are invisible for the entire inner loop -- candidate
scoring, selection and refit all happen inside the 110. Outer fold 0
contributes 28 predictions, not 3 or 55.

Per model, per outer fold: 55 candidates x 3 inner folds = **165 inner fits**,
then 1 refit. Across 50 outer folds and 2 model families: **16,500 inner fits
+ 100 refits** per task.

Totals for one reported number:

| | |
|---|---|
| outer fitted models | 5 folds x 10 repeats = **50** |
| final prediction values | 138 lines x 10 repeats = **1,380** |
| metrics computed | **10** (one pooled AUC/R2 per repeat) |
| reported | mean +/- SD across those 10 |

### How flat CV runs -- and one thing it does NOT do

**Flat CV does not select a config per fold and then compare folds.** There
is no "best config for fold 0". Each candidate is scored the same pooled way
the headline is:

```
for each of the 55 candidates:
    for each repeat (10):
        predict all 5 folds' held-out lines -> pool -> 138 predictions
        -> ONE score for that (candidate, repeat)
    -> 10 scores -> mean +/- SD across repeats
compare the 55 candidates by that mean; apply the one-SE rule once
```

So a flat row's `mean +/- SD` is across the **10 repeats**, exactly like a
nested row -- not across the 5 folds and not across the 28 lines in a fold.
The config is chosen **once**, from 55 means, not 5 times and then reduced.

The only structural difference from nested is *where the choice happens*:
nested chooses inside each outer fold (so the choice never sees the test
lines); flat chooses once, afterwards, over all the data. That is precisely
why the flat score is optimistically biased and is never reported.

### Best flat-CV configurations (donor_grouped)

Selection metric is MAE for regression and ROC-AUC for classification; R2 and
AUC shown for readability.

| model | raw argmax | one-SE choice | configs within 1 SE |
|---|---|---|---|
| ridge | k=5, alpha=0.01 (R2 0.6801) | **k=5, alpha=0.1** (R2 0.6804) | 2/55 |
| lasso | k=5, alpha=0.001 (R2 0.6838) | **k=5, alpha=0.001** (R2 0.6838) | 2/55 |
| logistic_l2 | k=6, C=0.01 (AUC 0.9638) | **k=4, C=0.01** (AUC 0.9629) | 2/55 |
| logistic_l1 | k=4, C=0.1 (AUC 0.9640) | **k=2, C=0.1** (AUC 0.9632) | 9/55 |

The one-SE rule only moves the answer where the band is genuinely wide:
`logistic_l1` has 9 statistically indistinguishable configs so it drops from
k=4 to k=2, and `logistic_l2` from k=6 to k=4. For ridge and lasso only 2
configs qualify, so it barely bites.

### Reconciling the two: they produce different things

The tension is real -- nested CV already does model selection, so why also
run flat? Because they answer different questions and yield different kinds
of output.

| | question | output | score trustworthy? |
|---|---|---|---|
| nested | how well does *the procedure* generalise? | a **number** | yes -- this is the reported performance |
| flat | which single configuration ranks best? | a **choice** | no -- never quoted |

Nested CV evaluates the whole pipeline *including its own tuning step*, which
is what makes it honest. That same property is why it **cannot hand you a
model**: it fits 50, each choosing its own config, and they disagree (the
modal choice wins only ~9-10 of 50). There is no single "the nested model" to
read a coefficient off.

So the division is:

> **Performance claims come from nested. Configuration choices come from
> flat. No number is ever taken from flat.**

Being straight about the residual nuance: the flat ranking is computed over
all 138 lines, so the chosen config has seen every outcome. For the
coefficients and SHAP that is fine -- they are descriptive summaries of a
full-data fit, not performance claims. For the LOCO and permutation deltas it
is slightly less clean: those are measured out-of-fold but at a config chosen
using all the data, a mild optimism. It largely cancels because the deltas
are within-fold with/without comparisons at the same config, but "largely
cancels" is not "does not exist". The alternative -- using the modal nested
config -- is also defensible and close (regression: one-SE k=5 vs modal k=4);
see "Configurations selected". The one-SE flat config was chosen because a
mode that wins 9-10 times out of 50 is a weak summary.

### What `donor_grouped` actually does

Two things it is **not**: it does not put the same donors in every fold, and
it does not average lines within a donor.

**One row is always one cell line.** 138 rows, 138 predictions. The only
averaging is pool-then-line on the features (mean per `(cell_line, pool)`,
then unweighted mean across that line's pools). Donors are never collapsed;
weighting is equal per line, not per donor.

What it does: **every line from a given donor lands in the same fold**, so a
donor is entirely in train or entirely in test, never split across the two.
Different folds hold out different donors.

`folds.py` passes `groups=lines["donor"]` to `StratifiedGroupKFold`. From the
committed assignments, repeat 0 (138 lines, 20 donors, 5 folds):

| fold | test lines | from donors |
|---|---|---|
| 0 | 28 | 2 |
| 1 | 28 | 5 |
| 2 | 27 | 3 |
| 3 | 27 | 5 |
| 4 | 28 | 5 |

Donors appearing in both train and test of the same fold: **0**.

Note fold 0 holds out 28 lines from only **2** donors -- that is the donor
with 18 lines. Folds are balanced by line count, not donor count.

**Why:** 136 of 138 lines share a donor with another line. Without grouping a
model can memorise a donor's genotype from one line and be scored on its
sibling. Grouping blocks that; it costs ~0.012 R2, which is the number
quantifying donor leakage.

Not to be confused with **LODO**, which is one fold per donor (20 folds,
each holding out that donor's entire set of lines).

## Leakage safety
- PCA/HVG refit **per fold**: exclude held-out lines' (or held-out donors' lines') cells from the fitting pool; project their cells into that fold's PCA space. Generalizes `005_d11_pca_features.py`'s `fit_mask` beyond just pool11.
- **pool11 rule frozen** from EDA (depth ~1,900 vs ~10-18K elsewhere) — always excluded from fitting when present in training, never tuned on results. Confirmed harmless: 4.2% of D11 cells (9,603/227,983); ~218K remain.
- **Pool correction**: residualize each feature on pool identity using that fold's *training* pool means, computed at the `(cell_line, pool)` level before line-averaging (slots into the existing `004`/`005` pool→line pipeline). Both corrected/uncorrected variants run.
- **No donor feature-correction** — donor is likely real signal, not nuisance; donor-*grouping* is the right tool, correcting it out risks deleting the signal we want.
- Known unfixable limitation: D11/D52 cell-type labels came from the paper's own global clustering, not refit per fold — diffuse, minor, unavoidable.

## Metrics
- Classification: ROC-AUC, PR-AUC, balanced accuracy, sensitivity/specificity @ 0.2, F1, confusion matrix, calibration (Brier); vs. majority-class + prevalence baselines. Keep the 0.2 *outcome* threshold separate from the probability *decision* threshold (default 0.5).
- Regression: MAE, RMSE, R²; vs. mean/median baseline; check [0,1] violations; residuals near 0.2.

## Reproducibility
- Fixed seeds; persist actual fold assignments (line IDs per scheme × repeat × fold) to a file, not just the seed.
- `models.SEED` pins `random_state` on both logistic factories. `saga` (used by `logistic_l1`) is stochastic; unseeded it drifted ~5.7e-4 between runs, which showed up as spurious 4th-decimal diffs in regenerated importance tables.
- Fold features are written incrementally per fold and a run can be resumed; see the D30 section for the provenance gap that makes resume unsafe across a settings change.

## Scope/cost warning
Full parallel grid (2 groupings × 2 tuning × 2 corrections × 2 model families, repeated K-fold, + LOCO/LODO) means potentially hundreds of fold-specific PCA refits (~90s each in `005`). Time one fold first; cache fold-specific fits before committing to the full run.

## Results distillation

**"Headline model" = the single configuration reported as *the* result**, as
opposed to the robustness grid of everything that was run. Without one there
are 16 numbers per task (4 schemes x 2 tuning modes x 2 model families) and
no answer -- and whichever gets quoted, a reader cannot tell whether it was
chosen before or after seeing the scores.

It is **pre-registered, not selected on performance**. Committed as constants
in `run_experiment.py` and enforced by `select_headline()`:

```python
HEADLINE       = {"scheme": "donor_grouped", "tuning": "nested", "pool_correction": False}
HEADLINE_MODEL = {"regression": "ridge", "classification": "logistic_l2"}
```

Each choice was made on design grounds, not on scores: `donor_grouped`
because it is the strictest realistic test; `nested` because flat CV tunes
and reports on the same rows; no pool-correction per the investigation below;
ridge / logistic_l2 as the L2 default, with L1 as the variant.

Note the L1 variants score *higher* (lasso R2 0.669 vs ridge 0.653;
logistic_l1 AUC 0.951 vs 0.947). A paired per-repeat comparison at matched
one-SE configs shows lasso's edge is real but tiny -- +0.0033 R2, winning
9/10 repeats, t-test p=0.0006. Most of the larger headline gap is nested CV's
per-fold reselection noise, not a model difference. They are reported as
labelled secondary results; quoting them as "the" headline is exactly the
cherry-picking the pre-registration exists to prevent, and that failure
actually occurred once -- see "Corrections from the audit" #3.

- Full grid still gets computed, but pre-register ONE headline config *before* looking at results, to avoid cherry-picking: **donor-grouped + nested CV + repeated stratified K-fold + Ridge/L2-logistic** (pool-correction dropped entirely -- see "Pool-correction dropped" above). Report this first, with mean ± SD.
- Everything else (plain grouping, flat CV, Lasso/L1, LOCO/LODO) becomes a secondary "robustness grid" — one table/heatmap (rows = config, cols = key metric) to scan for consistency, not N separate headline results.

## Architecture (Open-Closed + reusable for D30→D52)
- **Model registry, not conditionals**: a list of `ModelSpec(name, estimator_factory, param_grid)` entries (Ridge/Lasso/L2-logistic/L1-logistic to start). The CV harness only ever calls the sklearn estimator interface (`.fit`/`.predict`/`.predict_proba`) generically — adding a model later (e.g. random forest) = one new registry entry, zero changes to harness code.
- **Feature extraction (timepoint-parameterized) vs. CV harness (timepoint-agnostic)** — *aspirational; see "D30 readiness" above for what actually still hardcodes D11*: feature extraction takes a timepoint (D11 now, D30 later — same 138 qualifying lines, since the qualifying list already requires ≥10 cells at D11 *and* D30 *and* D52) and a fold's train/held-out line split, returns feature matrices. The CV/tuning/metrics harness takes `X_train, y_train, X_test, y_test` + the model registry and is completely agnostic to which timepoint produced the features. D30→D52 later = swap the feature-extraction call, harness untouched.

## Files
- `features.py` — feature-extraction, parameterized by timepoint (D11 now; reusable for D30).
- `folds.py` — fold construction (plain/donor-grouped repeated stratified K-fold, LOCO, LODO) — timepoint-agnostic.
- `models.py` — model registry (`ModelSpec` list).
- `harness.py` — CV/tuning/metrics harness using the registry + fold module.
- `run_regression.py` / `run_classification.py` — orchestration entry points.
- Summary table + comparison plots (headline config; plain vs. donor gap; flat vs. nested; corrected vs. uncorrected).

## Verify
- Per-fold fitting cell counts stay large; fold assignments reproducible; classification metrics only where both classes present; repeated-CV reported per-repeat not pooled; regression predictions checked against [0,1].

## Latest results

Authoritative current numbers. Regenerated 2026-09-09 in a single clean
run on one consistent code state (after the seed pin and the fixes in
"Second correctness review"), so no figure here is a mix of pre- and
post-fix runs. Everything below is the pre-registered configuration:
**donor_grouped grouping, nested CV, no pool-correction**.

### Headline (pre-registered model per task)

| task | model | baseline | qualifying-only |
|---|---|---|---|
| regression | **ridge** | R2 = 0.6531 +/- 0.0211 | **R2 = 0.6672 +/- 0.0124** |
| classification | **logistic_l2** | ROC-AUC = 0.9470 +/- 0.0070 | **ROC-AUC = 0.9483 +/- 0.0157** |

### Secondary (L1 variants -- NOT the headline; see "Results distillation")

| task | model | baseline | qualifying-only |
|---|---|---|---|
| regression | lasso | R2 = 0.6685 +/- 0.0150 | R2 = 0.6725 +/- 0.0108 |
| classification | logistic_l1 | ROC-AUC = 0.9512 +/- 0.0081 | ROC-AUC = 0.9500 +/- 0.0060 |

Restricting the PCA fit to the study population helps regression modestly
(+0.014 R2) and roughly halves its across-repeat SD; classification is a
wash (AUC +0.001, SD doubled). Full discussion in "Fitting population"
above, including why the two importance tables must NOT be compared by PC
label.

### Configurations selected

Nested CV re-selects per outer fold, so `modal` is the most common choice
across 50 folds, not "the" configuration -- the spread is wide (see
"Regenerated results" note below). The one-SE column is the single config
used for the feature-importance tables.

| run | model | modal nested config | one-SE config (importance tables) |
|---|---|---|---|
| baseline | ridge | k=4, alpha=0.01 | k=5, alpha=0.1 |
| baseline | lasso | k=5, alpha=0.001 | k=5, alpha=0.001 |
| baseline | logistic_l2 | k=2, C=100.0 | k=4, C=0.01 |
| baseline | logistic_l1 | k=2, C=1.0 | k=2, C=0.1 |
| qualifying-only | ridge | k=6, alpha=1.0 | k=6, alpha=0.1 |
| qualifying-only | lasso | k=6, alpha=0.001 | k=6, alpha=0.001 |
| qualifying-only | logistic_l2 | k=8, C=0.01 | k=8, C=0.01 |
| qualifying-only | logistic_l1 | k=0, C=0.1 | k=3, C=0.1 |

Note `logistic_l1` selects **k=0** most often on the qualifying-only basis
-- i.e. the modal nested model uses the cell-type proportions and no PCs
at all, and still reaches AUC 0.950. Consistent with the standing finding
that the proportions (`phat_NB` especially) carry most of the
classification signal.

### Hyperparameter selection is unstable -- read single-config tables with that in mind

Recovering the per-fold selections (they were previously computed and
dropped) shows no configuration dominates. Across 50 outer folds the modal
choice wins only ~9-10 times, spread over 13-19 distinct configs, with `C`
ranging the full four orders of magnitude. Full per-fold record in
`results/nested_selections_{task}{suffix}.csv` (516 rows each). This is why the
headline is pre-registered and why the one-SE rule is used -- and why any
single-config coefficient table is one draw from a wide distribution.

### Files backing this section

| | baseline | qualifying-only |
|---|---|---|
| grids | `results/results_{task}.csv` | `results/results_{task}_qualonly.csv` |
| per-fold selections | `results/nested_selections_{task}.csv` | `results/nested_selections_{task}_qualonly.csv` |
| importance | `results/feature_importance_table_{task}_{model}.csv` | `results/..._{model}_qualonly.csv` |
| fold features | `fold_data/fold_features_D11_full.csv` | `fold_data/fold_features_D11_qualonly.csv` |
| global PCA basis | `metadata_eda/pca/d11_pca_coords_per_line.csv` | `..._qualonly.csv` |
