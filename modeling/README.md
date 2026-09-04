# D11 -> D52 modeling: CV/modeling plan

Predicting D52 differentiation efficiency from D11 features (138 cell
lines, from `metadata_eda/qualifying_cell_line_pool_min10_per_timepoint.csv`
and `metadata_eda/d52_diff_efficiency_label.csv`). Designed to generalize
to a D30 -> D52 baseline model later (see Architecture).

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

| task | model | key metrics |
|---|---|---|
| regression | lasso | MAE=0.135, RMSE=0.172, **R2=0.668 ± 0.015** |
| regression | ridge | MAE=0.138, RMSE=0.176, R2=0.653 ± 0.021 |
| classification | logistic_l1 | **ROC-AUC=0.951 ± 0.008**, PR-AUC=0.970, balanced_acc=0.913 |
| classification | logistic_l2 | ROC-AUC=0.943 ± 0.008, PR-AUC=0.965, balanced_acc=0.906 |

D11 features are clearly predictive of D52 differentiation efficiency,
even under the strict donor-grouped (never-seen-donor) test. Going from
the 5-repeat first pass to the full 10-repeat + LOCO/LODO run barely
moved these numbers (R2 0.667->0.668, AUC 0.949->0.951) -- the estimates
were already stable.

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

**Caveat on the regression R2 (checked, see `pool_correction_investigation.md`
"Follow-up"): mostly reflects correctly separating success from failure,
not fine-grained precision.** The outcome is bimodal (42 failures at
0-0.185, 96 successes at 0.221-~0.92, real gap at the threshold).
Splitting predictions by true clump: within-success R2 is modest but real
(~0.17-0.22); within-failure R2 is strongly negative (~-14, i.e. far worse
than a trivial baseline -- absolute error ~0.135 is nearly the whole
failure clump's true range). Trust this model for success/failure
classification; don't trust the regression output as a precise efficiency
estimate, especially for predicted failures.

Visual confirmation in `plot_regression_diagnostics.png` /
`plot_classification_diagnostics.png` (predicted-vs-true and residual
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

`feature_importance.py` -> `feature_importance_table_{task}_{model}.csv`.
Design and caveats in `feature_importance_plan.md`. Results below are on
the full 258-fold data; they were essentially unchanged from the 5-repeat
first pass (selection frequencies firmed up slightly, LOCO deltas moved
by <0.005), so they look stable.

**Regression (lasso, one-SE-selected k=5):**

| Feature | Univariate ρ | Coef | Sel. freq | LOCO Δ | Perm Δ | SHAP | Technical covariate |
|---|---|---|---|---|---|---|---|
| phat_FPP | 0.19 | -0.045 | 0.88 | -0.0015 | 0.035 (grp) | 0.003 | High (η²=0.76) |
| **phat_NB** | **-0.73** | +0.067 | 1.00 | +0.0017 | 0.035 (grp) | 0.003 | **Low (η²=0.04)** |
| phat_P_FPP | 0.34 | reference | -- | -0.0007 (grp) | 0.035 (grp) | -- | High (0.65) |
| PC1 | 0.56 | 0.004 | 0.80 | -0.0024 | 0.005 | 0.005 | Moderate (0.49) |
| **PC2** | 0.67 | **+0.312** | 1.00 | **+0.0186** | **0.176** | **1.17** | **High (η²=0.76)** |
| PC3 | -0.22 | -0.251 | 1.00 | +0.0141 | 0.146 | 0.44 | Moderate (0.50) |
| PC4 | 0.53 | 0.044 | 1.00 | -0.0009 | 0.015 | 0.046 | Moderate (0.45) |
| PC5 | -0.02 | -0.052 | 1.00 | +0.0016 | 0.009 | 0.068 | High (0.79) |

**Classification (logistic_l1, one-SE-selected k=2** -- notably simpler
than the raw argmax's k=4, the one-SE rule working as intended):

| Feature | Univariate ρ | Coef | Sel. freq | LOCO Δ | Perm Δ | SHAP | Technical covariate |
|---|---|---|---|---|---|---|---|
| phat_FPP | 0.20 | 0.0 (regularized out) | 0.00 | 0.0 | 0.277 (grp) | 0.0 | High (0.76) |
| **phat_NB** | **-0.73** | **-1.535** | 0.90 | **+0.051** | **0.277 (grp)** | 0.068 | **Low (η²=0.04)** |
| phat_P_FPP | 0.41 | reference | -- | **+0.064** (grp) | 0.277 (grp) | -- | High (0.65) |
| PC1 | 0.58 | 0.0 (regularized out) | 0.10 | 0.0 | 0.000 | 0.0 | Moderate (0.49) |
| PC2 | 0.59 | +0.648 | 0.82 | +0.004 | 0.080 | 2.44 | High (0.76) |

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
2. **PC2 dominates regression but is heavily pool-confounded**
   (η²=0.76) -- largest coefficient, SHAP, and LOCO delta, but its
   predictive power may be substantially technical rather than
   biological. Treat PC2-based claims cautiously.
3. **Proportions matter greatly for classification, almost not at all
   for regression** -- dropping both costs +0.064 AUC (largest single
   effect in either table) but costs regression ~nothing (-0.0007).

Caveat: `phat_NB`'s univariate correlation is negative (-0.73) but its
*regression* coefficient is positive (+0.067) -- a suppression effect
from multivariate adjustment (PC2/PC3 absorb the shared signal). Don't
read that sign in isolation; the univariate direction and the
classification coefficient (-1.53) agree that more NB at D11 -> worse D52
outcome.

## Pool-correction dropped

Investigated (with an independent Codex review) why pool-correction cut
regression R2 by ~0.10-0.12 but barely touched classification AUC. No
single settled explanation emerged -- see `pool_correction_investigation.md`
for the concise summary (my reasoning, Codex's oracle-test refutation of
my leading theory, its alternative explanation, and two real
implementation issues found along the way). Decision: drop pool-correction
from the active pipeline rather than resolve this now -- it added
significant complexity for an unclear net benefit at this prototype
stage. `features.py`/`harness.py` still support it if revisited.

## Setup
- 138 lines, one row each. Regression (`diff_efficiency`) + classification (`>=0.2`), same features/folds for both.
- Features: 3 D11 proportions (always in) + up to 10 PCs, `k=0..10` tuned by truncating one PCA fit (no refit per k).
- Hyperparams: `k` (#PCs) and regularization strength (ridge/lasso `alpha`, logistic `C`).
- Models, full parallel arms: Ridge + Lasso (regression); L2 + L1 logistic (classification).

## CV — both variants run, neither "primary"
- **Plain split**: ordinary split over the 138 rows (grouping by cell_line is a no-op here — each row is already one whole line).
- **Donor-grouped split**: lines sharing a donor forced into the same fold (136/138 lines have a donor-sibling; one donor has 18 lines).
- **LOCO + LODO**: full leave-one-line-out and leave-one-donor-out runs, as comparability checks vs. the paper's (non-donor-grouped) LOOCV.
- **Flat CV + nested CV**: both run and reported (flat = grid search, report best mean CV score; nested = outer loop for honest estimate, inner loop tunes `k`/regularization). Inner tuning metric: PR-AUC (classification), MAE (regression).
- **Repeated K-fold**: report mean ± SD **per repeat** — never pool predictions across repeats (double-counts each line); pooling *within* one repeat is fine.
- **Weighting**: equal weight per line (not per donor).

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

## Scope/cost warning
Full parallel grid (2 groupings × 2 tuning × 2 corrections × 2 model families, repeated K-fold, + LOCO/LODO) means potentially hundreds of fold-specific PCA refits (~90s each in `005`). Time one fold first; cache fold-specific fits before committing to the full run.

## Results distillation
- Full grid still gets computed, but pre-register ONE headline config *before* looking at results, to avoid cherry-picking: **donor-grouped + nested CV + repeated stratified K-fold + Ridge/L2-logistic** (pool-correction dropped entirely -- see "Pool-correction dropped" above). Report this first, with mean ± SD.
- Everything else (plain grouping, flat CV, Lasso/L1, LOCO/LODO) becomes a secondary "robustness grid" — one table/heatmap (rows = config, cols = key metric) to scan for consistency, not N separate headline results.

## Architecture (Open-Closed + reusable for D30→D52)
- **Model registry, not conditionals**: a list of `ModelSpec(name, estimator_factory, param_grid)` entries (Ridge/Lasso/L2-logistic/L1-logistic to start). The CV harness only ever calls the sklearn estimator interface (`.fit`/`.predict`/`.predict_proba`) generically — adding a model later (e.g. random forest) = one new registry entry, zero changes to harness code.
- **Feature extraction (timepoint-parameterized) vs. CV harness (timepoint-agnostic)**: feature extraction takes a timepoint (D11 now, D30 later — same 138 qualifying lines, since the qualifying list already requires ≥10 cells at D11 *and* D30 *and* D52) and a fold's train/held-out line split, returns feature matrices. The CV/tuning/metrics harness takes `X_train, y_train, X_test, y_test` + the model registry and is completely agnostic to which timepoint produced the features. D30→D52 later = swap the feature-extraction call, harness untouched.

## Files
- `features.py` — feature-extraction, parameterized by timepoint (D11 now; reusable for D30).
- `folds.py` — fold construction (plain/donor-grouped repeated stratified K-fold, LOCO, LODO) — timepoint-agnostic.
- `models.py` — model registry (`ModelSpec` list).
- `harness.py` — CV/tuning/metrics harness using the registry + fold module.
- `run_regression.py` / `run_classification.py` — orchestration entry points.
- Summary table + comparison plots (headline config; plain vs. donor gap; flat vs. nested; corrected vs. uncorrected).

## Verify
- Per-fold fitting cell counts stay large; fold assignments reproducible; classification metrics only where both classes present; repeated-CV reported per-repeat not pooled; regression predictions checked against [0,1].
