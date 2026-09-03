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

## First-pass results (real data, ran successfully)

Headline (donor_grouped, nested CV, **no pool-correction** -- see
"Pool-correction dropped" below):

| task | model | key metrics |
|---|---|---|
| regression | lasso | MAE=0.135, RMSE=0.172, **R2=0.667** |
| regression | ridge | MAE=0.137, RMSE=0.176, R2=0.652 |
| classification | logistic_l1 | ROC-AUC=0.949, PR-AUC=0.970, balanced_acc=0.913 |
| classification | logistic_l2 | ROC-AUC=0.944, PR-AUC=0.966, balanced_acc=0.911 |

D11 features are clearly predictive of D52 differentiation efficiency,
even under the strict donor-grouped (never-seen-donor) test.

**Plain vs. donor_grouped gap is small** (regression R2 0.685->0.667,
classification barely moves at all) -- reassuring: the model isn't
leaning heavily on sibling-line/donor leakage, performance holds up close
to as well when donors are held out entirely.

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
