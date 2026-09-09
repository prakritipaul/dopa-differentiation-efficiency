# Feature importance: final table

One table per task (regression, classification):

| Feature | Univariate | Coefficient | Drop-one Δ | SHAP | Technical covariate |
|---|---|---|---|---|---|

- **Feature**: PC1..PC10 (up to winning k) + phat_FPP, phat_NB, phat_P_FPP.
- **Univariate**: Spearman corr with the outcome, no model.
- **Coefficient**: standardized coef from full-data fit at the winning `(k, params)`, picked via **one-SE rule** (not raw argmax). `phat_P_FPP` excluded from the fit (compositional — sums to 1.0 with FPP/NB; kept implicit as `1-FPP-NB`), shows "reference" in this column. Flag features **regularized to zero** (Lasso/L1) explicitly.
- **Drop-one Δ (LOCO)**: refit winning config minus one feature, same CV folds, report Δ in headline metric, computed as a **paired per-fold difference** (`metric_without − metric_with` within each fold, then averaged — not two independently-averaged CV scores subtracted at the end; reduces noise at n=138). PCs each tested individually. Proportions: **only FPP and NB are tested individually** (both are explicit features in the reference-coded baseline, so dropping either is well-defined); `phat_P_FPP`'s row instead shows the **grouped** delta (all 3 proportions dropped together), since P_FPP was never an explicit member of the baseline to begin with (it's the implicit reference) — there is no valid "with P_FPP" baseline to drop it from (Codex: a 3-proportion-explicit model is rank-deficient).
- **SHAP**: `shap.LinearExplainer` (exact for linear models, cheap) — `uv add shap`. Mean |SHAP| per feature.
- **Technical covariate**: eta² by pool, recomputed fresh on the *actual* model features (005's pool11-excluded PCA + proportions) — not reusing 008's numbers (different, uncorrected PCA basis). Bucketed High/Moderate/Low.

**Permutation importance** (separate check): out-of-fold, using the per-fold refit models already built, also as **paired per-fold deltas**. PCs permuted individually; the 3 proportions permuted **as one group** (shuffle which line's whole 3-vector goes to which row) — unlike LOCO, permuting a single proportion alone while its partners stay put creates impossible rows that don't sum to 1.

## Build order
1. Fix `feature_importance.py`: compositional handling (coefficient-fit only), one-SE rule, Lasso zero/selection-frequency reporting.
2. Univariate (Spearman).
3. LOCO: FPP/NB + all PCs individually, P_FPP as the grouped-proportions delta; paired per-fold deltas throughout.
4. Grouped permutation importance (proportions grouped, PCs individual), paired per-fold deltas.
5. SHAP (`uv add shap`, `LinearExplainer`).
6. Technical covariate assoc. (reuse `eta_squared_by_pool` from `008` via importlib, applied to actual model features).
7. Assemble + save both tables.

## Data
Build/validate against existing first-pass results (`fold_features_D11.csv` etc. — all inputs already exist). Label preliminary. Re-run against `fold_features_D11_full.csv` once the background 258-fold job finishes.

## Verify
- Compositional handling doesn't crash; P_FPP row reads "reference."
- Zeroed Lasso coefficients read clearly, not NaN.
- LOCO/permutation signs make sense (dropping a useful feature hurts performance).
- Technical-covariate numbers roughly sane vs. `008` (different basis, not exact match expected).

## Provenance
Design discussion: own reasoning + two rounds of independent Codex review
(first pass: overall methodology critique -- conditioning on one
data-selected config, fold-SD not a real standard error, Lasso needing
selection-frequency, the compositional proportions issue, recommended
permutation importance; second pass, after the user's own catch that
individual-LOCO-for-proportions seemed off: Codex confirmed FPP/NB can be
tested individually but P_FPP cannot, for a different reason -- it's the
implicit reference, never an explicit baseline feature, so there's no
valid "with P_FPP" model to drop it from; also recommended paired
per-fold deltas over independently-averaged CV scores).
