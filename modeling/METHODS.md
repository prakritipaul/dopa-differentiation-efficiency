# Methods

How day-52 dopaminergic yield is predicted from day-11 and day-30 single-cell
data, and why each choice was made.

Results are in [`../FINDINGS.md`](../FINDINGS.md); output files in
[`results/README.md`](results/README.md); the audit trail and the
label-variant machinery in [`README.md`](README.md).

---

## 1. Setup

### Cohort

| | |
|---|---|
| Lines | **136**, from 157 `(line, pool)` combinations |
| Donors | 20 — 6.8 lines per donor, largest 16; 134 of 136 lines share a donor |
| Pools | 10 differentiation runs |
| Inclusion | ≥ 10 untreated cells at D11 **and** D30 **and** D52 |

The 157 combinations are the same set at both timepoints, but a line's cell
count in a given pool is not: D11 and D30 per-combination counts correlate at
0.748. Cells are not tracked between timepoints — each is a separate harvest.

### Outcome

**`DA / all D52 cells`, untreated cells only**, binarised at 0.2 →
**61 success / 75 failure**.

Why not the published `(DA + Sert) / all D52 cells`: that metric counts
rotenone-treated cells, and rotenone preferentially damages the dopaminergic
neurons in its own numerator; and it sums two lineages that do not move together
in this cohort. The argument is in [`../FINDINGS.md`](../FINDINGS.md) §1 and is
not repeated here.

Both tasks are run on this outcome: **regression** on the continuous fraction,
**classification** on the binarised label.

### Features

Two blocks, built separately at D11 and D30.

| | D11 | D30 |
|---|---|---|
| Annotated cell types | 3 | 7 |
| Proportions entering a fit | 2 | 6 |
| Expression | PC1–PC10 from 2,000 highly variable genes | same |
| Depth-outlier pool excluded from fitting | `pool11` | `pool5` |

**Why cell-type proportions.** Composition is what a differentiation is scored
on, and each coordinate names a cell type a biologist can act on.

**Why principal components alongside them.** The annotation is a discrete
summary; the PCs carry graded state the labels discard. Ten is the number
computed; how many enter a fit is tuned (§2).

**Why one proportion is held out.** Proportions sum to exactly 1 per line, so
fitting all of them alongside an intercept makes the design matrix
rank-deficient. `phat_P_FPP` is the reference category at both timepoints and
never enters a fit — which is why D11 contributes two proportions and D30 six.
`harness.proportion_cols` derives this from the feature table and raises if the
reference column is missing, rather than trusting a per-timepoint flag.

### How features are built, per fold

- **HVG selection and PCA are fit on training-line cells only**, then *all*
  cells are projected into that basis. No held-out line influences the basis it
  is scored in.
- **Depth-outlier pools are excluded from fitting but still projected**, so no
  line is lost. The rule is per timepoint, not global: a pool is an outlier if
  its median raw UMI per cell falls below one third of the median of pool
  medians at that timepoint. That gives `pool11` at D11 (ratio 0.137, next
  lowest 0.747) and `pool5` at D30 (0.123, next 0.423). **It is not the same
  pool at both timepoints.**
- **PCs are standardised, then averaged pool-then-line** — mean per
  `(line, pool)`, then an unweighted mean of those pool means, so a line
  differentiated twice is not counted twice.

## 2. Model design

Four models, two per task, all on the same features and folds.

| model | task | regularisation | grid |
|---|---|---|---|
| `ridge` | regression | L2 | `alpha` ∈ {0.01, 0.1, 1, 10, 100} |
| `lasso` | regression | L1 | `alpha` ∈ {0.001, 0.01, 0.1, 1, 10} |
| `logistic_l2` | classification | L2 | `C` ∈ {0.01, 0.1, 1, 10, 100} |
| `logistic_l1` | classification | L1 | `C` ∈ {0.01, 0.1, 1, 10, 100} |

**Why regularised linear models.** 136 lines against 8–16 correlated features.
Regularisation is what keeps that from overfitting, and linear coefficients stay
readable — which matters because §5 reads them.

**Why both L1 and L2.** L2 for the reported fits; L1 as a sparsity check, since
`selection_frequency` — how often a feature survives across folds — is only
meaningful for a model that can zero a coefficient.

**Why lasso's grid sits a decade lower.** `C` is an inverse penalty and `alpha`
a direct one, and L1 on standardised features zeroes everything at `alpha` where
ridge is still fitting.

### Hyperparameters

Tuned jointly over two axes:

- **`k`** — number of leading PCs, 0 to 10. `k = 0` is the proportions-only
  model. PCA components are hierarchical, so this truncates columns from one
  fit per fold rather than needing a separate PCA per `k`.
- **Regularisation strength** — the five values above.

**11 × 5 = 55 configurations per model.**

Both logistic models use `class_weight="balanced"`; every estimator is pinned to
a fixed seed, because `saga` is stochastic and without it `logistic_l1`
coefficients drift between runs. Features are standardised on training folds
only.

## 3. How performance is measured

Hard class calls use p = 0.5 — a different threshold from the 0.2 that defines
the outcome.

*Classification*

| metric | why it is computed |
|---|---|
| **ROC-AUC** | threshold-free ranking — can lines be sorted by eventual yield |
| **PR-AUC** | ranking under 61/75 imbalance, where ROC-AUC is the more optimistic |
| **Balanced accuracy** | accuracy at p = 0.5, corrected for unequal class sizes |
| **Sensitivity** / **specificity** | the two error costs kept apart — discarding a line that would have worked, vs. carrying a failing one |
| **F1** | one number over precision and recall at p = 0.5 |
| **Brier** | calibration — whether the probability is usable, not just its order |

*Regression*

| metric | why it is computed |
|---|---|
| **R²** | share of across-line variance in the D52 dopaminergic fraction explained |
| **MAE** | error in DA-fraction units, robust to the long right tail |
| **RMSE** | same units, penalising the large misses MAE forgives |
| **Out-of-range** | fraction of predictions outside [0, 1] — a linear model can predict an impossible fraction |

`log_loss` and `accuracy` are computed in one place only: the paired test of the
full model against `phat_DA` alone, where ROC-AUC saturates near 0.95 and cannot
settle the question.

**The reported numbers are ROC-AUC for classification and R² for regression.**
Everything else is a check on them — chiefly Brier, which catches a model that
ranks well and is still badly calibrated.

**Reported is not the same as tuned.** The nested inner folds select on
**PR-AUC** (classification) and **negative MAE** (regression). MAE rather than
R² because R² is a ratio against the fold's own variance, which moves between
small folds. The flat one-SE configuration pick uses `roc_auc_mean` and
`mae_mean`.

**Aggregation.** Metrics are computed **per repeat**, over that repeat's
complete set of out-of-fold predictions, then averaged ± SD across repeats —
never pooled across repeats. Within one repeat the folds cover every line
exactly once, which is what makes the per-repeat metric well defined.

## 4. Cross-validation

Four schemes, run in full. Each answers a different question.

| scheme | folds | what it answers |
|---|---|---|
| **plain 5-fold** | 5 × 10 repeats | baseline. Sibling lines from one donor may land on both sides of a split |
| **donor-grouped 5-fold** | 5 × 10 repeats | no donor spans train and test. The strictest of the repeated schemes |
| **leave-one-line-out** | 136 | maximum training data per fit; the optimistic end |
| **leave-one-donor-out** | 20 | whole donors held out |

- Donor-grouped uses `StratifiedGroupKFold` on donor, stratified by the binary
  label; folds are uneven because donors contribute between 1 and 16 lines.
- LOCO and LODO are single passes by construction, so they carry no SD. A
  one-line fold cannot support a per-fold metric, so LOCO pools its predictions.
- Reporting all four is the point: agreement across schemes is what shows a
  result does not depend on one split.

**Tuning is nested throughout** — an inner 3-fold on training lines only, re-run
inside every outer fold, so no test line influences the configuration it is
scored under.

**Configuration is chosen by flat CV + one-SE:** the simplest configuration
within one standard error of the grid's best, on donor-grouped folds. Simplest
means fewer PCs first, then stronger regularisation. The tolerance is
SD / √n_repeats — a standard error, not a standard deviation.

**Performance is reported from nested CV, and both numbers are shown**, so the
gap between the score a configuration was selected on and its out-of-fold score
is visible rather than hidden.

The two do different jobs and are not interchangeable: nested CV produces *a
number*, flat CV produces *a choice*. Nested scores the whole selection
procedure — each outer fold re-runs the tuning and may land on a different `k`.

**Both timepoints use identical fold assignments**, verified to agree on every
fold's train/test membership, so the D11-vs-D30 comparison is paired.

### Leakage safety

- HVG selection and the PCA are refit inside every fold on training lines only.
- The depth-outlier exclusion is a frozen rule from EDA, never tuned on model
  performance.
- Fold assignments are persisted, so the same splits are reused across
  timepoints and reruns.
- One limitation that cannot be fixed here: the cell-type labels were assigned
  upstream, over all cells at once. Nothing in this pipeline can undo that.

## 5. How the most predictive features are picked

Seven measures, computed into one row-per-feature table
(`results/feature_importance_table_*.csv`). They answer different questions, and
a feature can rank high on one and near zero on another.

| measure | what it computes | what it tells us |
|---|---|---|
| **Univariate Spearman ρ** | feature against outcome, no model, all 136 lines | direction, and a model-free sanity check; the only column the reference proportion can appear in |
| **Standardised coefficient** | effect per 1 SD with the other features held fixed | a *conditional* pull — what the feature adds once its correlates are accounted for, which can oppose its standalone ρ |
| **Selection frequency** | fraction of folds where L1 keeps a non-zero coefficient | stability — whether the sparse model keeps choosing it |
| **LOCO Δ** (refit) | refit without the feature on the same folds, `score_with − score_without` | **necessity.** A feature that is real but redundant scores ≈ 0, because the refit recovers it from its correlates |
| **Grouped permutation Δ** (no refit) | shuffle the column in the test matrix of the already-fitted model | **reliance** — how hard the fitted model leans on it |
| **Mean \|SHAP\|** | mean \|coefficient × standardised value\| | **magnitude on a common scale** — how far the feature actually moves a prediction |
| **Pool η²** | share of the feature's across-line variance associated with which of the 10 runs the line went through | **run dependence, not importance.** High η² does not prove artefact; low η² does not prove portability |

**Features are ranked by mean |SHAP|.** It is the only measure that puts every
feature in the same units — how far this feature moves a prediction — so it
orders composition features and expression components on one scale. The others
are checks: LOCO Δ catches a feature that looks important but is redundant,
permutation Δ distinguishes reliance from necessity, and pool η² flags a feature
that may be reading the batch.

Two details that change the numbers:

- **Proportions are permuted and dropped as a block**, never individually.
  They are compositional; shuffling one alone implies a remainder that cannot
  exist.
- **Mean |SHAP| is computed in closed form**, `mean |coef × z|`, which is exact
  for a linear model and avoids a dependency on the `shap` package.

## 6. Reproducibility

- Every estimator and every split is seeded; fold assignments are persisted to
  `fold_data/` and reused rather than regenerated.
- The expression files are not in the repository. Point one environment
  variable at them — no source edits:

  ```bash
  export JERBER_DATA_DIR=/path/to/your/data   # holds day11.h5, day30.h5, day52.h5
  ```

- Run the pipeline:

  ```bash
  uv run python modeling/run_feature_extraction.py   # per-fold HVG + PCA
  uv run python modeling/run_regression.py
  uv run python modeling/run_classification.py
  uv run python modeling/feature_importance.py
  ```

- Feature extraction writes incrementally and resumes from what is already on
  disk; `verify_fold_features.py` is the gate that must pass before modelling.
- Tests: `uv run pytest modeling/ -m "not slow"` runs without the h5 files, on a
  miniature dataset built in the same on-disk shape. Run the full suite before
  treating any result as final.
