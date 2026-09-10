# pluricon-prototype

Predicting **day-52 dopaminergic neuron differentiation efficiency** from
**day-11** single-cell features, across 138 iPSC lines.

Practical question: 11 days into a midbrain differentiation, can you tell
whether the run will succeed — early enough to abort a bad one?

**Short answer: yes.** ROC-AUC 0.947, and 100% precision at 65% recall for
detecting failures, under never-seen-donor cross-validation.

Full write-ups: **[FINDINGS.md](FINDINGS.md)** (results + literature
context) · [`modeling/README.md`](modeling/README.md) (methods) ·
[`modeling/docs/PCA_interpretation.md`](modeling/docs/PCA_interpretation.md)
(what the PCs mean).

---

## Data

[Jerber et al. 2021](https://www.nature.com/articles/s41588-021-00801-6),
*Nat Genet* 53:304–312 — scRNA-seq of iPSC lines differentiated toward
midbrain dopaminergic fate, multiplexed into pools, sampled at D11, D30
and D52.

| | |
|---|---|
| Lines used | **138** (≥10 cells at D11 *and* D30 *and* D52), from 159 `(line, pool)` combos |
| Predictors | D11 only — 253,381 cells |
| Outcome | `diff_efficiency` = fraction of a line's D52 cells that are DA or Sert neurons |
| Binary label | success = efficiency ≥ 0.2 (the authors' threshold) → **96 success / 42 failure** |
| Outcome shape | bimodal: failures 0–0.185, successes 0.221–0.92, real gap at the threshold |

## Features (D11 only)

1. **Cell-type proportions** — `phat_FPP`, `phat_NB`, `phat_P_FPP`. These
   sum to 1, so `phat_P_FPP` is the implicit reference and only two enter
   any fit.
2. **PC1–PC10** — 2,000 HVGs → standardised → PCA. Per-cell coordinates
   averaged **pool-then-line** (mean per `(line, pool)`, then unweighted
   mean of pool means).

Both are computed **per fold**: HVG selection and the PCA are fit on
training-line cells only, then *all* cells are projected — so no held-out
line influences the basis it is scored in. `pool11` is excluded from every
fit (≈1,900 UMI/cell vs 10,000–18,000 elsewhere) but still projected, so
no line is lost.

## Cross-validation

| | |
|---|---|
| Primary scheme | **donor-grouped** 5-fold, 10 repeats, stratified by the binary label (136/138 lines share a donor with another line) |
| Also run | plain 5-fold, LOCO (138 folds), LODO (20 folds) |
| Tuning | **nested** — inner 3-fold on training lines only; inner metric PR-AUC (classification) / MAE (regression) |
| Metrics | computed **per repeat**, then averaged — never pooled across repeats |
| Total | 258 folds, PCA refit for each |

Model per task was **pre-registered before looking at results**:
donor-grouped + nested + ridge / logistic_l2.

## Models and hyperparameters

Tuned jointly over `k` (number of PCs) and regularisation strength — 55
configurations per model.

| model | task | grid |
|---|---|---|
| ridge *(headline)* | regression | `alpha` ∈ {0.01, 0.1, 1, 10, 100} |
| lasso | regression | `alpha` ∈ {0.001, 0.01, 0.1, 1, 10} |
| logistic_l2 *(headline)* | classification | `C` ∈ {0.01, 0.1, 1, 10, 100} |
| logistic_l1 | classification | `C` ∈ {0.01, 0.1, 1, 10, 100} |

`k` ∈ 0–10 for all. Both logistic models use `class_weight="balanced"`
and a fixed seed. Features standardised on training folds only.

---

## Results

Nested CV, mean ± SD across 10 repeats. `plain` ignores donors;
`donor-grouped` never lets a donor appear in both train and test.

| task | model | metric | plain 5-fold | donor-grouped 5-fold |
|---|---|---|---|---|
| classification | **logistic_l2** | ROC-AUC | 0.940 ± 0.012 | **0.947 ± 0.007** |
| classification | logistic_l1 | ROC-AUC | 0.948 ± 0.011 | 0.951 ± 0.008 |
| regression | **ridge** | R² | 0.676 ± 0.013 | **0.653 ± 0.021** |
| regression | lasso | R² | 0.680 ± 0.016 | 0.669 ± 0.015 |

Classification, headline model, donor-grouped: PR-AUC 0.968, balanced
accuracy 0.913, sensitivity 0.940, specificity 0.886, Brier 0.097.

**Best model: L2-logistic on the binary task.** The L1 variants edge it
but well within one SD, and were not pre-registered.

**Donor grouping costs regression a little and classification nothing** —
ridge −0.023 R², lasso −0.011, while both classification models are flat
or slightly better. LOCO/LODO agree (0.682/0.665 R², 0.942/0.963 AUC), so
donor leakage is real but small, about 0.012 R².

⚠️ **The regression R² mostly reflects separating the two outcome clumps,
not fine-grained accuracy.** Within-success R² is only ~0.17–0.22 and
within-failure R² is strongly negative. Use this for success/failure
calls, not as an efficiency estimate.

**The PCA fitting population barely matters.** Refitting the whole
pipeline on cells from the 138 study lines only shifts donor-grouped
results by +0.014 R² (ridge), +0.004 (lasso), ±0.001 AUC — all inside one
SD. It does roughly halve regression's across-repeat SD (0.021 → 0.012),
so it is somewhat more stable, but no conclusion changes.

## Which features matter

**`phat_NB`, the D11 neuroblast fraction, is the standout.** More
neuroblasts at D11 → *worse* D52 yield.

| evidence | value |
|---|---|
| Univariate Spearman ρ | **−0.73** — strongest of anything measured |
| LOCO Δ (classification) | +0.051 AUC, largest single feature |
| SHAP (classification) | 1.19 vs PC2's 0.51 |
| Pool η² | **0.04** — the only technically clean feature |

Dropping all proportions costs **+0.064 AUC** — the largest effect in
either table — but costs regression nothing. On the restricted PCA basis
the modal model selects **zero PCs** and still reaches AUC 0.950.

**PC1 is a proliferation axis** (ρ = +0.563; enriched for the Tirosh/Seurat
G2/M set at p = 2.3×10⁻²⁶; tracks proliferating-FPP fraction at ρ = +0.798).
More proliferative D11 cultures differentiate better.

## Technical confounders — read before interpreting any PC

Pool identity explains much of the variance in most PC features:

| feature | pool η² |
|---|---|
| `phat_NB` | 0.04 |
| PC1 / PC3 / PC4 | 0.45–0.50 |
| **PC2** | **0.76** — *and the most predictive PC* |
| PC6 / PC8 / PC9 | 0.91–0.96 — near-pure batch |

So **PC2 is model-predictive but technically entangled**; it may not
transfer to new pools. Separately, PC2 and PC3 are near-tied in variance
and **rotate between PCA fits** — `PCn` is not a stable label. PC1 is the
exception (replicates at r = 0.9995).

---

## Repository

```
FINDINGS.md      results + what is new vs. prior work
metadata_eda/    EDA outputs: cohort/ pca/ proportions/ qc/ technical/ plots/
modeling/        the pipeline (see modeling/README.md)
0NN_*.py         numbered EDA scripts, run in order
```

## Setup

```bash
uv sync
uv run pytest modeling/ -m "not slow"     # 51 tests
uv run python -m modeling.run_classification
```

Python 3.11, [uv](https://docs.astral.sh/uv/). Expression data (`day11.h5`,
`day30.h5`, `day52.h5`) is not in the repo; set the paths in
`modeling/features.py`.
