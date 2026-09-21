# dopa-differentiation-efficiency

**Predicting day-52 dopaminergic neuron yield from earlier single-cell
snapshots**, across 136 iPSC lines.

Practical question: partway through a midbrain differentiation, can you tell
whether it will produce dopaminergic neurons — early enough to act on it?

**Short answer: yes at D11 (ROC-AUC 0.906), and by D30 the outcome is largely
already determined.**

Two day-11 features carry that prediction, and they are not equally
interpretable:

- **More neuroblasts at day 11 → fewer dopaminergic neurons at day 52**
  (`phat_NB`, ρ = −0.554). Lines that have already started making neurons by
  day 11 do worse. This is the cleanest signal in the model: almost none of its
  variation across lines tracks which differentiation run a line went through
  (pool η² = 0.035), so it is the one most likely to mean the same thing in
  another lab.
- **`PC3` contributes the most of any single feature** (mean |SHAP| 1.40) and is
  an expression axis loaded with ribosome biogenesis, nucleotide and serine
  synthesis, and S-phase genes — a cell-growth programme. **But it is also the
  most run-associated feature in the model** (pool η² = 0.832), which is exactly
  what a growth-rate axis would look like if it were tracking culture conditions
  rather than the cell line. It predicts well within this study; it should not
  be read as biology. See [FINDINGS.md](FINDINGS.md#4-day-11--which-features-predict-and-what-they-suggest).

Full write-ups: **[FINDINGS.md](FINDINGS.md)** (results + literature context) ·
[`modeling/METHODS.md`](modeling/METHODS.md) (methods) ·
[`modeling/results/README.md`](modeling/results/README.md) (output files)

---

## Data

[Jerber et al. 2021](https://www.nature.com/articles/s41588-021-00801-6),
*Nat Genet* 53:304–312 — scRNA-seq of iPSC lines differentiated toward midbrain
dopaminergic fate, multiplexed into pools, sampled at D11, D30 and D52.

| | |
|---|---|
| Lines used | **136** (≥10 untreated cells at D11 *and* D30 *and* D52), from 157 `(line, pool)` combos, 20 donors |
| Predictors | **D11** (253,381 cells) and **D30** (250,923 cells), modelled separately |
| Outcome | `DA / all D52 cells` — the dopaminergic fraction at day 52 [†](#-on-the-outcome-definition) |
| Binary label | success = ≥ 0.2 → **61 success / 75 failure** |
| Outcome shape | unimodal, median 0.165, range 0.010–0.742 |

## Cell types and notation

`phat_X` denotes the **estimated proportion** of cell type *X* in a line — the
"hat" marks it as an estimate from a finite sample of cells, not a known
quantity. Proportions sum to 1 within a timepoint.

| label | name | timepoint |
|---|---|---|
| `FPP` | **floor-plate progenitor** — the on-target progenitor that gives rise to midbrain dopaminergic neurons | D11, D30, D52 |
| `P_FPP` | **proliferating floor-plate progenitor** — the same cell type, still cycling | D11, D30, D52 |
| `NB` | **neuroblast** — early post-mitotic neuronal precursor | D11 |
| `DA` | **dopaminergic neuron** — the target cell type | D30, D52 |
| `Sert` | **serotonergic neuron** — a floor-plate-derived neuron of different rostro-caudal identity | D30, D52 |
| `Epen1` | **ependymal-like 1** — ciliated / choroid-plexus-like off-target fate | D30, D52 |
| `U_Neur1`, `U_Neur2` | **unassigned neuron 1 / 2** — post-mitotic neurons with no lineage markers; off-target | D30 |
| `Astro` | **astrocyte** — off-target glial fate | D52 |

## Features

| | D11 | D30 |
|---|---|---|
| annotated cell types | 3 — floor-plate progenitor, neuroblast, proliferating FPP | 7 — adds dopaminergic, serotonergic, ependymal-like, 2 unassigned neuron types |
| free proportions in a fit | 2 | 6 |
| expression | PC1–PC10 from 2,000 HVGs | same |
| depth-outlier pool excluded from fitting | `pool11` | `pool5` |

- Proportions sum to 1, so `phat_P_FPP` is the implicit reference and does not
  enter a fit.
- PCs are standardised, then averaged **pool-then-line** — mean per
  `(line, pool)`, then unweighted mean of pool means.
- Both feature types are computed **per fold**. HVG selection and the PCA are
  fit on training-line cells only, then *all* cells are projected, so no
  held-out line influences the basis it is scored in.
- The depth-outlier pool is excluded from fitting but still projected, so no
  line is lost. **It is not the same pool at both timepoints.**

## Cross-validation

| | |
|---|---|
| Primary scheme | **donor-grouped** 5-fold, 10 repeats, stratified by the binary label (134/136 lines share a donor) |
| Also run | plain 5-fold, LOCO (136 folds), LODO (20 folds) |
| Tuning | **nested** — inner 3-fold on training lines only |
| Metrics | computed **per repeat**, then averaged — never pooled across repeats |
| Total | 256 folds, PCA refit for each, per timepoint |

- Configuration is chosen by **flat CV + one-SE** — the simplest config within
  one standard error of the grid's best, on donor-grouped folds.
- **Performance is reported from nested CV.** Both numbers are shown: the flat
  score is the selection criterion, not a result.
- Both timepoints use **identical fold assignments**, verified to agree on every
  fold's train/test membership, so the D11-vs-D30 comparison is paired.

## Models and hyperparameters

Tuned jointly over `k` (number of PCs) and regularisation — 55 configurations
per model.

| model | task | grid |
|---|---|---|
| ridge *(reported)* | regression | `alpha` ∈ {0.01, 0.1, 1, 10, 100} |
| lasso | regression | `alpha` ∈ {0.001, 0.01, 0.1, 1, 10} |
| logistic_l2 *(reported)* | classification | `C` ∈ {0.01, 0.1, 1, 10, 100} |
| logistic_l1 | classification | `C` ∈ {0.01, 0.1, 1, 10, 100} |

`k` ∈ 0–10 for all. Both logistic models use `class_weight="balanced"` and a
fixed seed. Features standardised on training folds only.

## Metrics

Every metric the harness computes. Hard class calls use p = 0.5, which is a
different threshold from the 0.2 that defines the *outcome*.

| metric | task | why it is reported |
|---|---|---|
| **ROC-AUC** | classification | threshold-free ranking — can lines be sorted by eventual yield |
| **PR-AUC** | classification | ranking under 61/75 class imbalance, where ROC-AUC is the more optimistic of the two |
| **Balanced accuracy** | classification | accuracy at p = 0.5, corrected for unequal class sizes |
| **Sensitivity** / **specificity** | classification | the two error costs kept apart — discarding a line that would have worked, vs. carrying a failing one for 41 more days |
| **F1** | classification | one number over precision and recall at p = 0.5 |
| **Brier** | classification | calibration — whether the probability itself is usable, not just its rank order |
| **R²** | regression | share of across-line variance in the D52 dopaminergic fraction explained |
| **MAE** | regression | error in DA-fraction units, robust to the long right tail |
| **RMSE** | regression | same units, penalising the large misses that MAE forgives |
| **Out-of-range** | regression | fraction of predictions outside [0, 1] — a linear model can predict an impossible fraction |
| **Log loss**, **accuracy** | classification, comparison only | used in the paired `phat_DA`-vs-model test, where ROC-AUC saturates near 0.95 and cannot settle the question |

Two things the table does not say:

- **Reported is not the same as tuned.** The nested inner folds select on
  **PR-AUC** (classification) and **negative MAE** (regression) — not on the
  ROC-AUC and R² quoted in the results. The flat-CV one-SE configuration pick
  uses `roc_auc_mean` and `mae_mean`.
- **Metrics are computed per repeat**, over that repeat's complete set of
  out-of-fold predictions, then averaged ± SD across repeats — never pooled
  across repeats. LOCO and LODO are single repeats by construction, so they
  carry no SD.

## Feature importance

Six measures, plus one association statistic that is not an importance measure
at all. They answer different questions — a feature can rank high on one and
near zero on another.

| measure | what it computes | what it tells us |
|---|---|---|
| **Univariate Spearman ρ** | feature against outcome, no model, all 136 lines | direction, and a model-free sanity check; the only column the reference proportion can appear in |
| **Standardised coefficient** | effect per 1 SD with the other features held fixed | a *conditional* pull, not a marginal one — what the feature adds once its correlates are already accounted for, which can point opposite to its standalone ρ |
| **Selection frequency** | fraction of folds where L1 keeps a non-zero coefficient | stability — whether the sparse model keeps choosing the feature, or it survived on one split |
| **LOCO Δ** (refit) | refit without the feature on the same folds, `score_with − score_without` | **necessity.** A feature that is real but redundant scores ≈ 0, because the refit recovers it from its correlates |
| **Grouped permutation Δ** (no refit) | shuffle the column in the test matrix of the already-fitted model | **reliance** — how hard the fitted model leans on it. Proportions are permuted as one block, since they are compositional and shuffling one alone implies an impossible remainder |
| **Mean \|SHAP\|** | mean \|coefficient × standardised value\| | **magnitude on a common scale** — how far the feature actually moves a prediction, comparable across features |
| **Pool η²** | share of the feature's across-line variance associated with which of the 10 differentiation runs the line went through | **run dependence, not importance.** High η² does not prove a technical artefact, and low η² does not prove portability — see [Caveats](#caveats) |

---

## Results

Donor-grouped CV, 10 repeats.

- **Nested** — the reported out-of-fold performance.
- **Flat** — the score on the grid the configuration was selected from, shown
  so the selection gap is visible rather than hidden.

*Selected configuration per task · donor-grouped, 10 repeats*

| Model | Configuration | Features | Nested | Flat |
|---|---|---|---|---|
| D11 classification | `logistic_l2`, C = 1.0, k = 4 | `phat_FPP`, `phat_NB` + PC1–PC4 | **0.906** | 0.918 |
| D11 regression | `ridge`, α = 10.0, k = 7 | same 2 proportions + PC1–PC7 | **0.503** | 0.506 |
| D30 classification | `logistic_l2`, C = 1.0, k = 6 | 6 proportions + PC1–PC6 | **0.945** | 0.956 |
| D30 regression | `ridge`, α = 10.0, **k = 0** | 6 proportions, **no PCs** | **0.802** | 0.812 |

- Classification is ROC-AUC, regression R².
- SD across the 10 repeats: 0.008, 0.015, 0.008, 0.009 — in table order.
- `phat_P_FPP` is the reference proportion and never enters a fit. Hence two
  proportions at D11, six at D30.
- The nested-to-flat gap is small (0.003–0.012) but **positive in all four
  cases** — what selection bias looks like, not what noise does.
- Full grids in [FINDINGS.md](FINDINGS.md#appendix-a--full-performance-tables).

### Day 11 — a forecast, 41 days before the readout

**The result that carries the project.**
ROC-AUC 0.906 ± 0.008 · R² 0.503 · Brier 0.125 · MAE 0.099.

1. **No dopaminergic neurons exist at D11.** The model forecasts a fate; it does
   not count one. The three annotated types are FPP, NB and P_FPP.
2. **Four cross-validation schemes agree.** 0.902 plain · 0.906 donor-grouped ·
   0.912 leave-one-line-out · 0.908 leave-one-donor-out. Nothing rests on one
   fortunate split.
3. **Premature neurogenesis predicts failure.** `phat_NB` is the strongest and
   most technically clean composition feature — ρ = −0.554, pool η² = 0.035,
   mean |SHAP| 1.05. More D11 neuroblasts, fewer D52 dopaminergic neurons.
   Lines that leave the cycling floor-plate pool early never build enough of it.
4. **Expression agrees, independently.** `PC2` (pool η² 0.114, ρ = −0.626) is a
   proneural axis — `NEUROD1 NHLH1 DLL3` against `HES1 OTX2` — read from genes
   alone, never from cell-type calls. Within a PCA basis it tracks the annotated
   neuroblast fraction at ρ ≈ +0.8. Counting neuroblasts and reading a
   neurogenic programme converge on one answer.

   `PC2` is a slot, not a fixed axis: it trades places with `PC3` in a minority
   of folds ([Caveats](#caveats)).
5. **One caveat kept in view.** The most influential D11 feature is `PC3` (mean
   |SHAP| 1.40, LOCO Δ +0.064) and also the most run-associated (pool η² 0.832).
   It improves prediction on held-out lines. Whether it ports to a
   differentiation run done elsewhere is untested.

### Day 30 — mostly already decided

**The higher score is not the better result.**
ROC-AUC 0.945 ± 0.008 · R² 0.802.

1. **D30's advantage is mostly definitional.** Its annotated cell types already
   *include* DA, so the outcome partly exists in the predictors — construct
   overlap, not leakage. The right reading is *how much of the D52 phenotype is
   already established by D30*: ~48% of cells have arrived, ~31% are still
   uncommitted progenitors, ~20% are off-target.
2. **Expression adds nothing beyond cell-type composition.** Three of four
   models select **zero PCs**. Grouped permutation Δ: 0.360 for the proportions
   against 0.018 for the best PC.
3. **One raw column out-ranks the full model — but is badly calibrated.**

   | D30 | `phat_DA` alone | full model |
   |---|---|---|
   | ROC-AUC | **0.956** | 0.939 |
   | log loss | 0.571 | **0.328** |

   To *rank* lines, one number suffices. For trustworthy *probabilities*, the
   model earns its keep. AUC alone would mislead either way.

## Repository

```
eda/              numbered 0NN_*.py scripts, run in order — the EDA phase
modeling/         the pipeline
  METHODS.md      features, models, CV, metrics, feature importance
  README.md       audit trail, decision records, label-variant machinery
  docs/           PCA interpretation, cell-type markers, audit ledger
  results/        every output table (see modeling/results/README.md)
  tests/          pytest suite, incl. a synthetic-data smoke test
metadata_eda/     EDA outputs: cohort/ pca/ proportions/ qc/ technical/ plots/
tools/            reference-PDF builder
```

## Setup

```bash
uv sync
uv run pytest modeling/ -q      # passes on a fresh clone, no data needed
```

The expression files (`day11.h5`, `day30.h5`, `day52.h5`) are multi-GB and not
in the repo.

- **You do not need them to verify the code.** The suite builds a miniature
  dataset in the same on-disk shape and drives the real pipeline over it in
  ~3 seconds.
- The two tests that genuinely require the real files are marked
  `requires_data` and skip themselves when the files are absent. A fresh clone
  gets a green run; the same tests still execute where the data is present.

To run on the real data, point one environment variable at the directory
holding `day11.h5`, `day30.h5` and `day52.h5` — no source edits:

```bash
export JERBER_DATA_DIR=/path/to/your/data
```

Data: [E-MTAB-10018](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-10018).
Pipeline commands in [`modeling/METHODS.md`](modeling/METHODS.md).

Python 3.11, [uv](https://docs.astral.sh/uv/).

---

## † On the outcome definition

The published metric is `(DA + Sert) / all D52 cells`. This project uses
**`DA / all D52 cells`, untreated cells only** — a deliberate divergence, for
two reasons:

1. **The published metric counts rotenone-treated cells.** D52 alone contains
   219,238 `ROT` against 303,856 `NONE`, interleaved so every qualifying combo
   has both. Rotenone inhibits mitochondrial complex I and preferentially
   damages dopaminergic neurons — the cells in its own numerator. Verified
   against the authors' notebook: no treatment filter anywhere.
2. **It sums two lineages that track independently.** D30 DA predicts D52 DA at
   ρ = +0.932; D30 Sert predicts D52 DA at Pearson **+0.057**.

- An earlier phase reproduced the authors' definition exactly — 138 lines, D11
  only, ROC-AUC 0.947 — which is what made comparison to their published result
  possible. Retained in
  [FINDINGS.md](FINDINGS.md#appendix-c--phase-1-the-authors-dasert-outcome).
- **The two are not comparable.** `DA+Sert` is bimodal with a real gap at 0.2;
  `DA/all` is unimodal with none, so the same threshold does a different job. A
  lower number here is a harder target, not a regression.

## Caveats

- **No external validation.** 20 donors, one protocol, one study.
- **PCs are pool-confounded** (η² up to 0.96), and D11's strongest PC under this
  outcome is the worst offender (η² 0.832). See
  `modeling/docs/PCA_interpretation.md`.
- **`PCn` is a slot, not a stable axis** — PC2/PC3 rotate between fits. Do not
  compare PC indices across tables.
- Open items tracked in
  [`modeling/docs/audit_ledger.md`](modeling/docs/audit_ledger.md).
