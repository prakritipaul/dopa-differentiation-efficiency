# pluricon-prototype

**Predicting day-52 dopaminergic neuron yield from earlier single-cell
snapshots**, across 136 iPSC lines.

Practical question: partway through a midbrain differentiation, can you tell
whether it will produce dopaminergic neurons — early enough to act on it?

**Short answer: yes at D11 (ROC-AUC 0.906), and by D30 the outcome is largely
already determined.**

Full write-ups: **[FINDINGS.md](FINDINGS.md)** (results + literature context) ·
[`modeling/README.md`](modeling/README.md) (methods) ·
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

Proportions sum to 1, so `phat_P_FPP` is the implicit reference and does not
enter a fit. PCs are standardised, then averaged **pool-then-line** (mean per
`(line, pool)`, then unweighted mean of pool means).

Both feature types are computed **per fold**: HVG selection and the PCA are fit
on training-line cells only, then *all* cells are projected — so no held-out
line influences the basis it is scored in. The depth-outlier pool is excluded
from fitting but still projected, so no line is lost. **It is not the same pool
at both timepoints.**

## Cross-validation

| | |
|---|---|
| Primary scheme | **donor-grouped** 5-fold, 10 repeats, stratified by the binary label (134/136 lines share a donor) |
| Also run | plain 5-fold, LOCO (136 folds), LODO (20 folds) |
| Tuning | **nested** — inner 3-fold on training lines only |
| Metrics | computed **per repeat**, then averaged — never pooled across repeats |
| Total | 256 folds, PCA refit for each, per timepoint |

Model per task **pre-registered before looking at results**: donor-grouped +
nested + ridge / logistic_l2. Both timepoints use **identical fold
assignments**, verified to agree on every fold's train/test membership, so the
D11-vs-D30 comparison is paired.

## Models and hyperparameters

Tuned jointly over `k` (number of PCs) and regularisation — 55 configurations
per model.

| model | task | grid |
|---|---|---|
| ridge *(headline)* | regression | `alpha` ∈ {0.01, 0.1, 1, 10, 100} |
| lasso | regression | `alpha` ∈ {0.001, 0.01, 0.1, 1, 10} |
| logistic_l2 *(headline)* | classification | `C` ∈ {0.01, 0.1, 1, 10, 100} |
| logistic_l1 | classification | `C` ∈ {0.01, 0.1, 1, 10, 100} |

`k` ∈ 0–10 for all. Both logistic models use `class_weight="balanced"` and a
fixed seed. Features standardised on training folds only.

---

## Results

Donor-grouped, nested, mean ± SD across 10 repeats.

| | **D11 → D52** | **D30 → D52** |
|---|---|---|
| ridge — R² | 0.503 ± 0.015 | 0.802 ± 0.009 |
| ridge — MAE | 0.099 | 0.061 |
| logistic_l2 — ROC-AUC | 0.906 ± 0.008 | 0.945 ± 0.008 |
| logistic_l2 — Brier | 0.125 | 0.090 |

### Three findings

**1. D30's advantage is mostly definitional.** Its annotated cell types already
*include* DA, so the outcome partly exists in the predictors — construct
overlap, not leakage. The right reading is *how much of the D52 phenotype is
already established by D30*: ~48% of D30 cells have arrived, ~31% are still
uncommitted progenitors, ~20% are off-target.

**2. At D30, expression adds nothing beyond cell-type composition.** Three of
four models select **zero PCs**. Permutation importance: 0.360 for the
proportions vs 0.018 for the best PC.

**3. One raw column out-ranks the full model — but is badly calibrated.**

| D30 | `phat_DA` alone | full model |
|---|---|---|
| ROC-AUC | **0.956** | 0.939 |
| log loss | 0.571 | **0.328** |

To *rank* lines, one number suffices. For trustworthy *probabilities*, the model
earns its keep. AUC alone would mislead either way.

**At D11**, where no such shortcut exists, the neuroblast fraction `phat_NB` is
the strongest and most technically clean predictor (ρ = −0.55, pool η² = 0.035)
— more D11 neuroblasts, worse D52 yield.

## Repository

```
eda/              numbered 0NN_*.py scripts, run in order — the EDA phase
modeling/         the pipeline (see modeling/README.md)
  docs/           methods, PCA interpretation, audit ledger
  results/        every output table (see modeling/results/README.md)
  tests/          pytest suite, incl. a synthetic-data smoke test
metadata_eda/     EDA outputs: cohort/ pca/ proportions/ qc/ technical/ plots/
tools/            reference-PDF builder
```

## Setup

```bash
uv sync
uv run pytest modeling/ -q      # full suite — runs without the expression data
```

The expression files (`day11.h5`, `day30.h5`, `day52.h5`) are multi-GB and not
in the repo. **You do not need them to verify the code**: the suite builds a
miniature dataset in the same on-disk shape and drives the real pipeline over it
in ~3 seconds. To run on real data, set the paths in `modeling/features.py` —
commands in [`modeling/README.md`](modeling/README.md).

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

An earlier phase of this project reproduced the authors' definition exactly
(138 lines, D11 only, ROC-AUC 0.947), which is what made comparison to their
published result possible. It is retained in
[FINDINGS.md](FINDINGS.md#appendix--phase-1-the-authors-dasert-outcome).

**The two are not comparable.** `DA+Sert` is bimodal with a real gap at 0.2;
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
