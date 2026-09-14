# pluricon-prototype

**Can you tell, early in a midbrain differentiation, whether it will produce
dopaminergic neurons?**

Using [Jerber et al. 2021](https://www.nature.com/articles/s41588-021-00801-6)
(*Nat Genet* 53:304–312) — scRNA-seq of iPSC lines differentiated toward
midbrain dopaminergic fate, sampled at D11, D30 and D52.

**Short answer: yes at D11, and the useful part is *which* signals carry it.**
Under never-seen-donor cross-validation, D11 features predict D52 dopaminergic
yield at ROC-AUC 0.906. D30 reaches 0.945 — but most of that is the outcome
already existing in the predictors, which is itself the interesting finding.

Full write-ups: **[FINDINGS.md](FINDINGS.md)** (results, literature
assessment) · [`modeling/README.md`](modeling/README.md) (methods) ·
[`modeling/results/README.md`](modeling/results/README.md) (every output file)

---

## Two phases — read this before comparing any numbers

| | **Phase 1 — exploratory** | **Phase 2 — headline** |
|---|---|---|
| outcome | `(DA + Sert) / all D52 cells` | **`DA / all D52 cells`** |
| D52 cells | all, incl. rotenone-treated | **untreated only** |
| cohort | 138 lines | **136 lines**, 157 combos, 20 donors |
| timepoints | D11 → D52 | **D11 → D52 and D30 → D52** |
| why | reproduce the authors' own definition so results are comparable to theirs | answer the question the protocol is actually for |

Phase 1 adopted the authors' exact outcome definition — including the choices
we later departed from. That is what made a check against a published result
possible, and it surfaced the two problems that motivated phase 2:

1. **The published metric counts rotenone-treated cells.** D52 alone contains
   219,238 `ROT` cells against 303,856 `NONE`, interleaved so every qualifying
   combo contains both. Rotenone inhibits mitochondrial complex I and
   preferentially damages dopaminergic neurons — the cells in its numerator.
   Verified against the authors' own notebook: no treatment filter anywhere.
2. **It sums two lineages that track independently.** D30 DA predicts D52 DA at
   ρ = +0.932; D30 Sert predicts D52 DA at Pearson **+0.057**.

Phase 1 is kept in full, not deleted — the literature comparison rests on it.

> ⚠️ **The two phases' numbers are not comparable.** Different outcome,
> different cells, different cohort. `DA+Sert` is bimodal with a real gap at the
> 0.2 threshold; `DA/all` is unimodal with none. **A lower phase-2 number is a
> harder target, not a regression.**

## Results — phase 2 (headline)

Donor-grouped 5-fold × 10 repeats, nested tuning, pre-registered model per task.

| | **D11 → D52** | **D30 → D52** |
|---|---|---|
| ridge — R² | 0.503 ± 0.015 | 0.802 ± 0.009 |
| logistic_l2 — ROC-AUC | 0.906 ± 0.008 | 0.945 ± 0.008 |
| — Brier | 0.125 | 0.090 |

Both timepoints were scored on fold tables verified to agree on **every** fold's
train/test membership, so the comparison is paired and attributable to the
timepoint alone.

### The three findings worth your attention

**1. At D30, the transcriptome adds nothing beyond cell-type composition.**
Three of four models select **k = 0 PCs**. Permutation importance: 0.360 for the
proportions against 0.018 for the best single PC — a 20× gap.

**2. One raw column beats the whole model at ranking — but not at calibration.**
`phat_DA` alone vs the full D30 model, identical folds and nested selection:

| | `phat_DA` alone | full model | full better in |
|---|---|---|---|
| ROC-AUC | **0.9563** | 0.9393 | 0 / 10 repeats |
| log loss | 0.5707 | **0.3276** | 10 / 10 repeats |
| Brier | 0.1897 | **0.0930** | 10 / 10 repeats |

To *rank* lines, one column suffices. For trustworthy *probabilities*, the model
earns its keep. AUC alone would have given the wrong answer either way.

**3. D30's strength is largely definitional.** Its annotation set already
contains DA cells, so the outcome partly exists in the predictors. This is
**construct overlap, not leakage** — D30 precedes D52, separate cells, no D52
information in the features. The defensible claim is *how much of the D52
phenotype is already established by D30*. Marker analysis over 250,923 D30 cells
splits the seven types into **arrived** (DA, Sert; 48%), **undecided**
(FPP, P_FPP; 31%), and **off-target** (Epen1, U_Neur1/2; 20%) — see
[`modeling/docs/D30_celltype_interpretation.md`](modeling/docs/D30_celltype_interpretation.md).

## Method

**Features.** Cell-type proportions (2 free at D11, 6 at D30 — they sum to 1, so
one is the implicit reference) plus PC1–10 from 2,000 HVGs. Both computed
**per fold**: HVG selection and PCA are fit on training-line cells only, then
*all* cells are projected, so no held-out line influences the basis it is scored
in. The depth-outlier pool is excluded from fitting but still projected —
**pool11 at D11, pool5 at D30**; it is not the same pool at both timepoints.

**Cross-validation.** Donor-grouped 5-fold × 10 repeats (136 of 138 lines share
a donor), plus plain, LOCO and LODO as robustness. Tuning is **nested** — inner
3-fold on training lines only. Metrics are computed **per repeat**, then
averaged, never pooled across repeats.

**Models.** ridge / lasso (regression), logistic_l2 / logistic_l1
(classification), tuned jointly over `k` ∈ 0–10 and regularisation — 55
configurations. The headline model per task was **pre-registered before looking
at results**.

## Repository

```
eda/              numbered 0NN_*.py scripts, run in order — the EDA phase
modeling/         the pipeline (see modeling/README.md)
  docs/           methods write-ups, PCA interpretation, audit ledger
  results/        every output table (see modeling/results/README.md)
  tests/          pytest suite, incl. a synthetic-data smoke test
metadata_eda/     EDA outputs: cohort/ pca/ proportions/ qc/ technical/ plots/
tools/            reference-PDF builder
FINDINGS.md       results + literature assessment for both phases
```

## Setup

```bash
uv sync
uv run pytest modeling/ -q          # full suite, incl. the synthetic smoke test
```

**The expression data is not in the repo** (`day11.h5`, `day30.h5`, `day52.h5`
are multi-GB). You do **not** need it to verify the code: the suite includes an
end-to-end smoke test that builds a miniature dataset in the same on-disk shape
and drives the real pipeline over it — cohort construction, per-fold PCA feature
extraction, and a CV fit — in about 3 seconds.

To run on the real data, set the paths in `modeling/features.py`, then:

```bash
uv run python eda/003_qualifying_cell_lines.py --untreated-only --suffix _untreated
uv run python eda/009_d52_outcome_label.py --untreated-only --da-only \
    --cohort-suffix _untreated --suffix _da_untreated
uv run python -m modeling.run_feature_extraction --timepoint D11 \
    --label-variant da_untreated --restrict-fit-to-qualifying --suffix _qualonly
uv run python -m modeling.run_variant --timepoint D11 --suffix _qualonly \
    --label-variant da_untreated --out-suffix _D11_qualonly
```

Feature extraction is ~3 h per timepoint and resumable. Python 3.11,
[uv](https://docs.astral.sh/uv/).

## Caveats

- **No external validation.** 20 donors, one protocol, one study. Everything is
  internal cross-validation.
- **PCs are pool-confounded** (η² up to 0.96 on some components), and D11's
  strongest PC under the DA-only outcome is its *worst* offender (η² 0.832).
  Read `modeling/docs/PCA_interpretation.md` before interpreting any PC.
- **`PCn` is a slot, not a stable axis.** PC2/PC3 are near-tied in variance and
  rotate between fits. Do not compare PC indices across tables.
- **Open items** are tracked in
  [`modeling/docs/audit_ledger.md`](modeling/docs/audit_ledger.md), including a
  possible pool5 annotation bias at D30 that excluding it from the PCA fit does
  not correct.
