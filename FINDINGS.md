# Findings

**Predicting day-52 dopaminergic neuron yield from earlier single-cell
snapshots** — 136 iPSC lines, 20 donors,
[Jerber et al. 2021](https://www.nature.com/articles/s41588-021-00801-6).

Midbrain dopaminergic differentiation, scored at day 52. This asks whether
single-cell data from days 11 and 30 already contain signal to answer that
question. What are the most informative features, and what do they say about the
biology?

| | day 11 — predictor | day 30 — predictor | day 52 — outcome |
|---|---|---|---|
| what is there | floor-plate progenitors, neuroblasts; **no DA cells yet** | 7 cell types; DA neurons have appeared | the dopaminergic fraction being predicted |
| headline | **ROC-AUC 0.906** | ROC-AUC 0.945 | median 0.165, range 0.010–0.742 |
| reading | a real forecast, 41 days out | higher, but largely definitional | — |

Outcome: `DA / all D52 cells`, untreated cells only. Success = ≥ 0.2 (61 / 75
lines). Donor-grouped 5-fold × 10 repeats. Methods:
[`modeling/README.md`](modeling/README.md); result files in
[`modeling/results/README.md`](modeling/results/README.md).

---

## 1. Is there a reliable predictive model?


**Yes — from day 11, 41 days ahead of the readout.** ROC-AUC `0.906 ± 0.008`, R² `0.503`, Brier `0.125`.

Hyperparameters come from nested inner folds that never see the test lines; every PCA is refit per fold on training lines only. And four independent cross-validation schemes agree — nothing hinges on one lucky split.

**0.902** Plain 5-fold · **0.906** Donor-grouped · **0.912** Leave-one-line-out · **0.908** Leave-one-donor-out

### The model that would be deployed

Configuration is chosen by **flat CV + one-SE**: the simplest config within one standard error of the grid's best, on donor-grouped folds.

*Shipped configuration per task · donor-grouped, 10 repeats*

| Model | Configuration | Features | Nested | Flat |
|---|---|---|---|---|
| D11 classification *(headline)* | logistic_l2, C=1.0, k=4 | phat_FPP, phat_NB + PC1–PC4 | 0.906 | 0.918 |
| D11 regression | ridge, α=10.0, k=7 | same 2 proportions + PC1–PC7 | 0.503 | 0.506 |
| D30 classification | logistic_l2, C=1.0, k=6 | 6 proportions + PC1–PC6 | 0.945 | 0.956 |
| D30 regression | ridge, α=10.0, **k=0** | 6 proportions, **no PCs** | 0.802 | 0.812 |

Classification is ROC-AUC, regression R². `phat_P_FPP` is the reference proportion and never enters a fit — which is why D11 contributes two proportions and D30 six.

**Why both numbers.** Flat CV is the *selection criterion* — the config was picked on those folds, so its score there is optimistic. Nested CV is the out-of-fold estimate. The gap is small, 0.003 to 0.012, but **positive in all four cases**, which is what selection bias looks like and what noise does not.

**What nested actually estimates.** It scores the selection *procedure*, not one fixed config — each outer fold re-runs the tuning and may land on a different `k`, so the figure describes “run this whole pipeline on new data”.

**D30 regression selects k=0** — proportions only, no expression component at all. The selection rule lands independently where the permutation importances already pointed.

> **One limit, stated up front.** At day 30 a single raw cell-type proportion — `phat_DA` — *ranks* lines better than the whole model (0.956 vs 0.939, in 10 of 10 repeats). What the model gives you instead is a probability you can act on: when it says a line has a 70% chance of succeeding, roughly 70% of such lines do. The raw proportion sorts lines correctly, but its numbers mean nothing on that scale — you could use it to pick the best ten lines, not to decide whether any one line is worth continuing.


## 2. Which features predict — and how much do they track the batch?

**Three claims that must not be blurred together:** whether a feature improves prediction on held-out lines, whether it would still work in someone else's batch, and whether its gene content means something. Only the first is measured here.

Features below are ranked by how hard the fitted model leans on them, with pool η² reported beside each one rather than used to sort them into bins. **Pool η² is the share of a feature's across-line variance associated with which of the ten differentiation runs a line went through.** It measures *potential run dependence* — not the fraction proven to be technical artefact, not the probability a feature is spurious, and not proof that a low-η² feature will generalise.

> **Importance against run association, day 11.** Horizontal: mean |SHAP|, how much the feature moves a prediction. Vertical: pool η². No threshold is drawn, because η² runs continuously from 0.035 to 0.832 and no line in it has been validated as a transfer test.

### The most important day-11 feature is also the most run-associated

`PC3` leads on both counts: the largest out-of-fold contribution of any feature (`+0.064`, 4.5× `phat_NB`'s) **and** η² = `0.832`. Both are true and neither cancels the other. The defensible statement is that **PC3 materially improves prediction among held-out lines, while whether it ports to a new run is unestablished and biological attribution is hazardous.**

Its loadings are worth reading as *hypothesis-generating*, not as mechanism. The negative pole is a coherent biosynthetic / growth programme — `PTMA RANBP1 SNRPB HSP90AA1 PA2G4 NASP TYMS PHGDH PSAT1` — ribosome biogenesis, chaperones, nucleotide and serine synthesis, with S-phase genes at median rank 57 of 2000 (p = 8.1 × 10⁻⁸). **A growth-rate axis is exactly what would differ between culture runs**, which makes η² = 0.832 unsurprising rather than mysterious — and is why this gene list should not anchor a biological conclusion.

### `phat_NB` is the least run-associated of the important features

η² = `0.035`, an order of magnitude below anything else that matters, with a clear direction: more day-11 neuroblasts → worse day-52 yield (ρ −0.554). That makes it the best portability *candidate*. But its out-of-fold contribution is modest (`+0.014`), so “it carries the signal” would oversell it.

### `PC2` corroborates it from expression alone

η² = `0.114` and the strongest univariate correlation in the table (ρ −0.626), so its loadings can be read with more confidence. They describe a **neurogenesis axis**: `NEUROD1 NHLH1 ELAVL3 DLL3 STMN2 MLLT11 ONECUT2` against `HES1 GPC3 FRZB BMP4 OTX2 CDH2`. Pan-neuronal genes sit at median rank 1958 of 2000 on the neurogenic pole (p = 7.1 × 10⁻⁷), and `NEUROD1`/`DLL3` against `HES1` is textbook Notch lateral inhibition.

**And it agrees with the composition measurement independently.** `PC2` tracks the neuroblast fraction at ρ = **+0.815** while being computed from genes alone, never from cell-type calls. That convergence — counting annotated neuroblasts, and reading a proneural expression programme — is what makes premature neurogenesis the most trustworthy day-11 finding, more than either number alone.

### At day 30 the picture is simpler

Composition carries the model and expression does not: grouped permutation Δ is `0.360` for the proportions against `0.018` for the best PC, and three of four model families select **zero** PCs. **There is no important low-η² expression component at day 30** — the two PCs with real contribution (`PC1`, `PC5`) are both at η² 0.518, and the two with low η² contribute essentially nothing, because they re-measure cell counts the annotation already provides.

> **The day-30 headline feature carries its own caveat.** `phat_DA` has η² = **0.315**, nearly ten times `phat_NB`'s. Least run-associated *among the day-30 proportions* is not the same as low. Full day-30 table and gene loadings in the appendix.


## 3. What does it say about the biology?

**Premature neurogenesis at day 11 predicts failure.** Lines that have already made neuroblasts by day 11 yield *fewer* dopaminergic neurons at day 52 — `phat_NB` correlates ρ = −0.554.

The leading expression axis says the same thing from the other side. **`PC1` at day 11** — 5.4% of the variance across 2,000 highly variable genes — is a cell-cycle axis: its top loadings are textbook G2/M (`HMGB2 PTTG1 NUSAP1 CENPF TOP2A MKI67`), and it tracks the *proliferating* progenitor fraction at `+0.83`. Keeping a cycling progenitor pool at day 11 is the good outcome; differentiating early is not.

**`PC1` at day 30** (9.8% of variance) **is that axis inverted.** The cell cycle moves to the opposite pole and the leading genes are pan-neuronal (`MLLT11 TUBB2B STMN2 GAP43 MAPT`). It now measures neuronal maturation — and it is driven by composition (R² 0.891) far more than by cell state (R² 0.484).

`PC1` is the only component read this way: it replicates across PCA fits at r = 0.9995, while the lower components are near-tied in variance. It is also *not* the component that carries the day-11 prediction — that is `PC3`. PCA sign is arbitrary, so each pole is named by its own genes.

> **There is a real window at day 30.** About 31% of day-30 cells have genuinely not decided — floor-plate progenitors, SOX2⁺/HES1⁺/VIM-high — against 48% already arrived and 20% off-target terminal fates. But looking *only* at those undecided cells, ignoring the neurons a line has already made, predicts the outcome at AUC 0.727 — better than a coin flip, far worse than simply counting the neurons. **By day 30 most of the answer is in what a line has already become, not in what is left to decide.**

## 4. Why this outcome, not the published one

The published metric is `(DA + Sert) / all D52 cells`. This project scores
**`DA / all D52 cells`, untreated cells only** — a deliberate divergence, not a
bug fix, and it changes what every number here means.

1. **The published metric counts rotenone-treated cells.** D52 alone holds
   219,238 `ROT` against 303,856 `NONE`, interleaved through every qualifying
   combination. Rotenone inhibits mitochondrial complex I and preferentially
   damages dopaminergic neurons — the cells in its own numerator. Verified
   against the authors' notebook: it applies a ≥10-cell threshold and no
   treatment filter anywhere.
2. **It sums two lineages that behave differently here.** At D30, `phat_DA`
   predicts the D52 dopaminergic fraction at ρ **+0.932**; `phat_Sert` predicts
   the same outcome at Pearson **+0.057**. That dopaminergic and serotonergic
   neurons arise from distinct rostro-caudal floor-plate domains is established
   biology, not a finding here — but summing them presumes they behave as one
   quantity for the purpose of scoring efficiency, and in this cohort they do not.

**Cohort cost:** dropping treated cells puts 2 lines below the ≥10-cell
threshold (both in `pool5`), giving **136 lines / 157 combinations / 20 donors**
against phase 1's 138. The donor count is unchanged, so donor-grouped CV keeps
its structure.

**The two outcome families are not comparable and must never be differenced.**
`DA+Sert` is bimodal with a real gap at 0.2; `DA/all` is unimodal with none, so
the same threshold does a different job. A lower number here is a harder target,
not a regression.


## 5. Performance, in full

Every fitted model, every CV scheme, nested tuning throughout — nothing filtered. Both timepoints were scored on fold tables verified to agree on every fold's train/test membership, so the comparison is paired.

*D11 → D52 · classification*

| Model | Scheme | ROC-AUC | PR-AUC | Bal. acc | Sens. | Spec. | F1 | Brier |
|---|---|---|---|---|---|---|---|---|
| `logistic_l2` | plain | 0.902 ± 0.010 | 0.864 ± 0.010 | 0.837 ± 0.015 | 0.856 ± 0.034 | 0.817 ± 0.019 | 0.822 ± 0.017 | 0.129 ± 0.007 |
| `logistic_l2` | donor_grouped | 0.906 ± 0.008 | 0.874 ± 0.012 | 0.834 ± 0.021 | 0.864 ± 0.031 | 0.804 ± 0.026 | 0.821 ± 0.022 | 0.125 ± 0.006 |
| `logistic_l2` | loco | 0.912 ★ | 0.883 | 0.861 ★ | 0.869 | 0.853 ★ | 0.848 ★ | 0.121 ★ |
| `logistic_l2` | lodo | 0.908 | 0.886 | 0.820 | 0.852 | 0.787 | 0.806 | 0.126 |
| `logistic_l1` | plain | 0.908 ± 0.016 | 0.873 ± 0.018 | 0.839 ± 0.025 | 0.844 ± 0.034 | 0.833 ± 0.022 | 0.824 ± 0.027 | 0.123 ± 0.012 |
| `logistic_l1` | donor_grouped | 0.905 ± 0.011 | 0.872 ± 0.012 | 0.837 ± 0.016 | 0.872 ± 0.017 ★ | 0.803 ± 0.025 | 0.825 ± 0.016 | 0.126 ± 0.009 |
| `logistic_l1` | loco | 0.907 | 0.877 | 0.861 ★ | 0.869 | 0.853 ★ | 0.848 ★ | 0.123 |
| `logistic_l1` | lodo | 0.906 | 0.886 ★ | 0.821 | 0.869 | 0.773 | 0.809 | 0.126 |
| trivial baseline | — | 0.500 | 0.449 | 0.500 | 1.000 | 0.000 | 0.619 | 0.247 |

*D30 → D52 · classification*

| Model | Scheme | ROC-AUC | PR-AUC | Bal. acc | Sens. | Spec. | F1 | Brier |
|---|---|---|---|---|---|---|---|---|
| `logistic_l2` | plain | 0.938 ± 0.015 | 0.910 ± 0.020 | 0.884 ± 0.012 | 0.885 ± 0.020 | 0.883 ± 0.016 | 0.872 ± 0.013 | 0.092 ± 0.007 |
| `logistic_l2` | donor_grouped | 0.945 ± 0.008 | 0.906 ± 0.019 | 0.877 ± 0.013 | 0.880 ± 0.022 | 0.875 ± 0.009 | 0.865 ± 0.015 | 0.090 ± 0.006 |
| `logistic_l2` | loco | 0.945 ★ | 0.886 | 0.876 | 0.885 | 0.867 | 0.864 | 0.087 ★ |
| `logistic_l2` | lodo | 0.938 | 0.878 | 0.884 ★ | 0.902 ★ | 0.867 | 0.873 ★ | 0.090 |
| `logistic_l1` | plain | 0.934 ± 0.014 | 0.908 ± 0.022 | 0.881 ± 0.011 | 0.872 ± 0.020 | 0.891 ± 0.014 ★ | 0.869 ± 0.013 | 0.094 ± 0.008 |
| `logistic_l1` | donor_grouped | 0.940 ± 0.009 | 0.899 ± 0.023 | 0.870 ± 0.013 | 0.861 ± 0.018 | 0.879 ± 0.015 | 0.856 ± 0.014 | 0.093 ± 0.005 |
| `logistic_l1` | loco | 0.943 | 0.914 ★ | 0.883 | 0.885 | 0.880 | 0.871 | 0.090 |
| `logistic_l1` | lodo | 0.939 | 0.873 | 0.876 | 0.885 | 0.867 | 0.864 | 0.091 |
| trivial baseline | — | 0.500 | 0.449 | 0.500 | 1.000 | 0.000 | 0.619 | 0.247 |

*D11 → D52 · regression*

| Model | Scheme | R² | MAE | RMSE | Out-of-range |
|---|---|---|---|---|---|
| `ridge` | plain | 0.515 ± 0.018 ★ | 0.100 ± 0.002 | 0.133 ± 0.002 ★ | 0.056 ± 0.010 |
| `ridge` | donor_grouped | 0.503 ± 0.015 | 0.099 ± 0.001 | 0.134 ± 0.002 | 0.057 ± 0.008 |
| `ridge` | loco | 0.501 | 0.100 | 0.134 | 0.066 |
| `ridge` | lodo | 0.505 | 0.099 ★ | 0.134 | 0.051 |
| `lasso` | plain | 0.514 ± 0.015 | 0.100 ± 0.002 | 0.133 ± 0.002 | 0.050 ± 0.019 ★ |
| `lasso` | donor_grouped | 0.502 ± 0.023 | 0.100 ± 0.001 | 0.134 ± 0.003 | 0.051 ± 0.011 |
| `lasso` | loco | 0.507 | 0.099 | 0.134 | 0.066 |
| `lasso` | lodo | 0.513 | 0.099 | 0.133 | 0.066 |
| trivial baseline | — | 0.000 | 0.158 | 0.190 | 0.000 |

*D30 → D52 · regression*

| Model | Scheme | R² | MAE | RMSE | Out-of-range |
|---|---|---|---|---|---|
| `ridge` | plain | 0.817 ± 0.014 | 0.059 ± 0.002 | 0.081 ± 0.003 | 0.037 ± 0.010 |
| `ridge` | donor_grouped | 0.802 ± 0.009 | 0.061 ± 0.001 | 0.085 ± 0.002 | 0.037 ± 0.006 |
| `ridge` | loco | 0.820 | 0.058 | 0.081 | 0.044 |
| `ridge` | lodo | 0.804 | 0.060 | 0.084 | 0.051 |
| `lasso` | plain | 0.828 ± 0.005 ★ | 0.056 ± 0.001 ★ | 0.079 ± 0.001 ★ | 0.004 ± 0.005 ★ |
| `lasso` | donor_grouped | 0.811 ± 0.006 | 0.058 ± 0.001 | 0.083 ± 0.001 | 0.010 ± 0.009 |
| `lasso` | loco | 0.821 | 0.058 | 0.081 | 0.029 |
| `lasso` | lodo | 0.819 | 0.058 | 0.081 | 0.015 |
| trivial baseline | — | 0.000 | 0.158 | 0.190 | 0.000 |

No row is privileged. The configuration that would ship is named above; it is not necessarily the best cell in any column here. ★ marks the best value per column. `loco` and `lodo` are single repeats by construction, so they carry no SD. Sensitivity 1.000 and F1 0.619 come free from always predicting “success” — which is why neither is ever quoted alone.

### Why day 30 scores higher — and why that is not a discovery

Day 30's annotation set already contains dopaminergic cells, so the outcome partly exists in the predictors. This is **construct overlap, not leakage**: day 30 and day 52 are separate cells from separate harvests, day 30 precedes day 52, and no day-52 information enters the features. The defensible claim is *how much of the day-52 phenotype is already established by day 30* — not that something was predicted.

*D30 · full model against one raw column, identical folds, paired per repeat*

| Metric | `phat_DA` alone | Full model | Full model better in |
|---|---|---|---|
| ROC-AUC | 0.9563 | 0.9393 | 0 / 10 repeats |
| Log loss | 0.5707 | 0.3276 | 10 / 10 repeats |
| Brier | 0.1897 | 0.0930 | 10 / 10 repeats |
| Accuracy | 0.8676 | 0.8750 | 6 / 10 repeats |

**The answer is split, and both halves are unanimous.** To rank lines, one raw column beats the entire model. To get trustworthy probabilities, the model earns its keep. This is exactly why AUC could not settle the question — at 0.95 it is saturated and speaks only to ordering.


**The fit-free reference lines**, for context — donor-level bootstrap over 20
donors, 2,000 resamples, in-sample so not like-for-like with the paired test
above:

| D30 benchmark (no model, no fitted parameter) | ρ | ROC-AUC |
|---|---|---|
| **`phat_DA` alone** | +0.932 | **0.960** [0.922–0.986] |
| `phat_DA + phat_Sert` | +0.799 | 0.886 [0.834–0.941] |
| `phat_FPP + phat_P_FPP` (undecided) | −0.788 | 0.883 |
| progenitor balance (share of non-target) | +0.452 | 0.727 |

**A mis-specified benchmark nearly hid this result.** The comparator was carried
over from the previous outcome and summed `phat_DA + phat_Sert` (0.886) rather
than `phat_DA` (0.960). Against the wrong reference the model appeared to add ~6
AUC points. **A benchmark whose numerator does not match the outcome is not a
weaker check — it is a misleading one.**

Applies to D30 only: D11 has no `phat_DA` (its three annotated types are FPP, NB,
P_FPP), so there is no analogous single-column comparator there.

**One related caution on reading coefficients.** `phat_Sert` has univariate
ρ **+0.215** but a fitted coefficient of **−0.288** — positive alone, negative
once the other types are held fixed. A conditional coefficient is not a simple
biological effect.



## 6. Literature assessment: recapitulated vs. new

The baseline comparison is worth restating: **Jerber et al. predicted differentiation
efficiency from *iPSC-stage bulk RNA-seq*, before differentiation began.**
Nothing here reproduces that. These are different predictors, a different
outcome, and in one case a different timepoint.

**Scope limit, stated up front.** "New" below means *not present in the authors'
analysis that we actually read* — their `Figure_2` notebook and the outcome
definition it encodes. We did not audit the full paper or its supplements. Treat
"new" as "not found where we looked", not as a priority claim.

| # | Finding | Status | Basis |
|---|---|---|---|
| 1 | D11 scRNA-seq predicts D52 **dopaminergic** yield: ROC-AUC 0.906, R² 0.503 | **New** | The authors predicted from iPSC-stage bulk, and their efficiency metric sums DA with Sert; a DA-only D11 → D52 predictor is not part of their analysis |
| 2 | D30 → D52 dopaminergic yield: ROC-AUC 0.945, R² 0.802 | **New, but largely definitional** | No D30 → D52 predictor in their work. The strength is mostly construct overlap — D30 already contains DA cells — not predictive discovery. See §5's construct-overlap note |
| 3 | **DA and Sert track independently**: D30 DA → D52 DA ρ = +0.932, while D30 Sert → D52 DA is Pearson **+0.057** | **New; and it undercuts the combined metric** | The authors' `diff_efficiency` sums the two, which presumes they behave as one quantity. They do not: `DA/(DA+Sert)` spreads near-uniformly 0–1 across lines, consistent with a line-intrinsic rostro-caudal identity fixed before D30 |
| 4 | `phat_NB` (neuroblast proportion) remains D11's strongest *and technically cleanest* composition feature (pool η² 0.035) | **Recapitulates phase 1, and the authors' indirect observation** | They noted a poor-differentiation cluster correlating with D11 neuroblast proportion. Phase 1 made it the top direct predictor; it survives the outcome change with the same sign |
| 5 | **At D30 the transcriptome adds nothing beyond cell-type composition** — 3 of 4 models select k = 0 PCs; permutation Δ 0.360 for proportions vs 0.018 for the best PC | **New** | A 20× gap. Not addressed by the authors, who did not build D30 predictors |


## 7. What this does not establish

- **No external validation.** 20 donors, one protocol, one study, no held-out batch. Everything here is internal cross-validation.
- **No causal claim** anywhere. All associations.
- **Day 11's strongest expression component is 83% pool-explained** — it may not transfer to a new run. The composition features are the robust ones.
- **Pool 5 annotation bias is unaddressed.** Its ~10× day-30 sequencing-depth deficit may have biased cell-type annotation itself; excluding it from PCA fitting does nothing to correct biased labels.
- **Whether day-30 and day-52 cell types were annotated independently is unknown.** Not leakage either way, but correlated labelling error would inflate the apparent continuity.
- **Batch and line are confounded by design**, so pool η² ranks the effects rather
  than partitioning the variance.
- **Not comparable to the phase-1 family.** Different outcome, different D52
  cells, different cohort — see §4.


---

# Appendix A — phase 1: the authors' DA+Sert outcome

*The bridge to the published result. Phase 1 reproduced Jerber et al.'s own
definition — `(DA+Sert) / all D52 cells`, **all** cells including rotenone-treated,
**138** lines (the 136 above plus `HPSI0115i-melw_1` and `HPSI0115i-qecv_2`,
which fall below the ≥10-cell threshold once treated cells are removed), D11
only. It is what makes comparison to their published result possible, and what
surfaced the two problems in §4. It used a **pre-registered** headline model,
declared before results were seen; that is recorded faithfully and not rewritten.
The current analysis supersedes it with the flat-CV + one-SE rule in §1.*

| phase 1, donor-grouped, nested | value |
|---|---|
| `logistic_l2` ROC-AUC | **0.947 ± 0.007** |
| `ridge` R² | **0.653 ± 0.021** |
| cohort | 138 lines, 96 success / 42 failure at the authors' 0.2 threshold |
| LOCO / LODO agreement | ROC-AUC 0.942 / R² 0.682 |

> **Not comparable to the results above**, for the reasons in §4. A lower number
> in the main analysis is a harder target, not a regression.

Full phase-1 tables, its own literature assessment and its technical-confounder
analysis are preserved in the repository's git history at commit `9a23404`, and
the artifacts it produced are unchanged on disk
(`modeling/results/results_{regression,classification}.csv`).



---

# Appendix B — feature importance and run association, in detail

*Supporting detail for [§2](#2-which-features-predict--and-how-much-do-they-track-the-batch).
**"Run-associated" throughout, not "contaminated"** — a high η² shows a feature
varies with the differentiation run, which is not the same as showing it is a
technical artefact. Nothing here demonstrates contamination.*

Supporting detail for answer 02. Kept out of the main narrative because it is a methods comparison, not a finding. **“Run-associated” throughout, not “contaminated”** — a high η² shows a feature varies with the differentiation run, which is not the same as showing it is a technical artefact. Nothing here demonstrates technical contamination.

> **Importance against run association, day 30.** Same axes and scale as the day-11 chart, so the two are directly comparable. The two PCs with real contribution, `PC1` and `PC5`, are both at η² 0.518; the two with low η² contribute almost nothing.

*Every day-30 feature · logistic_l2*

| Feature | SHAP | ρ | LOCO Δ | Pool η² |
|---|---|---|---|---|
| `phat_DA` | 1.292 | +0.792 | −0.002 | 0.315 |
| `PC1` | 0.653 | +0.680 | +0.001 | 0.518 |
| `phat_Epen1` | 0.652 | −0.614 | +0.000 | 0.133 |
| `PC5` | 0.610 | +0.557 | +0.002 | 0.518 |
| `PC4` | 0.530 | +0.182 | +0.001 | 0.376 |
| `PC6` | 0.499 | −0.201 | +0.002 | 0.406 |
| `phat_FPP` | 0.470 | −0.580 | −0.001 | 0.211 |
| `phat_U_Neur2` | 0.406 | −0.309 | +0.000 | 0.083 |
| `phat_U_Neur1` | 0.325 | −0.454 | −0.002 | 0.110 |
| `phat_Sert` | 0.246 | +0.214 | −0.001 | 0.057 |
| `PC2` | 0.151 | −0.213 | −0.001 | 0.076 |
| `PC3` | 0.138 | −0.301 | −0.001 | 0.143 |
| phat_P_FPP (reference) | — | −0.482 | +0.006 | 0.140 |

**All day-30 LOCO Δ are ≈ 0.** With six correlated proportions plus six PCs, dropping any single feature costs nothing measurable — no individual feature is *necessary*, which is why the grouped permutation Δ (0.360 for all proportions together) is the informative statistic at day 30, not the per-feature column.

### Why the low-η² day-30 components contribute nothing

Their gene content is perfectly legible — they simply re-measure the annotation. `PC2` (η² 0.076) is a motile-cilia programme (`C11orf88 FAM183A ZMYND10 ROPN1L RSPH1 TPPP3 PIFO`), the ependymal signature, and tracks `phat_Epen1` at +0.706 and `phat_Sert` at −0.861. `PC3` (η² 0.143) is a cell-cycle axis (G2/M at median rank 1922 of 2000, p = 6.6 × 10⁻¹⁸) and tracks `phat_FPP` at +0.815. Once the proportions are in the model, a component that restates them adds nothing. **They are redundant, not uninformative** — which means the day-30 annotation is capturing the same structure the transcriptome shows.

*The two timepoints side by side*

|  | Day 11 | Day 30 |
|---|---|---|
| Least run-associated important feature | phat_NB, η² 0.035 | phat_DA, η² 0.315 |
| Most important feature | PC3 (SHAP 1.398, η² 0.832) | phat_DA (SHAP 1.292, η² 0.315) |
| A low-η² PC that also predicts? | yes — PC2 (η² 0.114, ρ −0.626) | no |
| PCs in the shipped model | 4 of 10 | 3 of 4 families select zero |
| Composition vs expression | comparable | 0.360 vs 0.018 — 20× |

The contrast is real but should not be over-read. It says expression adds something at day 11 and little at day 30, which follows from day 30's annotation already containing the outcome's numerator. It does **not** say day 11's expression features would survive a new run — that remains untested at both timepoints.
