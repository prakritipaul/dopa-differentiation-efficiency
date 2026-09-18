# Findings

**Predicting day-52 dopaminergic neuron yield from earlier single-cell
snapshots** — 136 iPSC lines, 20 donors,
[Jerber et al. 2021](https://www.nature.com/articles/s41588-021-00801-6).

Outcome: `DA / all D52 cells`, untreated cells only. Success = ≥ 0.2 (61 / 75
lines). Donor-grouped 5-fold × 10 repeats. Configuration chosen by **flat-CV +
one-SE**; performance reported from **nested CV**, with both numbers shown.
Methods: [`modeling/README.md`](modeling/README.md).

---

# The three answers

## 1. Is there a reliable predictive model?

**Yes — 41 days ahead of the readout.** From a day-11 snapshot:

| D11 → D52 | |
|---|---|
| **ROC-AUC** | **0.906 ± 0.008** |
| R² | 0.503 ± 0.015 |
| Brier (calibration) | 0.125 |

Hyperparameters are chosen
by nested inner folds that never see the test lines; every PCA is refit per fold
on training lines only; and **four independent CV schemes agree** (plain 0.902,
donor-grouped 0.906, leave-one-line-out 0.912, leave-one-donor-out 0.908).
Nothing hinges on one lucky split.

A day-30 snapshot scores higher (**ROC-AUC 0.945, R² 0.802**) — but **that
number is largely definitional, not a better discovery.** D30's annotation
already contains DA cells, so the outcome partly exists in the predictors.
Detail in [§5](#5-what-d30-can-and-cannot-claim).

**One limit stated up front.** At D30, a single raw cell-type proportion *ranks*
lines better than the whole model. What the model gives you instead is a
probability you can act on: when it says a line has a 70% chance of succeeding,
roughly 70% of such lines do. The raw proportion sorts lines correctly, but its
numbers mean nothing on that scale — you could use it to pick the best ten lines,
not to decide whether any one line is worth continuing. Detail in
[§6](#6--the-d30-model-does-not-beat-a-single-raw-column).

### The model that would be deployed

Configuration is chosen by **flat CV + one-SE** (the simplest config within one standard error of the grid's
best, `modeling/feature_importance.py:219`) on donor-grouped folds:

| model | configuration | features | nested | flat |
|---|---|---|---|---|
| **D11 classification** *(primary)* | `logistic_l2`, C = 1.0, k = 4 | `phat_FPP`, `phat_NB` + PC1–PC4 | **0.906** ± 0.008 | 0.918 |
| **D11 regression** | `ridge`, α = 10.0, k = 7 | same 2 proportions + PC1–PC7 | **0.503** ± 0.015 | 0.506 |
| D30 classification | `logistic_l2`, C = 1.0, k = 6 | 6 proportions + PC1–PC6 | **0.945** ± 0.008 | 0.956 |
| D30 regression | `ridge`, α = 10.0, **k = 0** | 6 proportions, **no PCs** | **0.802** ± 0.009 | 0.812 |

Classification is ROC-AUC, regression R². `phat_P_FPP` is the reference
proportion and never enters a fit — which is why D11 contributes 2 proportions
and D30 contributes 6.

**Why both numbers are shown.** Flat CV is the *selection criterion*: the config
was picked on those folds, so its score there is optimistic. Nested CV is the
out-of-fold estimate. The gap is small — 0.003 to 0.012 — but **positive in all
four cases**, which is what selection bias looks like and what noise does not.

**What the nested number actually estimates.** It scores the *selection
procedure*, not one fixed config — each outer fold re-runs the tuning and may land
on a different `k`, so the figure describes "run this whole pipeline on new
data", which for one shipped config is a slightly conservative proxy.

**D30 regression selects `k = 0`** — proportions only, no expression component at
all. The selection rule reaches independently the same conclusion the permutation
importances did in [§7](#7-which-features-matter--full-tables): at D30 the
transcriptome adds essentially nothing beyond cell-type composition.

## 2. Which features predict, and how much do they track the batch?

**Three claims that must not be blurred together**, and this section separates
them deliberately:

1. **Predictive accuracy** — does the feature improve prediction on held-out
   lines? Measured by out-of-fold LOCO Δ.
2. **Generalisation to a new run** — would it still work in someone else's
   batch? **Not measured anywhere in this project.**
3. **Biological interpretation** — does its gene content mean something?

Features are ranked by how hard the fitted model leans on them (mean |SHAP|),
with pool η² reported honestly beside each one rather than used to sort them
into bins.

| D11 feature | SHAP | univariate ρ | out-of-fold LOCO Δ | pool η² |
|---|---|---|---|---|
| **`PC3`** | **1.398** | +0.614 | **+0.064** (largest) | **0.832** |
| `phat_NB` | 1.050 | −0.554 | +0.014 | **0.035** |
| `PC2` | 0.718 | **−0.626** | −0.002 | 0.114 |
| `PC1` | 0.466 | +0.452 | −0.004 | 0.472 |
| `PC4` | 0.341 | +0.242 | +0.004 | 0.506 |
| `phat_FPP` | 0.065 | +0.143 | −0.006 | 0.761 |

**What pool η² is.** The fraction of a feature's variance *across lines* that is
associated with which of the ten differentiation runs a line went through. It
measures **potential run dependence — not the fraction proven to be technical
artefact.** It is specifically **not** the probability the feature is artefactual,
not the share of predictive performance caused by batch, not evidence of a causal
batch effect, not proof that a low-η² feature will generalise, and not proof that
a high-η² feature contains no biology.

### The most important D11 feature is also the most run-associated

`PC3` leads on both counts: the largest out-of-fold contribution of any feature
(+0.064, 4.5× `phat_NB`'s) **and** η² = 0.832. Those two facts sit together and
neither cancels the other. The defensible statement to an industry reader is:
**`PC3` materially improves prediction among held-out lines; whether it ports to a
new run is unestablished, and biological attribution is hazardous.**

Its loadings are worth reading as *hypothesis-generating*, not as mechanism. The
negative pole is a coherent biosynthetic / growth program —
`PTMA RANBP1 SNRPB HSP90AA1 PA2G4 NASP TYMS PHGDH PSAT1 EIF4EBP1` — ribosome
biogenesis, chaperones, nucleotide and serine synthesis, with S-phase genes at
median rank 57 of 2000 (p = 8.1 × 10⁻⁸). The positive pole is membrane and
neuronal (`GPM6B TSPAN7 SYT4 NREP CDH2 EFEMP1`). **A growth-rate axis is exactly
the kind of thing that differs between culture runs**, which makes η² = 0.832
unsurprising rather than mysterious — and is why the gene list should not anchor a
biological conclusion.

### `phat_NB` is the least run-associated of the important features

η² = 0.035, an order of magnitude below anything else that matters, with a clear
direction: more D11 neuroblasts → worse D52 dopaminergic yield (ρ −0.554). That
makes it the best portability *candidate* in the table. But its out-of-fold
contribution is modest (+0.014), so "it carries the signal" would oversell it —
much of its information is recoverable from the other features.

### `PC2` corroborates it from expression alone

η² = 0.114 and the strongest univariate correlation in the table (ρ −**0.626**),
so its loadings can be read with more confidence. They describe a **neurogenesis
axis**:

| `PC2` pole (2.7% of HVG variance) | top loadings |
|---|---|
| neurogenic | `NEUROD1 NHLH1 ELAVL3 ELAVL4 DLL3 STMN2 TUBB2B MAP1B MLLT11 ONECUT2 ST18 INA` |
| progenitor | `HES1 GPC3 FRZB BMP4 CAPN6 ZNF503 TPM2 OTX2 CDH2 ARHGAP29` |

Pan-neuronal genes sit at median rank 1958 of 2000 on the neurogenic pole
(p = 7.1 × 10⁻⁷). `NEUROD1`/`DLL3` against `HES1` is textbook Notch lateral
inhibition — the switch deciding whether a progenitor stays a progenitor.

**And it agrees with the composition measurement independently.** `PC2` tracks the
neuroblast fraction at ρ = **+0.815** while being computed from genes alone, never
from cell-type calls. Two different measurements — counting annotated neuroblasts,
and reading a proneural expression program — point the same way. That
convergence, not either number alone, is what makes premature neurogenesis the
most trustworthy D11 finding.

### At D30 the picture is simpler, and worth stating plainly

Composition carries the model and expression does not: grouped permutation Δ is
**0.360** for the proportions against **0.018** for the best PC, and three of four
model families select **zero** PCs. **There is no important low-η² expression
component at D30** — the two PCs with real contribution (`PC1`, `PC5`) are both at
η² = 0.518, and the two with low η² (`PC2` 0.076, `PC3` 0.143) contribute
essentially nothing, because they are near-restatements of cell counts the
annotation already provides.

The D30 headline feature carries its own caveat: **`phat_DA` has η² = 0.315**,
nearly ten times `phat_NB`'s. Least run-associated *among the D30 proportions* is
not the same as low.

Full D30 table, the D30 gene loadings and the timepoint-by-timepoint comparison
are in [Appendix B](#appendix-b--feature-importance-and-run-association-in-detail).


## 3. What does this say about the biology?

**Premature neurogenesis at D11 predicts failure.** The `phat_NB` sign is
negative: lines that have already made neuroblasts by day 11 yield *fewer*
dopaminergic neurons at day 52. Consistent with this, **`PC1` at D11** (5.4% of
HVG variance, the leading expression axis) is a **cell-cycle axis** — its top
loadings are textbook G2/M (`HMGB2 PTTG1 NUSAP1 CENPF TOP2A MKI67`; G2/M set at
median rank 1972 of 2000) — and it tracks the *proliferating* progenitor
fraction at **+0.83**. Lines that keep a cycling progenitor pool at day 11 do
better; lines that differentiate early do worse.

**`PC1` at D30** (9.8% of HVG variance) **is the inverted axis: it measures
neuronal maturation, not proliferation.** The cell cycle moves to the opposite
pole and the leading genes are pan-neuronal
(`MLLT11 TUBB2B STMN2 GAP43 MAPT`). Decomposing it shows it is
**composition-driven** (R² 0.891 from composition vs 0.484 from within-type
state).

`PC1` is the only component read this way: it replicates across PCA fits at
r = 0.9995, whereas the lower components are near-tied in variance. Note it is
*not* the same component that carries the D11 prediction — that is `PC3`
([§7](#7-which-features-matter--full-tables)). **Sign is arbitrary**, so poles
are named by their own gene content.

**About 31% of D30 cells have genuinely not decided yet** — floor-plate
progenitors, SOX2⁺/HES1⁺/VIM-high — against 48% already-arrived neurons and 20%
off-target terminal fates. So a real window exists. But looking *only* at those
undecided cells, ignoring the neurons a line has already made, predicts the
outcome at AUC 0.727 — better than a coin flip, far worse than simply counting
the neurons. **By day 30 most of the answer is in what a line has already become,
not in what is left to decide.**

# Details

Result files carry the `_da_untreated` suffix, meaning **`DA / all D52 cells`,
untreated cells only** — see
[`modeling/results/README.md`](modeling/results/README.md) for a file-by-file
guide. The earlier DA+Sert analysis is in the [appendix](#appendix-a--phase-1-the-authors-dasert-outcome);
its files are unchanged and byte-identical.

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

## 4. Full performance tables

Identical design at both timepoints: donor-grouped 5-fold × 10 repeats, nested
tuning, same model grids, same selection rule. Both were scored on fold tables
verified to agree on **every** fold's train/test membership (34,816 rows
compared), so the comparison is paired and any difference is attributable to the
timepoint alone.

Every fitted model, every CV scheme, nested tuning throughout. Mean ± SD across
10 repeats for `plain` and `donor_grouped`; `loco` (136 folds) and `lodo` (20
folds) are single repeats by construction, so they have no SD. **★ marks the best
value in each column** — lowest for Brier, MAE, RMSE and out-of-range, highest
elsewhere.

No row is privileged. The configuration that would ship is named in
[§1](#the-model-that-would-be-deployed); it is not necessarily the best cell in
any column here, and these tables are laid out so you can see that for yourself.

The baseline row is what a trivial model scores: always predicting "success" for
classification, the training mean for regression. **Sensitivity 1.000 and F1
0.619 come free from that**, which is why neither is ever quoted alone.

### D11 → D52 — classification

| model | scheme | ROC-AUC | PR-AUC | balanced acc | sensitivity | specificity | F1 | Brier |
|---|---|---|---|---|---|---|---|---|
| `logistic_l2` | plain | 0.902 ± 0.010 | 0.864 ± 0.010 | 0.837 ± 0.015 | 0.856 ± 0.034 | 0.817 ± 0.019 | 0.822 ± 0.017 | 0.129 ± 0.007 |
| `logistic_l2` | donor_grouped | 0.906 ± 0.008 | 0.874 ± 0.012 | 0.834 ± 0.021 | 0.864 ± 0.031 | 0.804 ± 0.026 | 0.821 ± 0.022 | 0.125 ± 0.006 |
| `logistic_l2` | loco | 0.912 ★ | 0.883 | 0.861 ★ | 0.869 | 0.853 ★ | 0.848 ★ | 0.121 ★ |
| `logistic_l2` | lodo | 0.908 | 0.886 | 0.820 | 0.852 | 0.787 | 0.806 | 0.126 |
| `logistic_l1` | plain | 0.908 ± 0.016 | 0.873 ± 0.018 | 0.839 ± 0.025 | 0.844 ± 0.034 | 0.833 ± 0.022 | 0.824 ± 0.027 | 0.123 ± 0.012 |
| `logistic_l1` | donor_grouped | 0.905 ± 0.011 | 0.872 ± 0.012 | 0.837 ± 0.016 | 0.872 ± 0.017 ★ | 0.803 ± 0.025 | 0.825 ± 0.016 | 0.126 ± 0.009 |
| `logistic_l1` | loco | 0.907 | 0.877 | 0.861 ★ | 0.869 | 0.853 ★ | 0.848 ★ | 0.123 |
| `logistic_l1` | lodo | 0.906 | 0.886 ★ | 0.821 | 0.869 | 0.773 | 0.809 | 0.126 |
| *trivial baseline* | — | *0.500* | *0.449* | *0.500* | *1.000* | *0.000* | *0.619* | *0.247* |

### D30 → D52 — classification

| model | scheme | ROC-AUC | PR-AUC | balanced acc | sensitivity | specificity | F1 | Brier |
|---|---|---|---|---|---|---|---|---|
| `logistic_l2` | plain | 0.938 ± 0.015 | 0.910 ± 0.020 | 0.884 ± 0.012 | 0.885 ± 0.020 | 0.883 ± 0.016 | 0.872 ± 0.013 | 0.092 ± 0.007 |
| `logistic_l2` | donor_grouped | 0.945 ± 0.008 | 0.906 ± 0.019 | 0.877 ± 0.013 | 0.880 ± 0.022 | 0.875 ± 0.009 | 0.865 ± 0.015 | 0.090 ± 0.006 |
| `logistic_l2` | loco | 0.945 ★ | 0.886 | 0.876 | 0.885 | 0.867 | 0.864 | 0.087 ★ |
| `logistic_l2` | lodo | 0.938 | 0.878 | 0.884 ★ | 0.902 ★ | 0.867 | 0.873 ★ | 0.090 |
| `logistic_l1` | plain | 0.934 ± 0.014 | 0.908 ± 0.022 | 0.881 ± 0.011 | 0.872 ± 0.020 | 0.891 ± 0.014 ★ | 0.869 ± 0.013 | 0.094 ± 0.008 |
| `logistic_l1` | donor_grouped | 0.940 ± 0.009 | 0.899 ± 0.023 | 0.870 ± 0.013 | 0.861 ± 0.018 | 0.879 ± 0.015 | 0.856 ± 0.014 | 0.093 ± 0.005 |
| `logistic_l1` | loco | 0.943 | 0.914 ★ | 0.883 | 0.885 | 0.880 | 0.871 | 0.090 |
| `logistic_l1` | lodo | 0.939 | 0.873 | 0.876 | 0.885 | 0.867 | 0.864 | 0.091 |
| *trivial baseline* | — | *0.500* | *0.449* | *0.500* | *1.000* | *0.000* | *0.619* | *0.247* |

### D11 → D52 — regression

| model | scheme | R² | MAE | RMSE | out-of-range |
|---|---|---|---|---|---|
| `ridge` | plain | 0.515 ± 0.018 ★ | 0.100 ± 0.002 | 0.133 ± 0.002 ★ | 0.056 ± 0.010 |
| `ridge` | donor_grouped | 0.503 ± 0.015 | 0.099 ± 0.001 | 0.134 ± 0.002 | 0.057 ± 0.008 |
| `ridge` | loco | 0.501 | 0.100 | 0.134 | 0.066 |
| `ridge` | lodo | 0.505 | 0.099 ★ | 0.134 | 0.051 |
| `lasso` | plain | 0.514 ± 0.015 | 0.100 ± 0.002 | 0.133 ± 0.002 | 0.050 ± 0.019 ★ |
| `lasso` | donor_grouped | 0.502 ± 0.023 | 0.100 ± 0.001 | 0.134 ± 0.003 | 0.051 ± 0.011 |
| `lasso` | loco | 0.507 | 0.099 | 0.134 | 0.066 |
| `lasso` | lodo | 0.513 | 0.099 | 0.133 | 0.066 |
| *trivial baseline* | — | *0.000* | *0.158* | *0.190* | *0.000* |

### D30 → D52 — regression

| model | scheme | R² | MAE | RMSE | out-of-range |
|---|---|---|---|---|---|
| `ridge` | plain | 0.817 ± 0.014 | 0.059 ± 0.002 | 0.081 ± 0.003 | 0.037 ± 0.010 |
| `ridge` | donor_grouped | 0.802 ± 0.009 | 0.061 ± 0.001 | 0.085 ± 0.002 | 0.037 ± 0.006 |
| `ridge` | loco | 0.820 | 0.058 | 0.081 | 0.044 |
| `ridge` | lodo | 0.804 | 0.060 | 0.084 | 0.051 |
| `lasso` | plain | 0.828 ± 0.005 ★ | 0.056 ± 0.001 ★ | 0.079 ± 0.001 ★ | 0.004 ± 0.005 ★ |
| `lasso` | donor_grouped | 0.811 ± 0.006 | 0.058 ± 0.001 | 0.083 ± 0.001 | 0.010 ± 0.009 |
| `lasso` | loco | 0.821 | 0.058 | 0.081 | 0.029 |
| `lasso` | lodo | 0.819 | 0.058 | 0.081 | 0.015 |
| *trivial baseline* | — | *0.000* | *0.158* | *0.190* | *0.000* |

**The CV schemes agree.** Across all four tables the spread between schemes is
small and shows no systematic penalty for holding out whole donors or lines —
nothing here hinges on one splitting strategy.

**D11's DA-only numbers are lower than its published DA+Sert ones** (R² 0.503 vs
0.653, AUC 0.906 vs 0.947). That is **not a regression** — it is a harder
target. `DA+Sert` is bimodal with a genuine gap at 0.2, so the threshold
separates two clumps; `DA/all` is unimodal with no gap, so the same threshold
cuts through a dense region. Different question, different difficulty. The two
families are not comparable and should never be differenced.

## 5. What D30 can and cannot claim

D30's stronger numbers are **expected and largely definitional**. Its annotation
set already contains DA cells, so the outcome partly exists in the predictors.
This is **construct overlap, not leakage**: D30 and D52 are separate cells from
separate harvests, D30 precedes D52, and no D52 information enters the features.

The defensible claim is *how much of the D52 phenotype is already established by
D30* — not that something was predicted. Marker analysis over all 250,923 D30
cells splits the seven types three ways:

- **already arrived** — dopaminergic + serotonergic neurons (48% of cells)
- **still undecided** — floor-plate progenitors, cycling and not (31%):
  SOX2⁺/HES1⁺/VIM-high, not yet neuronal
- **off-target, terminal** — ependymal-like (ciliated/choroid-plexus, TTR 134.7)
  and unassigned neurons (post-mitotic, no lineage markers) (20%)

So ~31% of D30 cells genuinely have not decided — but the benchmark isolating
that signal reaches only AUC 0.727, far below the maturation-driven numbers.
Detail in
[`modeling/docs/D30_celltype_interpretation.md`](modeling/docs/D30_celltype_interpretation.md).

## 6. ⚠️ The D30 model does not beat a single raw column

The fit-free reference lines, donor-level bootstrap over 20 donors, 2,000
resamples:

| D30 benchmark (no model, no fitted parameter) | ρ | ROC-AUC |
|---|---|---|
| **`phat_DA` alone** | +0.932 | **0.960** [0.922–0.986] |
| `phat_DA + phat_Sert` | +0.799 | 0.886 [0.834–0.941] |
| `phat_FPP + phat_P_FPP` (undecided) | −0.788 | 0.883 |
| progenitor balance (share of non-target) | +0.452 | 0.727 |

**The full D30 model scores 0.945. One raw D30 column scores 0.960.**

Those two are not like-for-like — the benchmark is in-sample, the model
out-of-fold — so the gap above settles nothing on its own. **It has since been
tested properly** (`modeling/test_incremental_value.py`): both feature sets put
through identical donor-grouped folds, nested inner selection on training lines
only, differences paired per repeat.

| metric | `phat_DA` alone | full model | full better in |
|---|---|---|---|
| ROC-AUC | **0.9563** ± 0.0018 | 0.9393 ± 0.0080 | **0 / 10 repeats** |
| log loss | 0.5707 | **0.3276** ± 0.0311 | **10 / 10 repeats** |
| Brier | 0.1897 | **0.0930** ± 0.0052 | **10 / 10 repeats** |
| accuracy | 0.8676 | 0.8750 | 6 / 10 |

**The answer is split, and both halves are unanimous across repeats.**
`phat_DA` alone *ranks* lines better — by a small margin (−0.017 AUC) but in
every single repeat. The full model is *calibrated* far better — log loss
−0.243, Brier −0.097, again in every repeat.

So: **to rank lines, one raw D30 column beats the entire model. To get
trustworthy probabilities, the model earns its keep.** This is precisely why
AUC could not settle the question — at 0.95 it is saturated and speaks only to
ordering.

Applies to D30 only: D11 has no `phat_DA` (its three annotated types are FPP,
NB, P_FPP), so there is no analogous single-column comparator there.

**A mis-specified benchmark nearly hid this.** The comparator was carried over
from the previous outcome and summed `phat_DA + phat_Sert` (0.886) rather than
`phat_DA` (0.960). Against the wrong reference the model appeared to add ~6 AUC
points. **A benchmark whose numerator does not match the outcome is not a weaker
check — it is a misleading one.**

## 7. Which features matter — full tables

One-SE config selection, paired per-fold LOCO and permutation deltas, full-fit
coefficients/SHAP from a PCA fit once on all 136 lines, pool η² on the same
features. Tables in `modeling/results/`.

### Recapitulates the published analysis

**`phat_NB`, the neuroblast proportion, is still D11's standout composition
feature**, with the same sign and the same technical cleanliness:

| | published (DA+Sert) | DA-only |
|---|---|---|
| `phat_NB` univariate ρ | −0.733 | −0.554 (cls) / −0.610 (reg) |
| `phat_NB` SHAP | 0.243 | **1.050** (cls) |
| `phat_NB` pool η² | 0.037 | **0.035** |

More D11 neuroblasts → worse D52 dopaminergic yield, and it remains the only
strongly predictive feature that is *not* pool-associated. The weaker
correlation is expected: DA-only is a harder target.

**Pool confounding of the PCs persists**, and in the same severity range
(η² 0.11–0.83 here, 0.45–0.76 published). **Composition still dominates the
transcriptome** on permutation importance at both timepoints.

### New

**1. At D30 the transcriptome adds essentially nothing beyond composition.**
Three of four models select **k = 0 PCs** under the one-SE rule (lasso, ridge,
logistic_l1; only logistic_l2 keeps 6). Where PCs are kept, they are dwarfed:

| D30, logistic_l2 | permutation Δ |
|---|---|
| proportions (grouped) | **0.360** |
| best single PC (PC1) | 0.018 |

A **20×** gap. This is the feature-level counterpart of §6's benchmark result —
`phat_DA` alone nearly matches the full model — from the other direction. At D11
the PCs do carry weight (top LOCO Δ +0.064).

**2. `phat_DA` dominates D30**, as the construct-overlap framing predicts:
univariate ρ +0.792 (classification) and **+0.932** (regression), SHAP 1.292,
the largest of any feature at either timepoint.

**3. D11's strongest PC under this outcome is also its most pool-confounded
feature.** `PC3` carries the largest LOCO Δ (+0.064) and SHAP (1.398) in the
D11 classification table, and η² = **0.832** — more pool-associated than
anything in the published table. So D11's DA-only performance leans on a
technically entangled component, and may not transfer across pools. Treat with
more caution than the published D11 result, not less.

**4. `phat_Sert` flips sign conditionally at D30.** Univariate ρ is **+0.215**,
but its fitted coefficient is **−0.288**. Positive alone, negative once the
other types are held fixed — consistent with the lineage separation above: Sert
is a different lineage, and at fixed composition more Sert means fewer cells
left to become DA.

**5. What PC1 actually is, at each timepoint.** Full analysis in
[`modeling/docs/PCA_interpretation.md`](modeling/docs/PCA_interpretation.md):

| | D11 PC1 (5.4% var) | D30 PC1 (9.8% var) |
|---|---|---|
| one pole | **G2/M cell cycle** — `HMGB2 PTTG1 NUSAP1 CENPF TOP2A MKI67` | **pan-neuronal** — `MLLT11 TUBB2B STMN2 GAP43 MAP1B MAPT` |
| G2/M median rank (of 2000) | **1972** vs 980 background | **432** vs 1015 — *opposite pole* |
| tracks, at line level | proliferating FPP **+0.83** | dopaminergic **+0.85** |

**The cell cycle switches poles between timepoints**, so these are closer to
inverses than to one axis seen twice. The D11 result **replicates the published
D11 finding** (ρ +0.56 with the outcome here vs +0.563 published) on a different
cohort and outcome.

D30's apparent paradox — a *generic* neuronal gene list that tracks DA at +0.85
but serotonergic at **+0.05** — resolves by decomposing the line mean:

| | R² of line-level PC1 |
|---|---|
| **composition alone** | **0.891** |
| within-type state alone | 0.484 |

**Composition-driven**, and composition-only PC1 tracks the outcome *better*
than the real thing (+0.859 vs +0.762). Sert scores **−4.4** on PC1 against
DA's **+14.4** — matching the marker evidence that D30 Sert cells are
transcriptionally immature (STMN2 4.4 vs 27.9). So PC1 orders cells by neuronal
*maturity*, and Sert has not travelled along it yet. Its near-zero correlation
is informative, not an artefact of an invariant fraction: Sert ranges
0.000–0.685 across lines.

Two things this does **not** license: PC1 is ~50% pool-explained at both
timepoints (η² 0.47 / 0.52), so nothing here is shown to be intrinsic to a line;
and D30's +0.85 DA correlation is partly *built in*, since any line mean tracks
composition whenever cell types occupy different PC regions.

### Caveat on the k = 0 rows

The `k = 0` selections mean the grouped-proportions LOCO row compares against an
**intercept-only** model (no features at all). That is the correct comparison,
and it is computed rather than skipped, but it is a different kind of baseline
than the other rows.

## 8. Why this outcome, not the published one

Two reasons, in brief (detail in the
[README footnote](README.md#-on-the-outcome-definition)):

1. **The published metric counts rotenone-treated cells** — 219,238 `ROT`
   against 303,856 `NONE`, interleaved through every qualifying combo. Rotenone
   preferentially damages dopaminergic neurons, the cells in its own numerator.
2. **It sums two lineages that track independently** — the `phat_DA` +0.932 vs
   `phat_Sert` +0.057 split shown above.

A deliberate divergence, not a bug fix: the authors' notebook applies a ≥10-cell
threshold and no treatment filter. Cost: 2 lines fall below the threshold (both
pool5), giving **136 lines / 157 combos / 20 donors** — donor count unchanged,
so donor-grouped CV keeps its structure.

## 9. Literature assessment: recapitulated vs. new

Same framing as the appendix's §5, applied to these results. The baseline
comparison is worth restating: **Jerber et al. predicted differentiation
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
| 2 | D30 → D52 dopaminergic yield: ROC-AUC 0.945, R² 0.802 | **New, but largely definitional** | No D30 → D52 predictor in their work. The strength is mostly construct overlap — D30 already contains DA cells — not predictive discovery. See §5 |
| 3 | **DA and Sert track independently**: D30 DA → D52 DA ρ = +0.932, while D30 Sert → D52 DA is Pearson **+0.057** | **New; and it undercuts the combined metric** | The authors' `diff_efficiency` sums the two, which presumes they behave as one quantity. They do not: `DA/(DA+Sert)` spreads near-uniformly 0–1 across lines, consistent with a line-intrinsic rostro-caudal identity fixed before D30 |
| 4 | `phat_NB` (neuroblast proportion) remains D11's strongest *and technically cleanest* composition feature (pool η² 0.035) | **Recapitulates phase 1, and the authors' indirect observation** | They noted a poor-differentiation cluster correlating with D11 neuroblast proportion. Phase 1 made it the top direct predictor; it survives the outcome change with the same sign |
| 5 | **At D30 the transcriptome adds nothing beyond cell-type composition** — 3 of 4 models select k = 0 PCs; permutation Δ 0.360 for proportions vs 0.018 for the best PC | **New** | A 20× gap. Not addressed by the authors, who did not build D30 predictors |

## 10. Caveats

- **No external validation.** 20 donors, one protocol, one study, no held-out
  batch. Everything here is internal cross-validation.
- **Not comparable to the published family.** Different outcome, different D52
  cells, different cohort. Never difference the two.
- **Not a causal claim** anywhere. All associations.
- **D11's strongest PC is 83% pool-explained** — the expression signal may not
  transfer to a new run. The composition features are the robust ones.
- **`pool5` annotation bias is unaddressed.** Its ~10× D30 sequencing-depth
  deficit may have biased D30 *cell-type annotation* itself; excluding it from
  HVG/PCA fitting does nothing to correct biased labels. It holds 18 of 157
  combos.
- **Shared annotation machinery unverified.** Whether D30 and D52 types were
  annotated independently is unknown. Not leakage either way, but correlated
  labelling error would inflate the apparent D30→D52 continuity.
- **The 0.2 threshold is inherited**, and means something different here. Its
  fold-level class balance was not profiled.

---

# Appendix A — phase 1: the authors' DA+Sert outcome

*Phase 1 used a **pre-registered** headline model, declared before results were
seen. That is faithfully recorded below and not rewritten. The current analysis
supersedes it with the flat-CV + one-SE selection rule described in §1.*

*An earlier phase, reproducing Jerber et al.'s own definition
(`(DA+Sert)/all D52 cells`, all cells, 138 lines, D11 only). Retained
because the comparison to their published result rests on it, and because
it is what surfaced the two problems motivating the headline analysis.*

> **Not comparable to the results above.** Different outcome, different
> cells, different cohort. `DA+Sert` is bimodal with a real gap at the 0.2
> threshold; `DA/all` is unimodal with none, so the same cut does a
> different job. A lower number above is a harder target, not a regression.

*Reproduces Jerber et al.'s own definition of differentiation efficiency, so
that results can be set against theirs. Superseded as the headline, retained
because the literature comparison depends on it.*


Summary of results from this prototype, with an assessment of which
findings **recapitulate** Jerber et al. (2021) or established stem-cell
biology, and which appear to be **new**.

Data: Jerber et al. 2021 (*Nat Genet* 53:304–312), 138 iPSC lines with
≥10 cells at D11, D30 and D52. Outcome: `diff_efficiency` = fraction of a
line's D52 cells that are DA or Sert neurons; "success" = ≥0.2, the
authors' own threshold. 96 successes / 42 failures.

Method detail lives in `modeling/README.md`; PC biology in
`modeling/docs/PCA_interpretation.md`.

---

## 1. Model performance

Nested CV, mean ± SD across 10 repeats. Model per task was
**pre-registered before looking at results** to avoid cherry-picking.
`plain` ignores donors; `donor-grouped` never lets a line from a training
donor appear in test.

**Classification**

| model | scheme | ROC-AUC | PR-AUC | balanced acc | sensitivity | specificity | F1 | Brier |
|---|---|---|---|---|---|---|---|---|
| **logistic_l2** *(headline)* | plain | **0.940 ± 0.012** | **0.964 ± 0.011** | **0.914 ± 0.013** ★ | **0.931 ± 0.016** | **0.898 ± 0.023** ★ | **0.943 ± 0.009** | **0.102 ± 0.006** |
| **logistic_l2** *(headline)* | **donor-grouped** | **0.947 ± 0.007** | **0.968 ± 0.004** | **0.913 ± 0.013** | **0.940 ± 0.007** | **0.886 ± 0.022** | **0.945 ± 0.007** ★ | **0.097 ± 0.008** |
| logistic_l1 *(secondary)* | plain | 0.948 ± 0.011 | 0.968 ± 0.012 | 0.904 ± 0.019 | 0.933 ± 0.005 | 0.874 ± 0.037 | 0.939 ± 0.009 | 0.077 ± 0.007 |
| logistic_l1 *(secondary)* | donor-grouped | 0.951 ± 0.008 ★ | 0.970 ± 0.009 ★ | 0.912 ± 0.013 | 0.942 ± 0.011 ★ | 0.883 ± 0.021 | 0.945 ± 0.008 ★ | 0.075 ± 0.004 ★ |
| *trivial baseline* | — | *0.500* | *0.696* | *0.500* | *1.000* | *0.000* | *0.821* | *0.212* |

**Regression**

| model | scheme | R² | MAE | RMSE | out-of-range |
|---|---|---|---|---|---|
| **ridge** *(headline)* | plain | **0.676 ± 0.013** | **0.134 ± 0.003** ★ | **0.170 ± 0.003** | **0.050 ± 0.007** |
| **ridge** *(headline)* | **donor-grouped** | **0.653 ± 0.021** | **0.138 ± 0.004** | **0.175 ± 0.005** | **0.043 ± 0.005** ★ |
| lasso *(secondary)* | plain | 0.680 ± 0.016 ★ | 0.134 ± 0.003 ★ | 0.169 ± 0.004 ★ | 0.049 ± 0.008 |
| lasso *(secondary)* | donor-grouped | 0.668 ± 0.015 | 0.135 ± 0.003 | 0.172 ± 0.004 | 0.046 ± 0.004 |
| *predict the mean* | — | *0.000* | *0.266* | *0.298* | *0.000* |

★ = best in that column among the four fitted rows (baselines excluded).
**Bold** marks the pre-registered headline row, which is a separate thing
from being best. Higher is better everywhere except **Brier, MAE, RMSE and
out-of-range**, where lower is better — so the star follows the minimum in
those columns.

Baseline row = what a trivial model scores: always predicting "success"
for classification (or random ranking, for the two AUCs), and predicting
the training mean for regression. **Sensitivity 1.000 and F1 0.821 from
that baseline are why neither is quoted alone** — with 70% successes,
guessing "success" every time already looks respectable on both.

**The PCA fitting population barely matters.** Every figure here uses the
baseline basis. Refitting the whole pipeline on cells from the 138 study
lines only moves the donor-grouped results by +0.014 R² (ridge), +0.004
(lasso), +0.001 AUC (L2), −0.001 (L1) — all inside about one SD. See §3.

**Best model: L2-penalised logistic regression on the binary task.**
The L1 variants score marginally higher but the difference is well inside
one SD, and they were not pre-registered.

**Classification substantially outperforms regression, and this is a
property of the outcome, not the models.** `diff_efficiency` is bimodal —
42 failures in 0–0.185, 96 successes in 0.221–0.92, with a real gap at the
threshold. Splitting predictions by true clump: within-success R² is only
~0.17–0.22 and within-failure R² is strongly negative (≈ −14). **The
headline R² mostly reflects correctly separating the two clumps, not
fine-grained precision.** Use this model to call success/failure; do not
use its regression output as a precise efficiency estimate.

Robustness across CV schemes (nested, ridge / logistic_l2):

| scheme | regression R² | classification ROC-AUC |
|---|---|---|
| LOCO (leave-one-line-out) | 0.682 | 0.942 |
| plain 5-fold | 0.676 | 0.940 |
| LODO (leave-one-donor-out) | 0.665 | 0.963 |
| donor-grouped 5-fold | 0.653 | 0.947 |

Donor leakage is real but small (~0.012 R²), consistent across two
independent contrasts (plain→donor-grouped, LOCO→LODO). Classification is
essentially flat across schemes.

---

## 2. Which features are informative

Assessed five ways — univariate association, standardised coefficient,
leave-one-covariate-out, permutation, and SHAP — so that agreement or
disagreement between measures is itself informative.

**`phat_NB`, the day-11 neuroblast fraction, is the standout feature.**

- Strongest univariate association of anything measured: Spearman
  **ρ = −0.73** (both tasks). More neuroblasts at D11 → *worse* D52 yield.
- Largest single LOCO contribution for classification (+0.051 AUC).
- Largest SHAP value for classification (1.19 vs PC2's 0.51).
- **By far the cleanest technically**: pool η² = 0.04, the only "Low"
  feature in either table.

Dropping all cell-type proportions costs **+0.064 AUC** — the largest
single effect in either table — but costs regression essentially nothing
(−0.0007). On the restricted basis, the modal nested `logistic_l1` model
selects **k=0 PCs** — proportions alone — and still reaches AUC 0.950.

**PC1 is a proliferation axis and is genuinely predictive** (ρ = +0.563).
Quantitatively enriched for the Tirosh/Seurat G2/M set (30 of the top 50
positive loadings; Mann-Whitney p = 2.3 × 10⁻²⁶), and it tracks the
*proliferating* FPP fraction at ρ = +0.798. More proliferative D11
cultures differentiate better.

**LOCO and permutation disagree for the regression proportions, and the
disagreement is the finding**: permutation Δ = 0.035 (the fitted model
does use them) but LOCO Δ ≈ 0 (refitting without them costs nothing) —
the signature of information that is redundantly available in the PCs.
For classification both measures agree the proportions are essential.
Running only one measure would have misled in either direction.

---

## 3. Technical confounders

**Pool identity explains a large share of the variance in most PC
features** (η² of the line-level value on pool):

| feature | pool η² | |
|---|---|---|
| `phat_NB` | **0.04** | clean |
| PC4 | 0.45 | moderate |
| PC1 | 0.49 | moderate |
| PC3 | 0.50 | moderate |
| **PC2** | **0.76** | **high — and the single most predictive PC** |
| PC5 | 0.79 | high |
| PC6 / PC9 / PC8 | 0.91 / 0.95 / **0.96** | near-pure batch |

Two consequences:

1. **PC2 is model-predictive but technically entangled.** It has the
   largest coefficient, LOCO delta and permutation delta of any PC for
   regression, while ~76% of its line-level variance is pool-associated.
   The correct statement is *"PC2 improves prediction under the observed
   design"*, not *"PC2 is biologically predictive"*. Its contribution may
   not transfer to new pools or altered processing.
2. **The excluded PCs are excluded for good reason.** PC6/PC8/PC9 are
   nearly pure batch signal; the tuned `k` never reaches them.

Other technical decisions, frozen from EDA before any modelling:

- **pool11 excluded from every PCA/HVG fit** (≈1,900 mean UMI/cell vs
  10,000–18,000 elsewhere) but its cells are still *projected*, so no
  line is lost from train or test.
- **Fitting population matters.** Restricting the PCA fit to the study
  population (dropping ~25,400 cells from 39 lines outside the 138)
  improved regression R² by +0.014 and roughly **halved** its
  across-repeat SD (0.021 → 0.012). Classification was unaffected.

---

## 4. Additional findings worth recording

**Hyperparameter selection is unstable at n=138.** Across 50 outer folds
the modal nested configuration is chosen only 9–10 times, spread over
13–19 distinct configurations, with `C` spanning four orders of magnitude.
No single configuration is "the" model — which is why the headline is
pre-registered and a one-SE rule is used, and why any single-config
coefficient table is one draw from a wide distribution.

---

## 5. Literature assessment: recapitulated vs. new

The single most important comparison: **Jerber et al. predicted
differentiation efficiency from *iPSC-stage bulk RNA-seq*, before
differentiation began — not from day-11 scRNA-seq.** The D11 → D52
prediction here is therefore a different analysis on the same data, not a
reproduction of theirs.

| # | Finding | Status | Basis |
|---|---|---|---|
| 1 | D11 scRNA-seq features predict D52 efficiency, ROC-AUC 0.947 | **New** | Authors predicted from iPSC-stage bulk RNA-seq; they did not build a D11 → D52 predictor |
| 2 | **100% precision at 65% recall** for detecting poor differentiation | **Extends prior result** | Authors reported 100% precision at **35%** recall from iPSC bulk under LOOCV. Ours is ~2× the recall under a *stricter* (donor-grouped) scheme — see caveat below |
| 3 | D11 neuroblast fraction predicts *worse* D52 outcome (ρ = −0.73) | **Recapitulates, and sharpens** | Authors noted their poor-differentiation-associated "cluster 2" correlated with D11 neuroblast proportion — an indirect observation. Here `phat_NB` is the single strongest direct predictor and the most technically clean feature |
| 4 | Proliferation (PC1) associated with better differentiation | **Consistent with established biology; new in this system** | Cell-cycle state governs fate propensity in hPSCs (Pauklin & Vallier 2013); cell-cycle lengthening precedes neurogenesis. Not previously quantified as a D11 → D52 predictor here |
| 5 | PC3 = neuronal differentiation vs. S-phase/replication | **Recapitulates textbook biology** | Cell-cycle exit accompanying neurogenesis is long established. Best read as a *positive control* that the PCA captures real biology |
| 6 | Classification ≫ regression because the outcome is bimodal | **New (methodological)** | Follows from the authors' own 0.2 threshold but is not analysed by them |
| 7 | D11 expression PCs are strongly pool-confounded (η² up to 0.96) | **New; complements the authors** | They showed the *outcome* is not associated with lines-per-pool (R² = 0.04, P = 0.46) and that line effects dominate. That is about the outcome; this is about the *features*, where pool structure is large |
| 8 | Donor leakage inflates R² by only ~0.012 | **New (methodological)** | Authors used LOOCV without donor grouping; this quantifies what that costs |

### The important caveat on finding #2

It is **not** a like-for-like "we beat them" comparison, and the practical
implication runs the other way:

- Their predictor uses **iPSC-stage** expression, available *before*
  committing to a differentiation run — the practically useful screening
  point.
- Ours requires **11 days of culture** before it can say anything.

So theirs answers *"is this line worth differentiating?"*; ours answers
*"is this run going to work?"* — a different, later, cheaper-to-abort
question. The gain in accuracy is real but is bought with 11 days.

---

## 6. What would strengthen this most

1. **Leave-one-pool-out validation.** The decisive test of whether PC2's
   contribution is line biology or batch. Nothing here separates them.
2. **Within-pool permutation importance** instead of global — global
   permutation can inflate importance by creating PC×pool combinations
   absent from the data.
3. **The D30 → D52 baseline**, to establish how much of the D11 signal is
   simply "differentiation is already visibly on track".
4. **Direct per-cell program scoring** (S-phase, G2/M, neuronal) verified
   *within* pools, rather than enrichment among extreme loadings.

## 7. What this does not show

- No causal claim. Nothing here shows that increasing D11 proliferation,
  or suppressing neuroblast formation, would improve D52 yield.
- No claim that individual top-loading genes are drivers or targets.
- No demonstrated transfer to new pools, protocols or laboratories.
- The regression model does not predict fine-grained efficiency (§1).

---

---

## Sources

1. **Jerber J, Seaton DD, Cuomo ASE, et al. (2021)** "Population-scale
   single-cell RNA-seq profiling across dopaminergic neuron
   differentiation." *Nature Genetics* 53:304–312.
   <https://www.nature.com/articles/s41588-021-00801-6> ·
   <https://pmc.ncbi.nlm.nih.gov/articles/PMC7610897/> — source of the
   data, the 0.2 threshold, the iPSC-stage predictor (100% precision at
   35% recall), the cluster-2/neuroblast observation, and the
   lines-per-pool null result.
2. **Tirosh I, Izar B, Prakadan SM, et al. (2016)** "Dissecting the
   multicellular ecosystem of metastatic melanoma by single-cell RNA-seq."
   *Science* 352(6282):189–196.
   <https://www.science.org/doi/10.1126/science.aad0501> — G1/S and G2/M
   gene sets used for the PC1 enrichment test.
3. **Seurat `cc.genes`** <https://satijalab.org/seurat/reference/cc.genes>
   — the exact gene lists matched against.
4. **Whitfield ML, Sherlock G, Saldanha AJ, et al. (2002)**
   "Identification of genes periodically expressed in the human cell cycle
   and their expression in tumors." *Mol Biol Cell* 13(6):1977–2000.
   <https://pubmed.ncbi.nlm.nih.gov/12058064/> — foundational human
   cell-cycle periodicity map underlying the G2/M gene identities.
5. **Pauklin S, Vallier L (2013)** "The cell-cycle state of stem cells
   determines cell fate propensity." *Cell* 155(1):135–147.
   <https://pmc.ncbi.nlm.nih.gov/articles/PMC3898746/> — establishes that
   cell-cycle state governs differentiation propensity in hPSCs; the prior
   for finding #4.
6. **Hardwick LJA, Ali FR, Azzarelli R, Philpott A (2014)** "Cell cycle
   regulation of proliferation versus differentiation in the central
   nervous system." *Cell Tissue Research* 359(1):187–200.
   <https://pmc.ncbi.nlm.nih.gov/articles/PMC4284380/> — G1 lengthening
   and cell-cycle exit preceding neuronal differentiation; the prior for
   finding #5.
7. **Luecken MD, Theis FJ (2019)** "Current best practices in single-cell
   RNA-seq analysis: a tutorial." *Mol Syst Biol* 15(6):e8746.
   <https://doi.org/10.15252/msb.20188746> — QC covariates and technical
   confounders; basis for the §3 caveats.
8. **Huang Y, McCarthy DJ, Stegle O (2019)** "Vireo: Bayesian
   demultiplexing of pooled single-cell RNA-seq data without genotype
   reference." *Genome Biology* 20:273.
   <https://pmc.ncbi.nlm.nih.gov/articles/PMC6909514/> — pooled/multiplexed
   design context: pooling controls *within-pool* batch variation but
   leaves pool-to-pool structure, which is what §3 measures.

---

---

# Appendix B — feature importance and run association, in detail

*Supporting detail for [§2](#2-which-features-predict-and-how-much-do-they-track-the-batch).
Kept out of the main narrative because it is a methods comparison, not a finding.
Regenerate the loadings with `uv run python -m modeling.pc1_loadings_analysis`;
importance columns come from `modeling/results/feature_importance_table_*.csv`.*

**Terminology.** "Run-associated" (or pool-associated) throughout, not
"contaminated" — a high η² shows a feature varies with the differentiation run,
which is not the same as showing it is a technical artefact. Nothing here
demonstrates technical contamination.

## Every D30 feature

| D30 feature | SHAP | univariate ρ | LOCO Δ | pool η² |
|---|---|---|---|---|
| `phat_DA` | **1.292** | **+0.792** | −0.002 | 0.315 |
| `PC1` | 0.653 | +0.680 | +0.001 | **0.518** |
| `phat_Epen1` | 0.652 | −0.614 | +0.000 | 0.133 |
| `PC5` | 0.610 | +0.557 | +0.002 | **0.518** |
| `PC4` | 0.530 | +0.182 | +0.001 | 0.376 |
| `PC6` | 0.499 | −0.201 | +0.002 | 0.406 |
| `phat_FPP` | 0.470 | −0.580 | −0.001 | 0.211 |
| `phat_U_Neur2` | 0.406 | −0.309 | +0.000 | 0.083 |
| `phat_U_Neur1` | 0.325 | −0.454 | −0.002 | 0.110 |
| `phat_Sert` | 0.246 | +0.214 | −0.001 | 0.057 |
| `PC2` | 0.151 | −0.213 | −0.001 | 0.076 |
| `PC3` | 0.138 | −0.301 | −0.001 | 0.143 |
| `phat_P_FPP` | *reference* | −0.482 | +0.006 (grouped) | 0.140 |

**All D30 LOCO Δ are ≈ 0.** With six correlated proportions plus six PCs, dropping
any single feature costs nothing measurable — no individual feature is *necessary*,
which is why the grouped permutation Δ (0.360 for all proportions together) is the
informative statistic at D30, not the per-feature column.

## Why the low-η² D30 components contribute nothing

Their gene content is perfectly legible — they are simply re-measuring the
annotation:

| D30 component | η² | top loadings | tracks |
|---|---|---|---|
| `PC2` (5.7% var) | 0.076 | `C11orf88 FAM183A ZMYND10 ROPN1L RSPH1 TPPP3 PIFO CCDC146` — motile cilia, the ependymal signature | `phat_Epen1` +0.706, `phat_Sert` **−0.861** |
| `PC3` (3.6% var) | 0.143 | `HMGB2 BIRC5 NUSAP1 TUBA1B H2AFZ VIM SPARC` — cell cycle (G2/M median rank 1922 of 2000, p = 6.6 × 10⁻¹⁸) | `phat_FPP` +0.815, `phat_Sert` −0.883 |

Once the proportions are in the model, a component that restates them adds
nothing. **They are redundant, not uninformative** — a distinction worth keeping,
because it means the D30 annotation is capturing the same structure the
transcriptome shows.

## The two timepoints side by side

| | D11 | D30 |
|---|---|---|
| least run-associated important feature | `phat_NB`, η² 0.035 | `phat_DA`, η² 0.315 |
| most important feature | `PC3` (SHAP 1.398, η² 0.832) | `phat_DA` (SHAP 1.292, η² 0.315) |
| a low-η² PC that also predicts? | **yes** — `PC2` (η² 0.114, ρ −0.626) | **no** |
| PCs in the shipped model | 4 of 10 | 3 of 4 families select zero |
| composition vs expression | comparable | **0.360 vs 0.018** — 20× |

The contrast is real but should not be over-read. It says expression adds
something at D11 and little at D30, which follows from D30's annotation already
containing the outcome's numerator. It does **not** say D11's expression features
would survive a new run — that remains untested at both timepoints.
