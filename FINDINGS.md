# Findings

Two analyses of the same data (Jerber et al. 2021), run in sequence. **The
second is the headline result; the first is how we got there.**

| | **Phase 1 — exploratory** | **Phase 2 — headline** |
|---|---|---|
| outcome | `(DA + Sert) / all D52 cells` | **`DA / all D52 cells`** |
| D52 cells counted | all, incl. rotenone-treated | **untreated only** |
| cohort | 138 lines | **136 lines** |
| timepoints | D11 → D52 | **D11 → D52 and D30 → D52** |
| purpose | reproduce the authors' own outcome definition, so results are comparable to theirs | answer the question actually of interest: dopaminergic yield |

**Why phase 1 first.** Adopting the authors' exact definition — including the
choices we later departed from — was what made it possible to check this
pipeline against a published result at all. It established that the D11 → D52
prediction works, calibrated what "good" looks like, and surfaced the two
problems that motivated phase 2. It is kept in full rather than deleted,
because the comparison to the literature rests on it.

**Why phase 2 is the headline.** The published outcome (a) counts
rotenone-treated cells, a complex-I inhibitor that preferentially kills the
very neurons in its numerator, and (b) sums two distinct lineages that we show
track independently. Dopaminergic yield is the quantity the protocol is
actually for.

**Read phase 2 first** ([§P2](#phase-2--headline-dopaminergic-only-efficiency)).
Phase 1 follows as context.

> **Numbers from the two phases are not comparable.** Different outcome,
> different cells, different cohort. A lower phase-2 number is a harder target,
> not a regression. See §P2.2.

---

# Phase 1 — exploratory: the authors' DA+Sert outcome

*Reproduces Jerber et al.'s own definition of differentiation efficiency, so
that results can be set against theirs. Superseded as the headline by phase 2,
retained because the literature comparison depends on it.*


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

**The PC2–PC3 plane is stable; the individual axes are not.** PC2 and
PC3 are near-tied in explained variance (0.0265 / 0.0259) and *reorder*
between PCA fitting variants: baseline PC2 ↔ qualifying-only PC3
(|r| = 0.985), while qualifying-only PC2 is a different component
(|r| = 0.564). But the *plane* they span is essentially identical
(principal angles 6.0° and 1.9°). **"PCn" is not a stable label; align
components by correlation before interpreting.** PC1 is the exception —
it replicates at r = +0.9995 with 25/25 top-gene overlap.

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
| 9 | Hyperparameter selection unstable; PC2/PC3 rotate between fits | **New (methodological)** | Caution for anyone interpreting a single PCA solution from n≈138 |

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

# Phase 2 — headline: dopaminergic-only efficiency

Phase 1 above uses the published outcome — `(DA + Sert) / all D52 cells`,
counting every D52 cell. Phase 2 covers a second analysis with a **different target**, run alongside
rather than replacing it. The published
results are unchanged and every file backing them is byte-identical.

Results live in `modeling/results/*_da_untreated.csv`; that suffix means
**`DA / all D52 cells`, untreated cells only**. See
[`modeling/results/README.md`](modeling/results/README.md) for the file-by-file
guide.

## P2.1 Why a different outcome

**The published outcome counts rotenone-treated cells.** D52 — and only D52 —
contains 219,238 `ROT` cells against 303,856 `NONE`, interleaved so that *every*
qualifying `(cell_line, pool)` combo contains both. Rotenone inhibits
mitochondrial complex I and preferentially damages dopaminergic neurons, the
cells in the numerator. Checked against the authors' own notebook: they apply
only a ≥10-cell threshold and never filter on treatment. **This is therefore a
deliberate divergence from the published definition, not a bug fix.**

**It also merges two lineages.** DA and Sert are both floor-plate-derived but
differ by rostro-caudal position, and they track *separately*:

| D30 predictor | → D52 `DA/all` |
|---|---|
| `phat_DA` | ρ = **+0.932** |
| `phat_Sert` | Pearson **+0.057** |

D30 Sert says essentially nothing about D52 DA. Summing them averages over a
real biological axis. Consistent with this, `DA/(DA+Sert)` spreads almost
uniformly 0–1 across lines — rostro-caudal identity looks line-intrinsic and
fixed before D30.

**Cost of the change:** excluding ROT cells drops 2 lines below the ≥10-cell
threshold (`HPSI0115i-melw_1`, `HPSI0115i-qecv_2` — both pool5, both reduced to
6 and 8 untreated D52 cells), giving **136 lines, 157 combos, 20 donors**.
Donor count is unchanged, so donor-grouped CV keeps its structure.

## P2.2 Results

Identical design to the published work: donor-grouped 5-fold × 10 repeats,
nested tuning, same model grids, same pre-registered headline. Both timepoints
were scored on fold tables verified to agree on **every** fold's train/test
membership (34,816 rows compared), so the comparison is paired and any
difference is attributable to the timepoint alone.

| | **D11 → D52** | **D30 → D52** |
|---|---|---|
| ridge — R² | 0.503 ± 0.015 | **0.802 ± 0.009** |
| ridge — MAE / RMSE | 0.099 / 0.134 | 0.061 / 0.085 |
| logistic_l2 — ROC-AUC | 0.906 ± 0.008 | **0.945 ± 0.008** |
| logistic_l2 — PR-AUC | 0.874 | 0.906 |
| logistic_l2 — balanced acc | 0.834 | 0.877 |
| logistic_l2 — Brier | 0.125 | 0.090 |

**D11's DA-only numbers are lower than its published DA+Sert ones** (R² 0.503 vs
0.653, AUC 0.906 vs 0.947). That is **not a regression** — it is a harder
target. `DA+Sert` is bimodal with a genuine gap at 0.2, so the threshold
separates two clumps; `DA/all` is unimodal with no gap, so the same threshold
cuts through a dense region. Different question, different difficulty. The two
families are not comparable and should never be differenced.

## P2.3 ⚠️ The D30 model does not beat a single raw column

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

## P2.4 What D30 can and cannot claim

D30's stronger numbers are **expected and largely definitional**. Its annotation
set already contains DA cells, so the outcome partly exists in the predictors.
This is **construct overlap, not leakage**: D30 and D52 are separate cells from
separate harvests, D30 precedes D52, and no D52 information enters the features.

The defensible claim is *how much of the D52 phenotype is already established by
D30* — not that something was predicted. Marker analysis over all 250,923 D30
cells splits the seven types three ways:

- **already arrived** — DA, Sert (48% of cells)
- **still undecided** — FPP, P_FPP (31%): SOX2⁺/HES1⁺/VIM-high progenitors, not
  yet neuronal, P_FPP still cycling
- **off-target, terminal** — Epen1 (ciliated/choroid-plexus-like, TTR 134.7),
  U_Neur1/2 (post-mitotic, no lineage markers) (20%)

So ~31% of D30 cells genuinely have not decided, and what they become is real
prediction — but the progenitor-balance benchmark that isolates it reaches only
AUC 0.727, well below the maturation-driven numbers. Full detail in
[`modeling/docs/D30_celltype_interpretation.md`](modeling/docs/D30_celltype_interpretation.md).

## P2.5 Caveats specific to phase 2

- **Not comparable to the published family.** Different outcome, different D52
  cells, different cohort.
- **`pool5` annotation bias is unaddressed.** Its ~10× D30 sequencing-depth
  deficit may have biased D30 *cell-type annotation* itself; excluding it from
  HVG/PCA fitting does nothing to correct biased labels. It holds 18 of 157
  combos.
- **Shared annotation machinery unverified.** Whether D30 and D52 types were
  annotated independently is unknown. Not leakage either way, but correlated
  labelling error would inflate the apparent D30→D52 continuity.
- **The 0.2 threshold is inherited**, and means something different here. Its
  fold-level class balance was not profiled.
- **Feature importance was not computed** for this family — `feature_importance.py`
  raises for any non-published variant, because its full-fit columns come from a
  PCA basis fit under the published cohort and mixing bases is a known defect.
- **No external validation.** 20 donors, one protocol, one study.

## P2.6 Which features matter — and how it differs from the published analysis

Same method as §2: one-SE config selection, paired per-fold LOCO and
permutation deltas, full-fit coefficients/SHAP from a PCA fit once on all 136
qualifying lines (`005 --timepoint {D11,D30} --label-variant da_untreated`),
pool η² on those same features. Tables:
`modeling/results/feature_importance_table_*_{D11,D30}_qualonly_da_untreated.csv`.

### Recapitulates the published analysis

**`phat_NB` is still D11's standout composition feature**, with the same sign
and the same technical cleanliness:

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

A **20×** gap. This is the feature-level counterpart of the benchmark result in
§A3 — `phat_DA` alone nearly matches the full model — and it says the same
thing from the other direction. At D11 the PCs do carry weight (top LOCO Δ
+0.064).

**2. `phat_DA` dominates D30**, as the construct-overlap framing predicts:
univariate ρ +0.792 (classification) and **+0.932** (regression), SHAP 1.292,
the largest of any feature at either timepoint.

**3. PC indices rotated between cohorts — concretely, not just in principle.**
The published PC2 (ρ +0.593, η² 0.764) does **not** correspond to the DA-only
PC2 (ρ −0.626, η² 0.114). The heavily pool-associated component is now **PC3**
(η² 0.832). `modeling/README.md` documents that PC2/PC3 are near-tied in
variance and rotate between fits; this is that caveat materialising. **`PCn` is
a slot, not a stable biological entity** — do not compare PC indices across
these tables or against the published ones.

**4. D11's strongest PC under this outcome is also its most pool-confounded
feature.** `PC3` carries the largest LOCO Δ (+0.064) and SHAP (1.398) in the
D11 classification table, and η² = **0.832** — more pool-associated than
anything in the published table. So D11's DA-only performance leans on a
technically entangled component, and may not transfer across pools. Treat with
more caution than the published D11 result, not less.

**5. `phat_Sert` flips sign conditionally at D30.** Univariate ρ is **+0.215**,
but its fitted coefficient is **−0.288**. Positive alone, negative once the
other types are held fixed — consistent with §A1: Sert is a different lineage,
and at fixed composition more Sert means fewer cells left to become DA.

### Caveats

The `k = 0` selections mean the grouped-proportions LOCO row compares against an
**intercept-only** model (no features at all). That is the correct comparison,
and it is computed rather than skipped, but it is a different kind of baseline
than the other rows.

SHAP and coefficients come from a single full-fit basis; LOCO and permutation
average over 256 differently-rotated per-fold bases. For **PC1** that is fine
(it replicates at r = 0.9995). For the near-tied components the two column
families do not strictly refer to the same direction — the rows are kept
per-PC to match the published table's shape, but should be read as *a* PC of
roughly that rank, not a fixed axis.

## P2.7 Literature assessment: recapitulated vs. new

Same framing as §5, applied to the dopaminergic-only results. The baseline
comparison is unchanged and worth restating: **Jerber et al. predicted
differentiation efficiency from *iPSC-stage bulk RNA-seq*, before
differentiation began.** Nothing here reproduces that. These are different
predictors, a different outcome, and in one case a different timepoint.

**Scope limit, stated up front.** "New" below means *not present in the
authors' analysis that we actually read* — their `Figure_2` notebook and the
outcome definition it encodes. We did not audit the full paper or its
supplements. Treat "new" as "not found where we looked", not as a priority
claim.

| # | Finding | Status | Basis |
|---|---|---|---|
| 1 | D11 scRNA-seq predicts D52 **dopaminergic** yield: ROC-AUC 0.906, R² 0.503 | **New** | The authors predicted from iPSC-stage bulk, and their efficiency metric sums DA with Sert; a DA-only D11 → D52 predictor is not part of their analysis |
| 2 | D30 → D52 dopaminergic yield: ROC-AUC 0.945, R² 0.802 | **New, but largely definitional** | No D30 → D52 predictor in their work. The strength is mostly construct overlap — D30 already contains DA cells — not predictive discovery. See §P2.4 |
| 3 | **DA and Sert track independently**: D30 DA → D52 DA ρ = +0.932, while D30 Sert → D52 DA is Pearson **+0.057** | **New; and it undercuts the combined metric** | The authors' `diff_efficiency` sums the two, which presumes they behave as one quantity. They do not: `DA/(DA+Sert)` spreads near-uniformly 0–1 across lines, consistent with a line-intrinsic rostro-caudal identity fixed before D30 |
| 4 | The published efficiency metric **counts rotenone-treated cells** | **New (methodological)** | Verified directly against their notebook: a ≥10-cell threshold and a sum of two fractions, with no treatment filter. All 159 qualifying (line, pool) combos contain both treated and untreated cells, so every line is affected. Rotenone inhibits complex I and preferentially damages DA neurons — the numerator |
| 5 | `phat_NB` remains D11's strongest *and technically cleanest* composition feature (pool η² 0.035) | **Recapitulates phase 1, and the authors' indirect observation** | They noted a poor-differentiation cluster correlating with D11 neuroblast proportion. Phase 1 made it the top direct predictor; it survives the outcome change with the same sign |
| 6 | **At D30 the transcriptome adds nothing beyond cell-type composition** — 3 of 4 models select k = 0 PCs; permutation Δ 0.360 for proportions vs 0.018 for the best PC | **New** | A 20× gap. Not addressed by the authors, who did not build D30 predictors |
| 7 | One raw column (`phat_DA`) **ranks** better than the full model (AUC 0.956 vs 0.939, 10/10 repeats), while the full model is **calibrated** far better (log loss 0.328 vs 0.571, 10/10) | **New (methodological)** | Both through identical folds and nested selection. The split answer is the finding: discrimination and calibration disagree, and AUC alone would have reported the wrong conclusion in either direction |
| 8 | ~31% of D30 cells are still uncommitted progenitors (FPP/P_FPP: SOX2⁺/HES1⁺/VIM-high, P_FPP cycling); Epen1 is ciliated choroid-plexus-like (TTR 134.7) | **Recapitulates their annotation; new marker-level quantification** | They defined and named these types. What is added is the three-way split into *arrived / undecided / off-target* and the evidence for it, which is what makes the D30 result interpretable rather than circular |
| 9 | D11's strongest PC under this outcome is **also its most pool-confounded feature** (η² 0.832) | **New; a caution, not a result** | Phase 1 found PCs pool-confounded generally (§3); here the single most important PC is the worst offender, so the DA-only D11 result may not transfer across pools |
| 10 | PC indices **rotate between cohorts**, concretely: published PC2 (ρ +0.593, η² 0.764) ≠ DA-only PC2 (ρ −0.626, η² 0.114) | **Recapitulates phase 1's caution, now demonstrated** | §5 finding 9 warned that PC2/PC3 rotate between fits. This shows it happening across two real analyses, and is why PC indices must not be compared between tables |
| 11 | Excluding treated cells costs 2 of 138 lines but leaves donors unchanged at 20 | **New (methodological)** | Both dropped lines are pool5, the lowest-yield pool at D52 (median 609 cells/line vs 3,485 in pool1) |

### What phase 2 does *not* establish

- **Not a better predictor than phase 1.** The DA-only numbers are lower
  (AUC 0.906 vs 0.947 at D11) because the target is harder, not because
  anything regressed. The two are not comparable at all.
- **Not a validated D30 biomarker.** D30's strength is substantially the
  outcome already existing in the predictors.
- **Not externally validated.** 20 donors, one protocol, one study, no held-out
  batch. Everything here is internal cross-validation.
- **Not a causal claim** anywhere. All associations.
