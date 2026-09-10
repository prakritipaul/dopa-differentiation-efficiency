# Findings: predicting D52 dopaminergic differentiation efficiency from D11

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

| task | model | metric | plain 5-fold | donor-grouped 5-fold |
|---|---|---|---|---|
| classification | **logistic_l2** *(pre-registered)* | ROC-AUC | 0.940 ± 0.012 | **0.947 ± 0.007** |
| classification | logistic_l1 *(secondary)* | ROC-AUC | 0.948 ± 0.011 | 0.951 ± 0.008 |
| regression | **ridge** *(pre-registered)* | R² | 0.676 ± 0.013 | **0.653 ± 0.021** |
| regression | lasso *(secondary)* | R² | 0.680 ± 0.016 | 0.669 ± 0.015 |

Classification, headline model, donor-grouped, full metric set: ROC-AUC
0.947, PR-AUC 0.968, balanced accuracy 0.913, sensitivity 0.940,
specificity 0.886, Brier 0.097.

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
