# Findings

**Predicting day-52 dopaminergic yield from earlier single-cell snapshots**,
136 iPSC lines, [Jerber et al. 2021](https://www.nature.com/articles/s41588-021-00801-6).

Outcome: `DA / all D52 cells`, untreated cells only. Success = ≥ 0.2
(61 / 75). Donor-grouped 5-fold × 10 repeats, nested tuning, pre-registered
headline model. Methods: [`modeling/README.md`](modeling/README.md).

| | **D11 → D52** | **D30 → D52** |
|---|---|---|
| ridge — R² | 0.503 ± 0.015 | 0.802 ± 0.009 |
| logistic_l2 — ROC-AUC | 0.906 ± 0.008 | 0.945 ± 0.008 |
| logistic_l2 — Brier | 0.125 | 0.090 |

Identical fold assignments at both timepoints, so the comparison is paired.

**In one line:** D11 predicts dopaminergic yield well; D30 predicts it better
but mostly because the answer is already visible in the predictors; and at D30 a
single raw cell-type proportion out-ranks the entire model.

---

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

---

# Results

Result files carry the `_da_untreated` suffix, meaning **`DA / all D52 cells`,
untreated cells only** — see
[`modeling/results/README.md`](modeling/results/README.md) for a file-by-file
guide. The earlier DA+Sert analysis is in the [appendix](#appendix--phase-1-the-authors-dasert-outcome);
its files are unchanged and byte-identical.

## 1. Why this outcome, not the published one

Two reasons, in brief (detail in the [README footnote](README.md#-on-the-outcome-definition)):

1. **The published metric counts rotenone-treated cells** — 219,238 `ROT`
   against 303,856 `NONE`, interleaved through every qualifying combo. Rotenone
   preferentially damages dopaminergic neurons, the cells in its own numerator.
2. **It sums two lineages that track independently:**

| D30 predictor | → D52 `DA/all` |
|---|---|
| `phat_DA` | ρ = **+0.932** |
| `phat_Sert` | Pearson **+0.057** |

`DA/(DA+Sert)` spreads near-uniformly 0–1 across lines — rostro-caudal identity
looks line-intrinsic and fixed before D30.

A deliberate divergence, not a bug fix: the authors' notebook applies a ≥10-cell
threshold and no treatment filter. Cost: 2 lines fall below the threshold
(both pool5), giving **136 lines / 157 combos / 20 donors** — donor count
unchanged, so donor-grouped CV keeps its structure.

## 2. Results

Identical design to the published work: donor-grouped 5-fold × 10 repeats,
nested tuning, same model grids, same pre-registered headline. Both timepoints
were scored on fold tables verified to agree on **every** fold's train/test
membership (34,816 rows compared), so the comparison is paired and any
difference is attributable to the timepoint alone.

★ marks the best value in that column among the four fitted rows; **bold**
marks the pre-registered headline row, which is a separate thing from being
best. Higher is better everywhere except **Brier, MAE, RMSE and out-of-range**,
where the star follows the minimum.

The baseline row is what a trivial model scores: always predicting "success"
for classification, and the training mean for regression. **Sensitivity 1.000
and F1 0.619 come free from that**, which is why neither is quoted alone.

### D11 → D52 — classification

| model | scheme | ROC-AUC | PR-AUC | balanced acc | sensitivity | specificity | F1 | Brier |
|---|---|---|---|---|---|---|---|---|
| **logistic_l2** *(headline)* | plain | 0.902 ± 0.010 | 0.864 ± 0.010 | 0.837 ± 0.015 | 0.856 ± 0.034 | 0.817 ± 0.019 | 0.822 ± 0.017 | 0.129 ± 0.007 |
| **logistic_l2** *(headline)* | **donor_grouped** | 0.906 ± 0.008 | 0.874 ± 0.012 ★ | 0.834 ± 0.021 | 0.864 ± 0.031 | 0.804 ± 0.026 | 0.821 ± 0.022 | 0.125 ± 0.006 |
| logistic_l1 | plain | 0.908 ± 0.016 ★ | 0.873 ± 0.018 | 0.839 ± 0.025 ★ | 0.844 ± 0.034 | 0.833 ± 0.022 ★ | 0.824 ± 0.027 | 0.123 ± 0.012 ★ |
| logistic_l1 | **donor_grouped** | 0.905 ± 0.011 | 0.872 ± 0.012 | 0.837 ± 0.016 | 0.872 ± 0.017 ★ | 0.803 ± 0.025 | 0.825 ± 0.016 ★ | 0.126 ± 0.009 |
| *trivial baseline* | — | *0.500* | *0.449* | *0.500* | *1.000* | *0.000* | *0.619* | *0.247* |

### D30 → D52 — classification

| model | scheme | ROC-AUC | PR-AUC | balanced acc | sensitivity | specificity | F1 | Brier |
|---|---|---|---|---|---|---|---|---|
| **logistic_l2** *(headline)* | plain | 0.938 ± 0.015 | 0.910 ± 0.020 ★ | 0.884 ± 0.012 ★ | 0.885 ± 0.020 ★ | 0.883 ± 0.016 | 0.872 ± 0.013 ★ | 0.092 ± 0.007 |
| **logistic_l2** *(headline)* | **donor_grouped** | 0.945 ± 0.008 ★ | 0.906 ± 0.019 | 0.877 ± 0.013 | 0.880 ± 0.022 | 0.875 ± 0.009 | 0.865 ± 0.015 | 0.090 ± 0.006 ★ |
| logistic_l1 | plain | 0.934 ± 0.014 | 0.908 ± 0.022 | 0.881 ± 0.011 | 0.872 ± 0.020 | 0.891 ± 0.014 ★ | 0.869 ± 0.013 | 0.094 ± 0.008 |
| logistic_l1 | **donor_grouped** | 0.940 ± 0.009 | 0.899 ± 0.023 | 0.870 ± 0.013 | 0.861 ± 0.018 | 0.879 ± 0.015 | 0.856 ± 0.014 | 0.093 ± 0.005 |
| *trivial baseline* | — | *0.500* | *0.449* | *0.500* | *1.000* | *0.000* | *0.619* | *0.247* |

### D11 → D52 — regression

| model | scheme | R² | MAE | RMSE | out-of-range |
|---|---|---|---|---|---|
| **ridge** *(headline)* | plain | 0.515 ± 0.018 ★ | 0.100 ± 0.002 | 0.133 ± 0.002 ★ | 0.056 ± 0.010 |
| **ridge** *(headline)* | **donor_grouped** | 0.503 ± 0.015 | 0.099 ± 0.001 ★ | 0.134 ± 0.002 | 0.057 ± 0.008 |
| lasso | plain | 0.514 ± 0.015 | 0.100 ± 0.002 | 0.133 ± 0.002 | 0.050 ± 0.019 ★ |
| lasso | **donor_grouped** | 0.502 ± 0.023 | 0.100 ± 0.001 | 0.134 ± 0.003 | 0.051 ± 0.011 |
| *trivial baseline* | — | *0.000* | *0.158* | *0.190* | *0.000* |

### D30 → D52 — regression

| model | scheme | R² | MAE | RMSE | out-of-range |
|---|---|---|---|---|---|
| **ridge** *(headline)* | plain | 0.817 ± 0.014 | 0.059 ± 0.002 | 0.081 ± 0.003 | 0.037 ± 0.010 |
| **ridge** *(headline)* | **donor_grouped** | 0.802 ± 0.009 | 0.061 ± 0.001 | 0.085 ± 0.002 | 0.037 ± 0.006 |
| lasso | plain | 0.828 ± 0.005 ★ | 0.056 ± 0.001 ★ | 0.079 ± 0.001 ★ | 0.004 ± 0.005 ★ |
| lasso | **donor_grouped** | 0.811 ± 0.006 | 0.058 ± 0.001 | 0.083 ± 0.001 | 0.010 ± 0.009 |
| *trivial baseline* | — | *0.000* | *0.158* | *0.190* | *0.000* |

**LOCO and LODO agree with the headline**, so nothing hinges on the 5-fold
scheme: ROC-AUC 0.912 / 0.908 at D11 and 0.945 / 0.938 at D30; R² 0.501 / 0.505
and 0.820 / 0.804. Each is a single repeat by construction, so their SD is
undefined — see [`modeling/results/README.md`](modeling/results/README.md).

**Donor grouping costs almost nothing.** R² 0.515 → 0.503 at D11 and 0.817 →
0.802 at D30 against plain 5-fold; classification is flat or slightly better.
Donor leakage is real here but small.

**D11's DA-only numbers are lower than its published DA+Sert ones** (R² 0.503 vs
0.653, AUC 0.906 vs 0.947). That is **not a regression** — it is a harder
target. `DA+Sert` is bimodal with a genuine gap at 0.2, so the threshold
separates two clumps; `DA/all` is unimodal with no gap, so the same threshold
cuts through a dense region. Different question, different difficulty. The two
families are not comparable and should never be differenced.

## 3. ⚠️ The D30 model does not beat a single raw column

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

## 4. What D30 can and cannot claim

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

## 5. Caveats

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
- **No external validation.** 20 donors, one protocol, one study.

## 6. Which features matter

One-SE config selection, paired per-fold LOCO and permutation deltas, full-fit
coefficients/SHAP from a PCA fit once on all 136 lines, pool η² on the same
features. Tables in `modeling/results/`.

### Recapitulates the published analysis

**`phat_NB`, the neuroblast proportion, is still D11's standout composition feature**, with the same sign
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

A **20×** gap. This is the feature-level counterpart of §3's benchmark result — `phat_DA` alone nearly matches the full model —
from the other direction. At D11 the PCs do carry weight (top LOCO Δ
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
other types are held fixed — consistent with §1: Sert is a different lineage,
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

## 7. Literature assessment: recapitulated vs. new

Same framing as the appendix's §5, applied to these results. The baseline
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
| 4 | `phat_NB` (neuroblast proportion) remains D11's strongest *and technically cleanest* composition feature (pool η² 0.035) | **Recapitulates phase 1, and the authors' indirect observation** | They noted a poor-differentiation cluster correlating with D11 neuroblast proportion. Phase 1 made it the top direct predictor; it survives the outcome change with the same sign |
| 5 | **At D30 the transcriptome adds nothing beyond cell-type composition** — 3 of 4 models select k = 0 PCs; permutation Δ 0.360 for proportions vs 0.018 for the best PC | **New** | A 20× gap. Not addressed by the authors, who did not build D30 predictors |

### What this does *not* establish

- **Not a better predictor than phase 1.** The DA-only numbers are lower
  (AUC 0.906 vs 0.947 at D11) because the target is harder, not because
  anything regressed. The two are not comparable at all.
- **Not a validated D30 biomarker.** D30's strength is substantially the
  outcome already existing in the predictors.
- **Not externally validated.** 20 donors, one protocol, one study, no held-out
  batch. Everything here is internal cross-validation.
- **Not a causal claim** anywhere. All associations.

---

# Appendix — phase 1: the authors' DA+Sert outcome

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
