# What the D30 cell types are, and what that means for a D30 → D52 model

Companion to [`PCA_interpretation.md`](PCA_interpretation.md), which does the
same job for the PCs. This one is about the **annotated cell types** — because
at D30, unlike D11, the annotation set already contains the outcome.

Regenerate every number here with:

```bash
uv run python -m modeling.d30_celltype_markers          # marker table
uv run python -m modeling.d30_single_feature_benchmarks # benchmark table
```

---

## 1. The problem in one line

`diff_efficiency`, the outcome, is **the fraction of a line's D52 cells that
are DA or Sert**. D11's three annotated types are `FPP`, `P_FPP`, `NB` — all
progenitor/early, none of them DA or Sert. **D30's seven types include DA and
Sert.**

So a D30 model is not automatically doing what the D11 model does. Part of its
input *is* the outcome, measured 22 days earlier.

This is **not leakage**. D30 and D52 are separate cells from separate harvests,
D30 precedes D52, and no D52 information enters the D30 features. It is
**construct overlap** — predictor and outcome quantify the same lineage state
at two times. The distinction matters for how results are framed, not for
whether the CV is valid.

## 2. Which D30 types are upstream of the outcome?

Mean log-normalised expression per annotated type, all 250,923 D30 cells,
34/34 canonical markers found. Panels were fixed a priori, not chosen after
seeing which separated the clusters.

| type | n (%) | SOX2 | HES1 | VIM | MKI67 | TOP2A | STMN2 | identity markers |
|---|---|---|---|---|---|---|---|---|
| **FPP** | 58,109 (23%) | **2.16** | **1.82** | **34.7** | 0.02 | 0.06 | 0.43 | LMX1A 1.00 |
| **P_FPP** | 18,988 (8%) | **2.14** | 1.26 | **21.3** | **2.10** | **4.36** | 0.41 | LMX1A 0.88 |
| DA | 69,007 (27%) | 0.37 | 0.07 | 3.7 | 0.00 | 0.01 | **27.9** | TH 1.78, NR4A2 1.60, DDC 1.74 |
| Sert | 53,774 (21%) | 1.53 | 0.32 | 1.7 | 0.01 | 0.02 | 2.28 | **TPH1 14.79**, GATA2 0.49 |
| Epen1 | 28,145 (11%) | 1.76 | **2.39** | **30.0** | 0.01 | 0.19 | 0.29 | **TTR 134.65**, PIFO 9.35, RSPH1 8.92, FOXJ1 3.64 |
| U_Neur1 | 21,268 (8%) | 0.10 | 0.05 | 5.1 | 0.01 | 0.05 | **40.3** | none |
| U_Neur2 | 1,632 (0.6%) | 0.16 | 0.04 | 2.4 | 0.01 | 0.02 | **37.3** | none |

### The three-way split

**Already arrived — `DA`, `Sert` (48% of cells).** This *is* the outcome.
Post-mitotic, STMN2/DCX/MAP2 for DA. Note `SLC6A3` (DAT) is **0.00** in D30 DA
cells and `SLC6A4` is 0.00–0.01 in Sert cells: these are **immature**
DA/serotonergic neurons. The D52 label counts the same annotation, not mature
transporter-positive neurons.

**Still undecided — `FPP`, `P_FPP` (31%).** SOX2⁺/HES1⁺/VIM-high, essentially
no pan-neuronal signal (STMN2 ≈ 0.4 vs ~28–40 in the neuronal types). `P_FPP`
is the same population still cycling (MKI67 2.10, TOP2A 4.36, CCNB1 3.68).
**These cells have not committed.** What they become between D30 and D52 is
genuine prediction, in exactly the sense the whole D11 model is.

**Off-target, terminal — `Epen1`, `U_Neur1`, `U_Neur2` (20%).**
- `Epen1` is **not** a Sert relative. It is ciliated ependymal / choroid-plexus-like:
  TTR 134.65 (enormous), plus the motile-cilia module FOXJ1 / PIFO / RSPH1 /
  SPAG6. It *retains* SOX2/HES1/VIM, so it is neuroepithelium-derived — same
  broad origin as FPP — but committed down a divergent branch. Its ρ = −0.770
  against the outcome fits: Epen1 is where the culture goes when it fails.
- `U_Neur1`/`U_Neur2` are dead ends: STMN2 ≈ 37–40, zero proliferation, no
  progenitor markers, no DA or Sert markers. Terminally differentiated
  off-target neurons. They cannot become anything else.

### Caveat that limits all of the above

Marker expression establishes cell **identity**, not **fate**. That FPP/P_FPP
carry progenitor and cell-cycle programs while U_Neur1/2 carry neither is
strong evidence, but proving Epen1 cells cannot still convert would need
lineage tracing, which this cross-sectional dataset cannot do.

## 3. Fit-free benchmarks

Single columns, no model, no fitted parameter. 138 lines, 20 donors,
donor-level bootstrap (donors, not lines — 136 of 138 lines share a donor, so a
line-level bootstrap would be too narrow), 2,000 resamples, seed 0.

| benchmark | ρ vs efficiency | ROC-AUC (oriented) | 95% CI | corr. with DA+Sert |
|---|---|---|---|---|
| `DA+Sert` — already arrived | +0.896 | **0.990** | 0.978–1.000 | — |
| `FPP+P_FPP` — still undecided | −0.835 | 0.956 | 0.929–0.982 | **−0.943** |
| **`FPP+P_FPP` share of non-target — progenitor balance** | **+0.661** | **0.874** | 0.816–0.928 | **+0.665** |

Directions are **declared a priori in the source**, not flipped after seeing
the result — ROC-AUC is direction-sensitive, and choosing the sign from the
data is a data-driven aggregation choice that would have to happen inside a
training fold.

These are **in-sample over the 138 study lines**. They are reference lines, not
validated model performance. There is nothing to cross-validate — with no
fitted parameter, a line's prediction is its own feature value whether or not
it was "in training", so a CV wrapper would produce a CV-shaped number that is
really in-sample.

### Reading these three rows correctly

**Rows 1 and 2 are not two findings.** They correlate **−0.943** — the same
maturation axis with the sign flipped ("how far along is this culture"). The
raw `FPP+P_FPP` fraction scoring 0.956 is mostly a restatement of `DA+Sert`
scoring 0.990.

**Row 3 is the one that answers the real question.** Dividing by the undecided
mass — *of the cells that have not become DA/Sert, what share are still
competent progenitors rather than committed to Epen1/U_Neur?* — removes the
overall maturation magnitude. Its correlation with `DA+Sert` drops from −0.943
to +0.665, and it still reaches **ρ = +0.661, AUC 0.874**.

For scale, D11's standout single feature `phat_NB` has ρ = −0.73. **The D30
progenitor balance carries genuine, D11-comparable predictive signal from cells
that have not yet decided.** It is not merely the outcome read early — but it is
smaller than the 0.990 headline, and the maturation magnitude has to be removed
before it is visible.

## 4. Consequences for the modelling

1. **`DA+Sert` alone is the ceiling, and belongs in the results table, not a
   footnote.** At AUC 0.990 there is almost no headroom; a fitted D30 model
   scoring 0.99 has added nothing. The question is whether the PCs add
   information *conditional on* that fraction, which needs Brier / log-loss /
   continuous R², not AUC.
2. **Two pre-registered CV variants** (headline fixed before any results were
   seen): `D30-full` (all 6 free proportions + PC1–10, the direct analogue of
   the D11 model) is the headline; `D30-progenitor-balance` is the sensitivity
   analysis.
3. **The sensitivity table renormalises rather than deletes.** Simply dropping
   `phat_DA`/`phat_Sert` leaves five proportions summing to `1 − DA − Sert`, so
   their common scale still encodes exactly the removed quantity and a model can
   recover it from any two of them. See
   [`../make_d30_variant_tables.py`](../make_d30_variant_tables.py).
4. **The PCs are not independent of composition** — expression reflects which
   cells are present. The progenitor-balance variant is a sensitivity analysis,
   not a clean decomposition.
5. **`phat_P_FPP` stays the reference category** at both timepoints. It is
   already D11's reference (so D11 results are bit-identical) and it is the
   weakest-correlated D30 type (ρ = −0.133), so it privileges nothing.

## 5. Unrelated D30 trap, recorded here because it was found alongside

The depth-outlier pool excluded from HVG/PCA fitting **is not the same pool at
D30 as at D11**, and the rule had been a single global constant:

| | D11 | D30 |
|---|---|---|
| pool11 | **1,740** median UMI/cell ← outlier | 13,058 — normal, in fact the deepest |
| pool5 | 15,407 | **1,048** ← outlier |

Reusing D11's `{pool11}` at D30 would have been wrong twice: excluding the
deepest pool and keeping the shallowest. The rule is now per timepoint and
stated — *median UMI/cell below ⅓ of the median-of-pool-medians* — which
reproduces D11's frozen `{pool11}` (0.137, next-lowest 0.747) and gives
`{pool5}` at D30 (0.123, next-lowest 0.423). Wide gaps on both sides, so ⅓ is a
separator, not a tuned knob. See `modeling/features.py:DEPTH_OUTLIER_POOLS`.

Open item flagged by the independent review and **not yet addressed**: pool5's
~10× depth deficit may have biased D30 *cell-type annotation* itself, and
excluding it from the PCA fit does nothing to correct biased labels. It holds
18 of 159 qualifying (line, pool) combos.
