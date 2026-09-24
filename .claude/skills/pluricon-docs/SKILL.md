---
name: pluricon-docs
description: Write or revise this project's three reader-facing documents — README.md, FINDINGS.md and modeling/METHODS.md. Covers what each document is for and in what order, which vocabulary is retired, which section titles source code references by name, where every number must be traced from, and the verification sweep to run before calling a documentation change done. Use for any edit to those three files, or to modeling/README.md.
---

# Pluricon documents

**Craft rules live in the `scientific-writeup` skill — load that first.** Prose
form, measurement-then-interpretation, explaining statistics, evidence
discipline, restructuring safely, figures. This file only covers what is
specific to this repository, and does not restate them.

## What each document is for

| document | job | order settled on |
|---|---|---|
| `README.md` | front door | summary → the day-11 punchline → the science (data, features, CV, metrics, importance, results) → repository and setup → outcome-definition footnote → caveats → **"How this was built with AI" last** |
| `FINDINGS.md` | the results narrative | the four questions stated up front, then answered in order: outcome definition → model → donor → D11 → D30 → literature → conclusions → references → three short appendices |
| `modeling/METHODS.md` | how and why | setup → model design → metrics → cross-validation → feature importance → reproducibility. **No results numbers.** |
| `modeling/README.md` | provenance | audit corrections, decision records, label-variant machinery. **No performance numbers at all** |

- FINDINGS conclusions are **numbered questions with crisp sub-bullets**, one
  claim per bullet.
- Appendices stay short. The full performance grids are the reason Appendix A
  exists; it needs no commentary beyond the ★ note.
- METHODS is the one place where explaining *why* a statistical choice was made
  is the content, not prose debt.

## Retired vocabulary

**"Headline" and "shipped" do not appear in prose.** The configuration table is
captioned *Selected configuration per task*.

**Exception — these are live code symbols and stay verbatim:**
`select_headline()`, `HEADLINE`, `HEADLINE_MODEL`. They are imported by
`run_regression.py`, `run_classification.py` and `run_variant.py`, and asserted
on in `tests/test_run_experiment.py` and `tests/test_gates_2_3_blind.py`.
Renaming the concept in prose without renaming the code desynchronises both.

## Where numbers come from

Every figure traces to `modeling/results/`:

| claim | file |
|---|---|
| ROC-AUC, PR-AUC, Brier etc. | `results_classification_{D11,D30}_qualonly_da_untreated.csv` |
| R², MAE, RMSE, out-of-range | `results_regression_{D11,D30}_qualonly_da_untreated.csv` |
| SHAP, ρ, LOCO Δ, permutation Δ, pool η² | `feature_importance_table_classification_logistic_l2_{D11,D30}_qualonly_da_untreated.csv` |
| `phat_DA`-alone vs full model | `incremental_value_D30_da_untreated.csv` |
| fit-free benchmarks + bootstrap CIs | `d30_single_feature_benchmarks_da_untreated.csv` |
| PC gene loadings, gene-set ranks | `pc1_loadings_{D11,D30}_da_untreated.csv`, `pc2_loadings_D11_da_untreated.csv`, `modeling/docs/PCA_interpretation.md` |
| cell-type markers, the 48/31/20 split | `d30_celltype_markers.csv`, `modeling/docs/D30_celltype_interpretation.md` |
| donor and pool η² on the outcome | **not on disk** — `modeling/donor_batch_variance.py` prints and saves nothing; run it from the repo root |

Filter to the current analysis: `tuning == "nested"`, `scheme ==
"donor_grouped"`, `logistic_l2` / `ridge`. The `_da_untreated` suffix is the
current outcome; anything without it is phase 1.

> **Nothing in `modeling/README.md` is a current number.** Every figure it ever
> carried was phase 1 — 138 lines, `(DA+Sert)/all`, D11 only. Do not quote it.
> Two benchmark tables also exist for different cohorts; do not mix
> `D30_celltype_interpretation.md`'s 0.990 / 0.956 / 0.874 (138-line DA+Sert)
> with the current 0.960 / 0.886 / 0.883 / 0.727.

## Do not rename these section titles

Source code and tests reference them as strings:

- **"Results distillation"** — `run_experiment.py:7,42`
- **"Second correctness review"** — `run_experiment.py:124`,
  `tests/test_run_experiment.py:158`
- **"Do NOT compare the two importance tables row by row"** —
  `modeling/docs/PCA_interpretation.md:81`

## Cross-document contract

Every claim and number in `README.md` appears in `FINDINGS.md` unchanged.
FINDINGS may add detail; it may not introduce a competing figure. When the two
disagree, one of them is a defect.

## Verification sweep

Run before calling a documentation change done:

1. **Numbers traced** to the files above; and diff the number multisets against
   the previous commit — expect nothing added, and only intended removals.
2. **README ↔ FINDINGS** shared claims agree.
3. `grep -in "shipped\|headline"` — clean, except the code symbols above.
4. **Tables well-formed** — header, separator and body agree on column count;
   escaped `\|` inside cells does not count as a delimiter.
5. **Anchors and links resolve** — every `](#...)` matches a heading, every
   relative path exists. Renumbering sections breaks `README.md`'s two
   `FINDINGS.md#...` links.
6. **Figures** regenerated if their inputs moved, and *looked at*.
7. `uv run pytest modeling/ -m "not slow"` — 120 passing. Run the full 124
   before treating results as final.
