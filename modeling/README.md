# `modeling/` — provenance and infrastructure

**Methods live in [`METHODS.md`](METHODS.md)** — features, models,
hyperparameters, cross-validation, metrics and feature importance. Results live
in [`../FINDINGS.md`](../FINDINGS.md); output files in
[`results/README.md`](results/README.md).

This file is the record of how the pipeline got to its current state: what audit
review found and fixed, which decisions were taken and why, and how the
label-variant machinery works.

> **No performance numbers here.** Every results table this file used to carry
> was from phase 1 — 138 lines, the `(DA+Sert)/all` outcome, D11 only — and is
> superseded by the 136-line `DA/all` untreated analysis. The current grids are
> in [`../FINDINGS.md`](../FINDINGS.md) Appendix A, the phase-1 bridge in its
> Appendix C, and the originals in git history at `9a23404`.

## Layout

```
modeling/
  *.py           the package: folds, features, models, harness,
                 run_* entry points, feature_importance
  tests/         pytest suite (30 fast + 1 slow integration test)
  fold_data/     per-fold feature tables and fold assignments
                 (`_full` = baseline basis, `_qualonly` = restricted-fit variant)
  results/       CV grids, nested per-fold selections, feature-importance tables
  docs/          PCA_interpretation.md, feature_importance_plan.md,
                 pool_correction_investigation.md, plus two standalone
                 HTML explainers (open in a browser):
                 eta_squared_and_r2.html  -- what the pool eta^2 column means
                 roc_auc_and_pr_auc.html  -- ROC/PR, and the 100%-precision point
  plots/         diagnostic figures
  archive/       superseded first-pass artifacts, kept for provenance only --
                 nothing reads these
```

Writers create their own subdirectory, so a fresh clone or a redirected
`out_dir` works without manual `mkdir`.
## Do NOT compare the two importance tables row by row

The components REORDER between bases, so "PC2" does not mean the same
thing in each. Correlating per-line coordinates across the 138 lines:

| | closest match | \|r\| | |
|---|---|---|---|
| baseline PC1 | qualonly PC1 | 0.995 | same component |
| baseline PC2 | qualonly **PC3** | 0.985 | same component |
| qualonly PC2 | baseline PC2 | 0.564 | **a different component** |

Read positionally, the qualonly table looks like the fit restriction
cleaned the pool signal out of the strongest PC (eta2 by pool
0.764 -> 0.114). It did not. The pool-associated axis MOVED from position
2 to position 3 and became slightly *more* pool-associated
(eta2 0.764 -> 0.832). What rose into position 2 is a genuinely different
component that happens to be largely pool-independent (eta2 0.114) and
strongly outcome-associated (Spearman rho -0.745 vs binary `success`,
-0.809 vs continuous `diff_efficiency`).

Explained variance is near-identical across bases
(0.0539/0.0265/0.0259/0.0203 -> 0.0539/0.0268/0.0254/0.0203), i.e. PC2 and
PC3 are near-tied, which is exactly the regime where ordering is unstable.
**Align components by correlation before interpreting anything
biologically.**
## Corrections from the audit

After the analysis was complete, an independent Codex correctness audit
(implementation-vs-intent, not methodology) checked every stated design
intention across `folds.py`, `features.py`, `harness.py`,
`feature_importance.py`, and `run_experiment.py`. Most passed. Three real
defects were found and fixed; all three are recorded here rather than
quietly patched, since two of them changed reported numbers.

1. **SHAP mixed unit systems (changed a conclusion).** Coefficients come
   from a model fit on standardized features (units: per 1 SD), but they
   were multiplied by *raw* `(x - mean)` deviations. Every contribution
   was therefore scaled by that feature's raw SD, systematically
   inflating wide-scale features (PCs, raw SD ~1-9) and deflating narrow
   ones (proportions, raw SD ~0.05) -- exactly the comparison the column
   exists to make. Fixed to use standardized deviations. **This reversed
   the classification SHAP ranking**: `phat_NB` went 0.068 -> 1.190 and
   PC2 went 2.438 -> 0.514, so `phat_NB` now leads rather than trailing
   PC2 by ~36x. It strengthens the `phat_NB` finding and weakens the PC2
   one; the earlier README claim that PC2 dominated SHAP was wrong.
2. **The "one-SE rule" used one SD, not one SE.** ~3x too wide at 10
   repeats, so it selected simpler models than a real one-SE rule would.
   Fixed to SD/sqrt(n_repeats). Verified the tolerance genuinely narrowed
   (lasso: 5 configs within tolerance -> 1), but the *selections* happened
   not to change (k=5 regression, k=2 classification both still qualify),
   so no downstream numbers moved.
3. **`select_headline` didn't filter by model, and the reporting
   compounded it.** The pre-registration named ridge / logistic_l2, but
   the function returned all models and the better-scoring non-registered
   models (lasso R2=0.668, logistic_l1 AUC=0.951) were being quoted as
   *the* result -- the exact cherry-picking the pre-registration exists
   to prevent. Fixed to filter on the pre-registered model; the reported
   numbers became ridge R2=0.653 and logistic_l2 AUC=0.943, with the L1
   variants labelled secondary. **The lesson outlives the
   pre-registration**: whatever declares the model -- a pre-registration
   then, the one-SE rule now -- the reporting function must filter on it,
   or the best-scoring row wins by default and the declaration is
   decorative.

Feature-importance tables are now produced for **all four** models
(`results/feature_importance_table_{task}_{model}.csv`) rather than just the L1
pair, for the same reason. The L1 tables remain the more informative ones
(only L1 produces sparsity, so `regularized_to_zero` /
`selection_frequency` are meaningful), but both are reported.
## Second correctness review (variant plumbing + test pollution)

A second independent Codex review, run specifically on the qualifying-only
work, found two more defects. Recorded here because one of them means a
previous entry in this README overstated what had been fixed.

All of these share one root cause: **adding a `--suffix` variant required
threading "which variant am I?" through several layers, and it was
threaded through some but not others.** Every resulting failure produced
plausible, non-crashing, wrong output.

1. **Mixed-basis importance tables.** `--suffix` redirected the
   fold-features input and the output filenames, but
   `load_full_fit_features()` always read 005's baseline global PCA. So a
   single `_qualonly` table had fold-level columns (LOCO, permutation,
   selection frequency) from one basis and full-fit columns (coefficient,
   SHAP, univariate, technical covariate) from another. The tell was that
   the variant's `univariate_spearman` column was byte-identical to
   baseline. Fixed by parameterizing `005` (`--restrict-to-qualifying`,
   `--suffix`) and threading `pca_csv` through every full-fit consumer.
2. **That fix was incomplete, and was described here as complete.**
   `compute_shap_like` was still calling `load_full_fit_features()` with
   no argument, so `shap_mean_abs` stayed on the baseline basis. The edit
   had been applied with a string replacement whose pattern did not match,
   and string replacement fails silently. Re-applied with an assertion
   that the edit landed plus a check that no no-arg call survives.
3. **The test suite was writing into the repo.** `run_all` saves
   `results/nested_selections_{task}.csv` to the package directory, and
   `test_run_experiment.py`'s fixture patched its inputs but not that
   path. Every `pytest` run overwrote the committed
   `nested_selections_regression.csv` with 32 rows of synthetic fixture
   data (`k=1, alpha=1.0`, and `repeat=1` for LOCO/LODO, which is
   structurally impossible) — and that corrupted file had been committed.
   `run_all` now takes an injectable `out_dir`, the fixture redirects it,
   and a regression test asserts `run_all` never touches committed files.
4. **`_infer_n_repeats` read a hardcoded table.** It always loaded
   `fold_data/fold_features_D11_full.csv` (falling back silently to
   `fold_features_D11.csv`), not the table the run was using, so the
   one-SE tolerance could come from an unrelated experiment. Harmless in
   fact here — both runs share fold assignments, so `n_repeats` was 10
   either way — but nothing enforced that. Now required and passed
   explicitly; a missing file is a loud error.
5. **`logistic_l1` was not reproducible.** `saga` is a stochastic solver
   and no `random_state` was set, so coefficients drifted run to run
   (measured max|diff| ~5.7e-4), producing spurious 4th-decimal diffs on
   re-runs. `models.SEED` now pins both logistic factories.

A guard now prevents defect 1 from recurring silently: with `--suffix`
set and no explicit basis, the matching suffixed basis is REQUIRED, and
its absence is a hard error naming the command that generates it.
## D30 readiness — NOT ready (read before starting D30 -> D52)

> **Closed.** D30 → D52 is complete; every blocker below was addressed, and
> the audit gates are recorded in [`docs/audit_ledger.md`](docs/audit_ledger.md).
> Kept because it records what the risks looked like before the work, not
> because any of it is still outstanding.

The earlier claim in "Architecture" that feature extraction is simply
"timepoint-parameterized" is too optimistic. What is genuinely reusable:
`folds.py` (fold construction is timepoint-independent — same 138 lines,
same D52 labels), `models.py`, most of `harness.py`, and
`run_experiment.run_all` given a correct feature table.

What must change first, ranked by risk. The loud failures are fine; the
SILENT ones are the danger, and they are marked.

1. **Cell-type schema differs and is hardcoded.** D11 has 3 types
   (`FPP, NB, P_FPP`); **D30 has 7** (`DA, Epen1, FPP, P_FPP, Sert,
   U_Neur1, U_Neur2`). `harness.ALL_PROPORTION_COLS` /
   `MODEL_PROPORTION_COLS` hardcode the D11 three. Fails loudly on a
   missing column, but the compositional reference category must be
   chosen and documented for 7 types, not patched ad hoc.
2. **SILENT: extraction is hardcoded to D11.** `run_feature_extraction.py`
   passes `"D11"` to `load_cell_metadata` and
   `compute_pca_features_for_fold` and names outputs
   `fold_features_D11{suffix}.csv`. Changing only the output filename
   would produce a plausible file named D30 containing D11 features.
3. **SILENT: resume has no provenance.** The incremental-write resume keys
   only on `(scheme, repeat, fold)`. It does not record timepoint,
   fitting population, or PCA settings, so reusing a destination after
   changing any of those silently keeps stale rows. Write a provenance
   header/sidecar and refuse to resume across a settings change.
4. **SILENT: full-fit inputs are D11-specific.** `005` is hardwired to
   `DAY11_FILE` and `d11_*` output names; `feature_importance.py`
   hardcodes `D11_PCA_CSV`, `D11_PROPORTIONS_CSV` and derives the variant
   basis as `d11_pca_coords_per_line{suffix}.csv`. Passing D30 fold
   features *unsuffixed* would combine D30 fold-level values with D11
   global PCA and D11 proportions. (A suffixed call would fail loudly on
   the missing `d11_...{suffix}` file.) Make these explicit parameters.
5. **`run_variant.py` input is hardcoded** to `fold_features_D11{suffix}`.
6. **Decide whether excluding pool11 is still justified at D30.** It is
   currently inherited from a D11 sequencing-depth observation without
   re-examination.

Note the qualifying `(cell_line, pool)` combos ARE the same 159 across
138 lines at both timepoints (qualification already requires >=10 cells at
D11 and D30 and D52), but the per-combo cell COUNTS are not: zero of 159
combos has identical D11 and D30 counts (D11 range 25-14,640; D30 range
13-13,112; correlation 0.748). Sampling precision therefore differs
substantially, which is a further reason D11 proportions must never be
reused as D30 features.
## Pool-correction dropped

Investigated (with an independent Codex review) why pool-correction cut
regression R2 by ~0.10-0.12 but barely touched classification AUC. No
single settled explanation emerged -- see `docs/pool_correction_investigation.md`
for the concise summary (my reasoning, Codex's oracle-test refutation of
my leading theory, its alternative explanation, and two real
implementation issues found along the way). Decision: drop pool-correction
from the active pipeline rather than resolve this now -- it added
significant complexity for an unclear net benefit at this prototype
stage. `features.py`/`harness.py` still support it if revisited.
## Results distillation

*Historical. This documents the pre-registration regime phase 1 used. The
current analysis selects its configuration by flat CV + one-SE instead — see
[`METHODS.md`](METHODS.md) §4. The constants below still exist in the code, so
they are documented here rather than deleted.*

Phase 1 committed one configuration per task **before** looking at scores, to
stop the best of 16 numbers per task being quoted after the fact. The constants
live in `run_experiment.py` and are enforced by `select_headline()`:

```python
HEADLINE       = {"scheme": "donor_grouped", "tuning": "nested", "pool_correction": False}
HEADLINE_MODEL = {"regression": "ridge", "classification": "logistic_l2"}
```

Each was chosen on design grounds, not on scores: `donor_grouped` because it is
the strictest realistic test; `nested` because flat CV tunes and reports on the
same rows; no pool-correction per the investigation above; ridge / logistic_l2
as the L2 default with L1 as the variant.

The L1 variants scored slightly higher in phase 1. A paired per-repeat
comparison at matched one-SE configurations showed lasso's edge was real but
tiny — +0.0033 R², winning 9 of 10 repeats, t-test p = 0.0006 — and most of the
apparent gap was nested CV's per-fold reselection noise rather than a model
difference.

**The lesson outlives the mechanism.** Whatever declares the configuration — a
pre-registration then, the one-SE rule now — the reporting function has to
filter on it, or the declaration is decorative. That failure actually occurred
once; see "Corrections from the audit" #3.

## Architecture (Open-Closed + reusable for D30→D52)
- **Model registry, not conditionals**: a list of `ModelSpec(name, estimator_factory, param_grid)` entries (Ridge/Lasso/L2-logistic/L1-logistic to start). The CV harness only ever calls the sklearn estimator interface (`.fit`/`.predict`/`.predict_proba`) generically — adding a model later (e.g. random forest) = one new registry entry, zero changes to harness code.
- **Feature extraction (timepoint-parameterized) vs. CV harness (timepoint-agnostic)** — *aspirational; see "D30 readiness" above for what actually still hardcodes D11*: feature extraction takes a timepoint (D11 now, D30 later — same 138 qualifying lines, since the qualifying list already requires ≥10 cells at D11 *and* D30 *and* D52) and a fold's train/held-out line split, returns feature matrices. The CV/tuning/metrics harness takes `X_train, y_train, X_test, y_test` + the model registry and is completely agnostic to which timepoint produced the features. D30→D52 later = swap the feature-extraction call, harness untouched.
## Files
- `features.py` — feature-extraction, parameterized by timepoint (D11 now; reusable for D30).
- `folds.py` — fold construction (plain/donor-grouped repeated stratified K-fold, LOCO, LODO) — timepoint-agnostic.
- `models.py` — model registry (`ModelSpec` list).
- `harness.py` — CV/tuning/metrics harness using the registry + fold module.
- `run_regression.py` / `run_classification.py` — orchestration entry points.
- Summary table + comparison plots (selected config; plain vs. donor gap; flat vs. nested; corrected vs. uncorrected).
## Verify
- Per-fold fitting cell counts stay large; fold assignments reproducible; classification metrics only where both classes present; repeated-CV reported per-repeat not pooled; regression predictions checked against [0,1].
# Label variants: running an alternative outcome

Added when a second outcome — `DA / all D52 cells`, untreated cells only — was
run alongside the published `(DA + Sert) / all D52 cells`. Results:
`modeling/results/*_da_untreated.csv`, documented file-by-file in
[`results/README.md`](results/README.md). Findings: [`../FINDINGS.md`](../FINDINGS.md)
"Addendum".

## Why a registry rather than threaded paths

`folds.py` was the single place that resolved the cohort and label, and the only
link in the chain without a variant parameter. The failure mode it guards
against is **not** "someone forgot to pass an argument" — it is *one of four
coupled values coming from a different variant than the other three*. Cohort
file, label file, threshold and expected line count are only meaningful
together; a label computed on one cohort, thresholded at another's cut-off and
checked against a third's line count is a silently wrong analysis that raises
nothing.

```python
@dataclass(frozen=True)
class LabelVariant:
    suffix: str          # applied by CODE when naming artifacts, never typed
    cohort_csv: Path
    label_csv: Path
    threshold: float
    n_lines: int
    description: str

VARIANTS = {"published": ..., "da_untreated": ...}   # dict key IS the name
```

Bundling them makes the mismatched combination unconstructible. Three
consequences worth keeping:

**No internal default.** `load_lines_with_label(variant)` requires the argument.
A default would let a call site that was never updated fall back to the
published cohort while the rest of the run used another — and the existing tests
would still pass, because they only ever exercise the default. "Existing tests
still pass" would then mean "the new branch was never exercised", not "the
refactor is safe". CLI entry points default to `"published"`, where the choice
is visible. Making it required immediately surfaced all 7 call sites.

**Duplicated cohort constants deleted.** `run_feature_extraction.py` held a
*second* hardcoded copy of the cohort path, independent of `folds.py`'s, and it
selected which combos features were computed on. The untreated cohort is 157
combos against 159, so the two could have disagreed silently. Removing it makes
a missed call site fail at import rather than in a results table.

**The variant suffix is appended by code, never typed.** An operator cannot aim
a new-label artifact at a published path, because they never get the chance to
write the suffix.

## Provenance: `label_variant` and `fit_population`

Every fold-feature row, results row and nested-selection row records which
variant produced it, and resuming across a differing value raises.

**An absent value is treated as UNKNOWN, not compatible.** This is the stricter
of the two rules in the file, and deliberately so: a label change is invisible
in the fold-features *filename*, which is all the earlier guard recorded. A
table predating the column could have been built under either variant, so
resuming onto it can produce exactly the mixed basis the check exists to
prevent. Consequence: the committed published tables cannot be resumed. They are
complete and never need to be; `resume=False` is the documented escape.

## Adding a variant

1. Build the cohort if it differs — `003_qualifying_cell_lines.py --untreated-only --suffix …`
2. Build the label — `009_d52_outcome_label.py --da-only --cohort-suffix … --suffix …`.
   `--cohort-suffix` is separate from `--suffix` because the cohort and the label
   are **different axes**: the cohort depends on cell counts, not on which cell
   types define the outcome, so `_da_untreated` labels sit on the `_untreated`
   cohort. No single suffix names both.
3. Register it in `folds.VARIANTS`. The parametrized test in `tests/test_folds.py`
   then checks the new bundle is self-consistent automatically.
4. Extract per timepoint — `run_feature_extraction --label-variant … --timepoint …`
5. Verify — `python -m modeling.verify_fold_features D11 --paired-with D30`
6. Score — `run_variant --label-variant … --timepoint …`

## `verify_fold_features.py`

A hard gate between an extraction finishing and anything running on its output.
Extraction takes ~3 hours and its failure mode does not raise: it emits a
plausible table that is short a fold, missing a line, or carrying stale
provenance. Expectations are **derived from the variant's own cohort** rather
than hardcoded — the LOCO fold count comes from the cohort's line count, and the
excluded depth-outlier pool is read from `features.DEPTH_OUTLIER_POOLS` rather
than restated, since a second copy could drift from the real rule.

`--paired-with` compares two timepoints' **feature tables** on per-fold
train/test membership. That is stronger than comparing the fold-assignment file,
and much stronger than trusting the message printed during extraction: a log
line is not a verification artifact, and it would not catch features derived
from the right assignments landing in the wrong folds.

## Resume integrity

A fold counts as complete only on **exact key-set equality** against the
variant's cohort — one row per expected `(cell_line, pool)`, no extras, no nulls
or non-finite values. A row *count* is not the contract: a fold missing one
expected combo and carrying one unexpected one has the right number of distinct
keys and would pass a cardinality check while describing the wrong population.

Truncation is judged in one place. `_read_resumable` handles only the case that
makes a file unreadable (a final row with too *many* fields, which raises); a
row with too *few* fields is silently NaN-padded by the parser and is
deliberately **kept**, so that `complete_folds` rejects the whole fold. Dropping
null-bearing rows in the reader would also discard a genuinely missing value and
hide whatever produced it.

## What deliberately does NOT support variants

`feature_importance.py` raises `NotImplementedError` for anything other than
`published`, **before reading any data**. Its full-fit columns (coefficient,
SHAP, univariate, technical covariate) come from a global PCA basis fit under
the published cohort, while its fold-level columns (LOCO, permutation, selection
frequency) would come from the new variant — byte-for-byte the mixed-basis
defect already recorded in this file. Correct support requires re-running
`005_d11_pca_features.py` against the other cohort first. It is not partially
threaded: that would imply support it does not have.

`010_d52_label_technical_covariates.py` independently re-derives the label
instead of reading the CSV, so its outputs still describe the published outcome
only. Worth revisiting — dropping rotenone cells changes the per-combo depth
distribution, so "is the new label pool/depth-confounded?" is genuinely open.
