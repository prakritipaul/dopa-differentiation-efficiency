"""
Shared orchestration: run flat + nested CV across both grouping schemes
(plain, donor_grouped) and both pool-correction variants, for a given
task. Produces the full robustness grid; run_regression.py and
run_classification.py each pull out the pre-registered headline result
(donor_grouped + nested + pool-corrected + Ridge/L2-logistic, per
modeling/README.md "Results distillation") on top of this.

Task-agnostic by construction -- the CV harness (harness.py) takes a
`task` string and never special-cases regression vs. classification
beyond picking the right label column and metric functions, so this
orchestration layer is shared rather than duplicated per task. Also
timepoint-agnostic: point FOLD_FEATURES_CSV at a D30-derived table later
for the D30 -> D52 baseline model, nothing else here changes.
"""

from pathlib import Path

import pandas as pd

from modeling.folds import load_lines_with_label
from modeling.harness import run_flat_cv, run_nested_cv, summarize_across_repeats

OUT_DIR = Path(__file__).parent
FOLD_FEATURES_CSV = OUT_DIR / "fold_features_D11_full.csv"

# LOCO/LODO have a single "repeat" (repeat=0 for every fold), so
# summarize_across_repeats naturally pools all their out-of-fold
# predictions into one group before computing metrics -- which is exactly
# the required handling, since a single-line LOCO test fold can't support
# per-fold classification metrics (see modeling/README.md).
SCHEMES = ["plain", "donor_grouped", "loco", "lodo"]
# Pool-correction dropped from the active pipeline (see
# pool_correction_investigation.md) -- the mechanism behind its effect on
# results turned out ambiguous/contested even after independent
# investigation, not clearly a net improvement, and added a lot of
# complexity for an unclear benefit at this prototype stage.
# harness.py/features.py still support it (pool_correction=True) if this
# gets revisited.
CORRECTIONS = [False]

# Pre-registered before looking at results (modeling/README.md "Results
# distillation") -- decided in advance to avoid cherry-picking the
# best-looking config after the fact. The pre-registration named a
# specific MODEL per task (ridge / logistic_l2), so select_headline must
# filter on it: an earlier version didn't, which let the better-scoring
# non-pre-registered model (lasso) get quoted as "the" headline -- exactly
# the cherry-picking the pre-registration exists to prevent.
HEADLINE = {"scheme": "donor_grouped", "tuning": "nested", "pool_correction": False}
HEADLINE_MODEL = {"regression": "ridge", "classification": "logistic_l2"}


def summarize_nested_selections(nested_preds: pd.DataFrame) -> pd.DataFrame:
    """Per model: the modal (k, params) chosen across outer folds, plus the
    full spread. Nested CV re-selects per outer fold, so a single value is a
    summary, not "the" configuration -- hence carrying both."""
    rows = []
    for model, group in nested_preds.groupby("model"):
        per_fold = group[["repeat", "fold", "selected_k", "selected_param"]].drop_duplicates()
        combos = per_fold["selected_k"].astype(str) + "|" + per_fold["selected_param"]
        counts = combos.value_counts()
        modal_k, modal_param = counts.index[0].split("|", 1)
        rows.append(
            {
                "model": model,
                "k": int(modal_k),
                "param_str": modal_param,
                "n_outer_folds": len(per_fold),
                "selection_distribution": "; ".join(
                    f"k={c.split('|')[0]},{c.split('|', 1)[1]} x{n}" for c, n in counts.items()
                ),
            }
        )
    return pd.DataFrame(rows)


def run_all(
    task: str,
    fold_features_csv: Path = FOLD_FEATURES_CSV,
    out_suffix: str = "",
    out_dir: Path | None = None,
) -> pd.DataFrame:
    """out_dir: where nested_selections_{task}{suffix}.csv is written.
    Defaults to the package directory. Injectable because run_all has a
    side effect on disk, and the test suite was silently overwriting the
    committed modeling/nested_selections_regression.csv with synthetic
    fixture data on every run -- which is how corrupt content ended up
    committed."""
    fold_features = pd.read_csv(fold_features_csv)
    lines = load_lines_with_label()

    summaries = []
    all_selections: list[pd.DataFrame] = []
    for scheme in SCHEMES:
        for correction in CORRECTIONS:
            flat_preds = run_flat_cv(fold_features, lines, scheme, task, correction)
            flat_summary = summarize_across_repeats(flat_preds, task, ["model", "k", "param_str"])
            flat_summary["scheme"], flat_summary["pool_correction"], flat_summary["tuning"] = scheme, correction, "flat"
            summaries.append(flat_summary)

            nested_preds = run_nested_cv(fold_features, lines, scheme, task, correction)
            # Metrics MUST stay grouped on ["model"] alone: within one repeat
            # the folds together cover all 138 lines exactly once, and that is
            # what makes the per-repeat metric valid. Grouping by the selected
            # config as well would split a repeat into config-specific subsets
            # of lines and silently change the reported number.
            #
            # The selected configs are instead summarized alongside. Nested CV
            # picks per outer fold, so there is no single "the" config -- k and
            # param_str carry the modal choice (schema-compatible with flat
            # rows) and selection_distribution carries the full spread.
            # Previously these were computed per fold and then dropped
            # entirely, leaving nested rows with blank k/param_str.
            nested_summary = summarize_across_repeats(nested_preds, task, ["model"])
            nested_summary = nested_summary.merge(
                summarize_nested_selections(nested_preds), on="model", how="left"
            )
            all_selections.append(nested_preds.assign(scheme=scheme, pool_correction=correction))
            nested_summary["scheme"], nested_summary["pool_correction"], nested_summary["tuning"] = scheme, correction, "nested"
            summaries.append(nested_summary)

    # Per-outer-fold selections persisted in full, not just the modal summary,
    # so the whole selection record survives rather than being reduced away.
    if all_selections:
        selections = pd.concat(all_selections, ignore_index=True)
        cols = ["scheme", "pool_correction", "model", "repeat", "fold", "selected_k", "selected_param"]
        # Suffixed HERE rather than written unsuffixed and renamed afterwards:
        # the rename approach overwrote the committed baseline file before
        # moving it aside, which deleted the baseline selections from the repo.
        selections[cols].drop_duplicates().to_csv(
            (out_dir or OUT_DIR) / f"nested_selections_{task}{out_suffix}.csv", index=False
        )

    return pd.concat(summaries, ignore_index=True)


def select_headline(results: pd.DataFrame, task: str | None = None) -> pd.DataFrame:
    """The single pre-registered configuration. Pass `task` to also apply
    the pre-registered model filter (ridge / logistic_l2); without it the
    other model family is returned alongside, which is fine for a
    secondary comparison but must not be quoted as "the" headline."""
    selected = results[
        (results["scheme"] == HEADLINE["scheme"])
        & (results["tuning"] == HEADLINE["tuning"])
        & (results["pool_correction"] == HEADLINE["pool_correction"])
    ]
    if task is not None:
        selected = selected[selected["model"] == HEADLINE_MODEL[task]]
    return selected
