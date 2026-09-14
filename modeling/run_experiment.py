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
FOLD_FEATURES_CSV = OUT_DIR / "fold_data/fold_features_D11_full.csv"

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


# One (scheme, correction) pair is the resume unit. LOCO is 138 folds and
# LODO 20, so a scheme is both the natural checkpoint boundary and the
# expensive thing to lose. No ETA is printed: the four schemes differ in cost
# by more than an order of magnitude, so a rate averaged over them would be
# actively misleading.
RESUME_KEY = ["scheme", "pool_correction", "tuning"]
# Suffixing happens at the path, not by renaming afterwards: the rename
# approach overwrote the committed baseline before moving it aside, which
# deleted the baseline selections from the repo.
SEL_COLS = ["scheme", "pool_correction", "model", "repeat", "fold", "selected_k", "selected_param"]


def _append(path: Path, chunk: pd.DataFrame) -> None:
    """Append `chunk`, reindexed to the header already on disk. The first
    write fixes column order and every later chunk must match it -- flat and
    nested summaries do not carry identical columns, so relying on chunk
    ordering would silently misalign values under the wrong headers."""
    if path.exists():
        chunk = chunk.reindex(columns=pd.read_csv(path, nrows=0).columns.tolist())
        chunk.to_csv(path, mode="a", header=False, index=False)
    else:
        chunk.to_csv(path, index=False)


def run_all(
    task: str,
    fold_features_csv: Path = FOLD_FEATURES_CSV,
    out_suffix: str = "",
    out_dir: Path | None = None,
    resume: bool = True,
    label_variant: str = "published",
) -> pd.DataFrame:
    """out_dir: where results_{task}{suffix}.csv and
    nested_selections_{task}{suffix}.csv are written. Defaults to the package
    directory. Injectable because run_all has a side effect on disk, and the
    test suite was silently overwriting the committed
    modeling/nested_selections_regression.csv with synthetic fixture data on
    every run -- which is how corrupt content ended up committed.

    resume: each (scheme, correction) pair's summary rows are appended as soon
    as that scheme finishes, so an interrupted run restarts from what is
    already on disk instead of recomputing LOCO's 138 folds. Pass resume=False
    to discard existing output and recompute from scratch.

    Resume trusts what is on disk. Every row records the fold-features file it
    came from, and resuming against a different basis is refused rather than
    silently mixing two provenances -- that exact failure produced the
    mixed-basis importance tables in modeling/README.md "Second correctness
    review". A change that does NOT alter the basis filename -- different
    hyperparameters, a fix inside harness.py, regenerated features at the same
    path -- is undetectable here. Use resume=False after any of those."""
    fold_features = pd.read_csv(fold_features_csv)
    lines = load_lines_with_label(label_variant)

    # The features table records which variant produced it. If it disagrees with
    # the variant being scored, the outcome/cohort/threshold do not match the
    # features and every number would be quietly wrong -- the fold-features
    # FILENAME cannot catch this, since a label change need not rename the file.
    ff_seen = set(fold_features.get("label_variant", pd.Series(dtype=object)).dropna().unique())
    if ff_seen and ff_seen != {label_variant}:
        raise ValueError(
            f"{fold_features_csv.name} was extracted under label_variant {sorted(ff_seen)} "
            f"but is being scored as {label_variant!r}."
        )

    base_dir = out_dir or OUT_DIR
    res_path = base_dir / f"results/results_{task}{out_suffix}.csv"
    sel_path = base_dir / f"results/nested_selections_{task}{out_suffix}.csv"
    res_path.parent.mkdir(parents=True, exist_ok=True)
    basis = Path(fold_features_csv).name

    done: set[tuple] = set()
    if res_path.exists():
        if resume:
            prev = pd.read_csv(res_path)
            # Strict rule for the label, unlike the laxer one below for the
            # basis filename: a label change is INVISIBLE in that filename, so
            # an absent column means "provenance unknown", not "compatible".
            seen_lv = set(prev.get("label_variant", pd.Series(dtype=object)).dropna().unique())
            if seen_lv != {label_variant}:
                raise ValueError(
                    f"{res_path} records label_variant="
                    f"{sorted(seen_lv) if seen_lv else 'UNRECORDED (no label_variant column)'} "
                    f"but this run uses {label_variant!r}. Resuming would mix two outcomes in "
                    f"one results table. Rerun with resume=False, or use a new out_suffix."
                )
            seen = set(prev.get("fold_features", pd.Series(dtype=object)).dropna().unique())
            if seen and seen != {basis}:
                raise ValueError(
                    f"{res_path} holds rows from {sorted(seen)} but this run uses "
                    f"{basis}. Resuming would mix feature bases in one table. "
                    f"Rerun with resume=False, or point out_suffix at a new file."
                )
            done = set(map(tuple, prev[RESUME_KEY].drop_duplicates().itertuples(index=False)))
            print(f"resuming from {res_path}: {len(done)} summary groups on disk", flush=True)
        else:
            res_path.unlink()
            sel_path.unlink(missing_ok=True)

    for scheme in SCHEMES:
        for correction in CORRECTIONS:
            if {(scheme, correction, "flat"), (scheme, correction, "nested")} <= done:
                print(f"skip {scheme} (correction={correction}): already on disk", flush=True)
                continue
            flat_preds = run_flat_cv(fold_features, lines, scheme, task, correction)
            flat_summary = summarize_across_repeats(flat_preds, task, ["model", "k", "param_str"])
            flat_summary["scheme"], flat_summary["pool_correction"], flat_summary["tuning"] = scheme, correction, "flat"

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
            nested_summary["scheme"], nested_summary["pool_correction"], nested_summary["tuning"] = scheme, correction, "nested"

            # Appended the moment the scheme finishes -- that is the whole
            # point of resume. Per-outer-fold selections are persisted in full,
            # not just the modal summary, and are appended per scheme for the
            # same reason: writing them once at the end would silently drop
            # every scheme that was restored from disk rather than recomputed.
            chunk = pd.concat([flat_summary, nested_summary], ignore_index=True)
            chunk["fold_features"] = basis
            chunk["label_variant"] = label_variant
            _append(res_path, chunk)
            _append(
                sel_path,
                nested_preds.assign(scheme=scheme, pool_correction=correction)[SEL_COLS]
                .drop_duplicates()
                .assign(label_variant=label_variant),
            )
            print(f"wrote {scheme} (correction={correction})", flush=True)

    if not res_path.exists():
        return pd.DataFrame()
    # Read back rather than returning this run's rows: on a resumed run the
    # in-memory summaries hold only the schemes recomputed this time.
    return pd.read_csv(res_path)


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
