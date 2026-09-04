"""
Timepoint-agnostic CV/tuning/metrics harness.

Consumes an already-extracted fold-features table (from
run_feature_extraction.py: one row per scheme/repeat/fold/cell_line/pool,
with PC1..PC10, phat_* proportions, and a train/test `split` column) plus
the model registry (models.py), and produces flat-CV and nested-CV
results. Never touches expression data or PCA fitting itself -- that
leakage-safe work already happened in features.py/run_feature_extraction.py.
Reused unchanged for D30 -> D52 later: just point it at a D30-derived
fold-features table.

Metrics are computed PER REPEAT (pooling only within that repeat's folds,
which together cover all 138 lines exactly once) and then reported as
mean +/- SD ACROSS repeats -- never pooled across repeats, which would
double-count each line's several repeat-predictions as if independent
(see modeling/README.md).

Known first-pass approximation: nested CV's inner loop selects
hyperparameters (k, regularization strength) by further splitting each
outer fold's TRAINING lines and reusing that outer fold's already-
extracted PCA features, rather than refitting PCA per inner sub-fold
(which would multiply the already-considerable feature-extraction cost
several-fold further). Flagged here and in modeling/README.md rather than
silently assumed.
"""

from itertools import product

import numpy as np
import pandas as pd
from sklearn.metrics import (
    average_precision_score,
    balanced_accuracy_score,
    brier_score_loss,
    confusion_matrix,
    f1_score,
    mean_absolute_error,
    r2_score,
    roc_auc_score,
    root_mean_squared_error,
)
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler

from modeling.features import pool_correct, pool_then_line_average
from modeling.models import PC_COUNT_GRID, ModelSpec, models_for_task

# The three D11 cell-type proportions sum to EXACTLY 1.0 for every row, so
# they are compositional: including all three alongside an intercept makes
# the design matrix rank-deficient. ALL_PROPORTION_COLS is for carrying the
# data around (line-level averaging, pool correction, reporting);
# MODEL_PROPORTION_COLS is what may enter a simultaneous fit, with
# phat_P_FPP kept implicit as 1 - FPP - NB.
#
# These two names are deliberately mirrored in (and imported by)
# feature_importance.py. An earlier version defined the model feature set
# separately in each module, which let them drift: the harness scored a
# 3-proportion model while the importance tables described a 2-proportion
# one. Keep this as the single definition.
ALL_PROPORTION_COLS = ["phat_FPP", "phat_NB", "phat_P_FPP"]
MODEL_PROPORTION_COLS = ["phat_FPP", "phat_NB"]
ALL_PC_COLS = [f"PC{i}" for i in range(1, 11)]
DECISION_THRESHOLD = 0.5  # probability -> class-call threshold; kept separate from the 0.2 outcome threshold


def _param_str(params: dict) -> str:
    return ",".join(f"{k}={v}" for k, v in sorted(params.items()))


def _param_grid_combos(param_grid: dict) -> list[dict]:
    keys = list(param_grid.keys())
    return [dict(zip(keys, combo)) for combo in product(*param_grid.values())]


def build_feature_matrix(line_level: pd.DataFrame, k: int) -> np.ndarray:
    """2 proportions + first k PCs (k=0..10). PCA components are
    hierarchical, so this truncates columns from one PCA fit rather than
    needing a distinct fit per k.

    Only MODEL_PROPORTION_COLS enters the matrix -- phat_P_FPP is the
    implicit reference category (1 - FPP - NB). Including all three made
    the design matrix rank-deficient once an intercept was present."""
    cols = MODEL_PROPORTION_COLS + ALL_PC_COLS[:k]
    return line_level[cols].to_numpy()


def fit_predict(
    model_spec: ModelSpec, params: dict, X_train: np.ndarray, y_train: np.ndarray, X_test: np.ndarray
) -> tuple[np.ndarray, np.ndarray | None]:
    """Standardize (fit on train only), fit a fresh estimator instance,
    predict. Returns (predictions, probabilities-or-None for classifiers)."""
    scaler = StandardScaler()
    X_train_s = scaler.fit_transform(X_train)
    X_test_s = scaler.transform(X_test)

    est = model_spec.estimator_factory()
    est.set_params(**params)
    est.fit(X_train_s, y_train)

    preds = est.predict(X_test_s)
    probs = est.predict_proba(X_test_s)[:, 1] if hasattr(est, "predict_proba") else None
    return preds, probs


def build_line_level_for_fold(
    fold_features: pd.DataFrame, scheme: str, repeat: int, fold: int, pool_correction: bool
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Returns (train_line_level, test_line_level): cell_line + proportion
    + PC columns, aggregated via pool-then-line two-stage averaging
    (features.pool_then_line_average). If pool_correction, each feature is
    residualized on pool identity using only this fold's TRAINING pool
    means (features.pool_correct) before line-level averaging."""
    sub = fold_features[
        (fold_features["scheme"] == scheme) & (fold_features["repeat"] == repeat) & (fold_features["fold"] == fold)
    ].copy()
    # All three proportions are carried through averaging/correction (they
    # are data, and downstream reporting uses phat_P_FPP); only
    # MODEL_PROPORTION_COLS reaches a fit, via build_feature_matrix.
    value_cols = ALL_PROPORTION_COLS + ALL_PC_COLS

    if pool_correction:
        train_mask = (sub["split"] == "train").to_numpy()
        sub = pool_correct(sub, value_cols, train_mask)

    train_line = pool_then_line_average(sub[sub["split"] == "train"], value_cols)
    test_line = pool_then_line_average(sub[sub["split"] == "test"], value_cols)
    return train_line, test_line


def compute_classification_metrics(y_true: np.ndarray, y_prob: np.ndarray, decision_threshold: float = DECISION_THRESHOLD) -> dict:
    if len(np.unique(y_true)) < 2:
        # rank-based / confusion-matrix metrics are undefined with only one class present
        return {k: np.nan for k in ("roc_auc", "pr_auc", "balanced_accuracy", "sensitivity", "specificity", "f1", "brier")}
    y_pred = (y_prob >= decision_threshold).astype(int)
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred, labels=[0, 1]).ravel()
    return {
        "roc_auc": roc_auc_score(y_true, y_prob),
        "pr_auc": average_precision_score(y_true, y_prob),
        "balanced_accuracy": balanced_accuracy_score(y_true, y_pred),
        "sensitivity": tp / (tp + fn) if (tp + fn) > 0 else np.nan,
        "specificity": tn / (tn + fp) if (tn + fp) > 0 else np.nan,
        "f1": f1_score(y_true, y_pred, zero_division=0),
        "brier": brier_score_loss(y_true, y_prob),
    }


def compute_regression_metrics(y_true: np.ndarray, y_pred: np.ndarray) -> dict:
    return {
        "mae": mean_absolute_error(y_true, y_pred),
        "rmse": root_mean_squared_error(y_true, y_pred),
        "r2": r2_score(y_true, y_pred),
        "frac_out_of_range": float(np.mean((y_pred < 0) | (y_pred > 1))),
    }


def _score_per_repeat(predictions: pd.DataFrame, task: str, group_cols: list[str]) -> pd.DataFrame:
    """One metric row per (group_cols + repeat) -- valid, since within one
    repeat each line appears in exactly one fold's test set."""
    rows = []
    for combo_key, combo_df in predictions.groupby(group_cols + ["repeat"], dropna=False):
        y_true = combo_df["y_true"].to_numpy()
        metrics = (
            compute_classification_metrics(y_true, combo_df["y_prob"].to_numpy())
            if task == "classification"
            else compute_regression_metrics(y_true, combo_df["y_pred"].to_numpy())
        )
        row = dict(zip(group_cols + ["repeat"], combo_key))
        row.update(metrics)
        row["n_lines"] = len(combo_df)
        rows.append(row)
    return pd.DataFrame(rows)


def summarize_across_repeats(predictions: pd.DataFrame, task: str, group_cols: list[str]) -> pd.DataFrame:
    """Mean +/- SD of each metric across repeats, for every distinct
    combo of group_cols (e.g. model/k/param_str for flat CV)."""
    per_repeat = _score_per_repeat(predictions, task, group_cols)
    metric_cols = [c for c in per_repeat.columns if c not in group_cols + ["repeat", "n_lines"]]
    summary = per_repeat.groupby(group_cols)[metric_cols].agg(["mean", "std"])
    summary.columns = [f"{m}_{stat}" for m, stat in summary.columns]
    return summary.reset_index()


def run_flat_cv(fold_features: pd.DataFrame, lines: pd.DataFrame, scheme: str, task: str, pool_correction: bool) -> pd.DataFrame:
    """For every (repeat, fold, model, k, param combo), fit on that fold's
    training lines and predict its test lines. Long-format per-line
    predictions -- hyperparameter selection (picking the best model/k/param
    by mean CV score) and performance reporting happen on the SAME rows
    here, which is what makes this "flat" (non-nested) CV; do that
    selection downstream via summarize_across_repeats."""
    y_col = "diff_efficiency" if task == "regression" else "success"
    label_lookup = lines.set_index("cell_line")[y_col]

    fold_keys = fold_features.loc[fold_features["scheme"] == scheme, ["repeat", "fold"]].drop_duplicates()
    rows = []
    for _, key in fold_keys.iterrows():
        repeat, fold = int(key["repeat"]), int(key["fold"])
        train_line, test_line = build_line_level_for_fold(fold_features, scheme, repeat, fold, pool_correction)
        y_train = label_lookup.loc[train_line["cell_line"]].to_numpy()
        y_test = label_lookup.loc[test_line["cell_line"]].to_numpy()

        for model_spec in models_for_task(task):
            for k in PC_COUNT_GRID:
                X_train = build_feature_matrix(train_line, k)
                X_test = build_feature_matrix(test_line, k)
                for param_combo in _param_grid_combos(model_spec.param_grid):
                    preds, probs = fit_predict(model_spec, param_combo, X_train, y_train, X_test)
                    for i, cell_line in enumerate(test_line["cell_line"]):
                        rows.append(
                            {
                                "scheme": scheme, "repeat": repeat, "fold": fold,
                                "model": model_spec.name, "k": k, "param_str": _param_str(param_combo),
                                "cell_line": cell_line, "y_true": y_test[i], "y_pred": preds[i],
                                "y_prob": probs[i] if probs is not None else None,
                            }
                        )
    return pd.DataFrame(rows)


def run_nested_cv(
    fold_features: pd.DataFrame,
    lines: pd.DataFrame,
    scheme: str,
    task: str,
    pool_correction: bool,
    inner_splits: int = 3,
    seed: int = 0,
) -> pd.DataFrame:
    """Outer loop = the already-extracted `scheme` folds. Inner loop =
    stratified K-fold over each outer fold's TRAINING lines only (reusing
    that outer fold's PCA features -- see module docstring), selecting
    (k, param) PER MODEL FAMILY by the inner tuning objective (PR-AUC for
    classification, MAE for regression). The selected combo is refit on
    the full outer training set and evaluated once on the outer test set
    -- the only prediction that counts toward the reported metric."""
    y_col = "diff_efficiency" if task == "regression" else "success"
    label_lookup = lines.set_index("cell_line")[y_col]
    success_lookup = lines.set_index("cell_line")["success"]

    fold_keys = fold_features.loc[fold_features["scheme"] == scheme, ["repeat", "fold"]].drop_duplicates()
    rows = []
    for _, key in fold_keys.iterrows():
        repeat, fold = int(key["repeat"]), int(key["fold"])
        train_line, test_line = build_line_level_for_fold(fold_features, scheme, repeat, fold, pool_correction)
        y_train_full = label_lookup.loc[train_line["cell_line"]].to_numpy()
        y_test = label_lookup.loc[test_line["cell_line"]].to_numpy()
        strat_train = success_lookup.loc[train_line["cell_line"]].to_numpy()

        inner_skf = StratifiedKFold(n_splits=inner_splits, shuffle=True, random_state=seed)
        inner_idx_splits = list(inner_skf.split(train_line, strat_train))

        for model_spec in models_for_task(task):
            best_score, best_combo = -np.inf, None
            for k in PC_COUNT_GRID:
                X_train_full = build_feature_matrix(train_line, k)
                for param_combo in _param_grid_combos(model_spec.param_grid):
                    inner_scores = []
                    for inner_train_idx, inner_val_idx in inner_idx_splits:
                        X_it, X_iv = X_train_full[inner_train_idx], X_train_full[inner_val_idx]
                        y_it, y_iv = y_train_full[inner_train_idx], y_train_full[inner_val_idx]
                        preds, probs = fit_predict(model_spec, param_combo, X_it, y_it, X_iv)
                        if task == "classification":
                            if len(np.unique(y_iv)) < 2:
                                continue
                            inner_scores.append(average_precision_score(y_iv, probs))
                        else:
                            inner_scores.append(-mean_absolute_error(y_iv, preds))
                    if not inner_scores:
                        continue
                    mean_inner = float(np.mean(inner_scores))
                    if mean_inner > best_score:
                        best_score, best_combo = mean_inner, (k, param_combo)

            if best_combo is None:
                continue
            k, param_combo = best_combo
            X_train_full = build_feature_matrix(train_line, k)
            X_test = build_feature_matrix(test_line, k)
            preds, probs = fit_predict(model_spec, param_combo, X_train_full, y_train_full, X_test)
            for i, cell_line in enumerate(test_line["cell_line"]):
                rows.append(
                    {
                        "scheme": scheme, "repeat": repeat, "fold": fold,
                        "model": model_spec.name, "selected_k": k, "selected_param": _param_str(param_combo),
                        "cell_line": cell_line, "y_true": y_test[i], "y_pred": preds[i],
                        "y_prob": probs[i] if probs is not None else None,
                    }
                )
    return pd.DataFrame(rows)
