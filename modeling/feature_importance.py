"""
Which D11 features actually drive the predictions?

Produces one table per task (regression, classification):

  Feature | Univariate | Coefficient | Drop-one Delta (LOCO) | SHAP | Technical covariate

See modeling/feature_importance_plan.md for the full design rationale
(own reasoning + two rounds of independent Codex review). Key points:

- Compositional proportions: phat_FPP + phat_NB + phat_P_FPP sum to
  exactly 1.0 for every line (verified). MODEL_FEATURE_COLS therefore
  drops phat_P_FPP from every SIMULTANEOUS multi-feature fit (coefficient
  fitting, SHAP, LOCO/permutation baselines) -- it's kept implicit as
  `1 - FPP - NB`. Univariate association and technical-covariate
  association don't need a joint fit, so they still cover all 3.
- Config selection uses a one-SE rule (simplest config within 1 SE of the
  flat-CV grid's best score), not the raw noisy argmax across 55 combos.
- LOCO: every PC individually, phat_FPP and phat_NB individually (both
  are explicit MODEL_FEATURE_COLS members, so dropping either is a valid
  comparison); phat_P_FPP's row instead shows the GROUPED delta (both
  proportions dropped together), since P_FPP was never an explicit
  baseline feature to drop in the first place (Codex: a 3-proportion
  simultaneous fit is rank-deficient, so there's no valid "with P_FPP"
  baseline). All deltas are PAIRED per-fold differences (metric_without -
  metric_with computed within each fold, then averaged), not two
  independently-averaged CV scores subtracted at the end.
- Permutation importance: PCs individually; the 2 model proportions
  permuted AS ONE GROUP (shuffle which line's (FPP, NB) pair goes to
  which row) -- unlike LOCO, permuting one alone while the other stays
  fixed would imply an impossible third coordinate. Also paired per-fold.
- SHAP: computed directly (coef_j * (x_ij - mean_j)), the exact
  closed-form for a linear model -- mathematically identical to
  shap.LinearExplainer's default (independent-feature) output. The `shap`
  package's build failed in this environment (llvmlite/numba wouldn't
  compile for this Python/platform combo) so this avoids that dependency
  entirely rather than fighting the build.
- Technical covariate association: eta-squared by pool, recomputed fresh
  on the ACTUAL features the models use (005's pool11-excluded PCA +
  proportions), not reusing 008's numbers (which used the uncorrected,
  pool11-included PCA -- a different basis).
"""

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.preprocessing import StandardScaler

from modeling.folds import load_lines_with_label
from modeling.harness import (
    ALL_PC_COLS,
    build_feature_matrix,
    build_line_level_for_fold,
    compute_classification_metrics,
    compute_regression_metrics,
)
from modeling.models import ModelSpec, models_for_task

REPO_ROOT = Path(__file__).parent.parent
OUT_DIR = Path(__file__).parent
D11_PCA_CSV = REPO_ROOT / "metadata_eda" / "d11_pca_coords_per_line.csv"
D11_PROPORTIONS_CSV = REPO_ROOT / "metadata_eda" / "d11_celltype_proportions_with_se.csv"
QUALIFYING_COMBOS_CSV = REPO_ROOT / "metadata_eda" / "qualifying_cell_line_pool_min10_per_timepoint.csv"

ALL_PROPORTION_COLS = ["phat_FPP", "phat_NB", "phat_P_FPP"]
MODEL_PROPORTION_COLS = ["phat_FPP", "phat_NB"]  # phat_P_FPP excluded from joint fits (compositional)


def _load_module(name: str):
    path = REPO_ROOT / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


m008 = _load_module("008_technical_covariate_associations")


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------


def load_full_fit_features() -> pd.DataFrame:
    """Line-level D11 features (3 proportions + PC1..PC10) fit on ALL
    qualifying lines (no held-out split) -- reuses 005's global
    (non-pool11) PCA fit and 004's proportions directly."""
    qualifying_lines = pd.read_csv(QUALIFYING_COMBOS_CSV)["cell_line"].unique()

    pcs = pd.read_csv(D11_PCA_CSV)
    pcs = pcs[pcs["cell_line"].isin(qualifying_lines)]

    props_long = pd.read_csv(D11_PROPORTIONS_CSV)
    props_wide = props_long.pivot(index="cell_line", columns="celltype", values="phat").reset_index()
    props_wide.columns = ["cell_line"] + [f"phat_{c}" for c in props_wide.columns[1:]]

    combined = pcs.merge(props_wide, on="cell_line", how="inner")
    assert len(combined) == 138, f"expected 138 lines, got {len(combined)}"
    return combined


def _model_feature_names(k: int) -> list[str]:
    return MODEL_PROPORTION_COLS + ALL_PC_COLS[:k]


# ---------------------------------------------------------------------------
# Config selection (one-SE rule)
# ---------------------------------------------------------------------------


def _simplicity_key(model: str, k: int, params: dict) -> tuple:
    """Lower = simpler. Smaller k first; then stronger regularization
    (larger alpha for ridge/lasso, smaller C for logistic, since C is the
    inverse of regularization strength)."""
    if "C" in params:
        reg_strength_rank = params["C"]  # smaller C = more regularization = simpler
    else:
        reg_strength_rank = -params.get("alpha", 0)  # larger alpha = more regularization = simpler
    return (k, reg_strength_rank)


def select_winning_config_one_se(
    results_csv: Path,
    model: str,
    scheme: str,
    metric: str,
    minimize: bool,
    n_repeats: int | None = None,
) -> tuple[int, dict]:
    """Simplest (k, params) within one STANDARD ERROR of the flat-CV
    grid's best mean score, rather than the raw (noisy, at n=138 across
    55 candidates) argmax.

    SE = SD_across_repeats / sqrt(n_repeats). An earlier version used the
    raw SD as the tolerance, which is ~sqrt(n_repeats) times too wide
    (~3x at 10 repeats) and therefore selected simpler models than a real
    one-SE rule would. n_repeats is read from the fold-assignment data
    when not supplied."""
    df = pd.read_csv(results_csv)
    flat = df[(df["tuning"] == "flat") & (df["model"] == model) & (df["scheme"] == scheme)].copy()

    mean_col = metric if metric.endswith("_mean") else f"{metric}_mean"
    std_col = mean_col.replace("_mean", "_std")

    if n_repeats is None:
        n_repeats = _infer_n_repeats(scheme)

    best_row = flat.loc[flat[mean_col].idxmin() if minimize else flat[mean_col].idxmax()]
    tolerance = best_row[std_col] / np.sqrt(n_repeats)
    if minimize:
        within_tol = flat[flat[mean_col] <= best_row[mean_col] + tolerance]
    else:
        within_tol = flat[flat[mean_col] >= best_row[mean_col] - tolerance]

    within_tol = within_tol.copy()
    within_tol["_simplicity"] = within_tol.apply(
        lambda r: _simplicity_key(model, int(r["k"]), _parse_param_str(r["param_str"])), axis=1
    )
    simplest = within_tol.sort_values("_simplicity").iloc[0]
    return int(simplest["k"]), _parse_param_str(simplest["param_str"])


def _parse_param_str(s: str) -> dict:
    out = {}
    for part in s.split(","):
        key, val = part.split("=")
        out[key] = float(val)
    return out


def _infer_n_repeats(scheme: str, fold_features_csv: Path | None = None) -> int:
    """Number of distinct repeats for a scheme, needed to turn the
    across-repeat SD into a standard error. LOCO/LODO have a single
    repeat by construction."""
    path = fold_features_csv or (OUT_DIR / "fold_features_D11_full.csv")
    if not path.exists():
        path = OUT_DIR / "fold_features_D11.csv"
    ff = pd.read_csv(path, usecols=["scheme", "repeat"])
    return int(ff.loc[ff["scheme"] == scheme, "repeat"].nunique())


# ---------------------------------------------------------------------------
# Coefficients (full-fit + fold stability)
# ---------------------------------------------------------------------------


def _fit_coefficients(
    model_spec: ModelSpec, params: dict, X: np.ndarray, y: np.ndarray, feature_names: list[str]
) -> pd.Series:
    scaler = StandardScaler()
    X_s = scaler.fit_transform(X)
    est = model_spec.estimator_factory()
    est.set_params(**params)
    est.fit(X_s, y)
    coef = np.ravel(est.coef_)
    return pd.Series(coef, index=feature_names)


def fit_full_model_coefficients(model_spec: ModelSpec, task: str, k: int, params: dict) -> pd.Series:
    features = load_full_fit_features()
    lines = load_lines_with_label()
    features = features.merge(lines[["cell_line", "diff_efficiency", "success"]], on="cell_line")

    feature_names = _model_feature_names(k)
    X = features[feature_names].to_numpy()
    y = features["diff_efficiency" if task == "regression" else "success"].to_numpy()
    return _fit_coefficients(model_spec, params, X, y, feature_names)


def fold_coefficient_stability(
    fold_features: pd.DataFrame, lines: pd.DataFrame, scheme: str, task: str, model_spec: ModelSpec, k: int, params: dict
) -> pd.DataFrame:
    y_col = "diff_efficiency" if task == "regression" else "success"
    label_lookup = lines.set_index("cell_line")[y_col]
    feature_names = _model_feature_names(k)

    fold_keys = fold_features.loc[fold_features["scheme"] == scheme, ["repeat", "fold"]].drop_duplicates()
    rows = []
    for _, key in fold_keys.iterrows():
        repeat, fold = int(key["repeat"]), int(key["fold"])
        train_line, _test_line = build_line_level_for_fold(fold_features, scheme, repeat, fold, pool_correction=False)
        X_train = train_line[feature_names].to_numpy()
        y_train = label_lookup.loc[train_line["cell_line"]].to_numpy()
        coefs = _fit_coefficients(model_spec, params, X_train, y_train, feature_names)
        rows.append({"repeat": repeat, "fold": fold, **coefs.to_dict()})
    return pd.DataFrame(rows)


def summarize_coefficients(fold_coefs: pd.DataFrame, full_fit_coefs: pd.Series, feature_names: list[str]) -> pd.DataFrame:
    rows = []
    for feat in feature_names:
        vals = fold_coefs[feat]
        nonzero = vals.abs() > 1e-10
        rows.append(
            {
                "feature": feat,
                "full_fit_coef": full_fit_coefs[feat],
                "regularized_to_zero": bool(abs(full_fit_coefs[feat]) < 1e-10),
                "selection_frequency": float(nonzero.mean()),
                "fold_coef_mean_given_nonzero": vals[nonzero].mean() if nonzero.any() else np.nan,
                "fold_coef_std_given_nonzero": vals[nonzero].std() if nonzero.sum() > 1 else np.nan,
            }
        )
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Univariate association
# ---------------------------------------------------------------------------


def univariate_association(task: str) -> pd.Series:
    features = load_full_fit_features()
    lines = load_lines_with_label()
    df = features.merge(lines[["cell_line", "diff_efficiency", "success"]], on="cell_line")
    y = df["diff_efficiency" if task == "regression" else "success"]

    feature_cols = ALL_PROPORTION_COLS + ALL_PC_COLS
    return pd.Series({f: spearmanr(df[f], y)[0] for f in feature_cols})


# ---------------------------------------------------------------------------
# Paired LOCO and permutation importance
# ---------------------------------------------------------------------------


def _score(task: str, y_true: np.ndarray, y_pred: np.ndarray, y_prob: np.ndarray | None) -> float:
    """Higher = better, uniformly, so deltas mean the same thing in both tasks."""
    if task == "regression":
        return -compute_regression_metrics(y_true, y_pred)["mae"]
    metrics = compute_classification_metrics(y_true, y_prob)
    return metrics["roc_auc"]


def _fit_predict(model_spec: ModelSpec, params: dict, X_train, y_train, X_test):
    scaler = StandardScaler()
    X_train_s = scaler.fit_transform(X_train)
    X_test_s = scaler.transform(X_test)
    est = model_spec.estimator_factory()
    est.set_params(**params)
    est.fit(X_train_s, y_train)
    preds = est.predict(X_test_s)
    probs = est.predict_proba(X_test_s)[:, 1] if hasattr(est, "predict_proba") else None
    return preds, probs


def paired_loco_deltas(
    fold_features: pd.DataFrame, lines: pd.DataFrame, scheme: str, task: str, model_spec: ModelSpec, k: int, params: dict
) -> pd.DataFrame:
    """Per feature (PCs individually, FPP/NB individually, and a
    'grouped_proportions' entry dropping both FPP+NB together): paired
    per-fold delta = score_with - score_without (positive = feature helps,
    since removing it lowered the score), averaged across folds."""
    y_col = "diff_efficiency" if task == "regression" else "success"
    label_lookup = lines.set_index("cell_line")[y_col]
    baseline_features = _model_feature_names(k)

    drop_targets = {f: [f] for f in baseline_features}
    drop_targets["grouped_proportions"] = MODEL_PROPORTION_COLS

    fold_keys = fold_features.loc[fold_features["scheme"] == scheme, ["repeat", "fold"]].drop_duplicates()
    per_fold_deltas = {name: [] for name in drop_targets}

    for _, key in fold_keys.iterrows():
        repeat, fold = int(key["repeat"]), int(key["fold"])
        train_line, test_line = build_line_level_for_fold(fold_features, scheme, repeat, fold, pool_correction=False)
        y_train = label_lookup.loc[train_line["cell_line"]].to_numpy()
        y_test = label_lookup.loc[test_line["cell_line"]].to_numpy()

        X_train_full = train_line[baseline_features].to_numpy()
        X_test_full = test_line[baseline_features].to_numpy()
        preds, probs = _fit_predict(model_spec, params, X_train_full, y_train, X_test_full)
        score_with = _score(task, y_test, preds, probs)

        for name, drop_cols in drop_targets.items():
            remaining = [c for c in baseline_features if c not in drop_cols]
            X_train_d = train_line[remaining].to_numpy()
            X_test_d = test_line[remaining].to_numpy()
            preds_d, probs_d = _fit_predict(model_spec, params, X_train_d, y_train, X_test_d)
            score_without = _score(task, y_test, preds_d, probs_d)
            # positive = removing this feature hurt performance = feature helps
            per_fold_deltas[name].append(score_with - score_without)

    return pd.DataFrame(
        {
            "feature": list(per_fold_deltas.keys()),
            "loco_delta_mean": [np.mean(v) for v in per_fold_deltas.values()],
            "loco_delta_std": [np.std(v) for v in per_fold_deltas.values()],
        }
    )


def paired_permutation_deltas(
    fold_features: pd.DataFrame,
    lines: pd.DataFrame,
    scheme: str,
    task: str,
    model_spec: ModelSpec,
    k: int,
    params: dict,
    seed: int = 0,
) -> pd.DataFrame:
    """Per feature (PCs individually; FPP+NB permuted together as one
    group): paired per-fold delta = score_unpermuted - score_permuted
    (positive = feature helps), using a model fit fresh per fold (same
    fold boundaries, standard baseline feature set)."""
    y_col = "diff_efficiency" if task == "regression" else "success"
    label_lookup = lines.set_index("cell_line")[y_col]
    baseline_features = _model_feature_names(k)
    pc_cols = [c for c in baseline_features if c.startswith("PC")]

    permute_targets = {pc: [pc] for pc in pc_cols}
    permute_targets["proportions_group"] = MODEL_PROPORTION_COLS

    rng = np.random.default_rng(seed)
    fold_keys = fold_features.loc[fold_features["scheme"] == scheme, ["repeat", "fold"]].drop_duplicates()
    per_fold_deltas = {name: [] for name in permute_targets}

    for _, key in fold_keys.iterrows():
        repeat, fold = int(key["repeat"]), int(key["fold"])
        train_line, test_line = build_line_level_for_fold(fold_features, scheme, repeat, fold, pool_correction=False)
        y_train = label_lookup.loc[train_line["cell_line"]].to_numpy()
        y_test = label_lookup.loc[test_line["cell_line"]].to_numpy()

        X_train = train_line[baseline_features].to_numpy()
        X_test = test_line[baseline_features].to_numpy()

        scaler = StandardScaler()
        X_train_s = scaler.fit_transform(X_train)
        X_test_s = scaler.transform(X_test)
        est = model_spec.estimator_factory()
        est.set_params(**params)
        est.fit(X_train_s, y_train)

        preds = est.predict(X_test_s)
        probs = est.predict_proba(X_test_s)[:, 1] if hasattr(est, "predict_proba") else None
        score_unpermuted = _score(task, y_test, preds, probs)

        for name, cols in permute_targets.items():
            idx = [baseline_features.index(c) for c in cols]
            X_test_perm = X_test_s.copy()
            perm_order = rng.permutation(X_test_perm.shape[0])
            X_test_perm[:, idx] = X_test_perm[perm_order][:, idx]

            preds_p = est.predict(X_test_perm)
            probs_p = est.predict_proba(X_test_perm)[:, 1] if hasattr(est, "predict_proba") else None
            score_permuted = _score(task, y_test, preds_p, probs_p)
            per_fold_deltas[name].append(score_unpermuted - score_permuted)

    return pd.DataFrame(
        {
            "feature": list(per_fold_deltas.keys()),
            "perm_delta_mean": [np.mean(v) for v in per_fold_deltas.values()],
            "perm_delta_std": [np.std(v) for v in per_fold_deltas.values()],
        }
    )


# ---------------------------------------------------------------------------
# SHAP (closed-form for linear models -- see module docstring)
# ---------------------------------------------------------------------------


def compute_shap_like(full_fit_coefs: pd.Series, feature_names: list[str]) -> pd.Series:
    """mean |coef_j * z_ij| where z is the STANDARDIZED feature value.

    The coefficients being reported come from a model fit on
    StandardScaler output, so their units are "per 1 SD of the feature".
    They must therefore be multiplied by standardized deviations, not raw
    ones -- an earlier version used raw (x - mean), which scaled every
    contribution by that feature's raw SD and so systematically inflated
    wide-scale features (PCs, raw SD ~1-9) relative to narrow ones
    (proportions, raw SD ~0.05). Using z = (x - mean)/std puts every
    feature on the same footing, which is the whole point of the column."""
    features = load_full_fit_features()
    subset = features[feature_names]
    z = (subset - subset.mean()) / subset.std(ddof=0)
    contributions = z.mul(full_fit_coefs[feature_names], axis=1)
    return contributions.abs().mean()


# ---------------------------------------------------------------------------
# Technical covariate association (recomputed on actual model features)
# ---------------------------------------------------------------------------


def technical_covariate_association() -> pd.Series:
    qualifying = pd.read_csv(QUALIFYING_COMBOS_CSV)[["cell_line", "pool"]]
    features = load_full_fit_features()
    df = features.merge(qualifying, on="cell_line", how="inner")

    eta_sq = {}
    for f in ALL_PROPORTION_COLS + ALL_PC_COLS:
        eta_sq[f] = m008.eta_squared_by_pool(df[f], df["pool"])
    return pd.Series(eta_sq)


def bucket_eta_sq(x: float) -> str:
    if x >= 0.5:
        return "High (pool-associated)"
    if x >= 0.15:
        return "Moderate"
    return "Low"


# ---------------------------------------------------------------------------
# Assembly
# ---------------------------------------------------------------------------


def run_for_model(
    fold_features_csv: Path, results_csv: Path, model_name: str, task: str, scheme: str, metric: str, minimize: bool
) -> pd.DataFrame:
    model_spec = next(m for m in models_for_task(task) if m.name == model_name)
    k, params = select_winning_config_one_se(results_csv, model_name, scheme, metric, minimize)
    print(f"{model_name}: winning config (one-SE rule) k={k}, params={params}")

    feature_names = _model_feature_names(k)
    full_fit_coefs = fit_full_model_coefficients(model_spec, task, k, params)
    lines = load_lines_with_label()

    fold_features = pd.read_csv(fold_features_csv)
    fold_coefs = fold_coefficient_stability(fold_features, lines, scheme, task, model_spec, k, params)
    coef_summary = summarize_coefficients(fold_coefs, full_fit_coefs, feature_names)

    univariate = univariate_association(task)
    loco = paired_loco_deltas(fold_features, lines, scheme, task, model_spec, k, params)
    perm = paired_permutation_deltas(fold_features, lines, scheme, task, model_spec, k, params)
    shap_like = compute_shap_like(full_fit_coefs, feature_names)
    tech_cov = technical_covariate_association()

    all_features = ALL_PROPORTION_COLS + [f"PC{i}" for i in range(1, k + 1)]
    rows = []
    for feat in all_features:
        is_p_fpp = feat == "phat_P_FPP"
        is_proportion = feat in ALL_PROPORTION_COLS
        coef_row = coef_summary[coef_summary["feature"] == feat]
        loco_row = loco[loco["feature"] == ("grouped_proportions" if is_p_fpp else feat)]
        # Permutation always groups the proportions (permuting one alone
        # implies an impossible third coordinate), so all 3 proportion
        # rows share the group's value; PCs get their individual value.
        perm_row = perm[perm["feature"] == ("proportions_group" if is_proportion else feat)]

        rows.append(
            {
                "feature": feat,
                "univariate_spearman": round(float(univariate.get(feat, np.nan)), 4),
                "coefficient": "reference (=1-FPP-NB)" if is_p_fpp else round(float(coef_row["full_fit_coef"].iloc[0]), 4),
                "regularized_to_zero": None if is_p_fpp else bool(coef_row["regularized_to_zero"].iloc[0]),
                "selection_frequency": None if is_p_fpp else round(float(coef_row["selection_frequency"].iloc[0]), 2),
                "loco_delta": round(float(loco_row["loco_delta_mean"].iloc[0]), 4) if len(loco_row) else np.nan,
                "loco_note": "grouped (FPP+NB dropped together)" if is_p_fpp else "individual",
                "permutation_delta": round(float(perm_row["perm_delta_mean"].iloc[0]), 4) if len(perm_row) else np.nan,
                "permutation_note": "grouped (all proportions permuted together)" if is_proportion else "individual",
                "shap_mean_abs": round(float(shap_like.get(feat, np.nan)), 4) if not is_p_fpp else np.nan,
                "technical_covariate_eta2": round(float(tech_cov[feat]), 3),
                "technical_covariate_label": bucket_eta_sq(tech_cov[feat]),
            }
        )

    table = pd.DataFrame(rows)
    out_csv = OUT_DIR / f"feature_importance_table_{task}_{model_name}.csv"
    table.to_csv(out_csv, index=False)
    print(table.to_string(index=False))
    print(f"Saved {out_csv}\n")
    return table


def main(fold_features_csv: Path = OUT_DIR / "fold_features_D11_full.csv") -> None:
    """All four models, not just the L1 pair.

    The L1 variants (lasso, logistic_l1) are the more informative ones
    here -- only they produce sparsity, so `regularized_to_zero` and
    `selection_frequency` are meaningful; ridge/L2 never zero anything.
    But ridge/logistic_l2 are the PRE-REGISTERED performance models, so
    both are reported to avoid any appearance of picking whichever model
    told the nicer story."""
    for model_name in ("lasso", "ridge"):
        run_for_model(
            fold_features_csv, OUT_DIR / "results_regression.csv", model_name, "regression", "donor_grouped", "mae_mean", True
        )
    for model_name in ("logistic_l1", "logistic_l2"):
        run_for_model(
            fold_features_csv,
            OUT_DIR / "results_classification.csv",
            model_name,
            "classification",
            "donor_grouped",
            "roc_auc_mean",
            False,
        )


if __name__ == "__main__":
    main()
