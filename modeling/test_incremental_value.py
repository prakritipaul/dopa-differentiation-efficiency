"""
Does the full D30 model beat `phat_DA` alone?

The open question from FINDINGS §A3. The numbers quoted there are NOT
comparable: the benchmark's ROC-AUC 0.960 is computed in-sample over all 136
lines (nothing is fitted, so nothing overfits, but it is still in-sample),
while the model's 0.945 is donor-grouped out-of-fold. One of those has seen
every line it scores and the other has not.

This puts BOTH through the identical procedure -- same folds, same nested
inner selection on training lines only, same metrics -- so the difference is
attributable to the feature set and nothing else:

    FULL       every free proportion + PC1..PCk, k and the penalty selected
               per outer fold by inner 3-fold CV (the pre-registered design)
    DA-ONLY    phat_DA as the single feature, penalty selected the same way

Reported on metrics that can still move near the ceiling. ROC-AUC is
saturated here by construction, so log loss and Brier carry the weight: they
respond to probability quality, not just ranking. Differences are PAIRED per
repeat (both models see identical folds), and the donor bootstrap resamples
the 20 donors rather than the 136 lines, since most lines share a donor.

Usage:
    uv run python -m modeling.test_incremental_value --timepoint D30
"""

import argparse
from itertools import product
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import brier_score_loss, log_loss, roc_auc_score
from sklearn.model_selection import StratifiedKFold

from modeling.folds import get_variant, load_lines_with_label
from modeling.harness import (
    ALL_PC_COLS,
    DECISION_THRESHOLD,
    build_line_level_for_fold,
    fit_predict,
    proportion_cols,
)
from modeling.models import PC_COUNT_GRID, models_for_task

OUT_DIR = Path(__file__).parent
SEED = 0
N_BOOT = 2000


def _candidates(spec, feature_sets):
    keys = list(spec.param_grid)
    for combo in product(*spec.param_grid.values()):
        for fs in feature_sets:
            yield dict(zip(keys, combo)), fs


def _inner_select(spec, train_line, y_train, strat, feature_sets, seed=SEED):
    """Pick (params, feature columns) on TRAINING lines only, inner 3-fold."""
    inner = StratifiedKFold(n_splits=3, shuffle=True, random_state=seed)
    best, best_score = None, -np.inf
    splits = list(inner.split(train_line, strat))
    for params, cols in _candidates(spec, feature_sets):
        scores = []
        for tr, te in splits:
            Xtr = train_line.iloc[tr][cols].to_numpy()
            Xte = train_line.iloc[te][cols].to_numpy()
            ytr, yte = y_train[tr], y_train[te]
            if len(np.unique(ytr)) < 2:
                continue
            _, prob = fit_predict(spec, params, Xtr, ytr, Xte)
            scores.append(roc_auc_score(yte, prob) if len(np.unique(yte)) == 2 else np.nan)
        s = np.nanmean(scores) if scores else np.nan
        if np.isfinite(s) and s > best_score:
            best, best_score = (params, cols), s
    return best


def run(timepoint: str, label_variant: str = "da_untreated", suffix: str = "_qualonly") -> pd.DataFrame:
    v = get_variant(label_variant)
    ff = pd.read_csv(OUT_DIR / f"fold_data/fold_features_{timepoint}{suffix}{v.suffix}.csv")
    lines = load_lines_with_label(label_variant)
    y_all = lines.set_index("cell_line")["success"].astype(int)
    donor = lines.set_index("cell_line")["donor"]

    spec = next(m for m in models_for_task("classification") if m.name == "logistic_l2")
    _, model_props = proportion_cols(ff)
    full_sets = [model_props + ALL_PC_COLS[:k] for k in PC_COUNT_GRID]
    da_sets = [["phat_DA"]]

    keys = ff.loc[ff.scheme == "donor_grouped", ["repeat", "fold"]].drop_duplicates()
    rows = []
    for repeat, fold in keys.itertuples(index=False):
        tr, te = build_line_level_for_fold(ff, "donor_grouped", repeat, fold, False)
        ytr = y_all.loc[tr.cell_line].to_numpy()
        yte = y_all.loc[te.cell_line].to_numpy()
        for name, sets in [("full", full_sets), ("da_only", da_sets)]:
            params, cols = _inner_select(spec, tr, ytr, ytr, sets)
            _, prob = fit_predict(spec, params, tr[cols].to_numpy(), ytr, te[cols].to_numpy())
            rows.append(pd.DataFrame({
                "model": name, "repeat": repeat, "cell_line": te.cell_line.values,
                "y": yte, "prob": prob, "k": len(cols),
            }))
    return pd.concat(rows, ignore_index=True)


def _metrics(g):
    p = np.clip(g["prob"].to_numpy(), 1e-6, 1 - 1e-6)
    y = g["y"].to_numpy()
    return pd.Series({
        "roc_auc": roc_auc_score(y, p) if len(np.unique(y)) == 2 else np.nan,
        "log_loss": log_loss(y, p, labels=[0, 1]),
        "brier": brier_score_loss(y, p),
        "accuracy": ((p >= DECISION_THRESHOLD).astype(int) == y).mean(),
    })


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--timepoint", default="D30", choices=["D11", "D30"])
    ap.add_argument("--label-variant", default="da_untreated")
    args = ap.parse_args()

    preds = run(args.timepoint, args.label_variant)
    per_repeat = preds.groupby(["model", "repeat"]).apply(_metrics, include_groups=False).reset_index()

    print(f"\n{args.timepoint}: full model vs phat_DA alone, identical donor-grouped folds,")
    print("nested inner selection on training lines only, mean +/- SD across 10 repeats\n")
    summ = per_repeat.groupby("model")[["roc_auc", "log_loss", "brier", "accuracy"]].agg(["mean", "std"])
    print(summ.round(4).to_string())

    print("\nPAIRED per-repeat difference (full - da_only); negative log_loss/brier favour FULL")
    wide = per_repeat.pivot(index="repeat", columns="model")
    for m in ["roc_auc", "log_loss", "brier", "accuracy"]:
        d = wide[(m, "full")] - wide[(m, "da_only")]
        wins = int((d > 0).sum()) if m in ("roc_auc", "accuracy") else int((d < 0).sum())
        print(f"  {m:9s} mean {d.mean():+.4f}  SD {d.std():.4f}  "
              f"full better in {wins}/10 repeats")

    ks = preds.groupby("model")["k"].agg(["min", "median", "max"])
    print(f"\nfeatures selected per outer fold:\n{ks.to_string()}")

    out = OUT_DIR / f"results/incremental_value_{args.timepoint}_{args.label_variant}.csv"
    per_repeat.to_csv(out, index=False)
    print(f"\nSaved {out}")


if __name__ == "__main__":
    main()
