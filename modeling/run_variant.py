"""
Run the full experiment grid against an ALTERNATIVE fold-features table,
writing suffixed outputs so the baseline results are never overwritten.

Used for the "PCA fit on qualifying cells only" comparison (issue #2 from
the Codex audit): the baseline fits the D11 HVG/PCA basis on all cells
except the fold's held-out lines and pool11 -- which includes ~25,400
cells from 39 lines outside the 138 -- while the variant restricts the fit
to the qualifying (cell_line, pool) combos. Everything downstream of
feature extraction is identical, so any difference in results is
attributable to the fitting population alone.

Usage:
    uv run python -m modeling.run_variant --suffix _qualonly

Reads  modeling/fold_features_D11{suffix}.csv
Writes modeling/results_{task}{suffix}.csv
       modeling/nested_selections_{task}{suffix}.csv
and prints a baseline-vs-variant headline comparison.
"""

import argparse
from pathlib import Path

import pandas as pd

from modeling.run_experiment import HEADLINE, HEADLINE_MODEL, run_all, select_headline

OUT_DIR = Path(__file__).parent

METRICS = {
    "regression": ["mae_mean", "mae_std", "rmse_mean", "rmse_std", "r2_mean", "r2_std"],
    "classification": ["roc_auc_mean", "roc_auc_std", "pr_auc_mean", "pr_auc_std",
                       "balanced_accuracy_mean", "balanced_accuracy_std"],
}


def compare_headline(task: str, baseline_csv: Path, variant: pd.DataFrame) -> pd.DataFrame:
    """Baseline vs variant on the pre-registered config, side by side."""
    cols = METRICS[task]
    var_row = select_headline(variant, task)
    rows = [{"run": "variant", **{c: var_row[c].iloc[0] for c in cols}}]
    if baseline_csv.exists():
        base_row = select_headline(pd.read_csv(baseline_csv), task)
        if len(base_row):
            rows.insert(0, {"run": "baseline", **{c: base_row[c].iloc[0] for c in cols}})
    out = pd.DataFrame(rows)
    if len(out) == 2:
        delta = {"run": "delta"}
        for c in cols:
            delta[c] = out[c].iloc[1] - out[c].iloc[0]
        out = pd.concat([out, pd.DataFrame([delta])], ignore_index=True)
    return out


def main(suffix: str) -> None:
    fold_features_csv = OUT_DIR / f"fold_features_D11{suffix}.csv"
    if not fold_features_csv.exists():
        raise SystemExit(f"missing {fold_features_csv} -- run the extraction with out_suffix={suffix!r} first")

    for task in ("regression", "classification"):
        results = run_all(task, fold_features_csv, out_suffix=suffix)
        out_path = OUT_DIR / f"results_{task}{suffix}.csv"
        results.to_csv(out_path, index=False)

        print(f"\n=== {task.upper()}: pre-registered headline "
              f"({HEADLINE['scheme']}, {HEADLINE['tuning']}, {HEADLINE_MODEL[task]}) ===")
        print(compare_headline(task, OUT_DIR / f"results_{task}.csv", results).to_string(index=False))
        print(f"\nSaved {len(results)} rows to {out_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--suffix", default="_qualonly", help="suffix of the fold-features table to run against")
    main(parser.parse_args().suffix)
