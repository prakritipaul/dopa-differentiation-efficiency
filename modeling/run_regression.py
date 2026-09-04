"""
Regression: predict continuous D52 diff_efficiency from D11 features.
See modeling/README.md for the full CV/tuning/correction/grouping design.
"""

from pathlib import Path

from modeling.run_experiment import HEADLINE, run_all, select_headline

OUT_DIR = Path(__file__).parent


def main() -> None:
    results = run_all(task="regression")
    out_path = OUT_DIR / "results_regression.csv"
    results.to_csv(out_path, index=False)

    cols = ["model", "mae_mean", "mae_std", "rmse_mean", "rmse_std", "r2_mean", "r2_std"]
    headline = select_headline(results, task="regression")
    print(f"=== HEADLINE (pre-registered: {HEADLINE['scheme']}, {HEADLINE['tuning']}, ridge) ===")
    print(headline[cols].to_string(index=False))

    other = select_headline(results)
    other = other[~other.index.isin(headline.index)]
    print("\n--- same config, other model family (secondary, NOT the headline) ---")
    print(other[cols].to_string(index=False))

    print(f"\nSaved full robustness grid ({len(results)} rows) to {out_path}")


if __name__ == "__main__":
    main()
