"""
Classification: predict D52 success (diff_efficiency >= 0.2) from D11
features. See modeling/README.md for the full CV/tuning/correction/
grouping design.
"""

from pathlib import Path

from modeling.run_experiment import HEADLINE, run_all, select_headline

OUT_DIR = Path(__file__).parent


def main() -> None:
    results = run_all(task="classification")
    out_path = OUT_DIR / "results/results_classification.csv"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    results.to_csv(out_path, index=False)

    cols = [
        "model", "roc_auc_mean", "roc_auc_std", "pr_auc_mean", "pr_auc_std",
        "balanced_accuracy_mean", "balanced_accuracy_std",
    ]
    headline = select_headline(results, task="classification")
    print(f"=== HEADLINE (pre-registered: {HEADLINE['scheme']}, {HEADLINE['tuning']}, logistic_l2) ===")
    print(headline[cols].to_string(index=False))

    other = select_headline(results)
    other = other[~other.index.isin(headline.index)]
    print("\n--- same config, other model family (secondary, NOT the headline) ---")
    print(other[cols].to_string(index=False))

    print(f"\nSaved full robustness grid ({len(results)} rows) to {out_path}")


if __name__ == "__main__":
    main()
