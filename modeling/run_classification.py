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
    out_path = OUT_DIR / "results_classification.csv"
    results.to_csv(out_path, index=False)

    headline = select_headline(results)
    print(f"=== HEADLINE ({HEADLINE['scheme']}, {HEADLINE['tuning']}, pool_correction={HEADLINE['pool_correction']}) ===")
    print(
        headline[
            ["model", "roc_auc_mean", "roc_auc_std", "pr_auc_mean", "pr_auc_std", "balanced_accuracy_mean", "balanced_accuracy_std"]
        ].to_string(index=False)
    )

    print(f"\nSaved full robustness grid ({len(results)} rows) to {out_path}")


if __name__ == "__main__":
    main()
