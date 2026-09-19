"""
Importance against run association: mean |SHAP| vs pool eta^2, per timepoint.

One scatter per timepoint, from the feature-importance table that
`feature_importance.py` already writes. Horizontal axis is how much the
fitted model leans on a feature; vertical axis is how much of that
feature's across-line variance goes with which differentiation run the
line went through. The two are independent questions, which is the whole
point of plotting them against each other -- a feature can be top-right
(important and run-associated, e.g. D11 PC3) or bottom-right (important
and run-independent, e.g. D11 phat_NB).

No threshold line is drawn. eta^2 runs continuously and no value of it
has been validated as a transfer test, so a line across the plot would
assert a cutoff the data does not support.

The reference proportion has no SHAP value (it never enters a fit) and is
dropped rather than plotted at zero.

Usage:
    uv run python modeling/plot_feature_importance.py
"""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

RESULTS_DIR = Path(__file__).parent / "results"
OUT_DIR = Path(__file__).parent / "plots"

TIMEPOINTS = ("D11", "D30")
LABEL_VARIANT = "qualonly_da_untreated"

COMPOSITION_COLOR = "#6a9c78"
EXPRESSION_COLOR = "#2f6b5e"

# Below this, a feature's across-line variance is mostly not the run.
# Used only to bold a label, never to filter or bin.
LOW_ETA2 = 0.15


def load_importance(timepoint: str) -> pd.DataFrame:
    """Plottable rows of the classification importance table.

    Drops the reference proportion, which carries no SHAP value.
    """
    path = RESULTS_DIR / f"feature_importance_table_classification_logistic_l2_{timepoint}_{LABEL_VARIANT}.csv"
    df = pd.read_csv(path)
    return df.dropna(subset=["shap_mean_abs"]).reset_index(drop=True)


def _label_side(df: pd.DataFrame) -> list[bool]:
    """True where a label goes left of its point instead of right.

    A label sits to the right by default, which collides when another
    point sits just to the right at nearly the same height (D30 PC5/PC1).
    Flipping the left member of such a pair separates them.
    """
    x_span = df["shap_mean_abs"].max() - df["shap_mean_abs"].min()
    xs = df["shap_mean_abs"].to_numpy()
    ys = df["technical_covariate_eta2"].to_numpy()
    return [
        any(
            0 < xs[j] - xs[i] < 0.10 * x_span and abs(ys[j] - ys[i]) < 0.04
            for j in range(len(xs))
        )
        for i in range(len(xs))
    ]


def plot_timepoint(timepoint: str) -> Path:
    df = load_importance(timepoint)
    is_pc = df["feature"].str.startswith("PC")

    fig, ax = plt.subplots(figsize=(7, 5))
    for mask, color, label in (
        (~is_pc, COMPOSITION_COLOR, "Composition feature"),
        (is_pc, EXPRESSION_COLOR, "Expression component"),
    ):
        ax.scatter(
            df.loc[mask, "shap_mean_abs"],
            df.loc[mask, "technical_covariate_eta2"],
            s=70,
            color=color,
            label=label,
            zorder=3,
        )

    for (_, row), left in zip(df.iterrows(), _label_side(df)):
        ax.annotate(
            row["feature"],
            (row["shap_mean_abs"], row["technical_covariate_eta2"]),
            textcoords="offset points",
            xytext=(-9 if left else 9, -4),
            ha="right" if left else "left",
            fontsize=9,
            color=EXPRESSION_COLOR if row["feature"].startswith("PC") else COMPOSITION_COLOR,
            fontweight="bold" if row["technical_covariate_eta2"] < LOW_ETA2 else "normal",
        )

    ax.set_xlabel("Mean |SHAP| — contribution to the prediction")
    ax.set_ylabel("Pool $\\eta^2$ — association with the run")
    ax.set_ylim(-0.05, 1.0)
    # Headroom on the right so the rightmost label is not clipped.
    ax.set_xlim(0, df["shap_mean_abs"].max() * 1.25)
    ax.grid(axis="y", color="#cccccc", linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(0, 1.12), ncol=2)

    fig.tight_layout()
    out_path = OUT_DIR / f"plot_feature_importance_{timepoint}_da_untreated.png"
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path


if __name__ == "__main__":
    for tp in TIMEPOINTS:
        print(f"wrote {plot_timepoint(tp)}")
