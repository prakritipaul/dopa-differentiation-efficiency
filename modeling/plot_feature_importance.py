"""
Importance against run association: mean |SHAP| vs pool eta^2, per timepoint.

One scatter per timepoint, from the feature-importance table that
`feature_importance.py` already writes. Horizontal axis is how much the
fitted model leans on a feature; vertical axis is how much of that
feature's across-line variance goes with which differentiation run the
line went through. Plotting them against each other separates two
questions that are easily conflated -- a feature can be top-right
(important and run-associated, e.g. D11 PC3) or bottom-right (important
and run-independent, e.g. D11 phat_NB).

No threshold line is drawn: eta^2 runs continuously, so a cutoff would
assert a distinction the data does not support.

The reference proportion has no SHAP value (it never enters a fit) and is
dropped rather than plotted at zero.

Usage:
    uv run python modeling/plot_feature_importance.py
"""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator

RESULTS_DIR = Path(__file__).parent / "results"
OUT_DIR = Path(__file__).parent / "plots"

TIMEPOINTS = ("D11", "D30")
LABEL_VARIANT = "qualonly_da_untreated"

COMPOSITION_COLOR = "#6a9c78"
EXPRESSION_COLOR = "#2f6b5e"
INK = "#3d4441"
MUTED = "#8a8f8c"
GRID = "#d8dcd9"

FIG_SIZE = (11, 6.5)
LABEL_FONTSIZE = 14
DPI = 300

# Font stacks, not single names: matplotlib walks the list and falls back,
# so these render on a machine with none of the preferred faces installed.
# Optima and PT Sans are deliberately absent -- both lack Greek, and the
# axis labels carry a literal eta and rho.
TEXT_FONTS = ["Avenir Next", "Helvetica Neue", "Charter", "DejaVu Sans"]
TICK_FONTS = ["Menlo", "DejaVu Sans Mono", "monospace"]
# Drawable area of the axes in points, after tight_layout. Only used to
# estimate label boxes for collision avoidance, so approximate is fine.
AXES_W_PTS, AXES_H_PTS = 700, 355

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


def _place_labels(df: pd.DataFrame, xlim, ylim) -> list[tuple[int, int, str]]:
    """Offset and alignment for each label, chosen to avoid overlap.

    Labels default to the right of their point. Where that would collide
    with another label or run off the axes -- the D30 plot has a dense
    cluster at low |SHAP| and low eta^2 -- successive fallbacks are tried:
    left, then above and below on either side.

    Boxes are estimated from the character count rather than measured from
    a renderer. Approximate, but enough to separate the dozen labels these
    two plots carry.
    """
    x_pts = AXES_W_PTS / (xlim[1] - xlim[0])
    y_pts = AXES_H_PTS / (ylim[1] - ylim[0])
    candidates = [
        (14, -5, "left"),
        (-14, -5, "right"),
        (14, 13, "left"),
        (-14, 13, "right"),
        (14, -22, "left"),
        (-14, -22, "right"),
    ]

    # Markers are obstacles too: a label flipped to the left would
    # otherwise be free to land on its neighbour's dot.
    marker_r = 9
    placed: list[tuple[float, float, float, float]] = [
        (
            r["shap_mean_abs"] * x_pts - marker_r,
            r["technical_covariate_eta2"] * y_pts - marker_r,
            r["shap_mean_abs"] * x_pts + marker_r,
            r["technical_covariate_eta2"] * y_pts + marker_r,
        )
        for _, r in df.iterrows()
    ]
    chosen_by_index = {}
    # Left to right, so a crowded label claims its slot before its
    # right-hand neighbour does.
    for i in df["shap_mean_abs"].sort_values().index:
        row = df.loc[i]
        x = row["shap_mean_abs"] * x_pts
        y = row["technical_covariate_eta2"] * y_pts
        # Bold labels (the low-eta^2 ones) set wider than regular ones.
        char_w = 0.78 if row["technical_covariate_eta2"] < LOW_ETA2 else 0.70
        w = len(row["feature"]) * LABEL_FONTSIZE * char_w
        h = LABEL_FONTSIZE * 1.3
        chosen = candidates[0]
        for dx, dy, ha in candidates:
            left = x + dx - (w if ha == "right" else 0)
            box = (left, y + dy - h / 2, left + w, y + dy + h / 2)
            if box[0] < xlim[0] * x_pts or box[2] > xlim[1] * x_pts:
                continue
            if any(
                box[0] < p[2] and p[0] < box[2] and box[1] < p[3] and p[1] < box[3]
                for p in placed
            ):
                continue
            chosen = (dx, dy, ha)
            placed.append(box)
            break
        chosen_by_index[i] = chosen
    return [chosen_by_index[i] for i in df.index]


def shared_xmax() -> float:
    """Largest |SHAP| across both timepoints.

    Both panels are drawn on this one scale so horizontal positions can be
    compared between them, which is what FINDINGS.md claims of the pair.
    Per-panel limits silently broke that: D11 reached 1.40 and D30 1.29, so
    the same distance meant different things in each.
    """
    return max(load_importance(tp)["shap_mean_abs"].max() for tp in TIMEPOINTS)


def plot_timepoint(timepoint: str, xmax: float | None = None) -> Path:
    with plt.rc_context({"font.family": "sans-serif", "font.sans-serif": TEXT_FONTS}):
        return _plot(timepoint, xmax)


def _plot(timepoint: str, xmax: float | None) -> Path:
    df = load_importance(timepoint)
    xmax = shared_xmax() if xmax is None else xmax
    is_pc = df["feature"].str.startswith("PC")

    fig, ax = plt.subplots(figsize=FIG_SIZE)
    for mask, color in ((~is_pc, COMPOSITION_COLOR), (is_pc, EXPRESSION_COLOR)):
        ax.scatter(
            df.loc[mask, "shap_mean_abs"],
            df.loc[mask, "technical_covariate_eta2"],
            s=200,
            color=color,
            edgecolors="white",
            linewidths=1.4,
            zorder=3,
        )

    xlim = (-0.04 * xmax, xmax * 1.12)
    ylim = (-0.06, 1.04)
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)

    for (_, row), (dx, dy, ha) in zip(df.iterrows(), _place_labels(df, xlim, ylim)):
        ax.annotate(
            row["feature"],
            (row["shap_mean_abs"], row["technical_covariate_eta2"]),
            textcoords="offset points",
            xytext=(dx, dy),
            ha=ha,
            fontsize=LABEL_FONTSIZE,
            color=EXPRESSION_COLOR if row["feature"].startswith("PC") else COMPOSITION_COLOR,
            fontweight="bold" if row["technical_covariate_eta2"] < LOW_ETA2 else "normal",
        )

    ax.set_xlabel("Mean |SHAP| — contribution to the prediction", fontsize=15, color=INK, labelpad=14)
    ax.set_ylabel("Pool \u03b7\u00b2 — association with the run", fontsize=15, color=INK, labelpad=14)
    ax.set_yticks([0.0, 0.25, 0.5, 0.75, 1.0])
    ax.xaxis.set_major_locator(MaxNLocator(4, steps=[1, 5, 10]))
    ax.grid(axis="y", color=GRID, linewidth=1.1, zorder=0)
    ax.set_axisbelow(True)
    for side in ax.spines:
        ax.spines[side].set_visible(False)
    ax.tick_params(length=0, pad=10, labelsize=13, colors=MUTED)
    for tick in ax.get_xticklabels() + ax.get_yticklabels():
        tick.set_fontfamily(TICK_FONTS)

    # Short rules rather than dots in the legend, as in the reference figure.
    ax.legend(
        handles=[
            Line2D([], [], color=COMPOSITION_COLOR, lw=3.2, label="Composition feature"),
            Line2D([], [], color=EXPRESSION_COLOR, lw=3.2, label="Expression component"),
        ],
        frameon=False,
        loc="upper left",
        bbox_to_anchor=(-0.02, 1.16),
        ncol=2,
        fontsize=13,
        handlelength=1.2,
        handletextpad=0.6,
        columnspacing=2.4,
        labelcolor=INK,
    )

    fig.tight_layout()
    out_path = OUT_DIR / f"plot_feature_importance_{timepoint}_da_untreated.png"
    fig.savefig(out_path, dpi=DPI)
    plt.close(fig)
    return out_path


if __name__ == "__main__":
    xmax = shared_xmax()
    for tp in TIMEPOINTS:
        print(f"wrote {plot_timepoint(tp, xmax)}")
