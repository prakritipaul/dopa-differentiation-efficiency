"""
D52 outcome label: differentiation efficiency = fraction of a cell line's
D52 cells that are DA (dopaminergic) or Sert (serotonergic) neurons -- the
two mature, "successfully differentiated" neuron types, per the paper's
own analysis (plotting_notebooks/Figure_2/fig2b_and_heatmap_extended.ipynb
in https://github.com/single-cell-genetics/singlecell_neuroseq_paper:
`diff_efficiency = df2[['DA_D52','Sert_D52']].sum(axis=1, skipna=False)`).

Computed with the exact same methodology as the D11 features (004/005/006):
per (cell_line, pool) binomial proportion + SE, restricted to the 138
qualifying lines, then two-stage pool-then-line averaging with SE
propagation for lines with multiple qualifying pools. Reuses
004_d11_celltype_proportion_se.py's compute_pool_level_proportions_and_se
and collapse_to_line_level directly (via importlib, same reuse pattern as
008_technical_covariate_associations.py) on a re-labeled celltype column
(DA/Sert -> "differentiated", everything else -> "other"), so DA and Sert
counts are combined into a single binomial indicator *before* computing
SE -- correct, since DA and Sert are mutually exclusive outcomes of the
same per-cell draw, not independent proportions whose SEs should be
propagated as if summing two separate estimates.
"""

import importlib.util
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

DAY52_FILE = "/Users/prakritipaul/Documents/2021_jerber/day52.h5"
OUT_DIR = Path(__file__).parent / "metadata_eda"
QUALIFYING_COMBOS_CSV = OUT_DIR / "cohort/qualifying_cell_line_pool_min10_per_timepoint.csv"
DIFFERENTIATED_CELLTYPES = {"DA", "Sert"}


def _load_module(name: str):
    path = Path(__file__).parent / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


m004 = _load_module("004_d11_celltype_proportion_se")


def load_d52_obs(untreated_only: bool = False) -> pd.DataFrame:
    """untreated_only: drop rotenone-treated cells. See main()'s docstring --
    this is a deliberate divergence from the published definition, and the
    default reproduces it."""
    with h5py.File(DAY52_FILE, "r") as f:
        cell_line = m004.read_obs_categorical(f, "donor_id")
        pool = m004.read_obs_categorical(f, "pool_id")
        celltype = m004.read_obs_categorical(f, "celltype")
        treatment = m004.read_obs_categorical(f, "treatment")

    is_differentiated = pd.Series(celltype).isin(DIFFERENTIATED_CELLTYPES)
    label = pd.Categorical(np.where(is_differentiated, "differentiated", "other"))

    obs = pd.DataFrame({"cell_line": cell_line, "pool": pool, "celltype": label,
                        "treatment": treatment})
    if untreated_only:
        before = len(obs)
        obs = obs[obs["treatment"] == "NONE"]
        print(f"untreated_only: kept {len(obs)} of {before} D52 cells "
              f"({before - len(obs)} rotenone-treated dropped)")
    return obs.drop(columns="treatment")


def plot_label(label: pd.DataFrame, out_suffix: str = "") -> None:
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    sub = label.sort_values("diff_efficiency")
    x = np.arange(len(sub))
    axes[0].errorbar(
        x, sub["diff_efficiency"], yerr=sub["diff_efficiency_se"],
        fmt="o", ms=3, elinewidth=1, capsize=0, alpha=0.7,
    )
    axes[0].set_xlabel("cell_line, sorted by diff_efficiency")
    axes[0].set_ylabel("diff_efficiency (DA + Sert fraction at D52)")
    axes[0].set_title("D52 differentiation efficiency per cell line")

    axes[1].hist(label["diff_efficiency"], bins=20)
    axes[1].set_xlabel("diff_efficiency")
    axes[1].set_ylabel("number of cell lines")
    axes[1].set_title("Distribution of D52 diff_efficiency")

    fig.tight_layout()
    fig.savefig(OUT_DIR / f"plots/plot_d52_diff_efficiency{out_suffix}.png", dpi=150)
    plt.close(fig)


def main(untreated_only: bool = False, out_suffix: str = "") -> None:
    """untreated_only: compute the label from untreated D52 cells only.

    A DELIBERATE DIVERGENCE from the published definition. The authors'
    notebook computes `diff_efficiency = df2[['DA_D52','Sert_D52']].sum(...)`
    with no treatment filter at all, so the default (False) reproduces them.
    But D52 contains 219,238 rotenone-treated cells against 303,856
    untreated, interleaved such that ALL 159 qualifying (cell_line, pool)
    combos contain both -- so an unfiltered per-combo fraction necessarily
    blends them. Rotenone inhibits mitochondrial complex I and
    preferentially damages dopaminergic neurons, the numerator of this very
    metric, so the published label conflates differentiation success with
    toxin sensitivity. Per-line effect is a mean shift of 0.025 and a max of
    0.405.

    IMPORTANT: out_suffix selects the COHORT file that is read as well as the
    label file that is written. The two must match -- a label computed on
    untreated cells against the treated-inclusive cohort would score lines
    that the untreated cohort excludes, which is the half-threaded variant
    flag this repo has already been bitten by."""
    OUT_DIR.mkdir(exist_ok=True)

    combos_csv = OUT_DIR / f"cohort/qualifying_cell_line_pool_min10_per_timepoint{out_suffix}.csv"
    if not combos_csv.exists():
        raise SystemExit(f"missing {combos_csv} -- run 003_qualifying_cell_lines.py "
                         f"with a matching --suffix first")
    print(f"cohort: {combos_csv.name}")
    qualifying = pd.read_csv(combos_csv)[["cell_line", "pool"]]
    obs = load_d52_obs(untreated_only).merge(qualifying, on=["cell_line", "pool"], how="inner")

    pool_level = m004.compute_pool_level_proportions_and_se(obs)
    line_level = m004.collapse_to_line_level(pool_level)

    label = line_level.loc[line_level["celltype"] == "differentiated", ["cell_line", "n_pools", "phat", "se"]]
    label = label.rename(columns={"phat": "diff_efficiency", "se": "diff_efficiency_se"})
    label = label.sort_values("cell_line").reset_index(drop=True)

    label.to_csv(OUT_DIR / f"cohort/d52_diff_efficiency_label{out_suffix}.csv", index=False)

    print(f"{len(label)} cell lines with a D52 diff_efficiency label.")
    print(label["diff_efficiency"].describe().to_string())
    print("\nHighest diff_efficiency:")
    print(label.sort_values("diff_efficiency", ascending=False).head(5).to_string(index=False))
    print("\nLowest diff_efficiency:")
    print(label.sort_values("diff_efficiency").head(5).to_string(index=False))

    plot_label(label, out_suffix)
    print(f"\nSaved CSV and plot to {OUT_DIR}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Build the D52 diff_efficiency label.")
    parser.add_argument("--untreated-only", action="store_true",
                        help="compute the label from untreated D52 cells only")
    parser.add_argument("--suffix", default="",
                        help="suffix selecting the cohort read AND the label written; "
                             "required with --untreated-only")
    args = parser.parse_args()
    if args.untreated_only and not args.suffix:
        parser.error("--untreated-only changes the label, so it needs a --suffix "
                     "rather than overwriting the published label file")
    main(untreated_only=args.untreated_only, out_suffix=args.suffix)
