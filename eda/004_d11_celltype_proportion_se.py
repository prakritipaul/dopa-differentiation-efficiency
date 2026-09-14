"""
Sanity check for a candidate D11 -> D52 predictive feature: per-cell_line
cell type proportions at D11.

Each data point is a cell_line, restricted to the 138 lines in
metadata_eda/qualifying_cell_line_pool_min10_per_timepoint.csv (the
paper-matching set: >= 10 cells in the same pool at every timepoint). A
line qualifying through more than one pool has its D11 cell type
proportion averaged across those pools (the authors' approach), not pooled
by summing raw counts:

    p_i = mean(phat_k for each qualifying pool k of line i)

Averaging proportions this way -- rather than summing counts -- means the
standard error of p_i is the SE of a mean of independent estimates, not a
single bigger-sample binomial SE:

    SE_i = sqrt(sum_k SE_k^2) / n_pools_i,   SE_k = sqrt(phat_k*(1-phat_k)/n_k)

(for a line with only one qualifying pool, n_pools_i=1 and this reduces to
the plain binomial SE of that one pool). If the variance of p_i observed
across cell lines is no bigger than this expected sampling variance, that
proportion is not a useful feature: the apparent spread would be explained
by counting noise alone, not by real line-to-line differences.

Per celltype, this script computes:
  - s^2_observed = Var(p_1, p_2, ...) across cell lines (step 2)
  - SE_i per cell line (step 3)
  - variance_ratio = s^2_observed / mean(SE_i^2) (step 4), and the same
    comparison on the proportion scale: observed_std vs. rms_se
    (= sqrt(mean(SE_i^2)), not the plain mean of SE_i)
"""

import re
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

DAY11_FILE = "/Users/prakritipaul/Documents/2021_jerber/day11.h5"
DONOR_RE = re.compile(r"^(HPSI\d+i)-")
OUT_DIR = Path(__file__).parent.parent / "metadata_eda"
QUALIFYING_COMBOS_CSV = OUT_DIR / "cohort/qualifying_cell_line_pool_min10_per_timepoint.csv"


def read_obs_categorical(f: h5py.File, column: str) -> pd.Categorical:
    """Read a legacy-anndata categorical obs column (codes + __categories)."""
    codes = f[f"obs/{column}"][:]
    categories = [c.decode() for c in f[f"obs/__categories/{column}"][:]]
    return pd.Categorical.from_codes(codes, categories=categories)


def load_d11_obs() -> pd.DataFrame:
    with h5py.File(DAY11_FILE, "r") as f:
        cell_line = read_obs_categorical(f, "donor_id")
        pool = read_obs_categorical(f, "pool_id")
        celltype = read_obs_categorical(f, "celltype")

    return pd.DataFrame({"cell_line": cell_line, "pool": pool, "celltype": celltype})


def compute_pool_level_proportions_and_se(obs: pd.DataFrame) -> pd.DataFrame:
    """Per (cell_line, pool, celltype): n_cells, n_total, phat, se."""
    counts = (
        obs.groupby(["cell_line", "pool", "celltype"], observed=True)
        .size()
        .reset_index(name="n_cells")
    )
    n_total = obs.groupby(["cell_line", "pool"], observed=True).size().rename("n_total")

    # fill in explicit zeros for celltypes a combo simply didn't produce
    counts = (
        counts.set_index(["cell_line", "pool", "celltype"])["n_cells"]
        .unstack("celltype", fill_value=0)
        .stack()
        .rename("n_cells")
        .reset_index()
        .merge(n_total.reset_index(), on=["cell_line", "pool"])
    )

    counts["phat"] = counts["n_cells"] / counts["n_total"]
    counts["se"] = np.sqrt(counts["phat"] * (1 - counts["phat"]) / counts["n_total"])
    return counts


def collapse_to_line_level(pool_level: pd.DataFrame) -> pd.DataFrame:
    """Average a line's per-pool phat across its qualifying pools, with SE
    propagated as the standard error of a mean of independent estimates."""

    def _agg(g: pd.DataFrame) -> pd.Series:
        n_pools = len(g)
        return pd.Series(
            {
                "n_pools": n_pools,
                "phat": g["phat"].mean(),
                "se": np.sqrt((g["se"] ** 2).sum()) / n_pools,
            }
        )

    return (
        pool_level.groupby(["cell_line", "celltype"], observed=True)
        .apply(_agg, include_groups=False)
        .reset_index()
    )


def summarize_by_celltype(proportions: pd.DataFrame) -> pd.DataFrame:
    def _agg(g: pd.DataFrame) -> pd.Series:
        # Step 2: observed variance of p_i across cell lines.
        observed_var = g["phat"].var(ddof=1)
        # Step 3/4: mean(SE_i^2) is the expected sampling variance; its
        # sqrt (RMS of SE_i, not the plain mean of SE_i) is the quantity
        # directly comparable to SD_observed on the proportion scale.
        mean_sampling_var = (g["se"] ** 2).mean()
        return pd.Series(
            {
                "n_lines": len(g),
                "mean_phat": g["phat"].mean(),
                "observed_std": np.sqrt(observed_var),
                "rms_se": np.sqrt(mean_sampling_var),
                "median_se": g["se"].median(),
                "variance_ratio": observed_var / mean_sampling_var,
            }
        )

    summary = proportions.groupby("celltype", observed=True).apply(_agg, include_groups=False)
    return summary.sort_values("variance_ratio", ascending=False)


def plot_proportion_vs_se(proportions: pd.DataFrame) -> None:
    celltypes = sorted(proportions["celltype"].unique())
    fig, axes = plt.subplots(len(celltypes), 1, figsize=(10, 3.5 * len(celltypes)), sharex=False)
    if len(celltypes) == 1:
        axes = [axes]

    for ax, ct in zip(axes, celltypes):
        sub = proportions.loc[proportions["celltype"] == ct].sort_values("phat")
        x = np.arange(len(sub))
        ax.errorbar(
            x, sub["phat"], yerr=sub["se"], fmt="o", ms=3, elinewidth=1, capsize=0, alpha=0.7
        )
        ax.axhline(sub["phat"].mean(), color="gray", linestyle="--", linewidth=1)
        rms_se = np.sqrt((sub["se"] ** 2).mean())
        ax.set_title(
            f"{ct}: observed_std={sub['phat'].std(ddof=1):.4f}, "
            f"rms_se={rms_se:.4f}"
        )
        ax.set_ylabel("D11 proportion (p_i)")
        ax.set_xlabel("cell_line (pool-averaged), sorted by proportion")

    fig.tight_layout()
    fig.savefig(OUT_DIR / "plots/plot_d11_celltype_proportion_vs_se.png", dpi=150)
    plt.close(fig)


def main() -> None:
    OUT_DIR.mkdir(exist_ok=True)

    qualifying = pd.read_csv(QUALIFYING_COMBOS_CSV)[["cell_line", "pool"]]
    obs = load_d11_obs()
    obs = obs.merge(qualifying, on=["cell_line", "pool"], how="inner")

    pool_level = compute_pool_level_proportions_and_se(obs)
    proportions = collapse_to_line_level(pool_level)

    n_lines = proportions["cell_line"].nunique()
    n_multi_pool = (proportions.groupby("cell_line")["n_pools"].first() > 1).sum()
    print(
        f"{n_lines} cell lines (from {len(qualifying)} qualifying combos); "
        f"{n_multi_pool} of them average across more than one pool."
    )

    proportions.to_csv(OUT_DIR / "proportions/d11_celltype_proportions_with_se.csv", index=False)

    summary = summarize_by_celltype(proportions)
    summary.to_csv(OUT_DIR / "proportions/d11_celltype_proportion_variance_vs_se.csv")

    print("\nPer-celltype: observed variation across cell lines vs. typical sampling noise at D11")
    print(summary)
    print(
        "\nvariance_ratio = s^2_observed / mean(SE_i^2) (equivalently, "
        "observed_std vs. rms_se on the proportion scale). >> 1 means real "
        "line-to-line signal dominates counting noise; ~1 or below means "
        "the proportion is mostly noise for this celltype."
    )

    plot_proportion_vs_se(proportions)
    print(f"\nSaved CSVs and plot to {OUT_DIR}")


if __name__ == "__main__":
    main()
