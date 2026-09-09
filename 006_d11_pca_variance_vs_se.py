"""
Sanity check for the D11 PCA-coordinate feature (005_d11_pca_features.py):
is the variance of each PC observed across cell lines bigger than the
sampling noise expected from averaging a finite number of cells, or is the
apparent spread just noise?

This is the same check 004_d11_celltype_proportion_se.py does for D11 cell
type proportions, generalized from a binomial proportion to a continuous
per-cell value: a line's PC_k feature is the mean of many single-cell PC_k
values, so its standard error is the standard error of that mean --
SE_i = std(cell-level PC_k in line i) / sqrt(n_i) -- exactly analogous to
SE_i = sqrt(p_i*(1-p_i)/n_i) for a proportion (which is also the SE of a
mean, just of 0/1-valued cells).

Each data point is a cell_line, restricted to the 138 lines in
metadata_eda/qualifying_cell_line_pool_min10_per_timepoint.csv. A line
qualifying through more than one pool has its PC_k value averaged across
those pools (matching 004 and 005's approach), with SE propagated as the
standard error of a mean of independent estimates:

    p_i = mean(pool_mean_k for each qualifying pool k of line i)
    SE_i = sqrt(sum_k SE_k^2) / n_pools_i,   SE_k = std(cell-level PC_k in pool k) / sqrt(n_k)
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

OUT_DIR = Path(__file__).parent / "metadata_eda"
PER_CELL_PCA_CSV = OUT_DIR / "pca/d11_pca_coords_per_cell_qualifying.csv"
QUALIFYING_COMBOS_CSV = OUT_DIR / "cohort/qualifying_cell_line_pool_min10_per_timepoint.csv"
N_PCS = 10
PC_COLS = [f"PC{i}" for i in range(1, N_PCS + 1)]


def compute_pool_level_mean_and_se(per_cell: pd.DataFrame) -> pd.DataFrame:
    """Per (cell_line, pool): n_cells, and for each PC, mean and SE of that mean."""
    grouped = per_cell.groupby(["cell_line", "pool"], observed=True)
    n_cells = grouped.size().rename("n_cells")
    pool_mean = grouped[PC_COLS].mean()
    pool_std = grouped[PC_COLS].std(ddof=1)
    pool_se = pool_std.div(np.sqrt(n_cells), axis=0)

    pool_mean.columns = [f"{c}_mean" for c in pool_mean.columns]
    pool_se.columns = [f"{c}_se" for c in pool_se.columns]
    return pd.concat([n_cells, pool_mean, pool_se], axis=1).reset_index()


def collapse_to_line_level(pool_level: pd.DataFrame) -> pd.DataFrame:
    """Average a line's per-pool PC means across its qualifying pools, with
    SE propagated as the standard error of a mean of independent estimates."""

    def _agg(g: pd.DataFrame) -> pd.Series:
        n_pools = len(g)
        out = {"n_pools": n_pools}
        for pc in PC_COLS:
            out[f"{pc}_mean"] = g[f"{pc}_mean"].mean()
            out[f"{pc}_se"] = np.sqrt((g[f"{pc}_se"] ** 2).sum()) / n_pools
        return pd.Series(out)

    return pool_level.groupby("cell_line", observed=True).apply(_agg, include_groups=False).reset_index()


def _pc_stats(p: pd.Series, se: pd.Series) -> dict:
    observed_var = p.var(ddof=1)
    mean_sampling_var = (se**2).mean()
    variance_ratio = observed_var / mean_sampling_var
    # observed_var already includes averaged sampling noise
    # (observed_var ~= between_line_var + mean_sampling_var), so the
    # debiased between-line variance is observed_var - mean_sampling_var,
    # and ICC = between_line_var / (between_line_var + mean_sampling_var)
    # = 1 - 1/variance_ratio. Bounded 0 (no real signal) to 1 (dominant
    # signal), vs. variance_ratio's unbounded scale where 1 (not 0) is
    # the "no signal" floor.
    icc = 1 - 1 / variance_ratio
    return {
        "n_lines": len(p),
        "mean_value": p.mean(),
        "observed_std": np.sqrt(observed_var),
        "rms_se": np.sqrt(mean_sampling_var),
        "median_se": se.median(),
        "variance_ratio": variance_ratio,
        "icc": icc,
    }


def summarize_by_pc(line_level: pd.DataFrame, pool11_lines: set[str]) -> pd.DataFrame:
    """Per-PC summary, computed both over all qualifying lines and again
    excluding the pool11-derived ones -- pool11 is a severe sequencing-depth
    batch outlier (see 005_d11_pca_features.py), and even after fitting PCA
    on non-pool11 cells, a residual pool11 signature can still show up as
    an artificially high variance_ratio/icc on some PCs: pool11 lines
    differing from everyone else *consistently* (same batch shift on every
    pool11 cell) looks identical, by this test alone, to a PC that
    genuinely separates lines biologically. Comparing the two columns
    reveals how much of each PC's apparent reliability is pool11-driven.
    """
    is_pool11 = line_level["cell_line"].isin(pool11_lines)
    rows = []
    for pc in PC_COLS:
        all_stats = _pc_stats(line_level[f"{pc}_mean"], line_level[f"{pc}_se"])
        excl = line_level.loc[~is_pool11]
        excl_stats = _pc_stats(excl[f"{pc}_mean"], excl[f"{pc}_se"])
        rows.append(
            {
                "PC": pc,
                **all_stats,
                "variance_ratio_excl_pool11": excl_stats["variance_ratio"],
                "icc_excl_pool11": excl_stats["icc"],
            }
        )
    return pd.DataFrame(rows).sort_values("variance_ratio", ascending=False)


def plot_pc_vs_se(line_level: pd.DataFrame) -> None:
    fig, axes = plt.subplots(N_PCS, 1, figsize=(10, 3 * N_PCS), sharex=False)

    for ax, pc in zip(axes, PC_COLS):
        sub = line_level.sort_values(f"{pc}_mean")
        x = np.arange(len(sub))
        ax.errorbar(
            x, sub[f"{pc}_mean"], yerr=sub[f"{pc}_se"], fmt="o", ms=3, elinewidth=1, capsize=0, alpha=0.7
        )
        ax.axhline(sub[f"{pc}_mean"].mean(), color="gray", linestyle="--", linewidth=1)
        observed_var = sub[f"{pc}_mean"].var(ddof=1)
        mean_sampling_var = (sub[f"{pc}_se"] ** 2).mean()
        icc = 1 - mean_sampling_var / observed_var
        ax.set_title(
            f"{pc}: observed_std={np.sqrt(observed_var):.4f}, "
            f"rms_se={np.sqrt(mean_sampling_var):.4f}, icc={icc:.3f}"
        )
        ax.set_ylabel(pc)
        ax.set_xlabel("cell_line (pool-averaged), sorted by value")

    fig.tight_layout()
    fig.savefig(OUT_DIR / "plots/plot_d11_pca_variance_vs_se.png", dpi=150)
    plt.close(fig)


def main() -> None:
    OUT_DIR.mkdir(exist_ok=True)

    per_cell = pd.read_csv(PER_CELL_PCA_CSV)

    pool_level = compute_pool_level_mean_and_se(per_cell)
    line_level = collapse_to_line_level(pool_level)

    n_lines = line_level["cell_line"].nunique()
    n_multi_pool = (line_level["n_pools"] > 1).sum()
    print(f"{n_lines} cell lines; {n_multi_pool} of them average across more than one pool.")

    line_level.to_csv(OUT_DIR / "pca/d11_pca_line_level_with_se.csv", index=False)

    qualifying = pd.read_csv(QUALIFYING_COMBOS_CSV)
    pool11_lines = set(qualifying.loc[qualifying["pool"] == "pool11", "cell_line"])

    summary = summarize_by_pc(line_level, pool11_lines)
    summary.to_csv(OUT_DIR / "pca/d11_pca_variance_vs_se.csv", index=False)

    print("\nPer-PC: observed variation across cell lines vs. typical sampling noise at D11")
    print(summary.to_string(index=False))
    print(
        "\nvariance_ratio = Var(p_i) / mean(SE_i^2) (equivalently, observed_std vs. "
        "rms_se on the PC's own scale). >> 1 means real line-to-line signal "
        "dominates sampling noise; ~1 or below means that PC is mostly noise.\n"
        "icc = 1 - 1/variance_ratio (bounded 0-1 reframing of the same ratio): "
        "the share of observed variance that is real between-line signal "
        "rather than sampling noise. ~0 = mostly noise; -> 1 = strong signal.\n"
        "*_excl_pool11 = the same two stats recomputed with the 16 pool11-derived "
        "lines dropped, since pool11 is a severe sequencing-depth batch outlier "
        "(005_d11_pca_features.py) whose consistent shift can inflate these "
        "reliability metrics without reflecting real biological signal."
    )

    plot_pc_vs_se(line_level)
    print(f"\nSaved CSVs and plot to {OUT_DIR}")


if __name__ == "__main__":
    main()
