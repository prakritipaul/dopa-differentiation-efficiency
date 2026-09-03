"""
How many D11 cells does each cell line have? For lines profiled in more
than one pool, shows both the per-pool breakdown and the total.

Uses metadata_eda/cell_line_donor_pool_timepoint_n_cells.csv (built by
001_eda.py), which already has per-(cell_line, pool, timepoint) counts for
every cell line in the dataset -- not just the 138 qualifying for the
004-006 feature-validity checks. Flags which lines are among those 138 for
cross-reference (a line can have plenty of D11 cells and still not
qualify, since qualifying also requires >=10 cells in the same pool at
D30 and D52).
"""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

OUT_DIR = Path(__file__).parent / "metadata_eda"
ALL_COUNTS_CSV = OUT_DIR / "cell_line_donor_pool_timepoint_n_cells.csv"
QUALIFYING_COMBOS_CSV = OUT_DIR / "qualifying_cell_line_pool_min10_per_timepoint.csv"


def build_table() -> pd.DataFrame:
    counts = pd.read_csv(ALL_COUNTS_CSV)
    d11 = counts.loc[counts["timepoint"] == "D11"]

    def _pool_breakdown(g: pd.DataFrame) -> str:
        pairs = g.sort_values("pool")[["pool", "n_cells"]].itertuples(index=False)
        return ", ".join(f"{pool}:{n_cells}" for pool, n_cells in pairs)

    per_line = (
        d11.groupby("cell_line", observed=True)
        .apply(
            lambda g: pd.Series(
                {
                    "donor": g["donor"].iloc[0],
                    "n_pools_D11": g["pool"].nunique(),
                    "pool_breakdown_D11": _pool_breakdown(g),
                    "total_n_cells_D11": g["n_cells"].sum(),
                }
            ),
            include_groups=False,
        )
        .reset_index()
    )

    qualifying_lines = set(pd.read_csv(QUALIFYING_COMBOS_CSV)["cell_line"])
    per_line["is_qualifying_138"] = per_line["cell_line"].isin(qualifying_lines)

    return per_line.sort_values("total_n_cells_D11", ascending=False).reset_index(drop=True)


def plot_counts(per_line: pd.DataFrame) -> None:
    sub = per_line.sort_values("total_n_cells_D11")
    fig, ax = plt.subplots(figsize=(10, max(4, len(sub) * 0.05)))
    colors = sub["is_qualifying_138"].map({True: "tab:blue", False: "tab:gray"})
    ax.barh(range(len(sub)), sub["total_n_cells_D11"], color=colors)
    ax.set_yticks([])
    ax.set_xlabel("total_n_cells_D11")
    ax.set_title("D11 cells per cell line (blue = among the 138 qualifying lines)")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "plot_d11_cell_counts_per_line.png", dpi=150)
    plt.close(fig)


def plot_histogram(per_line: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.hist(
        [
            per_line.loc[per_line["is_qualifying_138"], "total_n_cells_D11"],
            per_line.loc[~per_line["is_qualifying_138"], "total_n_cells_D11"],
        ],
        bins=30,
        stacked=True,
        color=["tab:blue", "tab:gray"],
        label=["qualifying (138)", "not qualifying"],
    )
    ax.set_xlabel("total_n_cells_D11")
    ax.set_ylabel("number of cell lines")
    ax.set_title("Distribution of D11 cells per cell line")
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT_DIR / "plot_d11_cell_counts_histogram.png", dpi=150)
    plt.close(fig)


def main() -> None:
    OUT_DIR.mkdir(exist_ok=True)

    per_line = build_table()
    per_line.to_csv(OUT_DIR / "d11_cell_counts_per_line.csv", index=False)

    n_multi_pool = (per_line["n_pools_D11"] > 1).sum()
    print(f"{len(per_line)} distinct D11 cell lines; {n_multi_pool} span more than one pool at D11.")
    print(
        per_line["total_n_cells_D11"].describe()[["min", "25%", "50%", "75%", "max"]].to_string()
    )
    print(f"\n{per_line['is_qualifying_138'].sum()} of these are among the 138 qualifying lines.")
    print("\nTop 5 by total D11 cells:")
    print(per_line.head(5).to_string(index=False))
    print("\nExample multi-pool line(s):")
    print(per_line.loc[per_line["n_pools_D11"] > 1].head(5).to_string(index=False))

    plot_counts(per_line)
    plot_histogram(per_line)
    print(f"\nSaved CSV and plots to {OUT_DIR}")


if __name__ == "__main__":
    main()
