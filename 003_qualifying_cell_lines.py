"""
Identify cell lines with at least MIN_CELLS_PER_TIMEPOINT cells at every
timepoint (D11, D30, D52) in the Jerber et al. dataset
(https://pmc.ncbi.nlm.nih.gov/articles/PMC7610897/), and record which
pool each was profiled in.

Filtering matches the authors' own analysis (plotting_notebooks/Figure_2/
fig2b_and_heatmap_extended.ipynb in
https://github.com/single-cell-genetics/singlecell_neuroseq_paper): the
unit of analysis is a (cell_line, pool) combination, not a cell_line alone.
For each timepoint, count cells per (cell_line, pool) and keep combos with
>= MIN_CELLS_PER_TIMEPOINT cells; a combo only survives overall if the same
(cell_line, pool) pair clears the threshold at all three timepoints (their
code does this via `pd.concat(...).dropna()` on the three per-timepoint
frames indexed by (donor_id, pool_id)). Summing a cell line's cells across
*different* pools before filtering (as an earlier version of this script
did) is looser and over-counts: it found 140 lines vs. the paper's reported
138.
"""

import re
from pathlib import Path

import h5py
import pandas as pd

DATA_FILES = [
    "/Users/prakritipaul/Documents/2021_jerber/day11.h5",
    "/Users/prakritipaul/Documents/2021_jerber/day30.h5",
    "/Users/prakritipaul/Documents/2021_jerber/day52.h5",
]

DONOR_RE = re.compile(r"^(HPSI\d+i)-")
OUT_DIR = Path(__file__).parent / "metadata_eda"
TIMEPOINTS = ["D11", "D30", "D52"]
MIN_CELLS_PER_TIMEPOINT = 10


def read_obs_categorical(f: h5py.File, column: str) -> pd.Categorical:
    """Read a legacy-anndata categorical obs column (codes + __categories)."""
    codes = f[f"obs/{column}"][:]
    categories = [c.decode() for c in f[f"obs/__categories/{column}"][:]]
    return pd.Categorical.from_codes(codes, categories=categories)


def load_obs(path: str) -> pd.DataFrame:
    with h5py.File(path, "r") as f:
        cell_line = read_obs_categorical(f, "donor_id")
        pool = read_obs_categorical(f, "pool_id")
        timepoint = read_obs_categorical(f, "time_point")
        # Only D52 has more than one level ('NONE', 'ROT'); D11 and D30 are
        # entirely untreated. Carried for every timepoint anyway so the filter
        # in main() is one uniform rule rather than a D52 special case.
        treatment = read_obs_categorical(f, "treatment")

    donor = [DONOR_RE.match(cl).group(1) for cl in cell_line]

    return pd.DataFrame(
        {
            "cell_line": cell_line,
            "donor": donor,
            "pool": pool,
            "timepoint": timepoint,
            "treatment": treatment,
        }
    )


def main(untreated_only: bool = False, out_suffix: str = "") -> None:
    """untreated_only: drop rotenone-treated cells before counting.

    This is a DELIBERATE DIVERGENCE from the published analysis, not a bug
    fix. The authors' own notebook (plotting_notebooks/Figure_2/
    fig2b_and_heatmap_extended.ipynb) applies only a >=10 cell threshold and
    never filters on treatment, so the default (False) reproduces them
    exactly. But D52 -- and only D52 -- contains rotenone-treated cells
    (303,856 NONE / 219,238 ROT), interleaved so that every qualifying
    (cell_line, pool) combo contains both. Rotenone is a mitochondrial
    complex I inhibitor that preferentially damages dopaminergic neurons,
    i.e. exactly the cell type in the numerator of `diff_efficiency`, so
    including treated cells mixes differentiation success with toxin
    sensitivity.

    Excluding them costs two lines whose untreated D52 count falls below the
    threshold (HPSI0115i-melw_1 and HPSI0115i-qecv_2, both pool5, 14 cells
    -> 6 and 8): 138 lines becomes 136, 159 combos becomes 157, donors stay
    at 20.

    Always pair with a distinct out_suffix. The unsuffixed outputs are the
    provenance for every result already published from them, and must not be
    overwritten by a run under different rules."""
    obs = pd.concat([load_obs(path) for path in DATA_FILES], ignore_index=True)
    if untreated_only:
        before = len(obs)
        obs = obs[obs["treatment"] == "NONE"]
        print(f"untreated_only: kept {len(obs)} of {before} cells "
              f"({before - len(obs)} rotenone-treated dropped)")

    per_combo_counts = (
        obs.groupby(["cell_line", "donor", "pool", "timepoint"], observed=True)
        .size()
        .reset_index(name="n_cells")
    )
    counts_wide = per_combo_counts.pivot(
        index=["cell_line", "donor", "pool"], columns="timepoint", values="n_cells"
    )
    for tp in TIMEPOINTS:
        if tp not in counts_wide.columns:
            counts_wide[tp] = 0
    counts_wide = counts_wide[TIMEPOINTS]

    # A (cell_line, pool) combo must clear the threshold at every timepoint
    # (missing entirely from a timepoint -> NaN -> fails `>=`, matching the
    # authors' dropna() after concatenating the three per-timepoint frames).
    qualifying_mask = (counts_wide >= MIN_CELLS_PER_TIMEPOINT).all(axis=1)
    result = counts_wide.loc[qualifying_mask].reset_index()
    result = result.rename(columns={tp: f"n_cells_{tp}" for tp in TIMEPOINTS})
    result = result.sort_values(["cell_line", "pool"]).reset_index(drop=True)

    OUT_DIR.mkdir(exist_ok=True)
    out_path = OUT_DIR / f"cohort/qualifying_cell_line_pool_min10_per_timepoint{out_suffix}.csv"
    result.to_csv(out_path, index=False)

    n_distinct_lines = result["cell_line"].nunique()

    print(result)
    print(
        f"\n{len(result)} (cell_line, pool) combos have >= {MIN_CELLS_PER_TIMEPOINT} "
        f"cells in that same pool at every timepoint ({', '.join(TIMEPOINTS)})."
    )
    print(f"{n_distinct_lines} distinct cell lines among those combos.")
    print(f"Saved to {out_path}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Build the qualifying (cell_line, pool) cohort.")
    parser.add_argument("--untreated-only", action="store_true",
                        help="drop rotenone-treated cells before counting (D52 only has any)")
    parser.add_argument("--suffix", default="", help="suffix for the output CSV; required with --untreated-only")
    args = parser.parse_args()
    if args.untreated_only and not args.suffix:
        parser.error("--untreated-only changes the cohort, so it needs a --suffix "
                     "rather than overwriting the published cohort file")
    main(untreated_only=args.untreated_only, out_suffix=args.suffix)
