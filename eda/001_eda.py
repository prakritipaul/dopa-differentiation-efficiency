"""
Exploratory Data Analysis (EDA) for the Jerber et al. dopaminergic neuron
differentiation dataset (https://pmc.ncbi.nlm.nih.gov/articles/PMC7610897/).

Builds a cell_line x donor x pool x timepoint x n_cells summary table from
the per-timepoint AnnData .h5 files. Only `obs` metadata is read (the files
are 10s of GB, so the expression matrix is never loaded).

Sample naming: obs["donor_id"] holds the full iPSC line id, e.g.
"HPSI0114i-eipl_1". The "HPSI0114i" prefix identifies the donor; the suffix
("eipl_1") identifies the specific clonal line derived from that donor. A
single donor can contribute multiple lines, so cell_line and donor are
tracked as separate columns.
"""

import re
from pathlib import Path

import h5py
import pandas as pd

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from data_paths import data_file  # noqa: E402  (repo root, added above)

DATA_FILES = [
    data_file("D11"),
    data_file("D30"),
    data_file("D52"),
]

DONOR_RE = re.compile(r"^(HPSI\d+i)-")


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
        celltype = read_obs_categorical(f, "celltype")

    donor = [DONOR_RE.match(cl).group(1) for cl in cell_line]

    return pd.DataFrame(
        {
            "cell_line": cell_line,
            "donor": donor,
            "pool": pool,
            "timepoint": timepoint,
            "celltype": celltype,
        }
    )


def make_table(obs: pd.DataFrame, group_cols: list[str]) -> pd.DataFrame:
    return (
        obs.groupby(group_cols, observed=True)
        .size()
        .reset_index(name="n_cells")
        .sort_values(group_cols)
        .reset_index(drop=True)
    )


OUT_DIR = Path(__file__).parent.parent / "metadata_eda"


def main() -> None:
    obs = pd.concat([load_obs(path) for path in DATA_FILES], ignore_index=True)
    OUT_DIR.mkdir(exist_ok=True)

    table = make_table(obs, ["donor", "cell_line", "timepoint", "pool"])
    out_path = OUT_DIR / "cohort/cell_line_donor_pool_timepoint_n_cells.csv"
    table.to_csv(out_path, index=False)
    print(table)
    print(f"\nSaved {len(table)} rows to {out_path}")

    celltype_table = make_table(
        obs, ["donor", "cell_line", "timepoint", "pool", "celltype"]
    )
    celltype_out_path = (
        OUT_DIR / "cohort/cell_line_donor_pool_timepoint_celltype_n_cells.csv"
    )
    celltype_table.to_csv(celltype_out_path, index=False)
    print(celltype_table)
    print(f"\nSaved {len(celltype_table)} rows to {celltype_out_path}")


if __name__ == "__main__":
    main()
