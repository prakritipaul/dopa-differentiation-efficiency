"""
Is the D52 outcome label (diff_efficiency, 009_d52_outcome_label.py)
itself confounded by technical covariates -- the same check
008_technical_covariate_associations.py did for the D11 features?

Matters because if pool or sequencing depth correlates with *both* a D11
feature and the D52 label, a "predictive" model could be learning that
shared confound rather than real biology, independent of whether each
side looks confounded on its own.

Same covariates as 008, computed at D52 instead of D11, at the same
(cell_line, pool) granularity: mean_total_counts/mean_n_genes_detected
(002's compute_qc_metrics on raw/X), n_cells, and pool identity
(eta-squared). Reuses 002/004/008 via importlib (same pattern as
008/009) rather than duplicating their logic.
"""

import importlib.util
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from data_paths import data_file  # noqa: E402  (repo root, added above)

DAY52_FILE = data_file("D52")
OUT_DIR = Path(__file__).parent.parent / "metadata_eda"
QUALIFYING_COMBOS_CSV = OUT_DIR / "cohort/qualifying_cell_line_pool_min10_per_timepoint.csv"
DIFFERENTIATED_CELLTYPES = {"DA", "Sert"}
NUMERIC_COVARIATES = ["n_cells", "mean_total_counts", "mean_n_genes_detected"]


def _load_module(name: str):
    path = Path(__file__).parent / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


m002 = _load_module("002_metadata_eda")
m004 = _load_module("004_d11_celltype_proportion_se")
m008 = _load_module("008_technical_covariate_associations")


def build_combo_level_table() -> pd.DataFrame:
    qualifying = pd.read_csv(QUALIFYING_COMBOS_CSV)[["cell_line", "pool"]]

    with h5py.File(DAY52_FILE, "r") as f:
        cell_line = m004.read_obs_categorical(f, "donor_id")
        pool = m004.read_obs_categorical(f, "pool_id")
        celltype = m004.read_obs_categorical(f, "celltype")
    n_genes_detected, total_counts = m002.compute_qc_metrics(DAY52_FILE)

    is_differentiated = pd.Series(celltype).isin(DIFFERENTIATED_CELLTYPES)
    binary_celltype = pd.Categorical(np.where(is_differentiated, "differentiated", "other"))

    per_cell = pd.DataFrame(
        {
            "cell_line": cell_line,
            "pool": pool,
            "celltype": binary_celltype,
            "total_counts": total_counts,
            "n_genes_detected": n_genes_detected,
        }
    )
    per_cell = per_cell.merge(qualifying, on=["cell_line", "pool"], how="inner")

    combo_qc = per_cell.groupby(["cell_line", "pool"], observed=True).agg(
        n_cells=("total_counts", "size"),
        mean_total_counts=("total_counts", "mean"),
        mean_n_genes_detected=("n_genes_detected", "mean"),
    )

    pool_level = m004.compute_pool_level_proportions_and_se(per_cell[["cell_line", "pool", "celltype"]])
    diff_eff = pool_level.loc[pool_level["celltype"] == "differentiated"].set_index(["cell_line", "pool"])["phat"]
    diff_eff = diff_eff.rename("diff_efficiency")

    return combo_qc.join(diff_eff, how="inner").reset_index()


def main() -> None:
    OUT_DIR.mkdir(exist_ok=True)

    combo = build_combo_level_table()
    print(f"{len(combo)} qualifying (cell_line, pool) combos.")

    row = {"target": "diff_efficiency"}
    for cov in NUMERIC_COVARIATES:
        r, p = stats.pearsonr(combo["diff_efficiency"], combo[cov])
        row[f"r_{cov}"] = r
        row[f"p_{cov}"] = p
    row["eta_sq_pool"] = m008.eta_squared_by_pool(combo["diff_efficiency"], combo["pool"])
    assoc = pd.DataFrame([row])
    assoc.to_csv(OUT_DIR / "technical/technical_covariate_correlations_d52_label.csv", index=False)

    print(assoc.to_string(index=False))
    print(
        f"\ndiff_efficiency vs mean_total_counts: r={row['r_mean_total_counts']:.3f} "
        f"(p={row['p_mean_total_counts']:.2e}); eta^2(pool)={row['eta_sq_pool']:.3f}"
    )

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    axes[0].scatter(combo["mean_total_counts"], combo["diff_efficiency"], s=15, alpha=0.6)
    axes[0].set_xlabel("mean_total_counts (D52)")
    axes[0].set_ylabel("diff_efficiency")
    axes[0].set_title(f"r={row['r_mean_total_counts']:.2f}, p={row['p_mean_total_counts']:.1e}")

    pool_order = combo.groupby("pool")["diff_efficiency"].median().sort_values().index
    data = [combo.loc[combo["pool"] == p, "diff_efficiency"] for p in pool_order]
    axes[1].boxplot(data, tick_labels=pool_order)
    axes[1].set_xlabel("pool")
    axes[1].set_ylabel("diff_efficiency")
    axes[1].set_title(f"eta^2(pool)={row['eta_sq_pool']:.2f}")
    axes[1].tick_params(axis="x", labelsize=7, rotation=45)

    fig.tight_layout()
    fig.savefig(OUT_DIR / "plots/plot_d52_label_technical_covariates.png", dpi=150)
    plt.close(fig)

    print(f"\nSaved CSV and plot to {OUT_DIR}")


if __name__ == "__main__":
    main()
