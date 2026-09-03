"""
Are the D11 features (PCA coordinates, cell type proportions) associated
with technical covariates rather than (or in addition to) biology?

006_d11_pca_variance_vs_se.py found this informally for one covariate: the
*uncorrected* PCA (fit on all D11 cells, before 005's fix that excludes
pool11 from fitting) had PC1 dominated by pool11's sequencing-depth batch
effect -- i.e. a correlation between pool and PC1. This script formalizes
that into a systematic screen: all 10 uncorrected PCs *and* the 3 D11
cell-type proportions (004_d11_celltype_proportion_se.py) against several
technical covariates, at the same (cell_line, pool) granularity used
throughout 004-006 (the 159 qualifying combos), since a batch effect at
that granularity is what actually reaches the final per-line features.

Covariates checked:
  - mean_total_counts, mean_n_genes_detected (per-cell QC metrics, from
    002_metadata_eda.py's compute_qc_metrics on raw/X, averaged per combo)
  - n_cells (the combo's D11 cell count -- pool size)
  - pool identity (categorical; the direct generalization of "pool11 was
    an outlier" into "how much of this feature's variance is explained by
    which pool it's from", i.e. eta-squared / one-way ANOVA)
  - sample_id is intentionally NOT included: it isn't 1:1 with a
    (cell_line, pool) combo (a pool can span multiple 10x samples), so it
    doesn't aggregate as cleanly as pool does at this granularity.

Reuses (via importlib, since these are numbered modules a plain `import`
can't reach) substantial existing logic rather than duplicating it:
compute_qc_metrics (002), load_d11_obs / compute_pool_level_proportions_and_se
(004), load_line_pool / compute_gene_stats / select_hvgs / extract_hvg_matrix
/ run_pca (005). The PCA here is deliberately run "uncorrected" -- fit_mask
covering every D11 cell for both fitting and transforming -- to match the
original, pre-005-fix PCA whose pool11 confound this script is diagnosing.
"""

import importlib.util
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

OUT_DIR = Path(__file__).parent / "metadata_eda"
QUALIFYING_COMBOS_CSV = OUT_DIR / "qualifying_cell_line_pool_min10_per_timepoint.csv"


def _load_module(name: str):
    path = Path(__file__).parent / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


m002 = _load_module("002_metadata_eda")
m004 = _load_module("004_d11_celltype_proportion_se")
m005 = _load_module("005_d11_pca_features")

PC_COLS = [f"PC{i}" for i in range(1, m005.N_PCS + 1)]
NUMERIC_COVARIATES = ["n_cells", "mean_total_counts", "mean_n_genes_detected"]


def compute_uncorrected_pcs() -> np.ndarray:
    """Per-cell PC1..PCk fit on ALL D11 cells (no pool11 exclusion) --
    matches the pre-005-fix PCA. Reuses 005's pipeline functions."""
    mean, var, _ = m005.compute_gene_stats(m005.DAY11_FILE, fit_mask=None)
    hvg_idx = m005.select_hvgs(mean, var, n_top=m005.N_HVG)
    X_hvg = m005.extract_hvg_matrix(m005.DAY11_FILE, hvg_idx)
    n_cells = X_hvg.shape[0]
    all_fit = np.ones(n_cells, dtype=bool)
    pcs, _ = m005.run_pca(X_hvg, mean, var, hvg_idx, all_fit, n_components=m005.N_PCS)
    return pcs


def build_combo_level_table() -> pd.DataFrame:
    qualifying = pd.read_csv(QUALIFYING_COMBOS_CSV)[["cell_line", "pool"]]

    line_pool = m005.load_line_pool()
    n_genes_detected, total_counts = m002.compute_qc_metrics(m005.DAY11_FILE)
    pcs = compute_uncorrected_pcs()

    per_cell = pd.concat(
        [
            line_pool,
            pd.DataFrame(pcs, columns=PC_COLS),
            pd.DataFrame({"total_counts": total_counts, "n_genes_detected": n_genes_detected}),
        ],
        axis=1,
    )
    per_cell = per_cell.merge(qualifying, on=["cell_line", "pool"], how="inner")

    combo_qc = per_cell.groupby(["cell_line", "pool"], observed=True).agg(
        n_cells=("total_counts", "size"),
        mean_total_counts=("total_counts", "mean"),
        mean_n_genes_detected=("n_genes_detected", "mean"),
        **{pc: (pc, "mean") for pc in PC_COLS},
    )

    obs = m004.load_d11_obs().merge(qualifying, on=["cell_line", "pool"], how="inner")
    pool_props = m004.compute_pool_level_proportions_and_se(obs)
    phat_wide = pool_props.pivot(index=["cell_line", "pool"], columns="celltype", values="phat")
    phat_wide.columns = [f"phat_{c}" for c in phat_wide.columns]

    combo = combo_qc.join(phat_wide, how="inner").reset_index()
    return combo


def eta_squared_by_pool(values: pd.Series, pool: pd.Series) -> float:
    grand_mean = values.mean()
    ss_total = ((values - grand_mean) ** 2).sum()
    group_means = values.groupby(pool).transform("mean")
    ss_within = ((values - group_means) ** 2).sum()
    return (ss_total - ss_within) / ss_total


def compute_associations(combo: pd.DataFrame, targets: list[str]) -> pd.DataFrame:
    rows = []
    for target in targets:
        row = {"target": target}
        for cov in NUMERIC_COVARIATES:
            r, p = stats.pearsonr(combo[target], combo[cov])
            row[f"r_{cov}"] = r
            row[f"p_{cov}"] = p
        row["eta_sq_pool"] = eta_squared_by_pool(combo[target], combo["pool"])
        rows.append(row)
    return pd.DataFrame(rows)


def plot_associations(pca_assoc: pd.DataFrame, prop_assoc: pd.DataFrame) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    for ax, assoc, title in [
        (axes[0, 0], pca_assoc, "Uncorrected PCs vs. numeric covariates (Pearson r)"),
        (axes[0, 1], prop_assoc, "Celltype proportions vs. numeric covariates (Pearson r)"),
    ]:
        r_cols = [f"r_{c}" for c in NUMERIC_COVARIATES]
        mat = assoc.set_index("target")[r_cols].to_numpy()
        im = ax.imshow(mat, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
        ax.set_xticks(range(len(NUMERIC_COVARIATES)))
        ax.set_xticklabels(NUMERIC_COVARIATES, rotation=30, ha="right")
        ax.set_yticks(range(len(assoc)))
        ax.set_yticklabels(assoc["target"])
        for i in range(mat.shape[0]):
            for j in range(mat.shape[1]):
                ax.text(j, i, f"{mat[i, j]:.2f}", ha="center", va="center", fontsize=8)
        ax.set_title(title)
        fig.colorbar(im, ax=ax, fraction=0.046)

    for ax, assoc, title in [
        (axes[1, 0], pca_assoc, "Uncorrected PCs: variance explained by pool (eta^2)"),
        (axes[1, 1], prop_assoc, "Celltype proportions: variance explained by pool (eta^2)"),
    ]:
        ax.barh(assoc["target"], assoc["eta_sq_pool"], color="tab:orange")
        ax.set_xlim(0, 1)
        ax.set_xlabel("eta^2 (share of variance explained by pool)")
        ax.set_title(title)
        ax.invert_yaxis()

    fig.tight_layout()
    fig.savefig(OUT_DIR / "plot_technical_covariate_correlations.png", dpi=150)
    plt.close(fig)


def main() -> None:
    OUT_DIR.mkdir(exist_ok=True)

    print("Recomputing uncorrected PCA (fit on all D11 cells, including pool11)...")
    combo = build_combo_level_table()
    print(f"{len(combo)} qualifying (cell_line, pool) combos.")

    pca_assoc = compute_associations(combo, PC_COLS)
    pca_assoc.to_csv(OUT_DIR / "technical_covariate_correlations_pca_uncorrected.csv", index=False)

    prop_targets = [c for c in combo.columns if c.startswith("phat_")]
    prop_assoc = compute_associations(combo, prop_targets)
    prop_assoc.to_csv(OUT_DIR / "technical_covariate_correlations_celltype_proportions.csv", index=False)

    print("\n=== Uncorrected PCs vs. technical covariates ===")
    print(pca_assoc.sort_values("eta_sq_pool", ascending=False).to_string(index=False))
    print("\n=== Celltype proportions vs. technical covariates ===")
    print(prop_assoc.sort_values("eta_sq_pool", ascending=False).to_string(index=False))

    pc1_row = pca_assoc.loc[pca_assoc["target"] == "PC1"].iloc[0]
    print(
        f"\nPC1 vs mean_total_counts: r={pc1_row['r_mean_total_counts']:.3f} "
        f"(p={pc1_row['p_mean_total_counts']:.2e}); eta^2(pool)={pc1_row['eta_sq_pool']:.3f} "
        f"-- consistent with pool11 (low sequencing depth) driving uncorrected PC1."
    )
    max_prop_eta = prop_assoc["eta_sq_pool"].max()
    print(
        f"Largest celltype-proportion eta^2(pool) = {max_prop_eta:.3f} "
        f"(vs. up to {pca_assoc['eta_sq_pool'].max():.3f} for the uncorrected PCs)."
    )

    plot_associations(pca_assoc, prop_assoc)
    print(f"\nSaved CSVs and plot to {OUT_DIR}")


if __name__ == "__main__":
    main()
