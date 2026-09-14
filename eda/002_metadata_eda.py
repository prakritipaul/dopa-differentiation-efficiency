"""
Metadata exploration and brief EDA for the Jerber et al. dopaminergic neuron
differentiation dataset (https://pmc.ncbi.nlm.nih.gov/articles/PMC7610897/).

Companion to 001_eda.py: that script builds the cell_line x donor x pool x
timepoint summary tables. This script explores the rest of the `obs`
metadata (celltype, cluster_id, treatment, sample_id), checks `var`
(gene) consistency across files, and computes light per-cell QC metrics
(total UMI counts, genes detected) by streaming the raw counts layer in
bounded-size chunks -- the only step here that reads real expression data
(raw/X) rather than just metadata, since raw/X holds hundreds of millions
to ~1.6B nonzero entries per file (several GB) and is never loaded whole.
"""

import re
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

DATA_FILES = [
    "/Users/prakritipaul/Documents/2021_jerber/day11.h5",
    "/Users/prakritipaul/Documents/2021_jerber/day30.h5",
    "/Users/prakritipaul/Documents/2021_jerber/day52.h5",
]

DONOR_RE = re.compile(r"^(HPSI\d+i)-")
OUT_DIR = Path(__file__).parent.parent / "metadata_eda"
UMAP_PLOT_MAX_CELLS = 50_000  # downsample for scatter legibility/speed only


def read_obs_categorical(f: h5py.File, column: str) -> pd.Categorical:
    """Read a legacy-anndata categorical obs column (codes + __categories)."""
    codes = f[f"obs/{column}"][:]
    categories = [c.decode() for c in f[f"obs/__categories/{column}"][:]]
    return pd.Categorical.from_codes(codes, categories=categories)


def load_full_obs(path: str) -> pd.DataFrame:
    with h5py.File(path, "r") as f:
        cell_line = read_obs_categorical(f, "donor_id")
        pool = read_obs_categorical(f, "pool_id")
        timepoint = read_obs_categorical(f, "time_point")
        celltype = read_obs_categorical(f, "celltype")
        treatment = read_obs_categorical(f, "treatment")
        sample_id = read_obs_categorical(f, "sample_id")
        cluster_id = f["obs/cluster_id"][:]
        umap = f["obsm/X_umap"][:]

    donor = [DONOR_RE.match(cl).group(1) for cl in cell_line]

    return pd.DataFrame(
        {
            "cell_line": cell_line,
            "donor": donor,
            "pool": pool,
            "timepoint": timepoint,
            "celltype": celltype,
            "treatment": treatment,
            "sample_id": sample_id,
            "cluster_id": cluster_id,
            "umap_1": umap[:, 0],
            "umap_2": umap[:, 1],
        }
    )


def compute_qc_metrics(
    path: str, chunk_cells: int = 20_000
) -> tuple[np.ndarray, np.ndarray]:
    """Per-cell (n_genes_detected, total_counts) from raw/X, streamed in chunks.

    n_genes_detected is free (just a diff of the CSR indptr). total_counts
    requires reading raw/X/data, but only chunk_cells worth of rows at a
    time, so peak memory stays bounded regardless of file size.
    """
    with h5py.File(path, "r") as f:
        indptr = f["raw/X/indptr"][:]
        data_ds = f["raw/X/data"]
        n_cells = len(indptr) - 1
        n_genes_detected = np.diff(indptr)
        total_counts = np.empty(n_cells, dtype=np.float64)

        for start in range(0, n_cells, chunk_cells):
            end = min(start + chunk_cells, n_cells)
            lo, hi = int(indptr[start]), int(indptr[end])
            chunk_data = data_ds[lo:hi]
            offsets = indptr[start : end + 1] - lo
            cumsum = np.concatenate(([0.0], np.cumsum(chunk_data, dtype=np.float64)))
            total_counts[start:end] = cumsum[offsets[1:]] - cumsum[offsets[:-1]]

    return n_genes_detected, total_counts


def make_table(obs: pd.DataFrame, group_cols: list[str]) -> pd.DataFrame:
    return (
        obs.groupby(group_cols, observed=True)
        .size()
        .reset_index(name="n_cells")
        .sort_values(group_cols)
        .reset_index(drop=True)
    )


def print_metadata_report(obs: pd.DataFrame) -> None:
    print("=" * 70)
    print("METADATA REPORT")
    print("=" * 70)

    for col in ["treatment", "celltype", "pool", "sample_id", "donor", "cell_line"]:
        print(f"\n--- {col}: distinct values per timepoint ---")
        print(obs.groupby("timepoint", observed=True)[col].nunique())

    print("\n--- cluster_id (numeric, per-timepoint clustering) ---")
    print(obs.groupby("timepoint", observed=True)["cluster_id"].agg(["min", "max", "nunique"]))

    print("\n--- treatment counts by timepoint (rotenone condition is D52-only) ---")
    print(obs.groupby(["timepoint", "treatment"], observed=True).size())

    print("\n--- celltype counts by timepoint ---")
    print(obs.groupby(["timepoint", "celltype"], observed=True).size())


def check_var_consistency() -> None:
    print("\n" + "=" * 70)
    print("VAR (GENE) CONSISTENCY CHECK")
    print("=" * 70)

    gene_symbols_by_file = []
    for path in DATA_FILES:
        with h5py.File(path, "r") as f:
            gene_symbols = f["var/index"][:]
            gene_symbols_by_file.append(gene_symbols)

            gene_id_keys = sorted(k for k in f["var"].keys() if k.startswith("gene_ids"))
            if gene_id_keys:
                first = f[f"var/{gene_id_keys[0]}"][:]
                identical_within_file = all(
                    np.array_equal(first, f[f"var/{k}"][:]) for k in gene_id_keys[1:]
                )
                print(
                    f"{Path(path).name}: n_genes={len(gene_symbols)}, "
                    f"gene_ids columns={len(gene_id_keys)}, "
                    f"identical across those columns={identical_within_file}"
                )
            else:
                print(
                    f"{Path(path).name}: n_genes={len(gene_symbols)}, "
                    "no per-sample gene_ids columns (only var/index)"
                )

    identical_across_files = all(
        np.array_equal(gene_symbols_by_file[0], g) for g in gene_symbols_by_file[1:]
    )
    print(f"gene symbol order identical across all 3 files: {identical_across_files}")
    print(
        "\nX/data is float32 log-normalized expression; "
        "raw/X/data is float32-encoded raw UMI counts."
    )


def save_crosstabs(obs: pd.DataFrame) -> None:
    treatment_x_timepoint = make_table(obs, ["timepoint", "treatment"])
    treatment_x_timepoint.to_csv(OUT_DIR / "qc/metadata_crosstab_treatment.csv", index=False)

    celltype_x_timepoint = make_table(obs, ["timepoint", "celltype"])
    celltype_x_timepoint.to_csv(OUT_DIR / "qc/metadata_crosstab_celltype.csv", index=False)

    celltype_x_cluster = make_table(obs, ["timepoint", "cluster_id", "celltype"])
    celltype_x_cluster.to_csv(
        OUT_DIR / "qc/metadata_crosstab_celltype_by_cluster.csv", index=False
    )

    print("\nSaved crosstab CSVs: treatment, celltype, celltype_by_cluster")


def save_qc_summary(obs: pd.DataFrame) -> pd.DataFrame:
    qc_summary = (
        obs.groupby(["timepoint", "celltype"], observed=True)[
            ["total_counts", "n_genes_detected"]
        ]
        .describe()
    )
    qc_summary.to_csv(OUT_DIR / "qc/qc_summary_by_timepoint_celltype.csv")
    print("\n--- QC summary (total_counts, n_genes_detected) by timepoint x celltype ---")
    print(qc_summary)
    return qc_summary


def plot_celltype_and_pool_bars(obs: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    celltype_counts = obs.groupby(["celltype", "timepoint"], observed=True).size().unstack("timepoint")
    celltype_counts.plot(kind="bar", ax=axes[0])
    axes[0].set_title("n_cells by celltype x timepoint")
    axes[0].set_ylabel("n_cells")

    pool_counts = obs.groupby(["pool", "timepoint"], observed=True).size().unstack("timepoint")
    pool_counts.plot(kind="bar", ax=axes[1])
    axes[1].set_title("n_cells by pool x timepoint")
    axes[1].set_ylabel("n_cells")
    axes[1].tick_params(axis="x", labelsize=7)

    fig.tight_layout()
    fig.savefig(OUT_DIR / "plots/plot_celltype_pool_bars.png", dpi=150)
    plt.close(fig)


def plot_umap_by_celltype(obs: pd.DataFrame) -> None:
    timepoints = list(obs["timepoint"].cat.categories)
    fig, axes = plt.subplots(1, len(timepoints), figsize=(6 * len(timepoints), 5))

    for ax, tp in zip(axes, timepoints):
        sub = obs.loc[obs["timepoint"] == tp]
        if len(sub) > UMAP_PLOT_MAX_CELLS:
            sub = sub.sample(UMAP_PLOT_MAX_CELLS, random_state=0)

        for ct in sub["celltype"].cat.categories:
            pts = sub.loc[sub["celltype"] == ct]
            ax.scatter(pts["umap_1"], pts["umap_2"], s=2, alpha=0.4, label=ct)

        ax.set_title(f"{tp} UMAP by celltype (n={len(sub)} plotted)")
        ax.legend(markerscale=5, fontsize=7, loc="best")

    fig.tight_layout()
    fig.savefig(OUT_DIR / "plots/plot_umap_by_celltype.png", dpi=150)
    plt.close(fig)


def plot_qc_distributions(obs: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    for tp in obs["timepoint"].cat.categories:
        vals = obs.loc[obs["timepoint"] == tp, "total_counts"]
        axes[0].hist(np.log10(vals + 1), bins=80, alpha=0.5, label=tp)
    axes[0].set_title("log10(total_counts + 1) by timepoint")
    axes[0].legend()

    for tp in obs["timepoint"].cat.categories:
        vals = obs.loc[obs["timepoint"] == tp, "n_genes_detected"]
        axes[1].hist(vals, bins=80, alpha=0.5, label=tp)
    axes[1].set_title("n_genes_detected by timepoint")
    axes[1].legend()

    fig.tight_layout()
    fig.savefig(OUT_DIR / "plots/plot_qc_distributions.png", dpi=150)
    plt.close(fig)


def main() -> None:
    OUT_DIR.mkdir(exist_ok=True)
    check_var_consistency()

    per_file_obs = []
    for path in DATA_FILES:
        obs = load_full_obs(path)
        print(f"\nStreaming raw counts for QC metrics: {Path(path).name} ...")
        n_genes_detected, total_counts = compute_qc_metrics(path)
        obs["n_genes_detected"] = n_genes_detected
        obs["total_counts"] = total_counts
        per_file_obs.append(obs)

    obs = pd.concat(per_file_obs, ignore_index=True)
    obs["timepoint"] = obs["timepoint"].astype("category")
    obs["celltype"] = obs["celltype"].astype("category")

    print_metadata_report(obs)
    save_crosstabs(obs)
    save_qc_summary(obs)

    plot_celltype_and_pool_bars(obs)
    plot_umap_by_celltype(obs)
    plot_qc_distributions(obs)

    print(f"\nSaved plots and CSVs to {OUT_DIR}")


if __name__ == "__main__":
    main()
