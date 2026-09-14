"""
Second candidate D11 -> D52 predictive feature: per-cell_line PCA
coordinates from D11 single-cell expression.

Pipeline (per the user's spec):
  1. Take all D11 cells across all cell lines.
  2. Normalize/log-transform expression.
  3. Identify D11 highly variable genes (HVGs).
  4. Run PCA on those HVGs across all D11 cells -> every cell gets
     PC1..PCk.
  5. Per cell line: mean PCs per (cell_line, pool) first, then average
     those pool-means for the cell line (same two-stage averaging used for
     the cell-type-proportion feature in 004_d11_celltype_proportion_se.py).

Step 2 is already done in the source file: day11.h5's X is already
log-normalized float32 expression (raw/X is raw UMI counts; confirmed in
002_metadata_eda.py's check_var_consistency), so this script uses X as-is
rather than recomputing normalization from raw counts.

day11.h5's X has ~724M nonzero entries (~5.8 GB for data+indices) across
253,381 cells x 32,738 genes, so this avoids scanpy/anndata (not a project
dependency) and instead streams X in bounded cell chunks, same pattern as
compute_qc_metrics in 002_metadata_eda.py, extended to two passes:
  - Pass 1 (compute_gene_stats): per-gene mean/variance across all D11
    cells, accumulated chunk by chunk via sparse csr_matrix ops.
  - HVG selection (select_hvgs): genes binned by mean expression into 20
    bins; within each bin, z-score dispersion = log(var/mean); take the
    top N by that normalized dispersion. This is a simplified variant of
    Scanpy's flavor='seurat' method operating directly on already-log
    values (skipping its internal un-log/expm1 step) -- standard and
    defensible, but not a byte-for-byte reproduction of Scanpy's output.
  - Pass 2 (extract_hvg_matrix): stream X again, this time keeping only
    the selected HVG columns, into one dense (n_cells, n_hvg) array.

PCA is run on cells scaled to zero mean / unit variance per HVG (using the
pass-1 stats restricted to HVGs) and clipped at +/-10 (matches Scanpy's
sc.pp.scale default).
"""

import re
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp
from sklearn.decomposition import PCA

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from data_paths import data_file  # noqa: E402  (repo root, added above)

DAY11_FILE = data_file("D11")
DONOR_RE = re.compile(r"^(HPSI\d+i)-")
OUT_DIR = Path(__file__).parent.parent / "metadata_eda"
QUALIFYING_COMBOS_CSV = OUT_DIR / "cohort/qualifying_cell_line_pool_min10_per_timepoint.csv"

# This script is the FULL-FIT (nothing held out) counterpart of
# run_feature_extraction.py's per-fold PCA, and modeling/feature_importance.py
# puts columns from both into one table. So the two MUST agree on the h5 file,
# the cohort, and the excluded depth-outlier pool. Both therefore resolve from
# the same sources the fold extractor uses -- features.TIMEPOINT_FILES,
# features.DEPTH_OUTLIER_POOLS, folds.VARIANTS -- rather than being restated
# here, because a second copy is free to drift and nothing would raise.
TIMEPOINT_OUT_PREFIX = {"D11": "d11", "D30": "d30"}

N_HVG = 2000
N_HVG_BINS = 20
N_PCS = 10
CHUNK_CELLS = 20_000
SCALE_CLIP = 10.0


def read_obs_categorical(f: h5py.File, column: str) -> pd.Categorical:
    """Read a legacy-anndata categorical obs column (codes + __categories)."""
    codes = f[f"obs/{column}"][:]
    categories = [c.decode() for c in f[f"obs/__categories/{column}"][:]]
    return pd.Categorical.from_codes(codes, categories=categories)


def load_line_pool(h5_file: str = DAY11_FILE) -> pd.DataFrame:
    with h5py.File(h5_file, "r") as f:
        cell_line = read_obs_categorical(f, "donor_id")
        pool = read_obs_categorical(f, "pool_id")
    return pd.DataFrame({"cell_line": cell_line, "pool": pool})


def _iter_chunks(f: h5py.File, chunk_cells: int):
    """Yield (start, end, csr_matrix) for successive cell chunks of X."""
    indptr = f["X/indptr"][:]
    n_cells = len(indptr) - 1
    n_genes = f["var/index"].shape[0]
    data_ds = f["X/data"]
    indices_ds = f["X/indices"]

    for start in range(0, n_cells, chunk_cells):
        end = min(start + chunk_cells, n_cells)
        lo, hi = int(indptr[start]), int(indptr[end])
        chunk_data = data_ds[lo:hi]
        chunk_indices = indices_ds[lo:hi]
        chunk_indptr = indptr[start : end + 1] - lo
        chunk = sp.csr_matrix(
            (chunk_data, chunk_indices, chunk_indptr), shape=(end - start, n_genes)
        )
        yield start, end, chunk


def compute_gene_stats(
    path: str, chunk_cells: int = CHUNK_CELLS, fit_mask: np.ndarray | None = None
) -> tuple[np.ndarray, np.ndarray, int]:
    """Pass 1: per-gene mean and (population) variance across D11 cells.

    If fit_mask is given (boolean, aligned with the file's cell order),
    only cells where it's True contribute -- used to exclude pool11 (a
    severe sequencing-depth batch outlier, ~1,900 vs ~10,000-18,000 mean
    UMI/cell in every other pool) from the stats used to select HVGs and
    fit/scale the PCA, so that batch effect doesn't dominate the learned
    space. Excluded cells are still scored later via extract_hvg_matrix +
    PCA.transform, just not used to define the space itself.
    """
    with h5py.File(path, "r") as f:
        n_genes = f["var/index"].shape[0]
        gene_sum = np.zeros(n_genes, dtype=np.float64)
        gene_sumsq = np.zeros(n_genes, dtype=np.float64)
        n_fit_cells = 0

        for start, end, chunk in _iter_chunks(f, chunk_cells):
            if fit_mask is not None:
                chunk = chunk[fit_mask[start:end]]
            n_fit_cells += chunk.shape[0]
            gene_sum += np.asarray(chunk.sum(axis=0)).ravel()
            gene_sumsq += np.asarray(chunk.multiply(chunk).sum(axis=0)).ravel()

    mean = gene_sum / n_fit_cells
    var = gene_sumsq / n_fit_cells - mean**2
    var = np.clip(var, 0, None)  # guard tiny negative values from float roundoff
    return mean, var, n_fit_cells


def select_hvgs(mean: np.ndarray, var: np.ndarray, n_top: int = N_HVG, n_bins: int = N_HVG_BINS) -> np.ndarray:
    """Top-n_top genes by mean-expression-binned, z-scored dispersion."""
    with np.errstate(divide="ignore", invalid="ignore"):
        dispersion = np.log(var / mean)
    dispersion[~np.isfinite(dispersion)] = np.nan

    valid = ~np.isnan(dispersion)
    bins = pd.cut(mean[valid], bins=n_bins, duplicates="drop")

    norm_dispersion = np.full(mean.shape, np.nan)
    disp_valid = pd.Series(dispersion[valid])
    grouped = disp_valid.groupby(bins.codes if hasattr(bins, "codes") else bins)
    bin_mean = grouped.transform("mean")
    bin_std = grouped.transform("std").replace(0, np.nan)
    z = (disp_valid - bin_mean) / bin_std
    norm_dispersion[np.where(valid)[0]] = z.values

    order = np.argsort(np.nan_to_num(norm_dispersion, nan=-np.inf))[::-1]
    return np.sort(order[:n_top])


def extract_hvg_matrix(path: str, hvg_idx: np.ndarray, chunk_cells: int = CHUNK_CELLS) -> np.ndarray:
    """Pass 2: dense (n_cells, n_hvg) matrix of just the HVG columns."""
    with h5py.File(path, "r") as f:
        n_cells = f["X/indptr"].shape[0] - 1
        out = np.empty((n_cells, len(hvg_idx)), dtype=np.float32)

        for start, end, chunk in _iter_chunks(f, chunk_cells):
            out[start:end, :] = chunk[:, hvg_idx].toarray()

    return out


def run_pca(
    X_hvg: np.ndarray,
    mean: np.ndarray,
    var: np.ndarray,
    hvg_idx: np.ndarray,
    fit_mask: np.ndarray,
    n_components: int = N_PCS,
    return_model: bool = False,
):
    """Fit PCA (and the scaling it uses) on cells where fit_mask is True
    only, then transform every cell (fit_mask True or False) into that
    space -- so excluded (pool11) cells still get PC coordinates, without
    having influenced what the axes mean."""
    hvg_mean = mean[hvg_idx]
    hvg_std = np.sqrt(var[hvg_idx])
    hvg_std[hvg_std == 0] = 1.0  # avoid div-by-zero for zero-variance genes (shouldn't be selected, but be safe)

    X_scaled = (X_hvg - hvg_mean) / hvg_std
    np.clip(X_scaled, -SCALE_CLIP, SCALE_CLIP, out=X_scaled)

    pca = PCA(n_components=n_components, svd_solver="randomized", random_state=0)
    pca.fit(X_scaled[fit_mask])
    pcs = pca.transform(X_scaled)
    # Default stays a 2-tuple: modeling/features.py unpacks exactly two
    # values, so returning the model is opt-in.
    if return_model:
        return pcs, pca.explained_variance_ratio_, pca
    return pcs, pca.explained_variance_ratio_


def collapse_to_line_level(pc_df: pd.DataFrame, pc_cols: list[str], qualifying: pd.DataFrame) -> pd.DataFrame:
    pc_df = pc_df.merge(qualifying, on=["cell_line", "pool"], how="inner")
    pool_means = pc_df.groupby(["cell_line", "pool"], observed=True)[pc_cols].mean()
    line_means = pool_means.groupby("cell_line", observed=True)[pc_cols].mean()
    return line_means.reset_index()


def plot_scree(explained_variance_ratio: np.ndarray, out_suffix: str = "", prefix: str = "d11") -> None:
    fig, ax = plt.subplots(figsize=(6, 4))
    x = np.arange(1, len(explained_variance_ratio) + 1)
    ax.plot(x, explained_variance_ratio, marker="o")
    ax.set_xlabel("PC")
    ax.set_ylabel("Explained variance ratio")
    ax.set_title("D11 PCA scree plot")
    ax.set_xticks(x)
    fig.tight_layout()
    fig.savefig(OUT_DIR / f"plots/plot_{prefix}_pca_scree{out_suffix}.png", dpi=150)
    plt.close(fig)


def plot_scatter(line_level: pd.DataFrame, out_suffix: str = "", prefix: str = "d11") -> None:
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.scatter(line_level["PC1"], line_level["PC2"], s=15, alpha=0.7)
    ax.set_xlabel("PC1 (line-level mean)")
    ax.set_ylabel("PC2 (line-level mean)")
    ax.set_title("D11 PCA: per-cell-line PC1 vs PC2")
    fig.tight_layout()
    fig.savefig(OUT_DIR / f"plots/plot_{prefix}_pca_scatter{out_suffix}.png", dpi=150)
    plt.close(fig)


def main(restrict_fit_to_qualifying: bool = False, out_suffix: str = "",
         timepoint: str = "D11", label_variant: str | None = None) -> None:
    """restrict_fit_to_qualifying: when True the HVG/PCA fit additionally
    excludes cells outside the qualifying (cell_line, pool) combos -- the
    study population. Default False reproduces the original global fit,
    which also uses cells from the 39 lines outside the 138.

    This is the FULL-FIT (no held-out split) counterpart of
    run_feature_extraction.py's per-fold restriction. Both must agree:
    modeling/feature_importance.py mixes fold-level results (LOCO,
    permutation) with full-fit results (coefficients, SHAP, univariate),
    so if only one of the two is restricted, a single importance table
    ends up describing two different PCA bases.

    timepoint: which h5 to read and which depth-outlier pool to exclude. The
    exclusion is NOT the same pool at both timepoints -- pool11 at D11, pool5
    at D30 -- and hardcoding D11's rule here would silently fit the D30 basis
    on its shallowest pool while excluding its deepest.

    label_variant: whose cohort defines "qualifying". Defaults to the
    published cohort when None, preserving the original behaviour.

    Pass a distinct out_suffix so the variant's outputs sit beside the
    originals instead of overwriting them."""
    from modeling.features import DEPTH_OUTLIER_POOLS, TIMEPOINT_FILES
    from modeling.folds import get_variant

    OUT_DIR.mkdir(exist_ok=True)
    h5_file = TIMEPOINT_FILES[timepoint]
    (excluded_pool,) = DEPTH_OUTLIER_POOLS[timepoint]
    prefix = TIMEPOINT_OUT_PREFIX[timepoint]
    combos_csv = QUALIFYING_COMBOS_CSV if label_variant is None else get_variant(label_variant).cohort_csv

    line_pool = load_line_pool(h5_file)
    fit_mask = (line_pool["pool"] != excluded_pool).to_numpy()
    if restrict_fit_to_qualifying:
        qual = pd.read_csv(combos_csv)[["cell_line", "pool"]]
        in_combos = pd.MultiIndex.from_arrays(
            [line_pool["cell_line"], line_pool["pool"]]
        ).isin(pd.MultiIndex.from_frame(qual))
        fit_mask = fit_mask & np.asarray(in_combos)
    print(
        f"Fitting HVGs/PCA on {fit_mask.sum()} of {len(fit_mask)} {timepoint} cells "
        f"(excluding {excluded_pool} -- the sequencing-depth outlier at this timepoint"
        + (f", and restricted to {combos_csv.name}" if restrict_fit_to_qualifying else "")
        + "). Excluded cells are still projected into the resulting PCA space."
    )

    print("Pass 1/2: computing per-gene mean/variance (fit cells only)...")
    mean, var, n_fit_cells = compute_gene_stats(h5_file, fit_mask=fit_mask)
    print(f"{n_fit_cells} fit cells, {len(mean)} genes.")

    hvg_idx = select_hvgs(mean, var, n_top=N_HVG)
    print(f"Selected {len(hvg_idx)} HVGs.")

    with h5py.File(h5_file, "r") as f:
        gene_symbols = np.array([g.decode() for g in f["var/index"][:]])

    hvg_table = pd.DataFrame(
        {
            "gene": gene_symbols[hvg_idx],
            "gene_index": hvg_idx,
            "mean": mean[hvg_idx],
            "var": var[hvg_idx],
        }
    )
    hvg_table.to_csv(OUT_DIR / f"pca/{prefix}_hvg_genes{out_suffix}.csv", index=False)

    print("Pass 2/2: extracting HVG-only expression matrix (all D11 cells)...")
    X_hvg = extract_hvg_matrix(h5_file, hvg_idx)

    print(f"Running PCA (n_components={N_PCS}, fit on non-pool11 cells, transform all)...")
    pcs, explained_variance_ratio, pca_model = run_pca(
        X_hvg, mean, var, hvg_idx, fit_mask, n_components=N_PCS, return_model=True
    )

    # Gene loadings: pca.components_ is (n_components, n_HVGs) in the space
    # the PCA was fit on -- STANDARDIZED expression, so a loading is the
    # weight per 1 SD of that gene, comparable across genes of different
    # absolute expression. Saved long-format for all 10 PCs.
    #
    # Sign is arbitrary: PCA fixes an axis, not a direction, so only the
    # loading's sign RELATIVE to other genes on the same PC is meaningful,
    # and the direction can flip between bases or reruns of a different fit.
    loadings = pd.DataFrame(
        pca_model.components_.T,
        columns=[f"PC{i}" for i in range(1, N_PCS + 1)],
    )
    loadings.insert(0, "gene", gene_symbols[hvg_idx])
    loadings.insert(1, "gene_index", hvg_idx)
    loadings_long = loadings.melt(
        id_vars=["gene", "gene_index"], var_name="PC", value_name="loading"
    )
    loadings_long["abs_loading"] = loadings_long["loading"].abs()
    loadings_long = loadings_long.sort_values(
        ["PC", "abs_loading"], ascending=[True, False]
    ).reset_index(drop=True)
    loadings_long.to_csv(OUT_DIR / f"pca/{prefix}_pca_gene_loadings{out_suffix}.csv", index=False)
    print(f"Saved gene loadings for {N_PCS} PCs x {len(hvg_idx)} HVGs.")

    pc_cols = [f"PC{i}" for i in range(1, N_PCS + 1)]
    variance_table = pd.DataFrame({"PC": pc_cols, "explained_variance_ratio": explained_variance_ratio})
    variance_table.to_csv(OUT_DIR / f"pca/{prefix}_pca_variance_explained{out_suffix}.csv", index=False)
    print(variance_table)

    pc_df = pd.concat([line_pool, pd.DataFrame(pcs, columns=pc_cols)], axis=1)

    qualifying = pd.read_csv(combos_csv)[["cell_line", "pool"]]

    # Per-cell PCs, restricted to qualifying combos -- reused by
    # 006_d11_pca_variance_vs_se.py to check these features the same way
    # 004 checks cell type proportions (variance across lines vs. SE).
    pc_df_qualifying = pc_df.merge(qualifying, on=["cell_line", "pool"], how="inner")
    pc_df_qualifying.to_csv(OUT_DIR / f"pca/{prefix}_pca_coords_per_cell_qualifying{out_suffix}.csv", index=False)

    line_level = collapse_to_line_level(pc_df, pc_cols, qualifying)
    line_level.to_csv(OUT_DIR / f"pca/{prefix}_pca_coords_per_line{out_suffix}.csv", index=False)

    print(f"\n{line_level['cell_line'].nunique()} cell lines in the final PCA feature table.")
    print(line_level.head())

    plot_scree(explained_variance_ratio, out_suffix, prefix)
    plot_scatter(line_level, out_suffix, prefix)
    print(f"\nSaved CSVs and plots to {OUT_DIR}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--restrict-to-qualifying", action="store_true",
                        help="fit HVGs/PCA only on cells in qualifying (cell_line, pool) combos")
    parser.add_argument("--suffix", default="", help="suffix for all output filenames")
    parser.add_argument("--timepoint", default="D11", choices=["D11", "D30"],
                        help="which h5 to fit; also picks the depth-outlier pool to exclude "
                             "(pool11 at D11, pool5 at D30) and the d11_/d30_ output prefix")
    parser.add_argument("--label-variant", default=None,
                        help="whose cohort defines 'qualifying'; default is the published cohort")
    args = parser.parse_args()
    main(restrict_fit_to_qualifying=args.restrict_to_qualifying, out_suffix=args.suffix,
         timepoint=args.timepoint, label_variant=args.label_variant)
