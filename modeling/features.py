"""
Per-fold, leakage-safe feature extraction, parameterized by timepoint (D11
now; D30 later for the baseline comparison model -- same 138 qualifying
lines, since qualifying already requires >=10 cells at D11 *and* D30
*and* D52).

Reuses 004_d11_celltype_proportion_se.py and 005_d11_pca_features.py's
streaming HVG/PCA pipeline (imported via importlib, since they're
numbered modules) generalized from "exclude pool11 always" (005's
original behavior) to "exclude whichever cells belong to this fold's
held-out lines, AND pool11 (frozen rule from EDA)". Fitting stats/HVGs/
PCA components must never see a held-out line's cells; held-out lines'
cells are still projected into that fold-specific PCA space -- the same
fit/project pattern 005 already used for pool11, now applied at the CV
boundary too.

Also builds D11 cell-type-proportion features (004's pattern, reused
directly) and the pool-correction step: residualize a pool-level feature
on pool identity using only the fold's TRAINING pool means, applied
before the pool-then-line two-stage averaging.
"""

import importlib.util
from pathlib import Path

import h5py
import pandas as pd

REPO_ROOT = Path(__file__).parent.parent
QUALIFYING_COMBOS_CSV = REPO_ROOT / "metadata_eda" / "cohort/qualifying_cell_line_pool_min10_per_timepoint.csv"

TIMEPOINT_FILES = {
    "D11": "/Users/prakritipaul/Documents/2021_jerber/day11.h5",
    "D30": "/Users/prakritipaul/Documents/2021_jerber/day30.h5",
}

# Depth-outlier pools, excluded from HVG/PCA *fitting* whenever present in a
# fold's training cells (their cells are still projected, so no line is lost).
# Never tuned on downstream model performance.
#
# THE RULE IS PER TIMEPOINT, NOT GLOBAL. It is: a pool is a depth outlier if
# its median raw UMI/cell is below 1/3 of the median-of-pool-medians at that
# timepoint. Applied to each timepoint's own depth table it gives:
#
#   D11: pool11 ratio 0.137, next-lowest pool12 at 0.747  -> {pool11}
#   D30: pool5  ratio 0.123, next-lowest pool12 at 0.423  -> {pool5}
#
# Both sides of the threshold are a wide gap, so 1/3 is a separator and not a
# tuned knob; the D11 entry reproduces the value frozen by EDA (005/006/008)
# exactly. Carrying D11's {pool11} over to D30 unchanged -- which an earlier
# single global constant did -- would have been wrong twice at D30: pool11 is
# the DEEPEST pool there (ratio 1.536) while pool5 collapses to ~1,050 median
# UMI/cell. See modeling/README.md "D30 depth outlier".
DEPTH_OUTLIER_POOLS = {
    "D11": frozenset({"pool11"}),
    "D30": frozenset({"pool5"}),
}


def _load_module(name: str):
    path = REPO_ROOT / "eda" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


m004 = _load_module("004_d11_celltype_proportion_se")
m005 = _load_module("005_d11_pca_features")


def load_cell_metadata(timepoint: str) -> pd.DataFrame:
    """cell_line, pool, celltype for every cell at this timepoint (row
    order matches the h5 file's own cell order)."""
    path = TIMEPOINT_FILES[timepoint]
    with h5py.File(path, "r") as f:
        cell_line = m004.read_obs_categorical(f, "donor_id")
        pool = m004.read_obs_categorical(f, "pool_id")
        celltype = m004.read_obs_categorical(f, "celltype")
    return pd.DataFrame({"cell_line": cell_line, "pool": pool, "celltype": celltype})


def compute_pca_features_for_fold(
    timepoint: str,
    held_out_lines: set[str],
    meta: pd.DataFrame | None = None,
    n_pcs: int = 10,
    exclude_pools: frozenset[str] | None = None,
    restrict_to_combos: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Fit HVGs/PCA on cells excluding held_out_lines and exclude_pools;
    project ALL cells (including held-out lines') into that space.
    Returns a per-cell DataFrame: cell_line, pool, PC1..PCk.

    exclude_pools: None (the default) resolves to THIS TIMEPOINT's depth
    outliers via DEPTH_OUTLIER_POOLS -- {pool11} at D11, {pool5} at D30.
    Defaulting to a single global set instead is how the D11 rule would
    silently follow the pipeline to D30 and exclude the wrong pool; pass an
    explicit frozenset only to override the rule deliberately (e.g. a
    sensitivity run), never to restate the default.

    restrict_to_combos: optional (cell_line, pool) frame further limiting
    the FIT population to those combos. None (the default) reproduces the
    original behaviour, in which cells outside the qualifying combos --
    including ~25,400 from 39 lines that are not among the 138 -- also
    contribute to gene stats, HVG selection and the PCA fit. Passing the
    qualifying combos gives the stricter "fit only on the study
    population" variant. Projection is unaffected either way: every cell
    is still projected into the fitted space."""
    path = TIMEPOINT_FILES[timepoint]
    if meta is None:
        meta = load_cell_metadata(timepoint)
    if exclude_pools is None:
        exclude_pools = DEPTH_OUTLIER_POOLS[timepoint]

    fit_mask = (~meta["cell_line"].isin(held_out_lines)) & (~meta["pool"].isin(exclude_pools))
    if restrict_to_combos is not None:
        eligible = pd.MultiIndex.from_frame(restrict_to_combos[["cell_line", "pool"]])
        in_combos = pd.MultiIndex.from_arrays([meta["cell_line"], meta["pool"]]).isin(eligible)
        fit_mask = fit_mask & pd.Series(in_combos, index=meta.index)
    fit_mask = fit_mask.to_numpy()

    mean, var, _n_fit = m005.compute_gene_stats(path, fit_mask=fit_mask)
    hvg_idx = m005.select_hvgs(mean, var, n_top=m005.N_HVG)
    X_hvg = m005.extract_hvg_matrix(path, hvg_idx)
    pcs, _explained = m005.run_pca(X_hvg, mean, var, hvg_idx, fit_mask, n_components=n_pcs)

    pc_cols = [f"PC{i}" for i in range(1, n_pcs + 1)]
    return pd.concat([meta[["cell_line", "pool"]].reset_index(drop=True), pd.DataFrame(pcs, columns=pc_cols)], axis=1)


def compute_proportion_features(meta: pd.DataFrame) -> pd.DataFrame:
    """Per (cell_line, pool): proportion of each celltype at this
    timepoint. Reuses 004's binomial-proportion logic directly (its `se`
    output isn't needed for modeling features and is dropped here)."""
    long = m004.compute_pool_level_proportions_and_se(meta[["cell_line", "pool", "celltype"]])
    wide = long.pivot(index=["cell_line", "pool"], columns="celltype", values="phat").reset_index()
    wide.columns = ["cell_line", "pool"] + [f"phat_{c}" for c in wide.columns[2:]]
    return wide


def pool_then_line_average(per_cell_or_pool: pd.DataFrame, value_cols: list[str]) -> pd.DataFrame:
    """Two-stage averaging (004/005's established pattern): mean per
    (cell_line, pool) if not already at that granularity, then mean of
    those pool-level values per cell_line."""
    if {"cell_line", "pool"}.issubset(per_cell_or_pool.columns) and len(per_cell_or_pool) > per_cell_or_pool[
        ["cell_line", "pool"]
    ].drop_duplicates().shape[0]:
        pool_level = per_cell_or_pool.groupby(["cell_line", "pool"], observed=True)[value_cols].mean().reset_index()
    else:
        pool_level = per_cell_or_pool
    return pool_level.groupby("cell_line", observed=True)[value_cols].mean().reset_index()


def pool_correct(pool_level: pd.DataFrame, value_cols: list[str], train_mask) -> pd.DataFrame:
    """Residualize each value column on pool identity, using only
    TRAINING rows' pool means (train_mask: boolean aligned with
    pool_level's rows). Applied to every row (train and held-out) using
    the training-derived pool means -- a held-out row's own value never
    informs the correction applied to it."""
    corrected = pool_level.copy()
    train_pool_means = pool_level.loc[train_mask].groupby("pool", observed=True)[value_cols].mean()
    for col in value_cols:
        corrected[col] = corrected[col] - corrected["pool"].map(train_pool_means[col])
    return corrected
