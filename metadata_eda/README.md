# metadata_eda

Outputs from the EDA scripts (`001_eda.py`, `002_metadata_eda.py`,
`003_qualifying_cell_lines.py`, in the repo root) run against the Jerber et
al. dopaminergic neuron differentiation dataset
(https://pmc.ncbi.nlm.nih.gov/articles/PMC7610897/): `day11.h5`, `day30.h5`,
`day52.h5`.

## Source metadata fields (in the `.h5` files themselves)

`obs` columns, beyond `celltype`:

- **`donor_id`** — full iPSC line identifier (e.g. `HPSI0114i-eipl_1`);
  parsed into `donor` (the `HPSI####i` prefix) and the specific clonal line
- **`pool_id`** — which multiplexed differentiation pool (7-24 lines each)
  the cell came from
- **`sample_id` / `sample_index`** — the specific 10x sequencing sample
  within a pool (finer-grained than pool)
- **`time_point`** — D11 / D30 / D52
- **`treatment`** — `NONE` or `ROT` (rotenone oxidative-stress exposure,
  D52 only)
- **`cluster_id`** — a plain numeric per-timepoint clustering assignment
  (e.g. Leiden/Louvain cluster number), distinct from the annotated
  `celltype` label

Beyond `obs`:

- **`obsm/X_umap`** — precomputed 2D UMAP coordinates per cell, per
  timepoint
- **`var`** — gene symbols (`var/index`, 32,738 genes, identical set/order
  across all three files) and, for day11/day30 only, redundant per-sample
  Ensembl gene IDs (`var/gene_ids-0..N`)
- **`X`** / **`raw/X`** — the expression features themselves: `X` is
  log-normalized expression (float32), `raw/X` is raw UMI counts
  (float32-encoded integers), both genes x cells sparse matrices

No other per-cell covariates (e.g. sex, age, batch date, sequencing depth as
a stored field) are present in `obs` — total counts / genes detected had to
be computed from `raw/X` directly (see `qc_summary_by_timepoint_celltype.csv`
below).

## From `001_eda.py`

- **`cell_line_donor_pool_timepoint_n_cells.csv`** — n_cells grouped by
  `donor, cell_line, timepoint, pool`.
- **`cell_line_donor_pool_timepoint_celltype_n_cells.csv`** — same, with
  `celltype` added as an extra grouping column.

## From `002_metadata_eda.py`

- **`metadata_crosstab_treatment.csv`** — n_cells by `timepoint x treatment`.
  Surfaces that `treatment` is `NONE` for D11/D30 but `NONE`/`ROT` for D52
  (the paper's rotenone oxidative-stress condition).
- **`metadata_crosstab_celltype.csv`** — n_cells by `timepoint x celltype`.
- **`metadata_crosstab_celltype_by_cluster.csv`** — n_cells by
  `timepoint x cluster_id x celltype`, to see how the numeric Leiden/Louvain
  `cluster_id` maps onto annotated cell types.
- **`qc_summary_by_timepoint_celltype.csv`** — `describe()` (count, mean,
  std, min/25/50/75/max) of `total_counts` and `n_genes_detected` per cell,
  grouped by `timepoint x celltype`. These two per-cell metrics are computed
  from `raw/X` (raw UMI counts), streamed in 20k-cell chunks so the full
  multi-GB sparse arrays are never loaded at once.
- **`plot_celltype_pool_bars.png`** — bar charts of n_cells by celltype and
  by pool, grouped by timepoint.
- **`plot_umap_by_celltype.png`** — per-timepoint UMAP scatter (using each
  file's precomputed `obsm/X_umap`), colored by celltype. Downsampled to at
  most 50,000 cells per timepoint for plot legibility/speed.
- **`plot_qc_distributions.png`** — histograms of `total_counts` (log10) and
  `n_genes_detected` per cell, overlaid by timepoint.

## From `003_qualifying_cell_lines.py`

- **`qualifying_cell_line_pool_min10_per_timepoint.csv`** — every
  `(cell_line, pool)` combination with >= 10 cells in that *same* pool at
  all three timepoints (D11, D30, D52). 159 combos, spanning 138 distinct
  cell lines (some lines qualify through more than one pool).

  Filtering matches the method in the paper's own analysis notebook
  (`plotting_notebooks/Figure_2/fig2b_and_heatmap_extended.ipynb` in
  https://github.com/single-cell-genetics/singlecell_neuroseq_paper), which
  reports "138 lines": count cells per `(cell_line, pool)` at each
  timepoint, keep combos with >= 10 cells, then keep only combos that
  survive at all three timepoints. An earlier version of this script summed
  a cell line's cells across *different* pools before filtering, which is
  looser and over-counted (140 lines) — a cell line split as, e.g., 6 cells
  in one pool and 8 in another at a timepoint would incorrectly pass.

## From `004_d11_celltype_proportion_se.py`

Sanity check for a candidate feature (predicting D52 differentiation
efficiency from D11 state). **Each data point is a cell_line**, restricted
to the 138 lines in `qualifying_cell_line_pool_min10_per_timepoint.csv`
(the paper-matching set). A line qualifying through more than one pool has
its D11 cell type proportion **averaged across those pools** (the authors'
approach — 21 of the 138 lines average across 2+ pools):

```
p_i = mean(phat_k for each qualifying pool k of line i)
```

Averaging proportions (rather than pooling raw counts) means `SE_i` is the
standard error of a mean of independent estimates, not a single
bigger-sample binomial SE:

```
SE_i = sqrt(sum_k SE_k^2) / n_pools_i,   SE_k = sqrt(phat_k*(1-phat_k)/n_k)
```

(a line with only one qualifying pool has `n_pools_i=1`, reducing to the
plain binomial SE of that pool). The question: is the variance of a D11
cell type's proportion observed *across cell lines* bigger than the
sampling variance expected *within* each line, or is the apparent spread
just counting noise?

- **`d11_celltype_proportions_with_se.csv`** — per-line `n_pools`,
  `phat`/`se` (the `p_i`/`SE_i` above) for each D11 celltype (FPP, P_FPP,
  NB).
- **`d11_celltype_proportion_variance_vs_se.csv`** — per-celltype summary:
  `observed_std` (std of `p_i` across lines, i.e. `sqrt(Var(p_1, p_2,
  ...))`) vs. `rms_se` (`sqrt(mean(SE_i^2))` — the sampling-noise SD on the
  same scale; not the plain mean of `SE_i`), plus `variance_ratio` =
  `Var(p_i) / mean(SE_i^2)`. All three D11 celltypes have `observed_std`
  4-8x `rms_se` (variance_ratio 14-58x) — real line-to-line signal clearly
  dominates counting noise, so D11 cell type proportions look like a
  reasonable feature to use.
- **`plot_d11_celltype_proportion_vs_se.png`** — per celltype, `phat` per
  cell line (sorted, with `se` error bars) vs. a dashed line at the mean;
  visually, the spread across lines is far wider than any individual
  error bar.

## From `005_d11_pca_features.py`

Third candidate feature: per-cell_line PCA coordinates from D11 single-cell
expression. Pipeline: take all D11 cells -> use the file's existing
log-normalized `X` (already normalized; no renormalization done) -> select
2000 HVGs (genes binned into 20 bins by mean expression, top genes by
within-bin z-scored dispersion — a simplified variant of Scanpy's
`flavor='seurat'`, operating directly on log values) -> scale each HVG to
zero mean/unit variance and clip at +/-10 -> PCA to 10 components -> every
D11 cell gets `PC1..PC10` -> per cell line, mean PCs per `(cell_line,
pool)` first, then mean of those pool-means (matching `004`'s averaging
approach), restricted to the 138 lines in
`qualifying_cell_line_pool_min10_per_timepoint.csv`.

**Important caveat, found and corrected during this analysis:** `pool11`
is a severe sequencing-depth batch outlier (~1,900 mean UMI/cell vs.
~10,000-18,000 in every other pool). Fitting PCA on all D11 cells made PC1
almost entirely a "is this pool11" indicator (17% of variance, every
pool11 line pinned to one extreme, every other line near the mean) rather
than a differentiation-relevant axis — and since 15 of `pool11`'s 16
qualifying lines have no other qualifying pool, simply dropping `pool11`
would have shrunk the feature table from 138 to 123 lines and mismatched
the `004` feature table. Fix: HVG selection and PCA fitting (mean/variance
for scaling, and the PCA components themselves) use only the 243,623
non-`pool11` D11 cells; `pool11` cells are then *projected* into that
learned space (`PCA.transform`, not `.fit_transform`) so all 138 lines
still get PC coordinates. After the fix, PC1 explains only 5.4% of
variance (vs. 17% before) with a smooth scree decay, and the PC1-PC2
scatter of line-level means shows a continuous cloud, not a pool11
cluster. A residual, much smaller batch signature remains on PC8 (0.85%
of variance) where `pool11` lines still show higher spread — worth
keeping in mind if using the later PCs for modeling.

- **`d11_hvg_genes.csv`** — the 2000 selected HVGs (gene symbol, index,
  mean, variance — from non-`pool11` cells).
- **`d11_pca_variance_explained.csv`** — explained variance ratio per PC.
- **`d11_pca_coords_per_line.csv`** — the 138-line x `(PC1..PC10)` feature
  table.
- **`plot_d11_pca_scree.png`** — explained variance ratio per PC.
- **`plot_d11_pca_scatter.png`** — PC1 vs PC2 of the 138 line-level means.

## Regenerating

From the repo root:

```
uv run python 001_eda.py
uv run python 002_metadata_eda.py   # streams raw/X for QC metrics; slower (~1-2 min)
uv run python 003_qualifying_cell_lines.py
uv run python 004_d11_celltype_proportion_se.py   # depends on 003's output CSV
uv run python 005_d11_pca_features.py             # depends on 003's output CSV; slower (~1-2 min)
```

All five scripts write into this directory (creating it if needed) and
read the source `.h5` files from the paths hardcoded in each script's
`DATA_FILES`.
