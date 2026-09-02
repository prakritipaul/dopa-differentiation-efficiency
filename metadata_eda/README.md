# metadata_eda

Outputs from the EDA scripts (`001_eda.py`, `002_metadata_eda.py`,
`003_qualifying_cell_lines.py`, in the repo root) run against the Jerber et
al. dopaminergic neuron differentiation dataset
(https://pmc.ncbi.nlm.nih.gov/articles/PMC7610897/): `day11.h5`, `day30.h5`,
`day52.h5`.

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

## Regenerating

From the repo root:

```
uv run python 001_eda.py
uv run python 002_metadata_eda.py   # streams raw/X for QC metrics; slower (~1-2 min)
uv run python 003_qualifying_cell_lines.py
```

All three scripts write into this directory (creating it if needed) and
read the source `.h5` files from the paths hardcoded in each script's
`DATA_FILES`.
