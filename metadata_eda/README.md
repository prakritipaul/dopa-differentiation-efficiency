# metadata_eda

Outputs from the EDA scripts (`001_eda.py`, `002_metadata_eda.py`,
`003_qualifying_cell_lines.py`, in the repo root) run against the Jerber et
al. dopaminergic neuron differentiation dataset
(https://pmc.ncbi.nlm.nih.gov/articles/PMC7610897/): `day11.h5`, `day30.h5`,
`day52.h5`.

## Layout

```
metadata_eda/
  cohort/        who is in the study: qualifying (cell_line, pool) combos,
                 the D52 label, per-line/pool/timepoint cell counts
  pca/           D11 PCA outputs -- per-line and per-cell coordinates, gene
                 loadings, variance explained, HVG lists.
                 `*_qualonly` = the qualifying-cells-only fitting variant
  proportions/   D11 cell-type proportions + their sampling SE
  qc/            QC summaries and metadata crosstabs
  technical/     technical-covariate (pool/batch) association analyses
  plots/         all figures
```

`cohort/qualifying_cell_line_pool_min10_per_timepoint.csv`,
`cohort/d52_diff_efficiency_label.csv`,
`proportions/d11_celltype_proportions_with_se.csv` and
`pca/d11_pca_coords_per_line.csv` are the interface consumed by
`modeling/` -- changing their paths means updating the modeling package.

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
be computed from `raw/X` directly (see `qc/qc_summary_by_timepoint_celltype.csv`
below).

## From `001_eda.py`

- **`cohort/cell_line_donor_pool_timepoint_n_cells.csv`** — n_cells grouped by
  `donor, cell_line, timepoint, pool`.
- **`cohort/cell_line_donor_pool_timepoint_celltype_n_cells.csv`** — same, with
  `celltype` added as an extra grouping column.

## From `002_metadata_eda.py`

- **`qc/metadata_crosstab_treatment.csv`** — n_cells by `timepoint x treatment`.
  Surfaces that `treatment` is `NONE` for D11/D30 but `NONE`/`ROT` for D52
  (the paper's rotenone oxidative-stress condition).
- **`qc/metadata_crosstab_celltype.csv`** — n_cells by `timepoint x celltype`.
- **`qc/metadata_crosstab_celltype_by_cluster.csv`** — n_cells by
  `timepoint x cluster_id x celltype`, to see how the numeric Leiden/Louvain
  `cluster_id` maps onto annotated cell types.
- **`qc/qc_summary_by_timepoint_celltype.csv`** — `describe()` (count, mean,
  std, min/25/50/75/max) of `total_counts` and `n_genes_detected` per cell,
  grouped by `timepoint x celltype`. These two per-cell metrics are computed
  from `raw/X` (raw UMI counts), streamed in 20k-cell chunks so the full
  multi-GB sparse arrays are never loaded at once.
- **`plots/plot_celltype_pool_bars.png`** — bar charts of n_cells by celltype and
  by pool, grouped by timepoint.
- **`plots/plot_umap_by_celltype.png`** — per-timepoint UMAP scatter (using each
  file's precomputed `obsm/X_umap`), colored by celltype. Downsampled to at
  most 50,000 cells per timepoint for plot legibility/speed.
- **`plots/plot_qc_distributions.png`** — histograms of `total_counts` (log10) and
  `n_genes_detected` per cell, overlaid by timepoint.

## From `003_qualifying_cell_lines.py`

- **`cohort/qualifying_cell_line_pool_min10_per_timepoint.csv`** — every
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
to the 138 lines in `cohort/qualifying_cell_line_pool_min10_per_timepoint.csv`
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

- **`proportions/d11_celltype_proportions_with_se.csv`** — per-line `n_pools`,
  `phat`/`se` (the `p_i`/`SE_i` above) for each D11 celltype (FPP, P_FPP,
  NB).
- **`proportions/d11_celltype_proportion_variance_vs_se.csv`** — per-celltype summary:
  `observed_std` (std of `p_i` across lines, i.e. `sqrt(Var(p_1, p_2,
  ...))`) vs. `rms_se` (`sqrt(mean(SE_i^2))` — the sampling-noise SD on the
  same scale; not the plain mean of `SE_i`), plus `variance_ratio` =
  `Var(p_i) / mean(SE_i^2)`. All three D11 celltypes have `observed_std`
  4-8x `rms_se` (variance_ratio 14-58x) — real line-to-line signal clearly
  dominates counting noise, so D11 cell type proportions look like a
  reasonable feature to use.
- **`plots/plot_d11_celltype_proportion_vs_se.png`** — per celltype, `phat` per
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
`cohort/qualifying_cell_line_pool_min10_per_timepoint.csv`.

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

- **`pca/d11_hvg_genes.csv`** — the 2000 selected HVGs (gene symbol, index,
  mean, variance — from non-`pool11` cells).
- **`pca/d11_pca_variance_explained.csv`** — explained variance ratio per PC.
- **`pca/d11_pca_coords_per_line.csv`** — the 138-line x `(PC1..PC10)` feature
  table.
- **`plots/plot_d11_pca_scree.png`** — explained variance ratio per PC.
- **`plots/plot_d11_pca_scatter.png`** — PC1 vs PC2 of the 138 line-level means.

## From `006_d11_pca_variance_vs_se.py`

Sanity check for the PCA feature, same idea as `004`'s check for cell type
proportions but generalized from a binomial proportion to a continuous
per-cell value: a line's `PC_k` feature is the mean of many single-cell
`PC_k` values, so its standard error is the standard error of that mean,
`SE_i = std(cell-level PC_k in line i) / sqrt(n_i)` (exactly analogous to
`SE_i = sqrt(p_i*(1-p_i)/n_i)`, which is also an SE-of-a-mean, just of
0/1-valued cells). Same two-stage pool-then-line averaging and SE
propagation as `004`/`005`.

Also adds an `icc` column: `icc = 1 - 1/variance_ratio`, a bounded (0-1)
reframing of the same ratio (`observed_var` already includes averaged
sampling noise, so `variance_ratio`'s "no signal" floor is 1, not 0; `icc`
maps that to 0, and -> 1 for strong signal) — the same quantity as
`variance_ratio`, just easier to read.

**Result, and a caveat found while checking it:** every PC looks
extremely reliable at face value (`variance_ratio` 24-2038x, `icc`
0.96-0.9995) — but PC8 and PC9's numbers are almost entirely a residual
`pool11` artifact (the sequencing-depth batch outlier from `005`):
excluding the 16 pool11-derived lines drops their `variance_ratio` by
30-60x (2038 -> 32, 1047 -> 49; see `variance_ratio_excl_pool11` /
`icc_excl_pool11` columns and the sharp, isolated jump for pool11 lines in
`plots/plot_d11_pca_variance_vs_se.png`'s PC8/PC9 panels). This is expected: a
technical batch that shifts every one of its cells the same way is
*reliably different* by this test, indistinguishable from real biology
using variance-vs-noise alone. PC1/PC2 are more robust (24 -> 12, 663 ->
389 when excluding pool11) but still somewhat inflated. All PCs remain
`>> 1` even excluding pool11 (12-390x), so all 10 look like real,
reliable signal — just noting that the *size* of that reliability
shouldn't be read at face value for PC8/PC9 specifically, and that this
check only establishes reliability, not that a PC predicts D52 efficiency
(a separate question for later cross-validated regression).

- **`pca/d11_pca_line_level_with_se.csv`** — per-line `n_pools` and, for each
  PC, `{PC}_mean`/`{PC}_se` (the `p_i`/`SE_i` above).
- **`pca/d11_pca_variance_vs_se.csv`** — per-PC summary: `observed_std`,
  `rms_se`, `variance_ratio`, `icc`, plus `variance_ratio_excl_pool11` /
  `icc_excl_pool11`.
- **`plots/plot_d11_pca_variance_vs_se.png`** — one panel per PC, each line's
  value sorted with an SE error bar.

## From `007_d11_cell_counts_per_line.py`

How many D11 cells does each cell line have — with a per-pool breakdown
for lines profiled in more than one pool? Uses
`cohort/cell_line_donor_pool_timepoint_n_cells.csv` (from `001`), so covers
every D11 cell line (177), not just the 138 qualifying ones. Also flags
which lines are among those 138, for cross-reference — a line can have
plenty of D11 cells and still not qualify, since qualifying also requires
>= 10 cells in the *same pool* at D30 and D52.

- **`cohort/d11_cell_counts_per_line.csv`** — per line: `donor`, `n_pools_D11`,
  `pool_breakdown_D11` (e.g. `pool4:2125, pool5:8321`),
  `total_n_cells_D11`, `is_qualifying_138`. Sorted by
  `total_n_cells_D11` descending. 25 of the 177 lines span more than one
  pool at D11; totals range from 1 to 14,640 cells.
- **`plots/plot_d11_cell_counts_per_line.png`** — sorted bar chart, colored by
  qualifying status.

## From `008_technical_covariate_associations.py`

Formalizes the ad hoc pool11 finding (`006`) into a systematic
correlation/association screen. Recomputes the **uncorrected** PCA (fit on
*all* D11 cells, i.e. the pre-`005`-fix version that pool11 dominated) and
checks it — plus the `004` cell type proportions, which had never been
checked this way — against technical covariates at the 159-combo
granularity: `mean_total_counts`/`mean_n_genes_detected` (Pearson r, from
`002`'s QC metrics), `n_cells` (Pearson r), and `pool` identity
(eta-squared / one-way-ANOVA variance-explained — the direct
generalization of "PC1 grouped by pool").

**Headline result:** as expected, uncorrected PC1 is almost entirely
explained by pool (`eta_sq_pool=0.999`) and strongly anti-correlated with
sequencing depth (`r=-0.87` with `mean_total_counts`, `r=-0.91` with
`mean_n_genes_detected`) — confirms the `005`/`006` pool11 story
quantitatively. PC5 and PC7 also show strong pool associations
(`eta_sq_pool` 0.87 and 0.83) not previously flagged; PC5's tracks
sequencing depth the same way PC1's does (`r≈-0.88`), while PC7's doesn't
(`r≈-0.03`) — a pool effect from something other than depth.

**More important: the cell type proportions are substantially
technical-covariate-confounded too**, which hadn't been checked before.
`phat_FPP` and `phat_P_FPP` correlate strongly with sequencing depth
(`r=0.73`/`r=-0.72` with `mean_total_counts`; `r=0.79`/`r=-0.75` with
`mean_n_genes_detected`) and pool (`eta_sq_pool` 0.79 and 0.69) — likely
because cell type calls are sensitive to capture depth (fewer marker
genes detected per cell biases classification). `phat_NB` looks clean
(`eta_sq_pool=0.03`, no significant correlations). This means the `004`
proportions feature, while *reliable* (validated via ICC in `004`/`006`),
is not free of the same technical confound found in the PCA feature —
worth accounting for (e.g. as a covariate, or checking whether it
survives controlling for pool/depth) before treating FPP/P_FPP as clean
predictors of D52 efficiency.

- **`technical/technical_covariate_correlations_pca_uncorrected.csv`** — 10 PCs x
  `{r, p}` per numeric covariate + `eta_sq_pool`.
- **`technical/technical_covariate_correlations_celltype_proportions.csv`** — same,
  for the 3 D11 celltypes.
- **`plots/plot_technical_covariate_correlations.png`** — correlation heatmaps
  (PCs and proportions vs. numeric covariates) and eta-squared-by-pool bar
  charts.
- **`technical/technical_covariate_analysis_plan.md`** — the design notes/plan
  written before implementing this and `007`.

## From `009_d52_outcome_label.py`

The D52 outcome label: **differentiation efficiency** = fraction of a
cell line's D52 cells that are DA (dopaminergic) or Sert (serotonergic)
neurons, the two mature "successfully differentiated" types — matching
the paper's own definition (`fig2b_and_heatmap_extended.ipynb`:
`diff_efficiency = DA_D52 + Sert_D52`). DA and Sert are combined into one
binary indicator *before* computing proportions/SE (not two proportions
summed afterward), since they're mutually exclusive outcomes of the same
per-cell draw. Otherwise identical methodology to `004`: binomial
proportion + SE per `(cell_line, pool)`, restricted to the 138 qualifying
lines, two-stage pool-then-line averaging with SE propagation for
multi-pool lines. Reuses `004`'s `compute_pool_level_proportions_and_se`
and `collapse_to_line_level` directly (via `importlib`, same pattern as
`008`).

138 lines get a label; wide spread (0.012 to 0.922, mean 0.45), with a
notable cluster of ~20 lines near 0 (essentially fail to produce DA/Sert
at all) and a broad spread from there — plenty of variance for a
downstream model to explain, and SE error bars are small relative to that
spread by eye.

- **`cohort/d52_diff_efficiency_label.csv`** — `cell_line`, `n_pools`,
  `diff_efficiency`, `diff_efficiency_se`.
- **`plots/plot_d52_diff_efficiency.png`** — sorted per-line values with SE
  error bars, and a histogram of the distribution.

## From `010_d52_label_technical_covariates.py`

Same technical-covariate check `008` ran on the D11 features, applied to
the D52 `diff_efficiency` label — matters because a shared confound on
both sides (a D11 feature *and* the D52 label) would let a "predictive"
model learn the confound instead of real biology, regardless of whether
either side looks confounded alone.

**Good news: the label looks clean.** `diff_efficiency` has weak,
non-significant correlations with D52 sequencing depth (`r=0.04`,
`p=0.63`) and genes detected (`r=0.10`, `p=0.22`), and a much smaller
`eta_sq_pool=0.13` than the D11 features saw (0.69-1.0). Notably, pool11
(the severe D11 sequencing-depth outlier) is unremarkable here — its
`diff_efficiency` box (median ~0.68, IQR 0.53-0.85) sits in the middle of
the pack, not off on its own. Pool medians do vary somewhat (0.25-0.68
across pools) but the boxes overlap heavily, consistent with the modest
eta-squared rather than a dominant batch effect.

- **`technical/technical_covariate_correlations_d52_label.csv`** — the association
  numbers above.
- **`plots/plot_d52_label_technical_covariates.png`** — scatter of
  `diff_efficiency` vs. mean D52 sequencing depth, and a boxplot by pool.

## Regenerating

From the repo root:

```
uv run python 001_eda.py
uv run python 002_metadata_eda.py   # streams raw/X for QC metrics; slower (~1-2 min)
uv run python 003_qualifying_cell_lines.py
uv run python 004_d11_celltype_proportion_se.py   # depends on 003's output CSV
uv run python 005_d11_pca_features.py             # depends on 003's output CSV; slower (~1-2 min)
uv run python 006_d11_pca_variance_vs_se.py        # depends on 005's per-cell-PC CSV
uv run python 007_d11_cell_counts_per_line.py      # depends on 001 and 003's output CSVs
uv run python 008_technical_covariate_associations.py  # depends on 003; slower (~1-2 min, recomputes uncorrected PCA)
uv run python 009_d52_outcome_label.py             # depends on 003's output CSV
uv run python 010_d52_label_technical_covariates.py  # depends on 003's output CSV
```

All ten scripts write into this directory (creating it if needed) and
read the source `.h5` files from the paths hardcoded in each script's
`DATA_FILES`. `008`, `009`, and `010` also import functions directly from
`002`/`004`/`005`/`008` via `importlib` (numbered modules aren't
importable with a plain `import`) rather than duplicating their logic.
