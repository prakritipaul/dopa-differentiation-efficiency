# D11 EDA completion: cell counts per line + technical covariate associations

Notes/plan for two follow-up EDA passes, written out for reference before
implementation (`007_d11_cell_counts_per_line.py` and
`008_technical_covariate_associations.py`).

## Comment on the requested logic (technical covariate correlation, task 2/3)

The instinct here is good methodology. Pool11 was identified as a batch
outlier essentially by eyeballing "PC1 grouped by pool" and noticing every
pool11 line sat at one extreme. Formalizing that into an actual
correlation/association screen — all 10 uncorrected PCs against relevant
technical covariates (sequencing depth, genes detected, pool identity) —
is exactly the standard next step in single-cell QC, and it can catch
subtler confounds that wouldn't be as visually obvious as pool11's ~9x
depth gap (e.g. a modest but real correlation on some other PC that isn't
as blatant but still isn't biology).

The extension to the cell-type-proportion feature (task 3) is a good
complementary check: `004`/`006` validated that feature's *reliability*
(ICC — real signal, not sampling noise) but never checked whether it's
*itself* confounded by the same technical covariates. Reliability and
technical-confound-free are different properties; a feature can be highly
reliable (low sampling noise relative to between-line spread) while still
being driven substantially by pool/batch rather than biology — which is
in fact exactly what turned out to be true for the PCA feature's PC8/PC9
in `006`. Worth checking whether the proportions feature has the same
issue, since it hasn't been checked yet.

## Context

Two more EDA passes requested to round out the D11 feature work
(`004`-`006`) before building the D52 outcome label:

1. A per-cell-line D11 cell count table (individual pool counts + total),
   for completeness/reference.
2. Formalize the ad hoc finding that pool11 dominated the *uncorrected*
   PC1 (from before the `005` fix that excludes pool11 from PCA fitting)
   into a systematic correlation/association screen: all 10 uncorrected
   PCs against relevant technical covariates (not just pool identity).
3. The same association screen, but for the D11 cell-type-proportion
   feature (`004`) — checking whether that feature (already validated as
   *reliable* via ICC) is itself confounded by the same technical
   covariates.

## Task 1: `007_d11_cell_counts_per_line.py`

Cheap, no h5 access needed — `metadata_eda/cohort/cell_line_donor_pool_timepoint_n_cells.csv`
(from `001_eda.py`) already has per-`(cell_line, pool, timepoint)` counts
for *every* D11 cell line (not just the 138 qualifying ones).

- Filter to `timepoint == "D11"`.
- Per `cell_line`: a per-pool breakdown (pool -> n_cells) and
  `total_n_cells_D11` (sum across pools), `n_pools_D11`.
- Flag whether each line is among the 138 in
  `qualifying_cell_line_pool_min10_per_timepoint.csv` (cross-reference,
  since that list required >=10 cells at D11 *and* D30 *and* D52 in the
  same pool — a line can have plenty of D11 cells and still not qualify).
- Save `metadata_eda/cohort/d11_cell_counts_per_line.csv`, sorted by
  `total_n_cells_D11` descending; print summary stats (how many lines
  split across >1 pool at D11, min/median/max total).
- `metadata_eda/plots/plot_d11_cell_counts_per_line.png` — sorted bar chart.

## Task 2+3: `008_technical_covariate_associations.py`

Granularity: the 159 qualifying `(cell_line, pool)` combos — the same
unit used throughout `004`-`006`, since that's what actually feeds the
final features; a batch effect at this granularity is what matters
downstream.

**Reuse note (deliberate exception to this project's usual "duplicate the
one small helper" pattern in `004`/`005`/`006`):** this script needs
substantial existing logic (`compute_qc_metrics` from `002`,
`compute_gene_stats`/`select_hvgs`/`extract_hvg_matrix`/`run_pca` from
`005`, `compute_pool_level_proportions_and_se` from `004`) — duplicating
~150+ lines of nontrivial streaming/PCA code would be worse than importing
it. Since these are numbered modules (`002_...py` etc.), a plain `import`
doesn't work; will use `importlib.util.spec_from_file_location` (same
mechanism already used ad hoc in this session to check things like
per-pool QC stats).

Steps:
1. **Uncorrected PCA**: call `005`'s `compute_gene_stats(path,
   fit_mask=None)` (already supports no mask -> uses all cells),
   `select_hvgs`, `extract_hvg_matrix`, then `run_pca(..., fit_mask=<all
   True>)` so both fit and transform use every D11 cell, matching the
   original (pre-fix) PCA exactly. Produces per-cell `PC1..PC10`
   (uncorrected).
2. **QC covariates**: `002`'s `compute_qc_metrics` on `day11.h5` ->
   per-cell `total_counts`, `n_genes_detected` (already-proven streaming
   code, ~15s).
3. **Cell type proportions**: `004`'s `compute_pool_level_proportions_and_se`
   (or just its phat computation) for per-`(cell_line, pool)`
   `phat_FPP`/`phat_P_FPP`/`phat_NB`.
4. Restrict everything to the 159 qualifying combos (inner join against
   `qualifying_cell_line_pool_min10_per_timepoint.csv`); aggregate to one
   row per combo: `n_cells`, `mean_total_counts`, `mean_n_genes_detected`,
   `pool`, mean of each uncorrected `PC1..PC10`, and the 3 `phat_*`
   columns.
5. **Association tables**:
   - Numeric covariates (`n_cells`, `mean_total_counts`,
     `mean_n_genes_detected`) vs. each of the 10 uncorrected PCs and each
     of the 3 `phat_*`: Pearson r (and p-value).
   - Categorical `pool` vs. the same 13 targets: eta-squared (one-way
     ANOVA variance-explained), the natural generalization of "PC1 grouped
     by pool" to a number — directly quantifies what was eyeballed for
     pool11.
   - `sample_id` intentionally excluded as a covariate here: it's not
     cleanly 1:1 with a `(cell_line, pool)` combo (a pool can span
     multiple 10x samples), so it doesn't aggregate as cleanly as `pool`
     does at this granularity.
6. Outputs into `metadata_eda/`:
   - `technical_covariate_correlations_pca_uncorrected.csv` (10 PCs x 3
     numeric covariates, Pearson r + p) and matching eta-squared-by-pool
     column.
   - `technical_covariate_correlations_celltype_proportions.csv` (3
     celltypes x same covariates).
   - `plot_technical_covariate_correlations.png` — correlation heatmap(s)
     (PCs and proportions vs. numeric covariates) plus a bar chart of
     eta-squared-by-pool per PC/celltype, so the pool11 signature is
     visually obvious the same way it was in `006`'s plots.
7. Print a short "strongest associations" summary (e.g. confirm PC1 vs.
   `mean_total_counts` correlation and eta-squared-by-pool line up with
   the known pool11 story) and note whether the `004` proportions show any
   comparable technical confound (expected to be much weaker, but this is
   the first time it's actually been checked).

## Verification

- `007`: row count matches total distinct D11 cell lines in the source
  data (no cells silently dropped); spot-check one known multi-pool line
  (e.g. `HPSI0114i-kolf_2`) shows the right per-pool breakdown and total.
- `008`: 159-row intermediate table before aggregation; sanity check that
  PC1's correlation with `mean_total_counts` and its eta-squared-by-pool
  are both strongly elevated (expected, given the known pool11 story) and
  that these numbers are directionally consistent with `006`'s
  `variance_ratio_excl_pool11` finding (PC1 partially, PC8/9 heavily
  pool-driven).
- Update `metadata_eda/README.md` with sections for both new scripts and
  add them to the "Regenerating" list.
