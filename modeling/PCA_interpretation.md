# PCA interpretation: what the D11 principal components mean

Companion to `modeling/README.md`. Covers the biological reading of the
D11 HVG/PCA basis used as model features, focused on **PC1**.

Loadings live in `metadata_eda/d11_pca_gene_loadings.csv` (baseline basis)
and `metadata_eda/d11_pca_gene_loadings_qualonly.csv` (qualifying-only basis) — all 10 PCs x 2000 HVGs,
long format, sorted by |loading| within each PC. Produced by
`005_d11_pca_features.py`.

---

## How the loadings were computed

`005_d11_pca_features.py` streams the D11 matrix, selects 2000 HVGs by
binned dispersion, standardizes each HVG to mean 0 / SD 1 (clipped at
`SCALE_CLIP`), fits `PCA(n_components=10, svd_solver="randomized",
random_state=0)` on the fitting subset, then projects all cells.

A loading is therefore **the weight on standardized expression** — the
contribution per 1 SD of that gene — which makes loadings comparable
across genes of very different absolute expression. They are **not**
fold-changes. Loadings are unit-norm per component (sum of squares = 1).

`random_state=0` is set, so the decomposition is deterministic:
regenerating reproduced the committed coordinate files byte-for-byte,
which is why these loadings correspond exactly to the PC values used in
every model.

**Sign is arbitrary.** PCA fixes an axis, not a direction. PC1 is
presented below oriented so that **positive = higher D52 efficiency**
(Spearman rho = +0.563, n = 138 lines). Do not carry that convention to
other PCs without re-checking their orientation.

---

## PC1 is a cell-cycle / proliferation axis

| property | value |
|---|---|
| variance explained | 5.4% of HVG variance |
| Spearman rho vs D52 `diff_efficiency` | **+0.563** (n = 138) |
| Spearman rho vs `phat_P_FPP` (proliferating FPP) | **+0.798** |
| Spearman rho vs `phat_NB` (neuroblast) | −0.552 |
| pool eta^2 (technical association) | 0.49 — *Moderate* |
| top 25 genes' share of total weight | 29.6% (median \|loading\| = 0.0012) |

### Quantitative enrichment, not eyeballing

Compared against the **Tirosh et al. (2016) G2/M and G1/S sets** as
distributed with Seurat (`cc.genes`):

- 39 of the 46 G2/M genes and 9 of the 43 G1/S genes survive HVG selection.
- **30 of the top 50 positive loadings are G2/M genes.** Zero are G1/S.
- **0 of the top 50 negative loadings** are in either set.
- The 39 G2/M genes have a median PC1 loading rank of **28 out of 2000**,
  median loading **+0.0990** against an all-HVG median of ~0.0000.
- Mann-Whitney (G2/M loadings > non-cell-cycle loadings): **p = 2.3e-26**.

The enrichment is at the **G2/M** end specifically, not cell-cycle
generally — G1/S genes are not enriched. That is the signature of cells
actively in mitosis, not merely cycling.

### Reading

Lines whose **D11 cultures are more proliferative differentiate better by
D52**. This is independently corroborated without going through the PCA
at all: PC1 tracks the *proliferating* FPP fraction (`phat_P_FPP`) at
rho = +0.798, and anti-correlates with the neuroblast fraction
(`phat_NB`, rho = −0.552) that separately predicts *worse* outcome.

The negative pole is ribosomal proteins, translation-initiation factors
and housekeeping/OXPHOS genes — the usual complement of a
high-proliferation axis, since a larger housekeeping/biosynthetic
transcript fraction is what non-dividing cells look like.

### Replication across the two PCA bases

PC1 loadings correlate **r = +0.9995** between the baseline and
qualifying-only bases, with **25/25** top-gene overlap. Unlike PC2/PC3 —
which reorder between bases (see `modeling/README.md`, "Do NOT compare the
two importance tables row by row") — **PC1 is stable and safe to compare
across bases.**

---

## Top 40 positive loadings (higher in high-efficiency lines)

| rank | gene | loading | canonical set |
|---:|---|---:|---|
| 1 | `HMGB2` | +0.1188 | G2/M |
| 2 | `PTTG1` | +0.1163 | — |
| 3 | `NUSAP1` | +0.1153 | G2/M |
| 4 | `UBE2C` | +0.1138 | G2/M |
| 5 | `CENPF` | +0.1136 | G2/M |
| 6 | `CKS2` | +0.1127 | G2/M |
| 7 | `TOP2A` | +0.1122 | G2/M |
| 8 | `PLK1` | +0.1113 | — |
| 9 | `CCNB1` | +0.1111 | — |
| 10 | `TUBB4B` | +0.1108 | G2/M |
| 11 | `KPNA2` | +0.1104 | — |
| 12 | `CCNB2` | +0.1087 | G2/M |
| 13 | `CDC20` | +0.1082 | G2/M |
| 14 | `BIRC5` | +0.1075 | G2/M |
| 15 | `CKS1B` | +0.1074 | G2/M |
| 16 | `AURKB` | +0.1071 | G2/M |
| 17 | `CCNA2` | +0.1071 | — |
| 18 | `TPX2` | +0.1064 | G2/M |
| 19 | `CDK1` | +0.1063 | G2/M |
| 20 | `ARL6IP1` | +0.1063 | — |
| 21 | `AURKA` | +0.1017 | G2/M |
| 22 | `KIF2C` | +0.1013 | G2/M |
| 23 | `SMC4` | +0.1010 | G2/M |
| 24 | `CDKN3` | +0.1005 | — |
| 25 | `MKI67` | +0.1003 | G2/M |
| 26 | `FAM64A` | +0.0995 | G2/M |
| 27 | `NEK2` | +0.0991 | — |
| 28 | `CDCA8` | +0.0990 | G2/M |
| 29 | `DLGAP5` | +0.0990 | G2/M |
| 30 | `GTSE1` | +0.0988 | G2/M |
| 31 | `UBE2S` | +0.0985 | — |
| 32 | `CENPA` | +0.0976 | — |
| 33 | `MAD2L1` | +0.0975 | — |
| 34 | `H2AFZ` | +0.0962 | — |
| 35 | `NUF2` | +0.0940 | G2/M |
| 36 | `HMGB1` | +0.0937 | — |
| 37 | `PSRC1` | +0.0929 | G2/M |
| 38 | `CKAP2` | +0.0928 | G2/M |
| 39 | `SGOL2` | +0.0925 | — |
| 40 | `KNSTRN` | +0.0924 | — |

## Top 40 negative loadings (higher in low-efficiency lines)

| rank | gene | loading | canonical set |
|---:|---|---:|---|
| 1 | `RPL12` | -0.0488 | — |
| 2 | `SNHG8` | -0.0486 | — |
| 3 | `RPL10` | -0.0461 | — |
| 4 | `RPS3` | -0.0406 | — |
| 5 | `RPS12` | -0.0399 | — |
| 6 | `EPB41L4A-AS1` | -0.0390 | — |
| 7 | `EIF3E` | -0.0387 | — |
| 8 | `SH3BGRL3` | -0.0353 | — |
| 9 | `RPS28` | -0.0340 | — |
| 10 | `ETFB` | -0.0329 | — |
| 11 | `MIAT` | -0.0325 | — |
| 12 | `APOE` | -0.0320 | — |
| 13 | `GAPDH` | -0.0319 | — |
| 14 | `MFNG` | -0.0318 | — |
| 15 | `UQCRB` | -0.0314 | — |
| 16 | `CCND2` | -0.0309 | — |
| 17 | `RPL37A` | -0.0298 | — |
| 18 | `SLC2A1` | -0.0297 | — |
| 19 | `COX7C` | -0.0296 | — |
| 20 | `PFDN5` | -0.0289 | — |
| 21 | `PCBP4` | -0.0288 | — |
| 22 | `TOMM7` | -0.0280 | — |
| 23 | `EIF3L` | -0.0267 | — |
| 24 | `QPRT` | -0.0263 | — |
| 25 | `EIF4A2` | -0.0260 | — |
| 26 | `SRCAP` | -0.0260 | — |
| 27 | `FTL` | -0.0260 | — |
| 28 | `RPL39` | -0.0257 | — |
| 29 | `BST2` | -0.0255 | — |
| 30 | `GPI` | -0.0243 | — |
| 31 | `UQCRH` | -0.0242 | — |
| 32 | `ENO2` | -0.0237 | — |
| 33 | `PDLIM4` | -0.0231 | — |
| 34 | `RGS16` | -0.0228 | — |
| 35 | `RAB3A` | -0.0226 | — |
| 36 | `CRABP2` | -0.0226 | — |
| 37 | `FAM57B` | -0.0218 | — |
| 38 | `TNNT1` | -0.0218 | — |
| 39 | `DLL3` | -0.0213 | — |
| 40 | `UBE2L6` | -0.0208 | — |

("canonical set" = membership in the Tirosh/Seurat `cc.genes` G2/M or G1/S
list. Symbols are as they appear in the source matrix; a few Seurat-list
genes use older symbols — `FAM64A` = *PIMREG*, `HN1` = *JPT1*, `MLF1IP` =
*CENPU* — so set membership is matched on the symbols present in the data.)

---

## Caveats

1. **PC1 is moderately pool-associated (eta^2 = 0.49).** Roughly half its
   line-level variance is explained by which pool a line came from.
2. **A proliferation-vs-ribosomal contrast is also the classic shape of a
   technical axis.** Library size, capture efficiency and cell viability
   all push ribosomal/housekeeping fraction one way and complexity the
   other. The biological reading is plausible and supported by the
   `phat_P_FPP` correlation, but PC1 is **not** as technically clean as
   `phat_NB` (eta^2 = 0.04, the only "Low" feature in the importance
   tables). Treat the proliferation story as a strong hypothesis, not a
   settled result.
3. **Cell-cycle signal is frequently regressed out** in scRNA-seq
   pipelines precisely because it dominates variance. It was deliberately
   *not* regressed out here — which is why PC1 looks like this, and is
   worth stating as a pipeline choice rather than a discovery.
4. **Correlation, not causation.** Nothing here shows that increasing D11
   proliferation would improve D52 yield.
5. **Loadings describe the axis, not any single line.** A line's PC1 value
   is a pool-then-line average of per-cell projections.

---

## Other components

`metadata_eda/d11_pca_gene_loadings.csv` carries all 10 PCs. Two cautions before
interpreting them:

- **PC2 and PC3 are near-tied in explained variance** (0.0265 / 0.0259) and
  **reorder between the two bases**: baseline PC2 corresponds to
  qualifying-only PC3 (\|r\| = 0.985), while qualifying-only PC2 is a
  different component entirely (\|r\| = 0.564 vs baseline PC2). Align
  components by correlation before interpreting.
- **PC6, PC8 and PC9 are close to pure batch signal** (pool eta^2 = 0.91 /
  0.96 / 0.95). Their loadings should be read as technical structure.

---

## Sources

Gene-function calls above rest on these; all were checked while writing
this file.

1. **Tirosh I, Izar B, Prakadan SM, et al. (2016)** "Dissecting the
   multicellular ecosystem of metastatic melanoma by single-cell RNA-seq."
   *Science* 352(6282):189–196.
   <https://www.science.org/doi/10.1126/science.aad0501>
   — origin of the G1/S and G2/M gene sets used here for enrichment.
2. **Seurat, `cc.genes` reference.**
   <https://satijalab.org/seurat/reference/cc.genes> and the cell-cycle
   vignette <https://satijalab.org/seurat/archive/v3.0/cell_cycle_vignette.html>
   — the exact gene lists matched against, and the standard scoring
   workflow.
3. **Whitfield ML, Sherlock G, Saldanha AJ, et al. (2002)**
   "Identification of genes periodically expressed in the human cell cycle
   and their expression in tumors." *Molecular Biology of the Cell*
   13(6):1977–2000. doi:10.1091/mbc.02-02-0030
   <https://pubmed.ncbi.nlm.nih.gov/12058064/>
   — the foundational human cell-cycle periodicity map (>850 periodic
   transcripts) establishing the G2/M identity of *TOP2A, CCNB1, CDC20,
   PLK1, AURKA, CENPF, BIRC5* and others listed above.
4. **Luecken MD, Theis FJ (2019)** "Current best practices in single-cell
   RNA-seq analysis: a tutorial." *Molecular Systems Biology* 15(6):e8746.
   <https://doi.org/10.15252/msb.20188746>
   (tutorial code: <https://github.com/theislab/single-cell-tutorial>)
   — quality-control covariates and why ribosomal/mitochondrial fractions
   and cell-cycle signal are treated as technical confounders; basis for
   caveats 2 and 3.
5. **Jerber J, Seaton DD, Cuomo ASE, et al. (2021)**
   "Population-scale single-cell RNA-seq profiling across dopaminergic
   neuron differentiation." *Nature Genetics* 53:304–312.
   <https://pmc.ncbi.nlm.nih.gov/articles/PMC7610897/>
   — source of the data, the D11/D30/D52 design, the cell-type labels
   (`FPP`, `NB`, `P_FPP`) and the 0.2 efficiency threshold.

Individual gene identities (*MKI67* proliferation marker; *AURKA/AURKB*
Aurora kinases; *PLK1* polo-like kinase 1; *CCNB1/CCNB2* B-type cyclins;
*CDK1*; *CDC20* APC/C activator; *BIRC5* survivin; *TPX2* spindle
assembly; *TOP2A* topoisomerase II-alpha; *PTTG1* securin; *SMC4*
condensin; *CKS1B/CKS2* CDK regulatory subunits) are standard textbook
assignments consistent with sources 1 and 3.
