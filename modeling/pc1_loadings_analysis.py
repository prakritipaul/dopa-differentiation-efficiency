"""
What is PC1, at D11 and at D30?

PC1 is the one principal component worth interpreting biologically. The others
are not: PC2 and PC3 are near-tied in variance and ROTATE between PCA fits, so
"PC2" names a slot rather than a fixed axis. PC1 replicates across fits at
r = 0.9995, which is what makes a gene-level reading of it meaningful.

This does three things, for each timepoint, under the da_untreated variant:

  1. Top loadings at each end of the axis.
  2. Enrichment of the Tirosh/Seurat cell-cycle sets among those loadings,
     by rank -- the same test the published D11 interpretation used, so the
     two are comparable.
  3. What PC1 tracks at the LINE level: its correlation with each cell-type
     proportion. A gene list alone is suggestive; agreement with an
     independent annotation is what makes it an interpretation.

SIGN IS ARBITRARY. PCA fixes an axis, not a direction, so "positive loading"
has no meaning on its own. Everything below is reported as one pole against
the other, and the poles are labelled by their own gene content.

Usage:
    uv run python -m modeling.pc1_loadings_analysis
Writes modeling/results/pc1_loadings_{D11,D30}_da_untreated.csv
"""

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr

from modeling.features import load_cell_metadata
from modeling.folds import get_variant, load_lines_with_label
from modeling.harness import proportion_cols

REPO_ROOT = __file__.rsplit("/", 2)[0]
PCA_DIR = f"{REPO_ROOT}/metadata_eda/pca"
OUT_DIR = f"{REPO_ROOT}/modeling/results"
N_TOP = 25

# Tirosh et al. 2016 / Seurat cc.genes, the same lists the published D11
# interpretation matched against. Kept verbatim rather than trimmed to the
# HVGs, so the enrichment test is against the standard set.
G2M = """HMGB2 CDK1 NUSAP1 UBE2C BIRC5 TPX2 TOP2A NDC80 CKS2 NUF2 CKS1B MKI67
TMPO CENPF TACC3 SMC4 CCNB2 CKAP2L CKAP2 AURKB BUB1 KIF11 ANP32E TUBB4B GTSE1
KIF20B HJURP CDCA3 CDC20 TTK CDC25C KIF2C RANGAP1 NCAPD2 DLGAP5 CDCA2 CDCA8
ECT2 KIF23 HMMR AURKA PSRC1 ANLN LBR CKAP5 CENPE CTCF NEK2 G2E3 GAS2L3 CBX5
CENPA""".split()
S_PHASE = """MCM5 PCNA TYMS FEN1 MCM2 MCM4 RRM1 UNG GINS2 MCM6 CDCA7 DTL PRIM1
UHRF1 CENPU HELLS RFC2 RPA2 NASP RAD51AP1 GMNN WDR76 SLBP CCNE2 UBR7 POLD3 MSH2
ATAD2 RAD51 RRM2 CDC45 CDC6 EXO1 TIPIN DSCC1 BLM CASP8AP2 USP1 CLSPN POLA1
CHAF1B BRIP1 E2F8""".split()

# A priori neuronal / progenitor markers, for reading the other pole.
NEURONAL = "STMN2 DCX MAP2 ELAVL4 SYT1 NEFL NEFM TUBB3 GAP43 MLLT11 NNAT".split()
PROGENITOR = "SOX2 HES1 VIM NES FABP7 PTPRZ1 SLC1A3 TTYH1 CD9 ID4".split()


def _enrichment(loadings: pd.DataFrame, genes: list[str], label: str) -> dict:
    """Rank-based: where do these genes sit among all HVGs by PC1 loading?

    Reported on the SIGNED loading, so a set concentrated at one pole shows up
    as a rank shift rather than being washed out by taking absolute values."""
    present = [g for g in genes if g in set(loadings.gene)]
    if len(present) < 5:
        return {"set": label, "n_in_hvgs": len(present), "note": "too few to test"}
    ranks = loadings.set_index("gene")["signed_rank"]
    hit = ranks.loc[present]
    rest = ranks.drop(present)
    u = mannwhitneyu(hit, rest, alternative="two-sided")
    return {
        "set": label,
        "n_in_hvgs": len(present),
        "median_rank": float(hit.median()),
        "median_rank_all": float(rest.median()),
        "p": float(u.pvalue),
    }


def analyse(timepoint: str, variant: str = "da_untreated") -> pd.DataFrame:
    v = get_variant(variant)
    prefix = {"D11": "d11", "D30": "d30"}[timepoint]
    load = pd.read_csv(f"{PCA_DIR}/{prefix}_pca_gene_loadings_qualonly{v.suffix}.csv")
    pc1 = load[load.PC == "PC1"].copy()
    pc1["signed_rank"] = pc1["loading"].rank()

    var = pd.read_csv(f"{PCA_DIR}/{prefix}_pca_variance_explained_qualonly{v.suffix}.csv")
    pct = float(var.loc[var.PC == "PC1", "explained_variance_ratio"].iloc[0]) * 100

    pos = pc1.nlargest(N_TOP, "loading")
    neg = pc1.nsmallest(N_TOP, "loading")

    print(f"\n{'=' * 78}\n{timepoint} PC1 — {pct:.1f}% of HVG variance, {len(pc1)} genes\n{'=' * 78}")
    print(f"\n  POSITIVE pole (top {N_TOP}):\n    " +
          "  ".join(pos.gene.tolist()[:N_TOP]))
    print(f"\n  NEGATIVE pole (top {N_TOP}):\n    " +
          "  ".join(neg.gene.tolist()[:N_TOP]))

    print("\n  Gene-set position along the axis (Mann-Whitney on signed rank):")
    for genes, label in [(G2M, "cell cycle G2/M"), (S_PHASE, "cell cycle S"),
                         (NEURONAL, "pan-neuronal"), (PROGENITOR, "progenitor/glial")]:
        e = _enrichment(pc1, genes, label)
        if "p" in e:
            side = "positive" if e["median_rank"] > e["median_rank_all"] else "negative"
            print(f"    {label:18s} n={e['n_in_hvgs']:3d}  median rank {e['median_rank']:7.1f} "
                  f"vs {e['median_rank_all']:7.1f} ({side} pole)  p={e['p']:.2e}")
        else:
            print(f"    {label:18s} {e['note']} (n={e['n_in_hvgs']})")

    # What does PC1 track at the line level?
    coords = pd.read_csv(f"{PCA_DIR}/{prefix}_pca_coords_per_line_qualonly{v.suffix}.csv")
    meta = load_cell_metadata(timepoint)
    qual = pd.read_csv(v.cohort_csv)[["cell_line", "pool"]]
    from modeling.features import compute_proportion_features, pool_then_line_average
    props = compute_proportion_features(meta.merge(qual, on=["cell_line", "pool"], how="inner"))
    phat = [c for c in props.columns if c.startswith("phat_")]
    line_props = pool_then_line_average(props, phat)

    lines = load_lines_with_label(variant)
    df = coords.merge(line_props, on="cell_line").merge(
        lines[["cell_line", "diff_efficiency"]], on="cell_line")
    print("\n  Line-level PC1 vs cell-type proportions (Spearman):")
    for c in sorted(phat):
        print(f"    {c:16s} {spearmanr(df['PC1'], df[c]).statistic:+.3f}")
    print(f"    {'diff_efficiency':16s} {spearmanr(df['PC1'], df['diff_efficiency']).statistic:+.3f}")

    out = pd.concat([pos.assign(pole="positive"), neg.assign(pole="negative")])
    out.to_csv(f"{OUT_DIR}/pc1_loadings_{timepoint}_{variant}.csv", index=False)
    return out


def decompose(timepoint: str, variant: str = "da_untreated") -> None:
    """Is line-level PC1 driven by COMPOSITION or by within-type state?

    A line's mean PC1 is a weighted average over its cell types:
        line_PC1 ~= sum_k  p_k * mean_PC1_within_type_k
    Either term can vary. Replacing each type's per-line score with a COMMON
    reference mean isolates the composition term; holding composition fixed at
    the cohort mean isolates the within-type term.

    This matters because a gene list alone cannot distinguish "lines differ in
    how mature their neurons are" from "lines differ in how many neurons they
    have" -- and those support very different claims."""
    v = get_variant(variant)
    prefix = {"D11": "d11", "D30": "d30"}[timepoint]
    cells = pd.read_csv(
        f"{PCA_DIR}/{prefix}_pca_coords_per_cell_qualifying_qualonly{v.suffix}.csv",
        usecols=["cell_line", "pool", "PC1"])
    meta = load_cell_metadata(timepoint)
    qual = pd.read_csv(v.cohort_csv)[["cell_line", "pool"]]
    meta_q = meta.merge(qual, on=["cell_line", "pool"], how="inner").reset_index(drop=True)
    if len(meta_q) != len(cells):
        raise ValueError(f"cell rows disagree: {len(meta_q)} metadata vs {len(cells)} coords")
    cells = cells.assign(celltype=meta_q["celltype"].values)

    print(f"\n  PC1 WITHIN each cell type (mean over all cells):")
    within = cells.groupby("celltype", observed=True)["PC1"].agg(["count", "mean", "std"])
    print(within.round(2).sort_values("mean", ascending=False).to_string())

    frac = (cells.groupby(["cell_line", "celltype"], observed=True).size()
            .unstack(fill_value=0).pipe(lambda d: d.div(d.sum(axis=1), axis=0)))
    lt = cells.groupby(["cell_line", "celltype"], observed=True)["PC1"].mean().unstack()
    actual = cells.groupby("cell_line", observed=True)["PC1"].mean()
    ref = cells.groupby("celltype", observed=True)["PC1"].mean()

    comp_only = (frac * ref).sum(axis=1)
    types = [c for c in frac.columns if c in lt.columns]
    within_only = (frac.mean() * lt[types]).sum(axis=1)
    ok = within_only.notna()

    r2 = lambda a, b: float(np.corrcoef(a, b)[0, 1] ** 2)
    print(f"\n  composition alone explains  R² {r2(actual, comp_only):.3f} of line-level PC1")
    print(f"  within-type state alone      R² {r2(actual[ok], within_only[ok]):.3f}")

    lines = load_lines_with_label(variant).set_index("cell_line")["diff_efficiency"]
    j = pd.DataFrame({"actual": actual, "comp": comp_only}).join(lines).dropna()
    print(f"\n  actual PC1 vs outcome            {spearmanr(j.actual, j.diff_efficiency).statistic:+.3f}")
    print(f"  composition-only PC1 vs outcome  {spearmanr(j.comp, j.diff_efficiency).statistic:+.3f}")

    print(f"\n  cell-type fractions across lines (is a ~0 correlation informative?):")
    print(frac.agg(["mean", "std", "min", "max"]).round(3).to_string())


def main() -> None:
    for tp in ["D11", "D30"]:
        analyse(tp)
    print(f"\n{'=' * 78}\nD30: composition vs within-type state\n{'=' * 78}")
    decompose("D30")
    print(f"\nSaved top loadings to {OUT_DIR}/pc1_loadings_*.csv")


if __name__ == "__main__":
    main()
