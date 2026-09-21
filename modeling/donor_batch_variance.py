"""How much of D52 dopaminergic yield is the donor, and how much is the batch?

FINDINGS reports pool eta^2 for the PREDICTORS, but never for the OUTCOME
itself. That leaves the question a cell-therapy reader asks first unanswered:
if a line differentiates badly, is that the line, or the run?

Three measurements, in increasing order of how much they settle:

  1. eta^2 of D52 efficiency by donor and by pool. Both are reported against a
     PERMUTATION NULL, because eta^2 rises with group count for free -- 20
     donors over 136 lines scores well above 0 on shuffled labels, so the raw
     number alone means nothing.

  2. Technical replicates: 21 lines were differentiated in more than one pool.
     Their agreement across pools is a direct reproducibility measure, and it
     needs no model.

  3. Donor-grouped vs plain CV, already in the results tables -- what donor
     leakage is actually worth in predictive terms.

Usage:
    uv run python -m modeling.donor_batch_variance
"""

import h5py
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

import sys
sys.path.insert(0, __file__.rsplit("/", 2)[0])
from data_paths import data_file
from modeling.folds import get_variant, load_lines_with_label

sys.path.insert(0, __file__.rsplit("/", 2)[0] + "/eda")
import importlib
m004 = importlib.import_module("004_d11_celltype_proportion_se")

RNG = np.random.default_rng(0)
N_PERM = 2000


def eta_sq(values: np.ndarray, groups: np.ndarray) -> float:
    """Fraction of variance in `values` explained by `groups` (one-way)."""
    grand = values.mean()
    total = ((values - grand) ** 2).sum()
    if total == 0:
        return float("nan")
    between = sum(
        len(values[groups == g]) * (values[groups == g].mean() - grand) ** 2
        for g in np.unique(groups)
    )
    return float(between / total)


def eta_sq_vs_null(values, groups, label: str) -> None:
    """eta^2 with a permutation null, since eta^2 is inflated by group count."""
    obs = eta_sq(values, groups)
    null = np.array([eta_sq(values, RNG.permutation(groups)) for _ in range(N_PERM)])
    p = (null >= obs).mean()
    print(f"  {label:34s} eta^2 {obs:.3f}   null median {np.median(null):.3f} "
          f"(95th {np.quantile(null, 0.95):.3f})   p = {p:.4f}")


def d52_efficiency_per_combo() -> pd.DataFrame:
    """DA / all cells at D52, untreated only, per (cell_line, pool).

    The published label is aggregated to one value per LINE, which cannot see
    batch effects at all. This keeps the pool axis."""
    with h5py.File(data_file("D52"), "r") as f:
        obs = pd.DataFrame({
            "cell_line": m004.read_obs_categorical(f, "donor_id"),
            "pool": m004.read_obs_categorical(f, "pool_id"),
            "celltype": m004.read_obs_categorical(f, "celltype"),
            "treatment": m004.read_obs_categorical(f, "treatment"),
        })
    obs = obs[obs["treatment"] == "NONE"]
    g = obs.groupby(["cell_line", "pool"], observed=True)
    return pd.DataFrame({
        "n_cells": g.size(),
        "da_over_all": g["celltype"].apply(lambda s: (s == "DA").mean()),
    }).reset_index()


def main() -> None:
    v = get_variant("da_untreated")
    lines = load_lines_with_label("da_untreated")
    cohort = pd.read_csv(v.cohort_csv)

    print("=" * 78)
    print("1. Is D52 dopaminergic yield a property of the donor?")
    print("=" * 78)
    print(f"  {len(lines)} lines, {lines.donor.nunique()} donors, "
          f"{len(lines) / lines.donor.nunique():.1f} lines per donor\n")
    eta_sq_vs_null(lines.diff_efficiency.values, lines.donor.values, "D52 efficiency by DONOR")

    combo = d52_efficiency_per_combo().merge(
        cohort[["cell_line", "pool", "donor"]], on=["cell_line", "pool"], how="inner")
    print(f"\n  per-(line,pool): {len(combo)} combos, {combo.pool.nunique()} pools")
    eta_sq_vs_null(combo.da_over_all.values, combo.pool.values, "D52 efficiency by POOL (batch)")
    # Donor again, but on the SAME rows as pool. The line-level eta^2 above is
    # over 136 points / 20 groups and the pool one over 157 / 10; putting the
    # two numbers side by side is only fair on identical units.
    eta_sq_vs_null(combo.da_over_all.values, combo.donor.values, "D52 efficiency by DONOR (same rows)")
    eta_sq_vs_null(combo.da_over_all.values, combo.cell_line.values, "D52 efficiency by CELL LINE")

    print("\n" + "=" * 78)
    print("2. Technical replicates: the same line, differentiated twice")
    print("=" * 78)
    rep = combo[combo.duplicated("cell_line", keep=False)].sort_values(["cell_line", "pool"])
    n_lines = rep.cell_line.nunique()
    print(f"  {n_lines} lines appear in more than one pool ({len(rep)} combos)\n")
    first = rep.groupby("cell_line", observed=True).nth(0).set_index("cell_line")
    second = rep.groupby("cell_line", observed=True).nth(1).set_index("cell_line")
    both = first[["da_over_all"]].join(second[["da_over_all"]], lsuffix="_a", rsuffix="_b").dropna()
    r = spearmanr(both.da_over_all_a, both.da_over_all_b)
    diff = (both.da_over_all_a - both.da_over_all_b).abs()
    print(f"  Spearman across pools        {r.statistic:+.3f}  (p = {r.pvalue:.1e}, n = {len(both)})")
    print(f"  median |difference|          {diff.median():.3f}")
    print(f"  vs cohort SD of efficiency   {lines.diff_efficiency.std():.3f}")
    print(f"  worst disagreement           {diff.max():.3f}")

    # Is that disagreement a batch effect, or just small-n counting noise?
    # Each efficiency is a proportion over n cells, so it carries binomial
    # error on its own. Comparing the observed spread against that error is
    # what separates "the run mattered" from "we counted too few cells".
    na = first.loc[both.index, "n_cells"].values
    nb = second.loc[both.index, "n_cells"].values
    pa = both.da_over_all_a.values
    pb = both.da_over_all_b.values
    se_pair = np.sqrt(pa * (1 - pa) / na + pb * (1 - pb) / nb)
    z = (pa - pb) / se_pair
    print(f"\n  cells per replicate combo    median {np.median(np.r_[na, nb]):.0f}, "
          f"min {np.r_[na, nb].min():.0f}")
    print(f"  expected |diff| from counting noise alone   "
          f"{np.median(se_pair) * 0.8:.3f}   (observed {diff.median():.3f})")
    print(f"  pairs differing by > 2 binomial SE          "
          f"{(np.abs(z) > 2).sum()} / {len(z)}")

    # The binomial floor assumes cells within a run are independent draws, and
    # cells sharing a well are probably not. Rather than assume that away, ask
    # how much the SE would have to be inflated before the effect disappears.
    for k in (1, 2, 5, 10, 15):
        print(f"    SE inflated x{k:<3d} -> {(np.abs(z / k) > 2).sum():2d}/{len(z)} pairs still beyond 2 SE")

    well = (na >= 100) & (nb >= 100)
    if well.sum() >= 8:
        rw = spearmanr(pa[well], pb[well])
        print(f"  restricted to >=100 cells both pools: rho {rw.statistic:+.3f} "
              f"(p = {rw.pvalue:.2f}, n = {well.sum()})")

    print("\n" + "=" * 78)
    print("3. What donor leakage costs, predictively")
    print("=" * 78)
    for tp in ["D11", "D30"]:
        reg = pd.read_csv(f"modeling/results/results_regression_{tp}_qualonly_da_untreated.csv")
        cls = pd.read_csv(f"modeling/results/results_classification_{tp}_qualonly_da_untreated.csv")
        for frame, model, col, name in [(reg, "ridge", "r2_mean", "R2"),
                                        (cls, "logistic_l2", "roc_auc_mean", "ROC-AUC")]:
            # tuning == "nested" is the headline row. The table also holds every
            # flat-grid config, and those are NOT comparable -- a flat row is the
            # best k chosen with the test fold visible.
            nested = frame[(frame.model == model) & (frame.tuning == "nested")]
            row = lambda s: nested[nested.scheme == s][col].iloc[0]
            print(f"  {tp} {name:8s} plain {row('plain'):.3f} -> donor_grouped "
                  f"{row('donor_grouped'):.3f}  (lodo {row('lodo'):.3f})")


if __name__ == "__main__":
    main()
