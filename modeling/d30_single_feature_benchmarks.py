"""
Fit-free single-column D30 benchmarks, with donor-level uncertainty.

Two aggregates of the D30 cell-type composition are reference lines that any
fitted D30 model has to be judged against:

  DA + Sert    the "already arrived" fraction. These are the very cell types
               whose D52 fraction IS the outcome, so this is not a discovered
               predictor -- it is the outcome measured 22 days early. It sets
               the CEILING: if the full model cannot beat this, the modelling
               added nothing.

  FPP + P_FPP  the "still undecided" fraction -- SOX2+/HES1+/VIM-high
               progenitors that are not yet neuronal (P_FPP still cycling:
               MKI67, TOP2A). These cells have not committed, so this is a
               genuine forecast in the same sense the D11 model is, and it is
               the honest analogue of D11's phat_NB result.

Neither has a fitted parameter, so there is nothing to cross-validate: the
prediction for a line is just its own feature value, identical whether or not
that line was "in training". Cross-validation would therefore report a
CV-shaped number that is really in-sample, which is worse than not reporting
one. Uncertainty is instead a DONOR-level bootstrap: 136 of 138 lines share a
donor with another line, so resampling lines would treat correlated lines as
independent and give intervals that are too narrow. Donors are the
independent unit, matching the donor-grouped CV used for the fitted models.

These numbers are in-sample over the 138 study lines and are reported as
reference lines, not as validated model performance.

Usage:
    uv run python -m modeling.d30_single_feature_benchmarks
Writes modeling/results/d30_single_feature_benchmarks.csv
"""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score

from modeling.features import compute_proportion_features, load_cell_metadata, pool_then_line_average
from modeling.folds import VARIANTS, get_variant, load_lines_with_label
from modeling.make_d30_variant_tables import OFF_TARGET, PROGENITOR, TARGET_LIKE

OUT_DIR = Path(__file__).parent
N_BOOT = 2000
SEED = 0

# (columns, expected direction). The DIRECTION IS DECLARED HERE, A PRIORI, on
# biology -- not read off the result. More target-like cells already made by
# D30 means a better D52 yield (+1); more cells still sitting undecided in the
# progenitor pool means fewer of them have become DA/Sert (-1).
#
# This matters because ROC-AUC is direction-sensitive: an unoriented feature
# that perfectly ranks failures above successes scores near 0.0, not near 1.0.
# Flipping a sign after seeing which way the number came out would be choosing
# the aggregation from the data, which is exactly what must not happen outside
# a training fold. Fixing it in the source makes the choice auditable.
#
# Each entry is (numerator columns, denominator columns or None, direction).
# A denominator makes the benchmark a SHARE rather than a fraction-of-all-cells,
# which is the only way to ask a question that is not just "how far along is
# this culture?". The first two both measure position on that single
# maturation axis and are near-mirror images of each other on the simplex --
# their mutual correlation is reported below so that is not mistaken for two
# independent findings.
BENCHMARKS = {
    # MATCH THE NUMERATOR TO THE OUTCOME. "Already arrived" only means what it
    # says when the cells counted here are the cells counted in the outcome.
    # Against the published (DA+Sert)/all outcome that is phat_DA + phat_Sert;
    # against the DA/all outcome it is phat_DA ALONE, and using the pair
    # instead understates the ceiling -- D30 phat_DA tracks D52 DA/all at
    # Spearman +0.93 while phat_DA+phat_Sert manages only +0.80, because Sert
    # is a different lineage that does not become DA. Both are reported so the
    # mismatched one cannot be quoted by accident.
    "DA alone (already arrived, matches a DA-only outcome)": (["phat_DA"], None, +1),
    "DA+Sert (already arrived, matches a DA+Sert outcome)": (TARGET_LIKE, None, +1),
    "FPP+P_FPP (still undecided)": (PROGENITOR, None, -1),
    # The progenitor-balance question: OF THE CELLS THAT HAVE NOT BECOME
    # TARGET-LIKE, what share are still competent progenitors rather than
    # committed off-target fates? Dividing by the undecided mass removes the
    # overall maturation magnitude, so this is not just a restatement of
    # DA+Sert. Direction +1 a priori: a remaining pool that is still
    # progenitor can yet become DA/Sert; one that is Epen1/U_Neur cannot.
    "FPP+P_FPP share of non-target (progenitor balance)": (
        PROGENITOR, PROGENITOR + OFF_TARGET, +1,
    ),
}


def d30_line_level_proportions(variant: str = "published") -> pd.DataFrame:
    """One row per qualifying cell line: every D30 phat_* column, aggregated
    pool-then-line -- the same two-stage averaging the fitted pipeline uses,
    so these benchmarks and the models see identical feature values.

    Recomputed from the h5 rather than read off a fold-features table because
    proportions have no fitted parameter and so do not vary by fold; taking
    them from one arbitrary fold would imply a fold-dependence that does not
    exist."""
    meta = load_cell_metadata("D30")
    qualifying = pd.read_csv(get_variant(variant).cohort_csv)[["cell_line", "pool"]]
    meta_q = meta.merge(qualifying, on=["cell_line", "pool"], how="inner")
    props = compute_proportion_features(meta_q)
    value_cols = [c for c in props.columns if c.startswith("phat_")]
    return pool_then_line_average(props, value_cols)


def _metrics(score: np.ndarray, efficiency: np.ndarray, success: np.ndarray, direction: int) -> dict:
    """Spearman is reported on the RAW feature, so its sign still shows which
    way the biology runs. ROC-AUC is reported on the feature oriented by the
    a-priori `direction`, because an AUC below 0.5 measures discrimination in
    the opposite direction and is not comparable to a model's AUC otherwise.

    ROC-AUC needs both classes present; a bootstrap resample of 20 donors can
    occasionally draw only successes, and scoring that would raise. Such a
    resample is undefined rather than zero, so it yields NaN and is dropped
    from the percentile, with the count reported."""
    out = {"spearman": spearmanr(score, efficiency).statistic}
    oriented = direction * score
    out["roc_auc"] = roc_auc_score(success, oriented) if len(np.unique(success)) == 2 else np.nan
    return out


def bootstrap_by_donor(
    lines: pd.DataFrame, score_col: str, direction: int, n_boot: int = N_BOOT, seed: int = SEED
) -> dict:
    """Resample DONORS with replacement (not lines): lines from one donor are
    not independent, so a line-level bootstrap understates the interval."""
    rng = np.random.default_rng(seed)
    donors = lines["donor"].unique()
    by_donor = {d: g for d, g in lines.groupby("donor", observed=True)}

    stats = {"spearman": [], "roc_auc": []}
    for _ in range(n_boot):
        drawn = rng.choice(donors, size=len(donors), replace=True)
        sample = pd.concat([by_donor[d] for d in drawn], ignore_index=True)
        m = _metrics(sample[score_col].to_numpy(), sample["diff_efficiency"].to_numpy(),
                     sample["success"].to_numpy(), direction)
        for k, v in m.items():
            stats[k].append(v)

    out = {}
    for k, vals in stats.items():
        arr = np.asarray(vals, dtype=float)
        ok = arr[~np.isnan(arr)]
        out[f"{k}_lo"], out[f"{k}_hi"] = np.percentile(ok, [2.5, 97.5])
        out[f"{k}_n_boot_used"] = len(ok)
    return out


def main(label_variant: str = "published") -> None:
    v = get_variant(label_variant)
    props = d30_line_level_proportions(label_variant)
    lines = load_lines_with_label(label_variant).merge(
        props, on="cell_line", how="inner", validate="one_to_one")
    if len(lines) != v.n_lines:
        raise ValueError(f"expected {v.n_lines} lines, got {len(lines)}")

    rows = []
    for name, (cols, denom_cols, direction) in BENCHMARKS.items():
        missing = [c for c in cols + (denom_cols or []) if c not in lines.columns]
        if missing:
            raise ValueError(f"{name}: {missing} not in the D30 proportion table")
        value = lines[cols].sum(axis=1)
        if denom_cols is not None:
            denom = lines[denom_cols].sum(axis=1)
            if (denom <= 0).any():
                raise ValueError(f"{name}: {int((denom <= 0).sum())} lines have an empty denominator")
            value = value / denom
        lines[name] = value
        point = _metrics(lines[name].to_numpy(), lines["diff_efficiency"].to_numpy(),
                         lines["success"].to_numpy(), direction)
        ci = bootstrap_by_donor(lines, name, direction)
        rows.append({"benchmark": name, "columns": "+".join(cols),
                     "denominator": "+".join(denom_cols) if denom_cols else "",
                     "direction": direction, "mean_fraction": lines[name].mean(), **point, **ci})

    out = pd.DataFrame(rows)
    out_path = OUT_DIR / f"results/d30_single_feature_benchmarks{v.suffix}.csv"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, index=False)

    print(f"{len(lines)} lines, {lines['donor'].nunique()} donors, "
          f"{int(lines['success'].sum())} success / {int((~lines['success']).sum())} failure")
    print(f"donor-level bootstrap, {N_BOOT} resamples, seed {SEED}\n")
    for r in rows:
        sign = "higher = better" if r["direction"] > 0 else "higher = WORSE (declared a priori)"
        print(f"{r['benchmark']}   (mean fraction {r['mean_fraction']:.3f}; {sign})")
        print(f"   Spearman rho  {r['spearman']:+.3f}   95% CI [{r['spearman_lo']:+.3f}, {r['spearman_hi']:+.3f}]"
              f"   (raw feature, sign meaningful)")
        print(f"   ROC-AUC       {r['roc_auc']:.3f}    95% CI [{r['roc_auc_lo']:.3f}, {r['roc_auc_hi']:.3f}]"
              f"   (oriented; {r['roc_auc_n_boot_used']}/{N_BOOT} resamples usable)")
    # Two benchmarks scoring well is not two independent findings if they are
    # the same axis with the sign flipped. Printed so the reader cannot miss it.
    names = list(BENCHMARKS)
    print("\nmutual Spearman between benchmarks (shared maturation axis, not independent evidence):")
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            print(f"   {a}  vs  {b}:  {spearmanr(lines[a], lines[b]).statistic:+.3f}")
    print(f"\nSaved {out_path}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--label-variant", default="published", choices=sorted(VARIANTS))
    main(parser.parse_args().label_variant)
