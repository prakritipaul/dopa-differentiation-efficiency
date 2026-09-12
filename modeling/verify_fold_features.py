"""
Structural verification of a completed per-fold feature table.

A hard gate between an extraction finishing and anything being run on its
output. Extraction takes ~3 hours and its failure mode does not raise -- it
emits a plausible table that is short a fold, or missing a line, or carrying a
stale provenance value. Every expectation is derived from the variant's own
cohort file rather than hardcoded, so this checks the table that was asked for
rather than one particular historical shape.

Exits non-zero on failure, so it can gate a shell pipeline.

Usage:
    uv run python -m modeling.verify_fold_features D11 --label-variant da_untreated
"""

import argparse
import sys

import pandas as pd

from modeling.features import DEPTH_OUTLIER_POOLS
from modeling.folds import VARIANTS, get_variant

OUT_DIR = __file__.rsplit("/", 1)[0]


def verify(timepoint: str, label_variant: str = "da_untreated", suffix: str = "_qualonly") -> bool:
    v = get_variant(label_variant)
    path = f"{OUT_DIR}/fold_data/fold_features_{timepoint}{suffix}{v.suffix}.csv"
    cohort = pd.read_csv(v.cohort_csv)
    n_combos, n_lines = len(cohort), cohort["cell_line"].nunique()

    d = pd.read_csv(path)
    state = {"ok": True}

    def check(label, cond, detail=""):
        state["ok"] &= bool(cond)
        print(f"  [{'PASS' if cond else 'FAIL'}] {label}{('  ' + detail) if detail else ''}")

    print(f"{path}\n  {len(d)} rows, {len(d.columns)} cols "
          f"(cohort: {n_combos} combos, {n_lines} lines)\n")

    folds = d[["scheme", "repeat", "fold"]].drop_duplicates()
    counts = folds.groupby("scheme").size().to_dict()
    # LOCO is one fold per line, so the total follows from the cohort.
    expected = {"plain": 50, "donor_grouped": 50, "loco": n_lines, "lodo": 20}
    check(f"{sum(expected.values())} folds", len(folds) == sum(expected.values()), str(len(folds)))
    check("scheme counts", counts == expected, str(counts))
    check("label_variant uniform and correct", set(d.label_variant) == {label_variant},
          str(set(d.label_variant)))
    check("fit_population uniform", set(d.fit_population) == {"qualifying_only"},
          str(set(d.fit_population)))
    check("no NaN", not d.isna().any().any())

    sizes = d.groupby(["scheme", "repeat", "fold"]).size()
    check(f"{n_combos} combos in EVERY fold", sizes.eq(n_combos).all(), str(sorted(sizes.unique())))
    per_fold_lines = d.groupby(["scheme", "repeat", "fold"])["cell_line"].nunique()
    check(f"{n_lines} lines in EVERY fold", per_fold_lines.eq(n_lines).all(),
          str(sorted(per_fold_lines.unique())))
    check("no duplicate (fold, line, pool)",
          not d.duplicated(["scheme", "repeat", "fold", "cell_line", "pool"]).any())

    check("split only train/test", set(d.split) == {"train", "test"})
    dup = d.groupby(["scheme", "repeat", "fold", "cell_line"])["split"].nunique()
    check("no line both train and test in a fold", dup.eq(1).all(), f"max {dup.max()}")

    loco = d[d.scheme == "loco"]
    lt = loco[loco.split == "test"].groupby(["repeat", "fold"])["cell_line"].nunique()
    check("LOCO holds out exactly 1 line", lt.eq(1).all(), str(sorted(lt.unique())))
    check("LOCO covers every line once",
          loco[loco.split == "test"]["cell_line"].nunique() == n_lines)

    phat = sorted(c for c in d.columns if c.startswith("phat_"))
    rs = d[phat].sum(axis=1)
    check(f"{len(phat)} proportions sum to 1.0", rs.sub(1.0).abs().max() < 1e-9,
          f"max dev {rs.sub(1.0).abs().max():.2e}")
    check("PC1..PC10 present", all(f"PC{i}" in d.columns for i in range(1, 11)))

    # The excluded pool is read from the same constant the extraction used --
    # a second copy restated here could disagree with the real rule.
    (excluded,) = DEPTH_OUTLIER_POOLS[timepoint]
    check(f"{excluded} still PRESENT (projected, excluded only from fitting)",
          (d["pool"] == excluded).any(), f"{int((d['pool'] == excluded).sum())} rows")
    check("all cohort pools present", d["pool"].nunique() == cohort["pool"].nunique(),
          f"{d['pool'].nunique()} vs {cohort['pool'].nunique()}")

    print(f"\n{'ALL CHECKS PASSED' if state['ok'] else 'SOME CHECKS FAILED'}")
    return state["ok"]


def paired_with(other_timepoint: str, timepoint: str, label_variant: str = "da_untreated",
                suffix: str = "_qualonly") -> bool:
    """Do two timepoints' FEATURE TABLES agree on per-fold train/test membership?

    Stronger than comparing the fold-assignment file: this checks the artifacts
    the models actually consume. A message printed at extraction time is not a
    verification artifact."""
    v = get_variant(label_variant)
    key = ["scheme", "repeat", "fold", "cell_line", "split"]
    a = pd.read_csv(f"{OUT_DIR}/fold_data/fold_features_{timepoint}{suffix}{v.suffix}.csv",
                    usecols=key).drop_duplicates().sort_values(key).reset_index(drop=True)
    b = pd.read_csv(f"{OUT_DIR}/fold_data/fold_features_{other_timepoint}{suffix}{v.suffix}.csv",
                    usecols=key).drop_duplicates().sort_values(key).reset_index(drop=True)
    same = a.equals(b)
    print(f"  [{'PASS' if same else 'FAIL'}] {timepoint} and {other_timepoint} feature tables "
          f"agree on every fold's train/test membership ({len(a)} rows compared)")
    return same


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("timepoint", choices=["D11", "D30"])
    parser.add_argument("--label-variant", default="da_untreated", choices=sorted(VARIANTS))
    parser.add_argument("--suffix", default="_qualonly", help="fit-population suffix of the table")
    parser.add_argument("--paired-with", default=None, choices=["D11", "D30"],
                        help="also check both tables agree on fold membership")
    args = parser.parse_args()

    ok = verify(args.timepoint, args.label_variant, args.suffix)
    if args.paired_with:
        ok &= paired_with(args.paired_with, args.timepoint, args.label_variant, args.suffix)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
