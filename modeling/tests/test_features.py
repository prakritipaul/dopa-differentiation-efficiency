"""Tests for modeling/features.py.

test_compute_proportion_features_matches_known_values is a regression
test against already-validated output from the earlier EDA phase
(metadata_eda/d11_celltype_proportions_with_se.csv), not a fresh
derivation.

test_integration_small_extraction is a slow (~3 min) integration test
that streams real D11 single-cell data via
modeling.run_feature_extraction.main(); it writes to a distinct
"_agent_check" suffix so it never touches the live full-extraction run's
output files (modeling/fold_features_D11.csv,
modeling/fold_assignments_first_pass.csv), and deletes its own output
when done.
"""

import math
from pathlib import Path

import pandas as pd
import pytest

from modeling.features import compute_proportion_features, load_cell_metadata

REPO_ROOT = Path(__file__).parent.parent.parent


@pytest.mark.requires_data
def test_compute_proportion_features_matches_known_values():
    meta = load_cell_metadata("D11")
    qualifying = pd.read_csv(REPO_ROOT / "metadata_eda" / "cohort/qualifying_cell_line_pool_min10_per_timepoint.csv")[
        ["cell_line", "pool"]
    ]
    meta_q = meta.merge(qualifying, on=["cell_line", "pool"], how="inner")

    props = compute_proportion_features(meta_q)
    row = props[(props["cell_line"] == "HPSI0114i-eipl_1") & (props["pool"] == "pool1")]
    assert len(row) == 1
    row = row.iloc[0]

    # known-correct reference values from metadata_eda/d11_celltype_proportions_with_se.csv
    assert math.isclose(row["phat_FPP"], 0.536413, abs_tol=1e-5)
    assert math.isclose(row["phat_NB"], 0.073147, abs_tol=1e-5)
    assert math.isclose(row["phat_P_FPP"], 0.390440, abs_tol=1e-5)


@pytest.mark.slow
@pytest.mark.requires_data
def test_integration_small_extraction():
    from modeling.run_feature_extraction import main

    out_suffix = "_agent_check"
    out_dir = REPO_ROOT / "modeling"
    features_csv = out_dir / f"fold_data/fold_features_D11{out_suffix}.csv"
    # Must match run_feature_extraction's actual output name. This previously
    # read "fold_assignments_first_pass{suffix}" -- the name used when this
    # test was written -- so cleanup silently missed the real file after the
    # writer was renamed.
    fold_assignments_csv = out_dir / f"fold_data/fold_assignments{out_suffix}.csv"

    try:
        # include_loco_lodo=False is required for the row count below to mean
        # anything: it defaults to True, which would add 138 LOCO + 20 LODO
        # folds (162 total, not 4). This test asserted 636 rows while
        # requesting the default, so it had been failing since LOCO/LODO were
        # added -- unnoticed because it is marked slow and every run used
        # `-m "not slow"`.
        main(n_splits=2, n_repeats=1, out_suffix=out_suffix, include_loco_lodo=False)

        assert features_csv.exists()
        df = pd.read_csv(features_csv)

        # 2 schemes (plain, donor_grouped) x (2 splits * 1 repeat) folds x 159 qualifying combos
        assert df.shape[0] == 4 * 159, f"expected 636 rows, got {df.shape[0]}"

        assert df.isna().sum().sum() == 0, "unexpected NaNs in extracted features"

        pc_cols = [f"PC{i}" for i in range(1, 11)]
        for col in pc_cols:
            assert col in df.columns
        for col in ("phat_FPP", "phat_NB", "phat_P_FPP"):
            assert col in df.columns

        # Independent re-derivation of the donor-grouping constraint directly
        # from the persisted output (not reusing folds.py's own Fold objects).
        qualifying = pd.read_csv(
            REPO_ROOT / "metadata_eda" / "cohort/qualifying_cell_line_pool_min10_per_timepoint.csv"
        )[["cell_line", "donor"]].drop_duplicates()
        dg = df[df["scheme"] == "donor_grouped"].merge(qualifying, on="cell_line", how="left")
        assert dg["donor"].isna().sum() == 0, "cell_line -> donor mapping failed for some rows"

        violations = []
        for (repeat, fold), g in dg.groupby(["repeat", "fold"]):
            train_donors = set(g.loc[g["split"] == "train", "donor"])
            test_donors = set(g.loc[g["split"] == "test", "donor"])
            overlap = train_donors & test_donors
            if overlap:
                violations.append((repeat, fold, overlap))
        assert violations == [], f"donor(s) split across train/test in persisted output: {violations}"

    finally:
        features_csv.unlink(missing_ok=True)
        fold_assignments_csv.unlink(missing_ok=True)
