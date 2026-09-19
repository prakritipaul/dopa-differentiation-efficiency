"""The importance scatter must plot the table, not a transformation of it.

The figure is the first thing a reader of FINDINGS.md sees for each
timepoint, so a silent mismatch between a plotted point and its row would
misstate a result without raising. The only row the plot is allowed to
drop is the reference proportion, which has no SHAP value because it never
enters a fit.

No h5 data needed -- the importance tables are in the repository.
"""

import pandas as pd
import pytest

from modeling.plot_feature_importance import (
    LABEL_VARIANT,
    RESULTS_DIR,
    TIMEPOINTS,
    load_importance,
)


def _table(timepoint: str) -> pd.DataFrame:
    return pd.read_csv(
        RESULTS_DIR / f"feature_importance_table_classification_logistic_l2_{timepoint}_{LABEL_VARIANT}.csv"
    )


@pytest.mark.parametrize("timepoint", TIMEPOINTS)
def test_plotted_points_equal_their_table_row(timepoint):
    raw = _table(timepoint).set_index("feature")
    for _, row in load_importance(timepoint).iterrows():
        source = raw.loc[row["feature"]]
        assert row["shap_mean_abs"] == source["shap_mean_abs"]
        assert row["technical_covariate_eta2"] == source["technical_covariate_eta2"]


@pytest.mark.parametrize("timepoint", TIMEPOINTS)
def test_only_rows_without_shap_are_dropped(timepoint):
    raw = _table(timepoint)
    dropped = set(raw["feature"]) - set(load_importance(timepoint)["feature"])
    assert dropped, "the reference proportion should always be dropped"
    assert raw[raw["feature"].isin(dropped)]["shap_mean_abs"].isna().all()
