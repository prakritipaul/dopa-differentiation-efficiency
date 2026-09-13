"""Blind gate-2/3 audit tests, derived only from the September 2026 contract.

This file was written before reading harness.py, run_experiment.py, or
run_variant.py.  Keep implementation-informed follow-ups in a separate file.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from sklearn.base import BaseEstimator

from modeling import harness, run_experiment
from modeling.models import ModelSpec


RESULTS = Path(__file__).parents[1] / "results"


def _nested_fixture():
    lines = pd.DataFrame({
        "cell_line": [f"line{i}" for i in range(12)],
        "diff_efficiency": np.linspace(0.05, 0.95, 12),
        "success": [0, 1] * 6,
    })
    rows = []
    for fold in range(3):
        for i, line in lines.iterrows():
            rows.append({
                "scheme": "plain", "repeat": 0, "fold": fold,
                "cell_line": line.cell_line, "pool": f"pool{i % 2}",
                "split": "test" if i % 3 == fold else "train",
                "phat_FPP": 0.1 + i / 100,
                "phat_NB": 0.2 + i / 200,
                "phat_P_FPP": 0.7 - 3 * i / 200,
                **{f"PC{j}": i + j / 100 for j in range(1, 11)},
            })
    return lines, pd.DataFrame(rows)


def test_outer_test_outcome_cannot_change_nested_selection(monkeypatch):
    """Oracle: a line absent from every inner split cannot influence selection."""
    class MeanRegressor(BaseEstimator):
        def __init__(self, token=1):
            self.token = token

        def fit(self, X, y):
            self.mean_ = float(np.mean(y))
            return self

        def predict(self, X):
            return np.repeat(self.mean_, len(X))

    spec = ModelSpec("mean", "regression", MeanRegressor, {"token": [1, 2]})
    monkeypatch.setattr(harness, "models_for_task", lambda _task: [spec])
    monkeypatch.setattr(harness, "PC_COUNT_GRID", [0, 1])
    lines, features = _nested_fixture()

    baseline = harness.run_nested_cv(
        features, lines, "plain", "regression", False, inner_splits=2
    )
    held_out = set(features.loc[(features.fold == 0) & (features.split == "test"), "cell_line"])
    changed = lines.copy()
    changed.loc[changed.cell_line.isin(held_out), "diff_efficiency"] += 10_000
    perturbed = harness.run_nested_cv(
        features, changed, "plain", "regression", False, inner_splits=2
    )

    key = ["model", "fold", "selected_k", "selected_param"]
    left = baseline.loc[baseline.fold == 0, key].drop_duplicates().reset_index(drop=True)
    right = perturbed.loc[perturbed.fold == 0, key].drop_duplicates().reset_index(drop=True)
    pd.testing.assert_frame_equal(left, right)


def test_scaling_is_invariant_to_test_features_and_outcomes():
    fitted = []

    class Recorder(BaseEstimator):
        def fit(self, X, y):
            fitted.append(X.copy())
            return self

        def predict(self, X):
            return np.zeros(len(X))

    spec = ModelSpec("recorder", "regression", Recorder, {})
    train_x = np.array([[0.0, 10.0], [2.0, 14.0], [4.0, 18.0]])
    train_y = np.array([0.1, 0.2, 0.3])
    harness.fit_predict(spec, {}, train_x, train_y, np.array([[6.0, 22.0]]))
    harness.fit_predict(spec, {}, train_x, train_y, np.array([[60_000.0, -90_000.0]]))
    np.testing.assert_array_equal(fitted[0], fitted[1])
    np.testing.assert_allclose(fitted[0].mean(axis=0), 0.0, atol=1e-12)
    np.testing.assert_allclose(fitted[0].std(axis=0), 1.0, atol=1e-12)


def test_repeat_metrics_are_computed_then_averaged_not_pooled():
    # Repeat A has perfect ranking; repeat B has reversed ranking. Their mean
    # ROC-AUC is .5, while pooling creates invalid cross-repeat comparisons.
    predictions = pd.DataFrame({
        "model": "m",
        "repeat": [0, 0, 1, 1, 1, 1],
        "y_true": [0, 1, 0, 0, 1, 1],
        "y_prob": [0.8, 0.9, 0.4, 0.3, 0.2, 0.1],
    })
    got = harness.summarize_across_repeats(predictions, "classification", ["model"])
    assert got.roc_auc_mean.item() == pytest.approx(0.5)
    assert got.roc_auc_std.item() == pytest.approx(np.std([1.0, 0.0], ddof=1))
    assert got.roc_auc_mean.item() != pytest.approx(
        harness.compute_classification_metrics(predictions.y_true, predictions.y_prob)["roc_auc"]
    )


@pytest.mark.parametrize("scheme", ["loco", "lodo"])
def test_leave_one_out_predictions_are_pooled_before_classification_metrics(scheme):
    predictions = pd.DataFrame({
        "scheme": scheme, "model": "m", "repeat": 0,
        "fold": [0, 1, 2, 3], "y_true": [0, 0, 1, 1],
        "y_prob": [0.1, 0.4, 0.6, 0.9],
    })
    got = harness.summarize_across_repeats(predictions, "classification", ["model"])
    assert got.roc_auc_mean.item() == 1.0
    assert got.pr_auc_mean.item() == 1.0
    assert np.isnan(got.roc_auc_std.item())


def test_nested_summary_groups_by_model_not_selected_configuration():
    predictions = pd.DataFrame({
        "model": "ridge", "repeat": [0] * 4,
        "fold": [0, 1, 2, 3], "cell_line": list("abcd"),
        "y_true": [0.0, 0.2, 0.8, 1.0], "y_pred": [0.1, 0.3, 0.7, 0.9],
        "selected_k": [0, 1, 0, 1],
        "selected_param": ["alpha=1", "alpha=10", "alpha=1", "alpha=10"],
    })
    got = harness.summarize_across_repeats(predictions, "regression", ["model"])
    assert len(got) == 1
    assert got.mae_mean.item() == pytest.approx(0.1)


@pytest.mark.parametrize("task,expected", [("regression", "ridge"), ("classification", "logistic_l2")])
def test_headline_is_exactly_the_preregistered_model(task, expected):
    frame = pd.DataFrame({
        "scheme": ["donor_grouped"] * 2,
        "tuning": ["nested"] * 2,
        "pool_correction": [False] * 2,
        "model": [expected, {"ridge": "lasso", "logistic_l2": "logistic_l1"}[expected]],
    })
    got = run_experiment.select_headline(frame, task)
    assert got.model.tolist() == [expected]


def _patched_run(monkeypatch, tmp_path, label_variant="da_untreated"):
    monkeypatch.setattr(run_experiment, "OUT_DIR", tmp_path)
    monkeypatch.setattr(run_experiment, "SCHEMES", ["plain"])
    monkeypatch.setattr(run_experiment, "CORRECTIONS", [False])
    monkeypatch.setattr(run_experiment, "load_lines_with_label", lambda _variant: pd.DataFrame())
    csv = tmp_path / "fold_features.csv"
    pd.DataFrame({"scheme": ["plain"], "label_variant": [label_variant]}).to_csv(csv, index=False)
    pred = pd.DataFrame({
        "scheme": ["plain", "plain"], "repeat": [0, 0], "fold": [0, 1],
        "model": ["ridge", "ridge"], "cell_line": ["a", "b"],
        "y_true": [0.0, 1.0], "y_pred": [0.0, 1.0], "y_prob": [None, None],
        "k": [0, 0], "param_str": ["alpha=1", "alpha=1"],
    })
    monkeypatch.setattr(run_experiment, "run_flat_cv", lambda *args: pred)
    monkeypatch.setattr(run_experiment, "run_nested_cv", lambda *args: pred.rename(
        columns={"k": "selected_k", "param_str": "selected_param"}
    ))
    return csv


def test_feature_table_label_variant_must_match_claim(monkeypatch, tmp_path):
    csv = _patched_run(monkeypatch, tmp_path, "published")
    with pytest.raises(ValueError, match="label_variant|variant|provenance"):
        run_experiment.run_all("regression", csv, label_variant="da_untreated")


@pytest.mark.parametrize("recorded", [None, "published"])
def test_resume_rejects_missing_or_different_label_variant(monkeypatch, tmp_path, recorded):
    csv = _patched_run(monkeypatch, tmp_path)
    prior = pd.DataFrame({
        "scheme": ["plain"], "tuning": ["flat"], "model": ["ridge"],
        "fold_features": [csv.name],
    })
    if recorded is not None:
        prior["label_variant"] = recorded
    results_dir = tmp_path / "results"
    results_dir.mkdir()
    prior.to_csv(results_dir / "results_regression.csv", index=False)
    with pytest.raises(ValueError, match="label_variant|variant|provenance|UNRECORDED"):
        run_experiment.run_all("regression", csv, label_variant="da_untreated", resume=True)


@pytest.mark.parametrize("task", ["regression", "classification"])
@pytest.mark.parametrize("timepoint", ["D11", "D30"])
def test_committed_results_table_contract(task, timepoint):
    path = RESULTS / f"results_{task}_{timepoint}_qualonly_da_untreated.csv"
    table = pd.read_csv(path)
    assert len(table) == 448
    assert set(table.scheme) == {"plain", "donor_grouped", "loco", "lodo"}
    assert set(table.tuning) == {"flat", "nested"}
    expected_models = {"ridge", "lasso"} if task == "regression" else {"logistic_l2", "logistic_l1"}
    assert set(table.model) == expected_models
    assert set(table.label_variant) == {"da_untreated"}
    assert set(table.fold_features) == {f"fold_features_{timepoint}_qualonly_da_untreated.csv"}

    counts = table.groupby(["scheme", "tuning", "model"]).size()
    assert (counts[counts.index.get_level_values("tuning") == "flat"] == 55).all()
    assert (counts[counts.index.get_level_values("tuning") == "nested"] == 1).all()

    means_and_sds = [c for c in table if c.endswith("_mean") or c.endswith("_std")]
    if task == "classification":
        bounded = [c for c in means_and_sds if c.startswith(("roc_auc", "pr_auc", "brier"))]
        assert table[bounded].apply(lambda s: s.dropna().between(0, 1).all()).all()
    else:
        assert table.r2_mean.le(1).all()

    headline_model = "ridge" if task == "regression" else "logistic_l2"
    headline = table[
        (table.scheme == "donor_grouped")
        & (table.tuning == "nested")
        & (table.model == headline_model)
    ]
    assert len(headline) == 1


SINGLE_REPEAT_SCHEMES = {"loco", "lodo"}


@pytest.mark.parametrize("task", ["regression", "classification"])
@pytest.mark.parametrize("timepoint", ["D11", "D30"])
def test_committed_results_metric_columns_are_finite(task, timepoint):
    """Every metric MEAN is finite everywhere. Standard deviations are finite
    for the repeated schemes and NaN for the single-repeat ones.

    The original contract said simply "metric columns finite", which this
    audit showed to be wrong rather than merely unmet: loco and lodo are
    exhaustive enumerations with exactly one repeat, and the sample SD of one
    observation is undefined. NaN is the correct value, and a finite sentinel
    such as 0.0 would be worse than useless -- it asserts zero across-repeat
    variability, which is a precision claim nobody measured and which could be
    quoted as one.

    So the spec is corrected, not relaxed, and the assertion is now
    bidirectional: NaN is REQUIRED where there is one repeat, not merely
    tolerated. That is a stronger test than the one it replaces, which would
    have passed had a sentinel been written into those cells."""
    table = pd.read_csv(RESULTS / f"results_{task}_{timepoint}_qualonly_da_untreated.csv")
    means = [c for c in table if c.endswith("_mean")]
    sds = [c for c in table if c.endswith("_std")]

    assert np.isfinite(table[means].to_numpy()).all(), "every metric mean must be finite"

    repeated = table[~table.scheme.isin(SINGLE_REPEAT_SCHEMES)]
    single = table[table.scheme.isin(SINGLE_REPEAT_SCHEMES)]
    assert len(repeated) and len(single), "fixture must cover both cases"

    assert np.isfinite(repeated[sds].to_numpy()).all(), \
        "repeated schemes average over 10 repeats, so their SDs are defined"
    assert table[table.scheme.isin(SINGLE_REPEAT_SCHEMES)][sds].isna().all().all(), \
        "single-repeat schemes must report NaN SD, not a fabricated finite value"
