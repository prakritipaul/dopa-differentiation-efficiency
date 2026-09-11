"""Fast, synthetic tests for the modeling CV harness."""

import numpy as np
import pandas as pd
import pytest
from sklearn.base import BaseEstimator

from modeling import harness
from modeling.models import ModelSpec


@pytest.fixture
def lines_and_fold_features():
    lines = pd.DataFrame({"cell_line": [f"line{i:02d}" for i in range(30)]})
    lines["diff_efficiency"] = np.linspace(0.02, 0.60, len(lines))
    lines["success"] = (np.arange(len(lines)) % 2).astype(int)

    rows = []
    for scheme in ("plain", "donor_grouped"):
        for repeat in range(2):
            for fold in range(2):
                for i, line in enumerate(lines.itertuples(index=False)):
                    p = line.diff_efficiency
                    rows.append(
                        {
                            "scheme": scheme,
                            "repeat": repeat,
                            "fold": fold,
                            "cell_line": line.cell_line,
                            "pool": f"pool{i % 4}",
                            "split": "test" if (i + repeat) % 2 == fold else "train",
                            "phat_FPP": p,
                            "phat_NB": (1 - p) * 0.4,
                            "phat_P_FPP": (1 - p) * 0.6,
                            **{f"PC{j}": i * 100 + j for j in range(1, 11)},
                        }
                    )
    return lines, pd.DataFrame(rows)


def _identity_regression_spec():
    class IdentityRegressor(BaseEstimator):
        def fit(self, X, y):
            return self

        def predict(self, X):
            # Scaling is deliberately bypassed in alignment tests by replacing
            # fit_predict; this estimator merely supplies a valid ModelSpec.
            return X[:, 0]

    return ModelSpec("identity", "regression", IdentityRegressor, {"token": [1]})


def test_build_feature_matrix_identity_and_order():
    frame = pd.DataFrame(
        {c: [i, i + 100] for i, c in enumerate(harness.ALL_PROPORTION_COLS + harness.ALL_PC_COLS)}
    )
    for k in (0, 3, 10):
        expected_cols = harness.MODEL_PROPORTION_COLS + harness.ALL_PC_COLS[:k]
        np.testing.assert_array_equal(harness.build_feature_matrix(frame, k), frame[expected_cols].to_numpy())


def test_proportion_cols_reproduces_the_pre_refactor_d11_constants():
    # Oracle: the literal lists that were hardcoded in harness.py before the
    # columns became data-derived for D30. If the derivation ever stops
    # agreeing with them, every published D11 number silently changes basis.
    old_all = ["phat_FPP", "phat_NB", "phat_P_FPP"]
    old_model = ["phat_FPP", "phat_NB"]
    d11 = pd.DataFrame({c: [0.3, 0.4] for c in old_all})
    assert harness.proportion_cols(d11) == (old_all, old_model)


def test_proportion_cols_handles_the_seven_type_d30_composition():
    d30_types = ["DA", "Epen1", "FPP", "P_FPP", "Sert", "U_Neur1", "U_Neur2"]
    frame = pd.DataFrame({f"phat_{t}": [1 / 7, 1 / 7] for t in d30_types})
    all_cols, model_cols = harness.proportion_cols(frame)
    assert all_cols == sorted(all_cols), "column order must be deterministic across runs"
    assert len(all_cols) == 7 and len(model_cols) == 6
    assert harness.REFERENCE_PROPORTION not in model_cols
    # Six free coordinates + intercept must be full rank on a real composition.
    rng = np.random.default_rng(0)
    raw = rng.random((20, 7))
    comp = pd.DataFrame(raw / raw.sum(axis=1, keepdims=True), columns=[f"phat_{t}" for t in d30_types])
    assert comp.sum(axis=1).round(12).eq(1.0).all(), "fixture must be compositional"
    X = harness.build_feature_matrix(comp, 0)
    assert X.shape[1] == 6
    assert np.linalg.matrix_rank(np.column_stack([np.ones(len(X)), X])) == 7


def test_proportion_cols_rejects_a_table_missing_the_reference_category():
    # Dropping the reference leaves a set that no longer sums to a known
    # constant; fitting it silently would be rank-deficient-adjacent and wrong.
    frame = pd.DataFrame({"phat_DA": [0.5], "phat_Sert": [0.5]})
    with pytest.raises(ValueError, match="reference category"):
        harness.proportion_cols(frame)
    with pytest.raises(ValueError, match="no phat_"):
        harness.proportion_cols(pd.DataFrame({"PC1": [1.0]}))


def test_build_feature_matrix_omits_one_implicit_compositional_coordinate():
    # The retained coordinates must vary INDEPENDENTLY, not just sum to 1. If
    # phat_P_FPP were held constant, phat_FPP + phat_NB would also be constant
    # and the affine dependency would survive dropping a column -- making this
    # test unable to flip to a pass once the defect is fixed.
    frame = pd.DataFrame(
        {
            "phat_FPP": [0.2, 0.3, 0.4, 0.25],
            "phat_NB": [0.3, 0.2, 0.1, 0.35],
            "phat_P_FPP": [0.5, 0.5, 0.5, 0.40],
            **{f"PC{i}": [i, i + 1, i + 2, i + 3] for i in range(1, 11)},
        }
    )
    assert frame[["phat_FPP", "phat_NB", "phat_P_FPP"]].sum(axis=1).eq(1.0).all(), "fixture must be compositional"
    X = harness.build_feature_matrix(frame, 0)
    assert X.shape[1] == 2
    assert np.linalg.matrix_rank(np.column_stack([np.ones(len(X)), X])) == 3


def test_fit_predict_scales_from_train_only_and_uses_fresh_estimator():
    fitted = []

    class RecordingEstimator(BaseEstimator):
        def set_params(self, **params):
            return self

        def fit(self, X, y):
            self.fit_X = X.copy()
            fitted.append(self)
            return self

        def predict(self, X):
            self.test_X = X.copy()
            return X[:, 0]

    spec = ModelSpec("recording", "regression", RecordingEstimator, {})
    X_train = np.array([[0.0], [2.0]])
    X_test = np.array([[101.0]])
    harness.fit_predict(spec, {}, X_train, np.array([0.0, 1.0]), X_test)
    harness.fit_predict(spec, {}, X_train, np.array([0.0, 1.0]), X_test)
    assert fitted[0] is not fitted[1]
    np.testing.assert_allclose(fitted[0].fit_X[:, 0], [-1, 1])
    np.testing.assert_allclose(fitted[0].test_X[:, 0], [100])


def test_build_line_level_partition_and_unweighted_pool_mean(lines_and_fold_features):
    lines, ff = lines_and_fold_features
    train, test = harness.build_line_level_for_fold(ff, "plain", 0, 0, False)
    assert set(train.cell_line).isdisjoint(test.cell_line)
    assert set(train.cell_line) | set(test.cell_line) == set(lines.cell_line)

    base = ff[(ff.scheme == "plain") & (ff.repeat == 0) & (ff.fold == 0)].copy()
    target = base.iloc[0].copy()
    target["cell_line"], target["pool"], target["split"] = "multi", "poolA", "train"
    target["PC1"] = 0.0
    additions = [target.copy() for _ in range(10)]
    target["pool"], target["PC1"] = "poolB", 100.0
    additions.append(target)
    augmented = pd.concat([base, pd.DataFrame(additions)], ignore_index=True)
    train, _ = harness.build_line_level_for_fold(augmented, "plain", 0, 0, False)
    assert train.loc[train.cell_line == "multi", "PC1"].item() == 50.0
    assert train.loc[train.cell_line == "multi", "PC1"].item() != pytest.approx(100 / 11)


@pytest.mark.xfail(strict=True, reason="Unseen test pools currently map to NaN during pool correction")
def test_pool_correction_has_defined_values_for_pool_seen_only_in_test():
    cols = harness.ALL_PROPORTION_COLS + harness.ALL_PC_COLS
    rows = []
    for line, pool, split, value in (("train", "poolA", "train", 1.0), ("test", "poolB", "test", 2.0)):
        rows.append({"scheme": "plain", "repeat": 0, "fold": 0, "cell_line": line,
                     "pool": pool, "split": split, **{c: value for c in cols}})
    train, test = harness.build_line_level_for_fold(pd.DataFrame(rows), "plain", 0, 0, True)
    assert not train[cols].isna().any().any()
    assert not test[cols].isna().any().any()


def test_flat_cv_grid_count_and_y_x_alignment(monkeypatch, lines_and_fold_features):
    lines, ff = lines_and_fold_features
    spec = _identity_regression_spec()
    spec = ModelSpec(spec.name, spec.task, spec.estimator_factory, {"token": [1, 2]})
    monkeypatch.setattr(harness, "models_for_task", lambda task: [spec])
    monkeypatch.setattr(harness, "PC_COUNT_GRID", [0, 1])

    def aligned_predict(_spec, _params, X_train, y_train, X_test):
        np.testing.assert_allclose(X_train[:, 0], y_train)
        return X_test[:, 0], None

    monkeypatch.setattr(harness, "fit_predict", aligned_predict)
    out = harness.run_flat_cv(ff, lines, "plain", "regression", False)
    assert len(out) == 2 * 2 * 1 * 2 * 2 * 15
    lookup = lines.set_index("cell_line").diff_efficiency
    np.testing.assert_allclose(out.y_true, out.cell_line.map(lookup))
    np.testing.assert_allclose(out.y_pred, out.y_true)


def test_nested_cv_alignment_and_inner_loop_excludes_outer_test(monkeypatch, lines_and_fold_features):
    lines, ff = lines_and_fold_features
    monkeypatch.setattr(harness, "models_for_task", lambda task: [_identity_regression_spec()])
    monkeypatch.setattr(harness, "PC_COUNT_GRID", [0])

    # phat_FPP uniquely identifies every line, making membership and alignment observable.
    outer_test_values = {}
    for fold in (0, 1):
        _, test = harness.build_line_level_for_fold(ff, "plain", 0, fold, False)
        outer_test_values[fold] = set(test.phat_FPP)
    active_fold = iter([0, 1])
    fold_state = {"fold": 0, "calls": 0}

    original_builder = harness.build_line_level_for_fold
    def tracking_builder(*args, **kwargs):
        result = original_builder(*args, **kwargs)
        fold_state["fold"] = args[3]
        fold_state["calls"] = 0
        return result

    def spy_predict(_spec, _params, X_train, y_train, X_test):
        np.testing.assert_allclose(X_train[:, 0], y_train)
        # First two calls per outer fold are the two inner fits; neither side may
        # contain a feature value belonging to that outer fold's held-out lines.
        if fold_state["calls"] < 2:
            assert set(X_train[:, 0]).isdisjoint(outer_test_values[fold_state["fold"]])
            assert set(X_test[:, 0]).isdisjoint(outer_test_values[fold_state["fold"]])
        fold_state["calls"] += 1
        return X_test[:, 0], None

    monkeypatch.setattr(harness, "build_line_level_for_fold", tracking_builder)
    monkeypatch.setattr(harness, "fit_predict", spy_predict)
    one_repeat = ff[(ff.scheme == "plain") & (ff.repeat == 0)]
    out = harness.run_nested_cv(one_repeat, lines, "plain", "regression", False, inner_splits=2)
    lookup = lines.set_index("cell_line").diff_efficiency
    np.testing.assert_allclose(out.y_true, out.cell_line.map(lookup))
    np.testing.assert_allclose(out.y_pred, out.y_true)


def test_classification_metrics_known_example_threshold_and_single_class():
    result = harness.compute_classification_metrics(np.array([0, 0, 1, 1]), np.array([0.1, 0.5, 0.4, 0.9]))
    assert result == pytest.approx({"roc_auc": 0.75, "pr_auc": 5 / 6, "balanced_accuracy": 0.5,
                                    "sensitivity": 0.5, "specificity": 0.5, "f1": 0.5,
                                    "brier": 0.1575})
    assert all(np.isnan(v) for v in harness.compute_classification_metrics(np.ones(3), np.array([.1, .2, .3])).values())


def test_regression_metrics_known_example():
    result = harness.compute_regression_metrics(np.array([0.0, 0.5, 1.0]), np.array([-0.1, 0.7, 1.2]))
    assert result == pytest.approx({"mae": 1 / 6, "rmse": np.sqrt(0.03), "r2": 0.82,
                                    "frac_out_of_range": 2 / 3})


def test_metrics_are_averaged_per_repeat_and_single_repeat_sd_is_nan():
    # Repeat MAEs are 0 and 2, so their mean is 1. Pooled MAE happens to be
    # 2/3 because repeats contain different numbers of predictions.
    preds = pd.DataFrame({"model": "m", "repeat": [0, 0, 0, 0, 1, 1],
                          "y_true": np.zeros(6), "y_pred": [0, 0, 0, 0, 2, 2]})
    summary = harness.summarize_across_repeats(preds, "regression", ["model"])
    assert summary.mae_mean.item() == 1.0
    assert summary.mae_mean.item() != pytest.approx(np.mean(np.abs(preds.y_pred - preds.y_true)))
    one = harness.summarize_across_repeats(preds[preds.repeat == 0], "regression", ["model"])
    assert np.isnan(one.mae_std.item())
