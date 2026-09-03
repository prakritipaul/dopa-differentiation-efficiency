"""Tests for modeling/models.py -- the ModelSpec registry."""

import pytest
from sklearn.base import BaseEstimator
from sklearn.linear_model import Lasso, LogisticRegression, Ridge

from modeling.models import PC_COUNT_GRID, models_for_task


def test_pc_count_grid():
    assert PC_COUNT_GRID == list(range(0, 11))


def test_models_for_task_regression():
    specs = models_for_task("regression")
    names = {s.name for s in specs}
    assert names == {"ridge", "lasso"}
    for s in specs:
        assert s.task == "regression"
        assert "alpha" in s.param_grid
        assert len(s.param_grid["alpha"]) > 0


def test_models_for_task_classification():
    specs = models_for_task("classification")
    names = {s.name for s in specs}
    assert names == {"logistic_l2", "logistic_l1"}
    for s in specs:
        assert s.task == "classification"
        assert "C" in s.param_grid
        assert len(s.param_grid["C"]) > 0


def test_models_for_task_unknown_raises():
    with pytest.raises(ValueError):
        models_for_task("not_a_real_task")


def test_estimator_factories_produce_correct_types():
    expected = {"ridge": Ridge, "lasso": Lasso, "logistic_l2": LogisticRegression, "logistic_l1": LogisticRegression}
    for spec in models_for_task("regression") + models_for_task("classification"):
        est = spec.estimator_factory()
        assert isinstance(est, BaseEstimator)
        assert isinstance(est, expected[spec.name])


def test_estimator_factory_produces_fresh_unfitted_instances():
    """Fitting one instance from a factory must not affect another
    instance produced by the same factory (each call is a fresh object,
    not a shared mutated one)."""
    import numpy as np

    spec = models_for_task("regression")[0]  # ridge
    est_a = spec.estimator_factory()
    est_b = spec.estimator_factory()
    assert est_a is not est_b

    X = np.array([[1.0], [2.0], [3.0]])
    y = np.array([1.0, 2.0, 3.0])
    est_a.fit(X, y)
    assert hasattr(est_a, "coef_")
    assert not hasattr(est_b, "coef_"), "fitting est_a leaked fitted state into est_b"
