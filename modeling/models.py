"""
Model registry.

Adding a new model means adding one ModelSpec entry here -- the CV/tuning
harness (harness.py) never special-cases a model class; it only calls the
sklearn estimator interface (.fit/.predict/.predict_proba for classifiers)
generically. Open for extension (new ModelSpec entries), closed for
modification (harness.py never changes to add a model).
"""

from dataclasses import dataclass
from typing import Callable, Sequence

from sklearn.base import BaseEstimator
from sklearn.linear_model import Lasso, LogisticRegression, Ridge

# Number of D11 PCs to include (k=0 is the proportions-only baseline; PCA
# components are hierarchical, so this truncates columns from one PCA fit
# per fold rather than needing a distinct PCA refit per k). Not part of
# any single ModelSpec's param_grid since it controls feature-matrix
# construction (see features.py), not an estimator hyperparameter.
PC_COUNT_GRID: list[int] = list(range(0, 11))


@dataclass(frozen=True)
class ModelSpec:
    """A pluggable model: a zero-arg factory (so every fold/candidate gets
    a fresh, unfitted estimator instance -- never a shared, mutated one)
    plus its hyperparameter grid, keyed by the constructor's own parameter
    names so it can be applied via `estimator.set_params(**candidate)`."""

    name: str
    task: str  # "regression" or "classification"
    estimator_factory: Callable[[], BaseEstimator]
    param_grid: dict[str, Sequence]


REGRESSION_MODELS: list[ModelSpec] = [
    ModelSpec(
        name="ridge",
        task="regression",
        estimator_factory=lambda: Ridge(),
        param_grid={"alpha": [0.01, 0.1, 1.0, 10.0, 100.0]},
    ),
    ModelSpec(
        name="lasso",
        task="regression",
        estimator_factory=lambda: Lasso(max_iter=10_000),
        param_grid={"alpha": [0.001, 0.01, 0.1, 1.0, 10.0]},
    ),
]

CLASSIFICATION_MODELS: list[ModelSpec] = [
    ModelSpec(
        # sklearn >=1.8 deprecated `penalty=` in favor of `l1_ratio`
        # (l1_ratio=0 == pure L2, l1_ratio=1 == pure L1). L2 works with
        # the default 'lbfgs' solver; L1 requires 'saga' (lbfgs only
        # supports L2/None).
        name="logistic_l2",
        task="classification",
        estimator_factory=lambda: LogisticRegression(
            l1_ratio=0.0, solver="lbfgs", max_iter=5_000, class_weight="balanced"
        ),
        param_grid={"C": [0.01, 0.1, 1.0, 10.0, 100.0]},
    ),
    ModelSpec(
        name="logistic_l1",
        task="classification",
        estimator_factory=lambda: LogisticRegression(
            l1_ratio=1.0, solver="saga", max_iter=5_000, class_weight="balanced"
        ),
        param_grid={"C": [0.01, 0.1, 1.0, 10.0, 100.0]},
    ),
]

ALL_MODELS: list[ModelSpec] = REGRESSION_MODELS + CLASSIFICATION_MODELS


def models_for_task(task: str) -> list[ModelSpec]:
    """task: 'regression' or 'classification'."""
    matches = [m for m in ALL_MODELS if m.task == task]
    if not matches:
        raise ValueError(f"No models registered for task={task!r}")
    return matches
