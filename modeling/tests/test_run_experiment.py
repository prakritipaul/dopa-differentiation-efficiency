"""Fast orchestration and pre-registration tests for run_experiment."""

import pandas as pd
import pytest

from modeling import run_experiment


def _summary_predictions(scheme, nested=False):
    models = ("ridge", "lasso")
    rows = []
    for model in models:
        for repeat in (0, 1):
            for line_number, value in enumerate((0.2, 0.4)):
                row = {"scheme": scheme, "repeat": repeat, "fold": line_number, "model": model,
                       "cell_line": f"{scheme}-{repeat}-{line_number}", "y_true": value,
                       "y_pred": value, "y_prob": None}
                if nested:
                    row.update(selected_k=1, selected_param="alpha=1.0")
                else:
                    row.update(k=1, param_str="alpha=1.0")
                rows.append(row)
    return pd.DataFrame(rows)


@pytest.fixture
def patched_run_all(monkeypatch, tmp_path):
    csv = tmp_path / "fold_features.csv"
    pd.DataFrame({"scheme": ["placeholder"]}).to_csv(csv, index=False)
    monkeypatch.setattr(run_experiment, "SCHEMES", ["plain", "donor_grouped", "loco", "lodo"])
    monkeypatch.setattr(run_experiment, "CORRECTIONS", [False])
    monkeypatch.setattr(run_experiment, "load_lines_with_label", lambda: pd.DataFrame())
    monkeypatch.setattr(run_experiment, "run_flat_cv",
                        lambda _ff, _lines, scheme, _task, _correction: _summary_predictions(scheme))
    monkeypatch.setattr(run_experiment, "run_nested_cv",
                        lambda _ff, _lines, scheme, _task, _correction: _summary_predictions(scheme, nested=True))
    return csv


def test_run_all_has_flat_and_nested_rows_for_every_scheme(patched_run_all):
    out = run_experiment.run_all("regression", patched_run_all)
    counts = out.groupby(["scheme", "tuning"]).size()
    assert set(counts.index) == {(s, t) for s in run_experiment.SCHEMES for t in ("flat", "nested")}
    assert (counts == 2).all()  # two model families per scheme/tuning


def test_run_all_preserves_nested_selected_configuration(patched_run_all):
    out = run_experiment.run_all("regression", patched_run_all)
    nested = out[out.tuning == "nested"]
    assert nested["k"].notna().all()
    assert nested["param_str"].notna().all()
    assert nested["selection_distribution"].notna().all()
    # Metrics must still be grouped on model alone: one row per model per
    # scheme. Grouping by the selected config instead would split a repeat
    # into config-specific subsets of lines and change the reported metric.
    assert (nested.groupby(["scheme", "model"]).size() == 1).all()


def _headline_frame():
    return pd.DataFrame(
        [
            {"scheme": "donor_grouped", "tuning": "nested", "pool_correction": False, "model": model}
            for model in ("ridge", "lasso", "logistic_l2", "logistic_l1")
        ]
        + [{"scheme": "plain", "tuning": "flat", "pool_correction": False, "model": "ridge"}]
    )


@pytest.mark.parametrize("task,model", [("regression", "ridge"), ("classification", "logistic_l2")])
def test_select_headline_returns_exactly_preregistered_model(task, model):
    selected = run_experiment.select_headline(_headline_frame(), task)
    assert len(selected) == 1
    assert selected.model.item() == model


def test_select_headline_without_task_returns_both_task_model_families():
    selected = run_experiment.select_headline(_headline_frame())
    assert set(selected.model) == {"ridge", "lasso", "logistic_l2", "logistic_l1"}


def test_preregistered_headline_exists_in_run_all_output(patched_run_all):
    out = run_experiment.run_all("regression", patched_run_all)
    selected = run_experiment.select_headline(out, "regression")
    assert len(selected) == 1
    assert selected.model.item() == "ridge"
