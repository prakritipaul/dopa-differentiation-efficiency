"""Implementation-informed integration tests, separate from the blind suite."""

import numpy as np
import pandas as pd
import pytest
from sklearn.decomposition import PCA


def _install_synthetic_expression(monkeypatch, features, expression):
    """Replace h5 streaming only; preserve features.py's masking contract."""
    state = {"expression": np.asarray(expression, dtype=float)}
    seen_masks = []

    def compute_gene_stats(_path, fit_mask):
        seen_masks.append(fit_mask.copy())
        fit = state["expression"][fit_mask]
        return fit.mean(axis=0), fit.var(axis=0), len(fit)

    monkeypatch.setattr(features.m005, "compute_gene_stats", compute_gene_stats)
    monkeypatch.setattr(
        features.m005, "select_hvgs",
        lambda mean, var, n_top: np.arange(len(mean)),
    )
    monkeypatch.setattr(
        features.m005, "extract_hvg_matrix",
        lambda _path, hvg_idx: state["expression"][:, hvg_idx].copy(),
    )

    def run_pca(X, mean, var, hvg_idx, fit_mask, n_components):
        sd = np.sqrt(var[hvg_idx])
        sd[sd == 0] = 1
        scaled = (X[:, hvg_idx] - mean[hvg_idx]) / sd
        model = PCA(n_components=n_components, svd_solver="full").fit(scaled[fit_mask])
        return model.transform(scaled), model.explained_variance_ratio_

    monkeypatch.setattr(features.m005, "run_pca", run_pca)
    return state, seen_masks


def test_pca_held_out_values_cannot_change_training_coordinates_and_all_are_projected(monkeypatch):
    from modeling import features

    meta = pd.DataFrame(
        {
            "cell_line": ["train_a", "train_a", "train_b", "train_b", "held", "held"],
            "pool": ["pool1", "pool2", "pool1", "pool2", "pool1", "pool2"],
        }
    )
    expression = np.array(
        [[0, 1, 2], [1, 2, 0], [3, 1, 1], [4, 0, 3], [5, 5, 5], [6, 6, 6]],
        dtype=float,
    )
    state, masks = _install_synthetic_expression(monkeypatch, features, expression)
    first = features.compute_pca_features_for_fold(
        "D30", {"held"}, meta=meta, n_pcs=2, exclude_pools=frozenset()
    )
    state["expression"][-2:] = [[-10_000, 8_000, 4_000], [9_000, -7_000, 3_000]]
    second = features.compute_pca_features_for_fold(
        "D30", {"held"}, meta=meta, n_pcs=2, exclude_pools=frozenset()
    )

    assert len(first) == len(second) == len(meta)
    assert all(mask.tolist() == [True, True, True, True, False, False] for mask in masks)
    train = meta.cell_line != "held"
    np.testing.assert_allclose(first.loc[train, ["PC1", "PC2"]], second.loc[train, ["PC1", "PC2"]])


def test_default_d30_pool_exclusion_affects_fit_mask_but_not_projection(monkeypatch):
    from modeling import features

    meta = pd.DataFrame(
        {
            "cell_line": ["a", "b", "c", "d", "held"],
            "pool": ["pool5", "pool11", "pool1", "pool2", "pool3"],
        }
    )
    expression = np.arange(15, dtype=float).reshape(5, 3)
    _, masks = _install_synthetic_expression(monkeypatch, features, expression)
    result = features.compute_pca_features_for_fold("D30", {"held"}, meta=meta, n_pcs=1)
    assert masks[0].tolist() == [False, True, True, True, False]
    assert len(result) == 5
    assert set(result.cell_line) == set(meta.cell_line)


@pytest.mark.parametrize("recorded_population", ["all_cells", None])
def test_resume_rejects_mixed_or_unrecorded_fit_population_before_any_pca_work(
    monkeypatch, tmp_path, recorded_population
):
    from modeling import run_feature_extraction as extraction
    from modeling.folds import Fold

    monkeypatch.setattr(extraction, "OUT_DIR", tmp_path)
    lines = pd.DataFrame(
        {"cell_line": ["a", "b"], "donor": ["x", "y"],
         "diff_efficiency": [0.1, 0.3], "success": [False, True]}
    )
    fold = Fold("plain", 0, 0, ("a",), ("b",))
    monkeypatch.setattr(extraction, "load_lines_with_label", lambda: lines)
    monkeypatch.setattr(extraction, "plain_repeated_kfold", lambda *a, **k: [fold])
    monkeypatch.setattr(extraction, "donor_grouped_repeated_kfold", lambda *a, **k: [])
    monkeypatch.setattr(
        extraction, "load_cell_metadata",
        lambda timepoint: pd.DataFrame(
            {"cell_line": ["a", "b"], "pool": ["p", "p"], "celltype": ["FPP", "FPP"]}
        ),
    )
    monkeypatch.setattr(
        extraction, "compute_proportion_features",
        lambda meta: pd.DataFrame(
            {"cell_line": ["a", "b"], "pool": ["p", "p"], "phat_FPP": [1.0, 1.0]}
        ),
    )
    qualifying = tmp_path / "qualifying.csv"
    pd.DataFrame({"cell_line": ["a", "b"], "pool": ["p", "p"]}).to_csv(qualifying, index=False)
    monkeypatch.setattr(extraction, "QUALIFYING_COMBOS_CSV", qualifying)
    monkeypatch.setattr(
        extraction, "compute_pca_features_for_fold",
        lambda *a, **k: pytest.fail("PCA must not run before provenance validation"),
    )

    out = tmp_path / "fold_data" / "fold_features_D30_audit.csv"
    out.parent.mkdir()
    previous = {"scheme": ["plain"], "repeat": [0], "fold": [0]}
    if recorded_population is not None:
        previous["fit_population"] = [recorded_population]
    pd.DataFrame(previous).to_csv(out, index=False)
    with pytest.raises(ValueError, match="mix two feature bases|fit_population|provenance"):
        extraction.main(
            timepoint="D30", n_splits=2, n_repeats=1, out_suffix="_audit",
            include_loco_lodo=False, restrict_fit_to_qualifying=True, resume=True,
        )


def test_declared_direction_controls_auc_and_is_not_inferred_from_data():
    from modeling.d30_single_feature_benchmarks import _metrics

    y = np.array([0, 0, 1, 1])
    score = np.array([4.0, 3.0, 2.0, 1.0])
    efficiency = y.astype(float)
    assert _metrics(score, efficiency, y, -1)["roc_auc"] == pytest.approx(1.0)
    assert _metrics(score, efficiency, y, +1)["roc_auc"] == pytest.approx(0.0)


def test_donor_bootstrap_reproducible_and_drops_single_class_resamples():
    from modeling.d30_single_feature_benchmarks import bootstrap_by_donor

    lines = pd.DataFrame(
        {"donor": ["fail", "fail", "success", "success"],
         "success": [0, 0, 1, 1], "diff_efficiency": [0.0, 0.1, 0.9, 1.0],
         "score": [4.0, 3.0, 2.0, 1.0]}
    )
    a = bootstrap_by_donor(lines, "score", -1, n_boot=200, seed=77)
    b = bootstrap_by_donor(lines, "score", -1, n_boot=200, seed=77)
    assert a == b
    assert 0 < a["roc_auc_n_boot_used"] < 200
    assert a["roc_auc_lo"] == pytest.approx(1.0)
    assert a["roc_auc_hi"] == pytest.approx(1.0)
