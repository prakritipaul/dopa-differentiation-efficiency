"""Blind, contract-derived tests for the D30 extension.

Written before inspecting features.py, harness.py, run_feature_extraction.py,
make_d30_variant_tables.py, or d30_single_feature_benchmarks.py.  Keep this
file independent of implementation-informed tests.
"""

from dataclasses import asdict

import numpy as np
import pandas as pd
import pytest

from modeling.folds import all_folds, load_lines_with_label


D11_PHATS = ["phat_FPP", "phat_NB", "phat_P_FPP"]
D30_PHATS = [
    "phat_DA",
    "phat_Epen1",
    "phat_FPP",
    "phat_P_FPP",
    "phat_Sert",
    "phat_U_Neur1",
    "phat_U_Neur2",
]
PROGENITOR_PHATS = [
    "phat_Epen1",
    "phat_FPP",
    "phat_P_FPP",
    "phat_U_Neur1",
    "phat_U_Neur2",
]


def _fold_bytes(folds_by_scheme):
    rows = []
    for scheme in sorted(folds_by_scheme):
        for fold in folds_by_scheme[scheme]:
            rows.append(asdict(fold))
    return repr(rows).encode()


def test_depth_outlier_pools_are_timepoint_specific_contract_values():
    from modeling.features import DEPTH_OUTLIER_POOLS

    assert DEPTH_OUTLIER_POOLS["D11"] == frozenset({"pool11"})
    assert DEPTH_OUTLIER_POOLS["D30"] == frozenset({"pool5"})
    assert "pool11" not in DEPTH_OUTLIER_POOLS["D30"]


def test_proportion_columns_and_feature_matrix_are_exact_and_ordered():
    from modeling.harness import build_feature_matrix, proportion_cols

    d11 = pd.DataFrame({name: [0.2, 0.3] for name in reversed(D11_PHATS)})
    assert proportion_cols(d11) == (D11_PHATS, D11_PHATS[:-1])

    rng = np.random.default_rng(9182)
    raw = rng.gamma(2, 1, size=(30, 7))
    comp = raw / raw.sum(axis=1, keepdims=True)
    d30 = pd.DataFrame(comp, columns=list(reversed(D30_PHATS)))
    all_cols, model_cols = proportion_cols(d30)
    assert all_cols == D30_PHATS
    assert model_cols == [c for c in D30_PHATS if c != "phat_P_FPP"]
    assert np.linalg.matrix_rank(
        np.column_stack((np.ones(len(d30)), d30[model_cols].to_numpy()))
    ) == 7

    for i in range(1, 11):
        d30[f"PC{i}"] = np.arange(len(d30)) + 100 * i
    for k in range(11):
        expected = d30[model_cols + [f"PC{i}" for i in range(1, k + 1)]].to_numpy()
        np.testing.assert_array_equal(build_feature_matrix(d30, k), expected)


@pytest.mark.parametrize(
    "frame, message",
    [
        (pd.DataFrame({"PC1": [1.0]}), "phat"),
        (pd.DataFrame({"phat_DA": [0.5], "phat_Sert": [0.5]}), "reference"),
    ],
)
def test_proportion_columns_reject_invalid_compositions(frame, message):
    from modeling.harness import proportion_cols

    with pytest.raises(ValueError, match=message):
        proportion_cols(frame)


def test_folds_are_deterministic_timepoint_free_and_obey_grouping_contract():
    lines = load_lines_with_label("published")
    first = all_folds(lines, n_splits=5, n_repeats=2, seed=1729)
    second = all_folds(lines.copy(), n_splits=5, n_repeats=2, seed=1729)
    # There is deliberately no timepoint argument: either early-timepoint data
    # set must induce these same bytes from the shared line metadata.
    assert _fold_bytes(first) == _fold_bytes(second)

    universe = set(lines.cell_line)
    donor_of = lines.set_index("cell_line").donor.to_dict()
    assert len(first["loco"]) == 138
    assert all(len(f.test_lines) == 1 for f in first["loco"])
    assert {f.test_lines[0] for f in first["loco"]} == universe

    assert len(first["lodo"]) == lines.donor.nunique() == 20
    for f in first["lodo"]:
        held_donors = {donor_of[x] for x in f.test_lines}
        assert len(held_donors) == 1
        held = next(iter(held_donors))
        assert set(f.test_lines) == set(lines.loc[lines.donor == held, "cell_line"])

    for f in first["donor_grouped"]:
        assert {donor_of[x] for x in f.train_lines}.isdisjoint(
            {donor_of[x] for x in f.test_lines}
        )
    for scheme in ("plain", "donor_grouped"):
        for repeat in range(2):
            test_lines = [
                x for f in first[scheme] if f.repeat == repeat for x in f.test_lines
            ]
            assert len(test_lines) == 138
            assert set(test_lines) == universe


def _variant_input():
    # Same conditional balance among the five progenitor types (1:2:3:4:5),
    # but non-target masses 0.9 and 0.2 respectively.
    relative = np.arange(1, 6, dtype=float)
    relative /= relative.sum()
    rows = []
    for line, kept_mass, da, sert, pc1 in (
        ("a", 0.9, 0.05, 0.05, 17.0),
        ("b", 0.2, 0.50, 0.30, -4.0),
    ):
        row = {"cell_line": line, "pool": "p", "PC1": pc1}
        row.update(dict(zip(PROGENITOR_PHATS, relative * kept_mass)))
        row.update(phat_DA=da, phat_Sert=sert)
        rows.append(row)
    return pd.DataFrame(rows)


def test_progenitor_balance_removes_target_magnitude_oracle():
    from modeling.make_d30_variant_tables import build_progenitor_balance_table

    source = _variant_input()
    result = build_progenitor_balance_table(source)
    assert len(result) == len(source)
    assert "phat_DA" not in result and "phat_Sert" not in result
    assert set(c for c in result if c.startswith("phat_")) == set(PROGENITOR_PHATS)
    np.testing.assert_allclose(result[PROGENITOR_PHATS].sum(axis=1), 1.0, atol=1e-9)
    np.testing.assert_allclose(
        result.loc[0, PROGENITOR_PHATS].to_numpy(float),
        result.loc[1, PROGENITOR_PHATS].to_numpy(float),
        atol=1e-12,
    )
    np.testing.assert_array_equal(result["PC1"], source["PC1"])


def test_progenitor_balance_rejects_zero_mass_and_unknown_cell_type():
    from modeling.make_d30_variant_tables import build_progenitor_balance_table

    zero = _variant_input().iloc[[0]].copy()
    zero[PROGENITOR_PHATS] = 0.0
    zero["phat_DA"], zero["phat_Sert"] = 0.7, 0.3
    with pytest.raises(ValueError):
        build_progenitor_balance_table(zero)

    unknown = _variant_input()
    unknown["phat_Mystery"] = 0.0
    with pytest.raises(ValueError, match="Mystery|unclassified|unknown"):
        build_progenitor_balance_table(unknown)


