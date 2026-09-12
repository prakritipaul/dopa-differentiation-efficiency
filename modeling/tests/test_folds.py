"""Tests for modeling/folds.py -- fold construction and leakage-safety
constraints (esp. donor-grouping never splitting a donor's lines)."""

import pandas as pd
import pytest

from modeling import folds as folds_mod
from modeling.folds import (
    donor_grouped_repeated_kfold,
    leave_one_donor_out,
    leave_one_line_out,
    load_lines_with_label,
    plain_repeated_kfold,
)


@pytest.fixture(scope="module")
def lines():
    return load_lines_with_label("published")


def test_load_lines_with_label_has_138_rows(lines):
    assert len(lines) == 138
    assert lines["cell_line"].nunique() == 138
    assert set(lines.columns) >= {"cell_line", "donor", "diff_efficiency", "success"}


def test_load_lines_with_label_donor_count(lines):
    # confirmed elsewhere in this project: 20 distinct donors among the 138 lines
    assert lines["donor"].nunique() == 20


def _assert_disjoint_and_complete(fold, all_lines: set[str]):
    train = set(fold.train_lines)
    test = set(fold.test_lines)
    assert train & test == set(), f"train/test overlap in fold {fold}"
    assert train | test == all_lines, f"fold {fold} doesn't cover all 138 lines"


def test_plain_repeated_kfold_fold_count_and_coverage(lines):
    all_lines = set(lines["cell_line"])
    folds = plain_repeated_kfold(lines, n_splits=5, n_repeats=3, seed=0)
    assert len(folds) == 5 * 3
    for f in folds:
        _assert_disjoint_and_complete(f, all_lines)


def test_donor_grouped_repeated_kfold_fold_count_and_coverage(lines):
    all_lines = set(lines["cell_line"])
    folds = donor_grouped_repeated_kfold(lines, n_splits=5, n_repeats=3, seed=0)
    assert len(folds) == 5 * 3
    for f in folds:
        _assert_disjoint_and_complete(f, all_lines)


def test_donor_grouped_repeated_kfold_never_splits_a_donor(lines):
    """Critical leakage-safety property: for every donor_grouped fold, no
    donor may have lines in both train and test."""
    donor_of = dict(zip(lines["cell_line"], lines["donor"]))
    folds = donor_grouped_repeated_kfold(lines, n_splits=5, n_repeats=10, seed=0)
    assert len(folds) == 50

    violations = []
    for f in folds:
        train_donors = {donor_of[c] for c in f.train_lines}
        test_donors = {donor_of[c] for c in f.test_lines}
        overlap = train_donors & test_donors
        if overlap:
            violations.append((f.repeat, f.fold, overlap))

    assert violations == [], f"donor(s) split across train/test: {violations}"


def test_leave_one_line_out_produces_138_singleton_folds(lines):
    all_lines = set(lines["cell_line"])
    folds = leave_one_line_out(lines)
    assert len(folds) == 138
    seen_test_lines = set()
    for f in folds:
        assert len(f.test_lines) == 1
        assert len(f.train_lines) == 137
        _assert_disjoint_and_complete(f, all_lines)
        seen_test_lines.add(f.test_lines[0])
    # every line held out exactly once across all folds
    assert seen_test_lines == all_lines


def test_leave_one_donor_out_produces_20_folds(lines):
    all_lines = set(lines["cell_line"])
    folds = leave_one_donor_out(lines)
    assert len(folds) == 20
    donor_of = dict(zip(lines["cell_line"], lines["donor"]))
    seen_test_donors = set()
    for f in folds:
        _assert_disjoint_and_complete(f, all_lines)
        test_donors = {donor_of[c] for c in f.test_lines}
        assert len(test_donors) == 1, "lodo fold's test set should be exactly one donor's lines"
        seen_test_donors |= test_donors
    assert seen_test_donors == set(lines["donor"])


# ---------------------------------------------------------------------------
# Label-variant registry and join contracts
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("name", sorted(folds_mod.VARIANTS))
def test_every_registered_variant_satisfies_its_own_contract(name):
    """Each variant declares cohort, label, threshold and expected line count as
    one bundle; this checks the bundle is self-consistent against the real
    files. It is what gives 136 the same force the literal 138 has elsewhere,
    and stops a future variant being registered with a wrong count."""
    v = folds_mod.get_variant(name)
    assert v.cohort_csv.exists(), f"{name}: missing {v.cohort_csv}"
    assert v.label_csv.exists(), f"{name}: missing {v.label_csv}"

    lines = folds_mod.load_lines_with_label(name)
    assert len(lines) == v.n_lines
    assert lines["cell_line"].is_unique
    assert lines["donor"].nunique() == 20
    assert (lines["success"] == (lines["diff_efficiency"] >= v.threshold)).all()


def test_variant_suffixes_are_distinct_so_two_variants_cannot_share_a_path():
    suffixes = [v.suffix for v in folds_mod.VARIANTS.values()]
    assert len(suffixes) == len(set(suffixes)), f"duplicate variant suffixes: {suffixes}"


def test_unknown_variant_name_raises_rather_than_falling_back():
    # A typo must not quietly resolve to the published default and produce a
    # different analysis than the operator asked for.
    with pytest.raises(KeyError, match="unknown label variant"):
        folds_mod.load_lines_with_label("da_untreeted")


def _variant_with(tmp_path, cohort, label, n_lines):
    c, lab = tmp_path / "cohort.csv", tmp_path / "label.csv"
    cohort.to_csv(c, index=False)
    label.to_csv(lab, index=False)
    return folds_mod.LabelVariant(suffix="_t", cohort_csv=c, label_csv=lab,
                                  threshold=0.2, n_lines=n_lines, description="fixture")


def _cohort(names):
    return pd.DataFrame({"cell_line": names, "pool": ["p"] * len(names),
                         "donor": [f"d{n}" for n in names]})


def _label(names, effs=None):
    return pd.DataFrame({"cell_line": names,
                         "diff_efficiency": effs if effs is not None else [0.5] * len(names)})


def test_duplicate_label_row_raises(tmp_path):
    v = _variant_with(tmp_path, _cohort(["a", "b"]), _label(["a", "b", "b"]), 2)
    with pytest.raises(ValueError, match="more than one row"):
        folds_mod.load_lines_with_label(v)


def test_missing_label_row_raises(tmp_path):
    v = _variant_with(tmp_path, _cohort(["a", "b"]), _label(["a"]), 2)
    with pytest.raises(ValueError, match="no label"):
        folds_mod.load_lines_with_label(v)


def test_a_duplicate_and_a_missing_row_that_cancel_on_count_still_raise(tmp_path):
    """Why a count check alone is not enough: one line duplicated and another
    absent leaves the row count exactly right, so an `n_lines` assertion passes
    while the modelling population is silently wrong."""
    v = _variant_with(tmp_path, _cohort(["a", "b"]), _label(["a", "a"]), 2)
    with pytest.raises(ValueError, match="more than one row|no label"):
        folds_mod.load_lines_with_label(v)


def test_out_of_range_efficiency_raises(tmp_path):
    v = _variant_with(tmp_path, _cohort(["a", "b"]), _label(["a", "b"], [0.5, 1.7]), 2)
    with pytest.raises(ValueError, match="out-of-range"):
        folds_mod.load_lines_with_label(v)


def test_a_line_mapped_to_two_donors_raises(tmp_path):
    cohort = pd.DataFrame({"cell_line": ["a", "a"], "pool": ["p1", "p2"], "donor": ["d1", "d2"]})
    v = _variant_with(tmp_path, cohort, _label(["a"]), 1)
    with pytest.raises(ValueError, match="more than one donor"):
        folds_mod.load_lines_with_label(v)
