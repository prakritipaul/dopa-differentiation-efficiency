"""Tests for modeling/folds.py -- fold construction and leakage-safety
constraints (esp. donor-grouping never splitting a donor's lines)."""

import pytest

from modeling.folds import (
    donor_grouped_repeated_kfold,
    leave_one_donor_out,
    leave_one_line_out,
    load_lines_with_label,
    plain_repeated_kfold,
)


@pytest.fixture(scope="module")
def lines():
    return load_lines_with_label()


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
