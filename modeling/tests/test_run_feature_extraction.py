"""Fast tests for the per-fold extraction runner's provenance guards.

These cover the second half of an audit finding: refusing to RESUME onto a
table of unknown provenance is not enough on its own, because the append path
reindexes each chunk to the header already on disk. Selecting those columns
silently drops any field the header predates, so appending to a legacy table
would strip `fit_population` off the new rows and leave a file with no
provenance anywhere and no error raised.
"""

import pandas as pd
import pytest

from modeling import folds
from modeling import run_feature_extraction as extraction


def _stub_inputs(monkeypatch, tmp_path):
    """Make main() runnable without touching multi-GB h5 files. PCA is
    stubbed to fail loudly: these tests assert the guards fire BEFORE any
    expensive fitting work begins."""
    monkeypatch.setattr(extraction, "OUT_DIR", tmp_path)
    monkeypatch.setattr(extraction, "load_lines_with_label", lambda _variant: pd.DataFrame(
        {"cell_line": ["a", "b"], "donor": ["d1", "d2"],
         "diff_efficiency": [0.1, 0.5], "success": [False, True]}
    ))
    monkeypatch.setattr(extraction, "plain_repeated_kfold",
                        lambda *a, **k: [extraction.Fold("plain", 0, 0, ("a",), ("b",))]
                        if hasattr(extraction, "Fold") else [])
    monkeypatch.setattr(extraction, "donor_grouped_repeated_kfold", lambda *a, **k: [])
    monkeypatch.setattr(extraction, "load_cell_metadata", lambda timepoint: pd.DataFrame(
        {"cell_line": ["a", "b"], "pool": ["p", "p"], "celltype": ["FPP", "FPP"]}
    ))
    monkeypatch.setattr(extraction, "compute_proportion_features", lambda meta: pd.DataFrame(
        {"cell_line": ["a", "b"], "pool": ["p", "p"], "phat_FPP": [1.0, 1.0]}
    ))
    # The cohort path now comes from the registry, not a module constant, so
    # register a throwaway variant rather than patching a global. This also
    # exercises the real resolution path instead of bypassing it.
    qualifying = tmp_path / "qualifying.csv"
    pd.DataFrame({"cell_line": ["a", "b"], "pool": ["p", "p"]}).to_csv(qualifying, index=False)
    monkeypatch.setitem(folds.VARIANTS, "test_tmp", folds.LabelVariant(
        suffix="", cohort_csv=qualifying, label_csv=qualifying,
        threshold=0.2, n_lines=2, description="throwaway fixture variant",
    ))
    return "test_tmp"


def test_resume_refuses_a_table_with_no_recorded_fit_population(monkeypatch, tmp_path):
    # "No fit_population column" means provenance UNKNOWN, not provenance
    # compatible. A table predating the column could have been fit on either
    # population, so resuming onto it can produce exactly the mixed basis the
    # check exists to prevent.
    variant = _stub_inputs(monkeypatch, tmp_path)
    monkeypatch.setattr(extraction, "compute_pca_features_for_fold",
                        lambda *a, **k: pytest.fail("must reject before any PCA work"))

    out = tmp_path / "fold_data" / "fold_features_D30_t.csv"
    out.parent.mkdir()
    pd.DataFrame({"scheme": ["plain"], "repeat": [0], "fold": [0]}).to_csv(out, index=False)

    with pytest.raises(ValueError, match="fit_population|label_variant|UNRECORDED"):
        extraction.main(timepoint="D30", label_variant=variant, n_splits=2, n_repeats=1,
                        out_suffix="_t", include_loco_lodo=False,
                        restrict_fit_to_qualifying=True, resume=True)


def _fold_rows(scheme, repeat, fold, n_combos):
    return pd.DataFrame({
        "scheme": [scheme] * n_combos, "repeat": [repeat] * n_combos, "fold": [fold] * n_combos,
        "cell_line": [f"line{i}" for i in range(n_combos)], "pool": ["p"] * n_combos,
    })


def test_a_partially_written_fold_is_not_counted_as_done():
    """The defect this guards. Resume used to treat the PRESENCE of a
    (scheme, repeat, fold) key as proof the whole fold was written. A process
    killed mid-append leaves some of that fold's rows on disk; it was then
    skipped forever, and the finished table was silently short a few rows of
    one fold -- plausible, non-crashing, and wrong."""
    prev = pd.concat([_fold_rows("plain", 0, 0, 5), _fold_rows("plain", 0, 1, 3)])
    complete, partial = extraction.complete_folds(prev, expected_combos=5)
    assert complete == {("plain", 0, 0)}
    assert partial == {("plain", 0, 1)}, "the short fold must be recomputed, not skipped"


def test_a_fold_with_duplicate_combos_is_treated_as_partial():
    # Right row count, wrong content: a retry that appended some rows twice.
    dup = pd.concat([_fold_rows("plain", 0, 0, 4), _fold_rows("plain", 0, 0, 1)])
    complete, partial = extraction.complete_folds(dup, expected_combos=5)
    assert complete == set() and partial == {("plain", 0, 0)}


def test_complete_folds_accepts_a_fully_written_fold():
    # The check must be able to FAIL in the other direction too -- a guard that
    # rejects everything would pass the two tests above while blocking all work.
    prev = _fold_rows("loco", 0, 7, 5)
    complete, partial = extraction.complete_folds(prev, expected_combos=5)
    assert complete == {("loco", 0, 7)} and partial == set()


def test_a_truncated_final_line_is_dropped_rather_than_failing_the_read(tmp_path):
    out = tmp_path / "t.csv"
    _fold_rows("plain", 0, 0, 3).to_csv(out, index=False)
    with out.open("a") as fh:
        fh.write("plain,0,0,line3")  # cut off mid-record, as a kill would leave it
    got = extraction._read_resumable(out)
    assert len(got) == 3, "the complete rows must survive; the partial one is dropped"


def test_appending_a_chunk_with_an_extra_column_raises_instead_of_dropping_it(tmp_path):
    """The root-cause guard. A legacy table lacking fit_population must not
    silently swallow it off the appended rows -- the old code selected the
    header's columns, which discarded the field with no error and produced a
    file with no provenance anywhere."""
    out = tmp_path / "table.csv"
    pd.DataFrame({"a": [1], "b": [2]}).to_csv(out, index=False)
    chunk = pd.DataFrame({"a": [3], "b": [4], "fit_population": ["qualifying_only"]})

    with pytest.raises(ValueError, match="fit_population"):
        extraction.append_fold_rows(out, chunk)
    assert len(pd.read_csv(out)) == 1, "nothing may be written when the schema disagrees"


def test_appending_a_chunk_missing_a_column_raises_instead_of_blanking_it(tmp_path):
    out = tmp_path / "table.csv"
    pd.DataFrame({"a": [1], "fit_population": ["all_cells"]}).to_csv(out, index=False)

    with pytest.raises(ValueError, match="fit_population"):
        extraction.append_fold_rows(out, pd.DataFrame({"a": [2]}))


def test_append_reorders_to_the_on_disk_header_and_preserves_values(tmp_path):
    # Reordering IS allowed -- only the column set is fixed. Values must land
    # under their own headers, not be positionally misaligned.
    out = tmp_path / "table.csv"
    pd.DataFrame({"a": [1], "b": [10]}).to_csv(out, index=False)
    extraction.append_fold_rows(out, pd.DataFrame({"b": [20], "a": [2]}))

    got = pd.read_csv(out)
    assert got["a"].tolist() == [1, 2]
    assert got["b"].tolist() == [10, 20]
