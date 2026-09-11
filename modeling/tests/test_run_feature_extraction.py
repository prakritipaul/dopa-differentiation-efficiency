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

from modeling import run_feature_extraction as extraction


def _stub_inputs(monkeypatch, tmp_path):
    """Make main() runnable without touching multi-GB h5 files. PCA is
    stubbed to fail loudly: these tests assert the guards fire BEFORE any
    expensive fitting work begins."""
    monkeypatch.setattr(extraction, "OUT_DIR", tmp_path)
    monkeypatch.setattr(extraction, "load_lines_with_label", lambda: pd.DataFrame(
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
    qualifying = tmp_path / "qualifying.csv"
    pd.DataFrame({"cell_line": ["a", "b"], "pool": ["p", "p"]}).to_csv(qualifying, index=False)
    monkeypatch.setattr(extraction, "QUALIFYING_COMBOS_CSV", qualifying)


def test_resume_refuses_a_table_with_no_recorded_fit_population(monkeypatch, tmp_path):
    # "No fit_population column" means provenance UNKNOWN, not provenance
    # compatible. A table predating the column could have been fit on either
    # population, so resuming onto it can produce exactly the mixed basis the
    # check exists to prevent.
    _stub_inputs(monkeypatch, tmp_path)
    monkeypatch.setattr(extraction, "compute_pca_features_for_fold",
                        lambda *a, **k: pytest.fail("must reject before any PCA work"))

    out = tmp_path / "fold_data" / "fold_features_D30_t.csv"
    out.parent.mkdir()
    pd.DataFrame({"scheme": ["plain"], "repeat": [0], "fold": [0]}).to_csv(out, index=False)

    with pytest.raises(ValueError, match="fit_population|UNRECORDED"):
        extraction.main(timepoint="D30", n_splits=2, n_repeats=1, out_suffix="_t",
                        include_loco_lodo=False, restrict_fit_to_qualifying=True, resume=True)


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
