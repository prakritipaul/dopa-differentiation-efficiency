"""Shared pytest configuration.

Provides the `requires_data` marker. The three h5 files are multi-GB and are
not in the repository, so a handful of tests cannot run on a fresh clone.
Those tests are marked, and skipped automatically when the files are absent.

The point is that `uv run pytest modeling/` should pass for someone who has
just cloned this and has no data -- which it did not before. An independent
runtime check caught it: with JERBER_DATA_DIR pointed at an empty directory,
two tests failed on a missing day11.h5 while the README claimed the suite
runs without the expression data.

The skip is conditioned on the files actually being MISSING, not on the
marker alone. On a machine that has the data these tests still run, so the
coverage is not quietly lost on the one machine where it matters.
"""

import pytest

from data_paths import DATA_DIR, ENV_VAR, TIMEPOINT_FILENAMES


def _missing() -> list[str]:
    return [n for n in TIMEPOINT_FILENAMES.values() if not (DATA_DIR / n).exists()]


def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "requires_data: needs the real h5 files; skipped when they are absent "
        f"(looked for under ${ENV_VAR}, currently {DATA_DIR})",
    )


def pytest_collection_modifyitems(config, items):
    missing = _missing()
    if not missing:
        return
    skip = pytest.mark.skip(
        reason=f"needs the expression data: {', '.join(missing)} not found in {DATA_DIR}. "
               f"Set ${ENV_VAR} to the directory holding them."
    )
    for item in items:
        if "requires_data" in item.keywords:
            item.add_marker(skip)
