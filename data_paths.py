"""
Where the expression data lives.

The three h5 files (day11, day30, day52) are multi-GB and are NOT in this
repository. Their location was previously hardcoded to one machine's home
directory in eight separate files, which meant nobody else could run anything
without editing source, and meant the path could drift between files.

Set the directory once, via the environment:

    export JERBER_DATA_DIR=/path/to/your/data
    uv run python eda/003_qualifying_cell_lines.py

Defaults to ~/Documents/2021_jerber if unset. The files are expected to be
named day11.h5, day30.h5, day52.h5 -- the names they ship with from
https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-10018 (Jerber
et al. 2021).

Nothing here imports from this project, so both the numbered eda/ scripts
and the modeling/ package can rely on it without an import cycle.
"""

import os
from pathlib import Path

ENV_VAR = "JERBER_DATA_DIR"
DEFAULT_DIR = "~/Documents/2021_jerber"

DATA_DIR = Path(os.environ.get(ENV_VAR, DEFAULT_DIR)).expanduser()

TIMEPOINT_FILENAMES = {"D11": "day11.h5", "D30": "day30.h5", "D52": "day52.h5"}


def data_file(timepoint: str) -> str:
    """Absolute path to one timepoint's h5.

    Returns a str rather than a Path because h5py.File and the streaming
    helpers throughout this project take strings, and several of them put the
    value straight into an f-string for logging.

    Deliberately does NOT check the file exists: the numbered scripts are
    imported by the test suite and by modeling/features.py for their helper
    functions alone, and raising at import time on a machine without the data
    would break both. Missing files surface at open time, where the error
    names the path."""
    try:
        return str(DATA_DIR / TIMEPOINT_FILENAMES[timepoint])
    except KeyError:
        raise KeyError(
            f"unknown timepoint {timepoint!r}; expected one of {sorted(TIMEPOINT_FILENAMES)}"
        ) from None


def describe() -> str:
    """One line for a script's startup output, so a run records where its
    inputs came from rather than leaving it to be inferred."""
    found = sum((DATA_DIR / n).exists() for n in TIMEPOINT_FILENAMES.values())
    return (f"data dir: {DATA_DIR} ({found}/{len(TIMEPOINT_FILENAMES)} files present; "
            f"override with {ENV_VAR})")
