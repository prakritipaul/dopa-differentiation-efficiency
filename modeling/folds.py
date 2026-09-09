"""
Fold construction for the D11 -> D52 modeling pipeline.

All fold schemes are timepoint-agnostic (they operate on cell_line/donor/
label metadata only, never touch expression data) and are used for BOTH
the regression and classification tasks, so both tasks share identical
fold assignments (per modeling/README.md).

Only `plain` and `donor_grouped` are STRATIFIED by the binary
success/failure label (diff_efficiency >= 0.2). `loco` and `lodo` are
exhaustive enumerations -- one fold per line, one fold per donor -- so
their folds are fully determined and there is no freedom left to
stratify with. (An earlier version of this docstring claimed all four
were stratified, which was false.)

Four schemes:
  - plain: ordinary repeated stratified K-fold over the 138 lines
    (grouping by cell_line here is a no-op -- each row already is one
    whole line, so this is just RepeatedStratifiedKFold).
  - donor_grouped: repeated stratified K-fold with lines sharing a donor
    forced into the same fold (StratifiedGroupKFold; sklearn has no
    "repeated" variant of this, so repeats are done manually across
    different random_state values).
  - loco: leave-one-cell-line-out (138 folds, singleton test sets).
  - lodo: leave-one-donor-out (one fold per donor, test set = every line
    from that donor).
"""

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.model_selection import RepeatedStratifiedKFold, StratifiedGroupKFold

REPO_ROOT = Path(__file__).parent.parent
QUALIFYING_COMBOS_CSV = REPO_ROOT / "metadata_eda" / "cohort/qualifying_cell_line_pool_min10_per_timepoint.csv"
LABEL_CSV = REPO_ROOT / "metadata_eda" / "cohort/d52_diff_efficiency_label.csv"
SUCCESS_THRESHOLD = 0.2


@dataclass(frozen=True)
class Fold:
    scheme: str  # "plain" | "donor_grouped" | "loco" | "lodo"
    repeat: int  # 0 for loco/lodo (no repeats -- deterministic by construction)
    fold: int
    train_lines: tuple[str, ...]
    test_lines: tuple[str, ...]


def load_lines_with_label() -> pd.DataFrame:
    """One row per qualifying cell line: cell_line, donor, diff_efficiency, success."""
    combos = pd.read_csv(QUALIFYING_COMBOS_CSV)[["cell_line", "donor"]].drop_duplicates()
    label = pd.read_csv(LABEL_CSV)[["cell_line", "diff_efficiency"]]
    lines = combos.merge(label, on="cell_line", how="inner")
    lines["success"] = lines["diff_efficiency"] >= SUCCESS_THRESHOLD
    lines = lines.sort_values("cell_line").reset_index(drop=True)
    assert len(lines) == 138, f"expected 138 qualifying lines with a label, got {len(lines)}"
    return lines


def plain_repeated_kfold(
    lines: pd.DataFrame, n_splits: int = 5, n_repeats: int = 10, seed: int = 0
) -> list[Fold]:
    rskf = RepeatedStratifiedKFold(n_splits=n_splits, n_repeats=n_repeats, random_state=seed)
    cell_line = lines["cell_line"].to_numpy()
    y = lines["success"].to_numpy()
    folds = []
    for i, (train_idx, test_idx) in enumerate(rskf.split(cell_line, y)):
        repeat, fold = divmod(i, n_splits)
        folds.append(Fold("plain", repeat, fold, tuple(cell_line[train_idx]), tuple(cell_line[test_idx])))
    return folds


def donor_grouped_repeated_kfold(
    lines: pd.DataFrame, n_splits: int = 5, n_repeats: int = 10, seed: int = 0
) -> list[Fold]:
    cell_line = lines["cell_line"].to_numpy()
    y = lines["success"].to_numpy()
    groups = lines["donor"].to_numpy()
    rng = np.random.default_rng(seed)
    folds = []
    for repeat in range(n_repeats):
        sgkf = StratifiedGroupKFold(
            n_splits=n_splits, shuffle=True, random_state=int(rng.integers(0, 2**31 - 1))
        )
        for fold, (train_idx, test_idx) in enumerate(sgkf.split(cell_line, y, groups)):
            folds.append(
                Fold("donor_grouped", repeat, fold, tuple(cell_line[train_idx]), tuple(cell_line[test_idx]))
            )
    return folds


def leave_one_line_out(lines: pd.DataFrame) -> list[Fold]:
    cell_line = lines["cell_line"].to_numpy()
    return [
        Fold("loco", 0, i, tuple(c for c in cell_line if c != held_out), (held_out,))
        for i, held_out in enumerate(cell_line)
    ]


def leave_one_donor_out(lines: pd.DataFrame) -> list[Fold]:
    donors = sorted(lines["donor"].unique())
    folds = []
    for i, donor in enumerate(donors):
        test_ids = tuple(lines.loc[lines["donor"] == donor, "cell_line"])
        train_ids = tuple(lines.loc[lines["donor"] != donor, "cell_line"])
        folds.append(Fold("lodo", 0, i, train_ids, test_ids))
    return folds


def all_folds(
    lines: pd.DataFrame, n_splits: int = 5, n_repeats: int = 10, seed: int = 0
) -> dict[str, list[Fold]]:
    return {
        "plain": plain_repeated_kfold(lines, n_splits, n_repeats, seed),
        "donor_grouped": donor_grouped_repeated_kfold(lines, n_splits, n_repeats, seed),
        "loco": leave_one_line_out(lines),
        "lodo": leave_one_donor_out(lines),
    }


def persist_folds(folds_by_scheme: dict[str, list[Fold]], out_path: Path) -> pd.DataFrame:
    """One row per (scheme, repeat, fold, cell_line, split) -- the actual
    line-to-fold assignment, not just the seed (protects against
    sklearn/library-version drift changing what a seed reproduces)."""
    rows = []
    for folds in folds_by_scheme.values():
        for f in folds:
            for cl in f.train_lines:
                rows.append(
                    {"scheme": f.scheme, "repeat": f.repeat, "fold": f.fold, "cell_line": cl, "split": "train"}
                )
            for cl in f.test_lines:
                rows.append(
                    {"scheme": f.scheme, "repeat": f.repeat, "fold": f.fold, "cell_line": cl, "split": "test"}
                )
    df = pd.DataFrame(rows)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, index=False)
    return df
