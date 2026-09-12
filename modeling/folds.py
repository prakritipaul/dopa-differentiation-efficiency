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
COHORT_DIR = REPO_ROOT / "metadata_eda" / "cohort"


@dataclass(frozen=True)
class LabelVariant:
    """Cohort file, label file, threshold and expected line count as ONE object.

    These four values are coupled: a label computed on one cohort, scored at
    another's threshold, checked against a third's line count, is a silently
    wrong analysis. The failure this prevents is not "someone forgot to pass an
    argument" -- it is "three of the four came from variant A and one from
    variant B". Bundling them makes that combination unconstructible.

    `suffix` is applied by CODE when naming derived artifacts, never typed by an
    operator, so a new-label run cannot land on a published path.
    """

    suffix: str
    cohort_csv: Path
    label_csv: Path
    threshold: float
    n_lines: int
    description: str


# Keyed by name; the key IS the name, so the two cannot disagree.
VARIANTS: dict[str, LabelVariant] = {
    "published": LabelVariant(
        suffix="",
        cohort_csv=COHORT_DIR / "qualifying_cell_line_pool_min10_per_timepoint.csv",
        label_csv=COHORT_DIR / "d52_diff_efficiency_label.csv",
        threshold=0.2,
        n_lines=138,
        description="(DA+Sert)/all D52 cells, ALL cells incl. rotenone-treated -- "
                    "reproduces Jerber et al.'s own definition exactly",
    ),
    "da_untreated": LabelVariant(
        suffix="_da_untreated",
        cohort_csv=COHORT_DIR / "qualifying_cell_line_pool_min10_per_timepoint_untreated.csv",
        label_csv=COHORT_DIR / "d52_diff_efficiency_label_da_untreated.csv",
        threshold=0.2,
        n_lines=136,
        description="DA/all D52 cells, UNTREATED cells only -- a deliberate "
                    "divergence from the published definition, not a bug fix",
    ),
}


def get_variant(variant: str | LabelVariant) -> LabelVariant:
    """Resolve a variant name. Unknown names raise and list what is valid --
    a typo must not fall back to the published default and quietly produce the
    wrong analysis."""
    if isinstance(variant, LabelVariant):
        return variant
    try:
        return VARIANTS[variant]
    except KeyError:
        raise KeyError(
            f"unknown label variant {variant!r}; valid names are {sorted(VARIANTS)}"
        ) from None


@dataclass(frozen=True)
class Fold:
    scheme: str  # "plain" | "donor_grouped" | "loco" | "lodo"
    repeat: int  # 0 for loco/lodo (no repeats -- deterministic by construction)
    fold: int
    train_lines: tuple[str, ...]
    test_lines: tuple[str, ...]


def load_lines_with_label(variant: str | LabelVariant) -> pd.DataFrame:
    """One row per qualifying cell line: cell_line, donor, diff_efficiency, success.

    `variant` is REQUIRED -- there is deliberately no internal default. A
    default here would let a call site that was never updated fall back to the
    published cohort while the rest of the run used another, and the existing
    tests would still pass because they only ever exercise the default. Entry
    points default to "published" at the CLI layer instead, where the choice is
    visible.

    The joins below are validated rather than assumed. An unchecked inner merge
    can lose or duplicate rows, and a duplicate row plus a missing row CANCEL
    OUT in a count check -- so `n_lines` alone does not prove the population is
    right."""
    v = get_variant(variant)

    combos = pd.read_csv(v.cohort_csv)
    if combos.duplicated(["cell_line", "pool"]).any():
        raise ValueError(f"{v.cohort_csv.name} has duplicate (cell_line, pool) keys")

    per_line = combos[["cell_line", "donor"]].drop_duplicates()
    if per_line["cell_line"].duplicated().any():
        clash = per_line.loc[per_line["cell_line"].duplicated(keep=False), "cell_line"].unique()
        raise ValueError(f"cell lines mapped to more than one donor in {v.cohort_csv.name}: {sorted(clash)}")

    label = pd.read_csv(v.label_csv)[["cell_line", "diff_efficiency"]]
    if label["cell_line"].duplicated().any():
        dup = label.loc[label["cell_line"].duplicated(keep=False), "cell_line"].unique()
        raise ValueError(f"{v.label_csv.name} has more than one row for {sorted(dup)}")

    # Set EQUALITY, not the subset an inner join would silently tolerate.
    missing = set(per_line["cell_line"]) - set(label["cell_line"])
    extra = set(label["cell_line"]) - set(per_line["cell_line"])
    if missing or extra:
        raise ValueError(
            f"cohort and label disagree for variant {v.suffix or 'published'!r}: "
            f"{len(missing)} cohort lines have no label ({sorted(missing)[:5]}...), "
            f"{len(extra)} labelled lines are not in the cohort ({sorted(extra)[:5]}...)"
        )

    lines = per_line.merge(label, on="cell_line", how="inner", validate="one_to_one")

    bad = ~lines["diff_efficiency"].between(0.0, 1.0) | lines["diff_efficiency"].isna()
    if bad.any():
        raise ValueError(
            f"{int(bad.sum())} lines have a non-finite or out-of-range diff_efficiency: "
            f"{lines.loc[bad, 'cell_line'].tolist()[:5]}"
        )

    lines["success"] = lines["diff_efficiency"] >= v.threshold
    lines = lines.sort_values("cell_line").reset_index(drop=True)
    # raise, not assert: an assert is stripped under -O, and this is a contract.
    if len(lines) != v.n_lines:
        raise ValueError(f"expected {v.n_lines} qualifying lines with a label, got {len(lines)}")
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
