"""
Extract per-fold D11 PCA + proportion features for all four CV schemes:
plain, donor_grouped (repeated stratified K-fold), plus LOCO and LODO
(full leave-one-line-out / leave-one-donor-out). The first-pass prototype
run (5 repeats, plain+donor_grouped only) validated the pipeline; this is
the full-scope run per modeling/README.md's original design.

For each fold: PCA is refit excluding that fold's held-out lines (and
pool11, frozen rule) and projected onto every cell (leakage safety);
proportions don't need per-fold refitting (no fitted parameter, just
tabulating existing celltype labels) so they're computed once and reused
across all folds.

FITTING POPULATION -- a deliberate deviation, stated explicitly because
an earlier version of this docstring described it wrongly. The PCA is fit
on all D11 cells except (a) the fold's held-out lines and (b) pool11 --
which INCLUDES ~26,400 cells outside the qualifying (cell_line, pool)
combos, ~25,400 of them from 39 cell lines that are not among the 138 and
therefore never appear in any train or test set. Full metadata is passed
to the PCA step below and the qualifying filter is applied only
afterwards, to the resulting coordinates.

This is not leakage: the held-out mask is applied by cell_line, so every
cell of a held-out line is excluded regardless of pool, and the extra
cells belong to lines that are never predicted on. Using additional
unlabeled cells to fit a basis is a legitimate (transductive) choice and
gives the basis more data. It is recorded here rather than silently
assumed.

Both variants are now produced and kept side by side for comparison:
`restrict_fit_to_qualifying=True` limits the fit to the qualifying combos
(226,874 cells vs 253,381 total), writing to a suffixed output. Note the
qualifying filter is a property of the (cell_line, pool) PAIR across all
three timepoints, so the strict variant also drops 1,109 D11 cells from 4
combos whose LINE is among the 138 but whose pool fails on D30/D52
availability -- e.g. HPSI0714i-kute_5/pool13 has 497 good D11 cells but
zero at D30.

Output: one row per (scheme, repeat, fold, cell_line, pool) combo, with
PC1..PC10 (pool-level means), phat_FPP/phat_NB/phat_P_FPP, and a `split`
column (train/test) -- saved to modeling/fold_features_D11.csv. This is
the leakage-safe per-fold feature table harness.py will consume
(pool-correction, if applied, happens downstream using the `split`
column to compute training-only pool means; dropped from the active
pipeline per modeling/pool_correction_investigation.md, but still
supported).
"""

import os
import time
from io import StringIO
from pathlib import Path

import numpy as np
import pandas as pd

from modeling.features import (
    TIMEPOINT_FILES,
    compute_pca_features_for_fold,
    compute_proportion_features,
    load_cell_metadata,
)
from modeling.folds import (
    VARIANTS,
    donor_grouped_repeated_kfold,
    get_variant,
    leave_one_donor_out,
    leave_one_line_out,
    load_lines_with_label,
    persist_folds,
    plain_repeated_kfold,
)

# There is deliberately no QUALIFYING_COMBOS_CSV here any more. It used to be a
# SECOND hardcoded copy of the cohort path, independent of folds.py's -- so a run
# could take its labels from one cohort and the combos its features are computed
# on from another, with nothing raising. The cohort now comes from the resolved
# LabelVariant, leaving nowhere for the two to disagree.
OUT_DIR = Path(__file__).parent
N_SPLITS = 5
N_REPEATS = 10
N_PCS = 10
SEED = 0


def append_fold_rows(out_path: Path, chunk: pd.DataFrame) -> None:
    """Append one fold's rows, creating the file on first write.

    Column ORDER is fixed by the first write, so later chunks are reindexed to
    the header already on disk rather than trusting dict ordering to stay
    stable. But reindexing must never change the column SET in either
    direction: selecting the header's columns alone silently drops any field
    the header predates, which is how appending to a legacy table strips
    `fit_population` off the new rows and leaves a file carrying no provenance
    at all, with no error raised anywhere. A schema difference is a bug in the
    caller, not something to paper over mid-run."""
    if not out_path.exists():
        chunk.to_csv(out_path, index=False)
        return

    header_cols = pd.read_csv(out_path, nrows=0).columns.tolist()
    extra = [c for c in chunk.columns if c not in header_cols]
    absent = [c for c in header_cols if c not in chunk.columns]
    if extra or absent:
        raise ValueError(
            f"{out_path} has a different schema than the rows being appended "
            f"(only in new rows: {extra}; only in file: {absent}). Appending would "
            f"silently drop or blank columns. Rerun with --no-resume to rebuild."
        )
    chunk[header_cols].to_csv(out_path, mode="a", header=False, index=False)


def complete_folds(
    prev: pd.DataFrame, expected_keys: set[tuple[str, str]]
) -> tuple[set[tuple], set[tuple]]:
    """Split the folds already on disk into (complete, partial).

    Resume previously treated the mere PRESENCE of a (scheme, repeat, fold)
    key as proof the whole fold had been written. It is not: a process killed
    part-way through appending leaves some of that fold's rows on disk, and on
    restart the fold is skipped forever. Nothing raises, and the finished table
    is quietly short a few rows -- this repo's characteristic failure.

    The contract is EXACT KEY-SET EQUALITY, not a row count. Counting is not
    enough, and the difference is not hypothetical: a fold missing one expected
    (cell_line, pool) and carrying one unexpected one has the right number of
    distinct keys and would pass a cardinality check while describing the wrong
    population. Independent review caught that in the first version of this
    function.

    A fold is complete only when it has exactly one row per expected key, no
    extras, and no nulls or non-finite values anywhere in it. Anything else is
    returned as partial so the caller can drop and recompute it -- malformed
    rows are never silently dropped on their own, because doing so would repair
    the symptom and hide whatever produced it."""
    keys = ["scheme", "repeat", "fold"]
    complete, partial = set(), set()
    for key, group in prev.groupby(keys, observed=True):
        key = (str(key[0]), int(key[1]), int(key[2]))
        got = set(map(tuple, group[["cell_line", "pool"]].itertuples(index=False)))
        numeric = group.select_dtypes(include="number")
        intact = (
            got == expected_keys
            and len(group) == len(expected_keys)      # one row per key, no duplicates
            and not group.isna().any().any()
            and bool(np.isfinite(numeric.to_numpy(dtype=float)).all())
        )
        (complete if intact else partial).add(key)
    return complete, partial


def _read_resumable(out_path: Path) -> pd.DataFrame:
    """Read a partially written table, tolerating a truncated final line.

    A crash mid-write can leave the last row cut off with TOO MANY fields for
    the parser, which raises and makes the whole file unreadable. Dropping that
    one line is safe -- the fold it belonged to is incomplete either way and
    complete_folds() will send it back for recomputation.

    The complementary case, a row with too FEW fields, is deliberately NOT
    handled here. pandas silently pads it with NaN, and an earlier version of
    this function dropped every null-bearing row on that basis. That was wrong
    twice over: it would also have discarded a genuinely missing value, hiding
    an upstream defect rather than surfacing it, and it is not a reliable
    truncation test anyway (a truncated numeric token can stay valid -- 1.2345
    cut to 1.23 is still a number). Nulls and non-finite values are instead
    judged by complete_folds(), which recomputes the whole fold.

    ponytail: a value truncated inside the FINAL column while staying
    syntactically valid is undetectable here. In practice the final column is
    `label_variant`, so a truncated one ("da_untre") fails the resume
    provenance check and raises. Per-fold atomic writes would close it
    properly; not worth it until this actually bites."""
    try:
        return pd.read_csv(out_path)
    except pd.errors.ParserError:
        text = out_path.read_text().splitlines()
        print(f"{out_path.name}: final line is unparseable, dropping it before resume")
        return pd.read_csv(StringIO("\n".join(text[:-1]) + "\n"))


def main(
    timepoint: str = "D11",
    label_variant: str = "published",
    n_splits: int = N_SPLITS,
    n_repeats: int = N_REPEATS,
    out_suffix: str = "",
    include_loco_lodo: bool = True,
    restrict_fit_to_qualifying: bool = False,
    resume: bool = True,
) -> None:
    """timepoint: "D11" or "D30" -- selects the h5 file, the cell-type
    annotation set (D11 has 3 types, D30 has 7), and the depth-outlier pool
    excluded from the HVG/PCA fit (features.DEPTH_OUTLIER_POOLS: pool11 at
    D11, pool5 at D30). It also names the output, so the two timepoints
    cannot overwrite each other. Fold assignments are NOT timepoint-specific
    -- they depend only on the 138 lines' donor/label metadata, so both
    timepoints get byte-identical folds and the D11-vs-D30 comparison is
    paired on the same held-out donors in every repeat.

    restrict_fit_to_qualifying: when True, the HVG/PCA fit is limited to
    cells in the qualifying (cell_line, pool) combos -- the study
    population the features are actually computed on. When False (the
    default, preserving the original run) cells outside those combos also
    contribute to the fit; see this module's docstring. Use a distinct
    out_suffix so the two variants sit side by side rather than one
    overwriting the other.

    resume: fold rows are appended to the output as each fold completes, so
    an interrupted run can be restarted and will skip folds already on
    disk. Pass resume=False to discard any existing output and start
    clean -- required if the fitting population changed, since otherwise
    rows computed under the old settings would be silently kept."""
    v = get_variant(label_variant)
    print(f"label variant: {label_variant} -- {v.description}")
    lines = load_lines_with_label(v)
    folds = {
        "plain": plain_repeated_kfold(lines, n_splits, n_repeats, SEED),
        "donor_grouped": donor_grouped_repeated_kfold(lines, n_splits, n_repeats, SEED),
    }
    if include_loco_lodo:
        folds["loco"] = leave_one_line_out(lines)
        folds["lodo"] = leave_one_donor_out(lines)

    # Fold assignments are shared across timepoints ON PURPOSE -- see the
    # docstring. Writing to the same un-suffixed path and refusing to change
    # it is what makes "D11 and D30 are evaluated on the same held-out
    # donors" a checked fact rather than an assumption: if anything ever
    # perturbs fold construction, the D30 run fails here instead of quietly
    # producing an unpaired comparison that still looks plausible.
    # The variant suffix is appended HERE, by code. An operator never types it,
    # so a new-label artifact cannot be aimed at a published path by mistake.
    out_suffix = f"{out_suffix}{v.suffix}"
    folds_path = OUT_DIR / f"fold_data/fold_assignments{out_suffix}.csv"
    if folds_path.exists():
        # Compare via a scratch file and NEVER overwrite the committed one:
        # the existing assignments are provenance for results already
        # published from them, so a mismatch must fail loudly with the
        # original still intact, not after it has been clobbered.
        scratch = folds_path.with_suffix(".regenerated.tmp.csv")
        persist_folds(folds, scratch)
        same = pd.read_csv(scratch).equals(pd.read_csv(folds_path))
        scratch.unlink()
        if not same:
            raise ValueError(
                f"regenerating folds for timepoint={timepoint} does not reproduce "
                f"{folds_path}. Fold assignments must be identical across "
                "timepoints for the D11-vs-D30 comparison to be paired. The "
                "existing file was left untouched; investigate before rerunning."
            )
        print(f"fold assignments reproduce {folds_path.name} exactly (paired with D11)")
    else:
        persist_folds(folds, folds_path)
    total_folds = sum(len(fs) for fs in folds.values())
    print(f"{total_folds} folds to process ({', '.join(f'{k}={len(v)}' for k, v in folds.items())})")

    meta = load_cell_metadata(timepoint)
    qualifying = pd.read_csv(v.cohort_csv)[["cell_line", "pool"]]
    meta_q = meta.merge(qualifying, on=["cell_line", "pool"], how="inner")

    # Proportions: no per-fold refit needed (pre-existing labels, no fitting step).
    props = compute_proportion_features(meta_q)

    out_path = OUT_DIR / f"fold_data/fold_features_{timepoint}{out_suffix}.csv"
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # Each fold's rows are appended as soon as they are computed, so a crash
    # at fold 250 of 258 costs one fold rather than the whole multi-hour run.
    # On restart, folds already present in the file are skipped.
    # Which cells were allowed to shape the HVG/PCA basis. Recorded on every
    # row so a resume cannot silently splice together folds computed under two
    # different fitting populations: the rows would look identical in schema
    # and the resulting table would be a plausible, non-crashing, wrong
    # feature basis -- the exact failure run_experiment.run_all already guards
    # against for its own inputs.
    fit_population = "qualifying_only" if restrict_fit_to_qualifying else "all_cells"

    done_keys: set[tuple] = set()
    if resume and out_path.exists():
        prev = _read_resumable(out_path)
        seen_label = set(prev.get("label_variant", pd.Series(dtype=object)).dropna().unique())
        if seen_label != {label_variant}:
            raise ValueError(
                f"{out_path} records label_variant="
                f"{sorted(seen_label) if seen_label else 'UNRECORDED (no label_variant column)'} "
                f"but this run uses {label_variant!r}. The outcome, cohort and threshold all "
                f"differ between variants, so resuming would mix two analyses in one table. "
                f"Rerun with --no-resume, or use a different --label-variant."
            )
        seen = set(prev.get("fit_population", pd.Series(dtype=object)).dropna().unique())
        # An EMPTY `seen` is not "compatible", it is "provenance unknown", and
        # those must not be conflated: a table predating this column could have
        # been fit on either population, so resuming onto it would produce
        # exactly the mixed basis this check exists to prevent. Legacy tables
        # (the committed D11 ones) are complete and never need resuming, so
        # refusing is the safe reading -- rerun with --no-resume to rebuild.
        if seen != {fit_population}:
            raise ValueError(
                f"{out_path} records fit_population="
                f"{sorted(seen) if seen else 'UNRECORDED (no fit_population column)'} "
                f"but this run fits on {fit_population!r}. Resuming would mix two feature "
                f"bases in one table. Rerun with --no-resume, or use a different --suffix."
            )
        expected_keys = set(map(tuple, qualifying[["cell_line", "pool"]].itertuples(index=False)))
        done_keys, partial = complete_folds(prev, expected_keys)
        if partial:
            # Rewrite without them rather than leaving them to be skipped. The
            # rewrite goes via a temp file and os.replace so an interruption
            # here cannot destroy the good folds too.
            print(f"discarding {len(partial)} partially written fold(s): {sorted(partial)}")
            keep = prev.set_index(["scheme", "repeat", "fold"]).index.isin(done_keys)
            tmp = out_path.with_suffix(".repair.tmp.csv")
            prev.loc[keep].to_csv(tmp, index=False)
            os.replace(tmp, out_path)
        print(f"resuming: {len(done_keys)} of {total_folds} folds already in {out_path.name} "
              f"(fit population {fit_population})")
    elif not resume and out_path.exists():
        out_path.unlink()

    rows_written = 0
    computed = 0  # folds actually computed this run; skipped ones must not
                  # count toward the rate, or the ETA is wildly wrong on resume
    t_start = time.time()
    fold_count = 0
    for scheme, fold_list in folds.items():
        for f in fold_list:
            fold_count += 1
            if (scheme, f.repeat, f.fold) in done_keys:
                continue
            t0 = time.time()
            held_out = set(f.test_lines)
            pcs = compute_pca_features_for_fold(
                timepoint, held_out, meta=meta, n_pcs=N_PCS,
                restrict_to_combos=qualifying if restrict_fit_to_qualifying else None,
            )
            pcs_q = pcs.merge(qualifying, on=["cell_line", "pool"], how="inner")

            pc_cols = [c for c in pcs_q.columns if c.startswith("PC")]
            pool_level_pcs = pcs_q.groupby(["cell_line", "pool"], observed=True)[pc_cols].mean().reset_index()

            combo = pool_level_pcs.merge(props, on=["cell_line", "pool"], how="inner")
            combo["scheme"] = scheme
            combo["repeat"] = f.repeat
            combo["fold"] = f.fold
            combo["split"] = combo["cell_line"].apply(lambda c: "test" if c in held_out else "train")
            combo["fit_population"] = fit_population
            combo["label_variant"] = label_variant

            # Append immediately; header only when creating the file. Column
            # order is fixed by the first write, so reindex every later chunk
            # to it rather than trusting dict ordering to stay stable.
            append_fold_rows(out_path, combo)
            rows_written += len(combo)
            computed += 1

            elapsed = time.time() - t0
            total_elapsed = time.time() - t_start
            eta = total_elapsed / computed * (total_folds - fold_count)
            print(
                f"[{fold_count}/{total_folds}] {scheme} repeat={f.repeat} fold={f.fold} "
                f"took {elapsed:.1f}s (elapsed {total_elapsed / 60:.1f}min, ETA {eta / 60:.1f}min)",
                flush=True,
            )

    total_rows = len(pd.read_csv(out_path, usecols=["scheme"]))
    print(f"\nSaved {total_rows} rows to {out_path} ({rows_written} written this run)")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--timepoint", default="D11", choices=sorted(TIMEPOINT_FILES))
    parser.add_argument("--label-variant", default="published", choices=sorted(VARIANTS),
                        help="which cohort+label+threshold bundle to use (see folds.VARIANTS)")
    parser.add_argument("--suffix", default="", help="suffix for the output table")
    parser.add_argument(
        "--restrict-fit-to-qualifying", action="store_true",
        help="fit HVGs/PCA only on cells in the qualifying (cell_line, pool) combos -- the study "
             "population itself. Without it, cells from lines outside the 138 also shape the basis. "
             "Changing this changes the features, so pair it with a distinct --suffix (_qualonly).",
    )
    parser.add_argument("--no-resume", action="store_true", help="discard existing output and recompute")
    args = parser.parse_args()
    main(
        timepoint=args.timepoint,
        label_variant=args.label_variant,
        out_suffix=args.suffix,
        restrict_fit_to_qualifying=args.restrict_fit_to_qualifying,
        resume=not args.no_resume,
    )
