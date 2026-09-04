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
assumed. Restricting the fit to qualifying combos would require
re-extracting all 258 folds and would change every PC-derived number; not
done, since there is no correctness argument for it.

Output: one row per (scheme, repeat, fold, cell_line, pool) combo, with
PC1..PC10 (pool-level means), phat_FPP/phat_NB/phat_P_FPP, and a `split`
column (train/test) -- saved to modeling/fold_features_D11.csv. This is
the leakage-safe per-fold feature table harness.py will consume
(pool-correction, if applied, happens downstream using the `split`
column to compute training-only pool means; dropped from the active
pipeline per modeling/pool_correction_investigation.md, but still
supported).
"""

import time
from pathlib import Path

import pandas as pd

from modeling.features import (
    compute_pca_features_for_fold,
    compute_proportion_features,
    load_cell_metadata,
)
from modeling.folds import (
    donor_grouped_repeated_kfold,
    leave_one_donor_out,
    leave_one_line_out,
    load_lines_with_label,
    persist_folds,
    plain_repeated_kfold,
)

OUT_DIR = Path(__file__).parent
QUALIFYING_COMBOS_CSV = Path(__file__).parent.parent / "metadata_eda" / "qualifying_cell_line_pool_min10_per_timepoint.csv"
N_SPLITS = 5
N_REPEATS = 10
N_PCS = 10
SEED = 0


def main(
    n_splits: int = N_SPLITS,
    n_repeats: int = N_REPEATS,
    out_suffix: str = "",
    include_loco_lodo: bool = True,
) -> None:
    lines = load_lines_with_label()
    folds = {
        "plain": plain_repeated_kfold(lines, n_splits, n_repeats, SEED),
        "donor_grouped": donor_grouped_repeated_kfold(lines, n_splits, n_repeats, SEED),
    }
    if include_loco_lodo:
        folds["loco"] = leave_one_line_out(lines)
        folds["lodo"] = leave_one_donor_out(lines)

    persist_folds(folds, OUT_DIR / f"fold_assignments{out_suffix}.csv")
    total_folds = sum(len(fs) for fs in folds.values())
    print(f"{total_folds} folds to process ({', '.join(f'{k}={len(v)}' for k, v in folds.items())})")

    meta = load_cell_metadata("D11")
    qualifying = pd.read_csv(QUALIFYING_COMBOS_CSV)[["cell_line", "pool"]]
    meta_q = meta.merge(qualifying, on=["cell_line", "pool"], how="inner")

    # Proportions: no per-fold refit needed (pre-existing labels, no fitting step).
    props = compute_proportion_features(meta_q)

    all_rows = []
    t_start = time.time()
    fold_count = 0
    for scheme, fold_list in folds.items():
        for f in fold_list:
            fold_count += 1
            t0 = time.time()
            held_out = set(f.test_lines)
            pcs = compute_pca_features_for_fold("D11", held_out, meta=meta, n_pcs=N_PCS)
            pcs_q = pcs.merge(qualifying, on=["cell_line", "pool"], how="inner")

            pc_cols = [c for c in pcs_q.columns if c.startswith("PC")]
            pool_level_pcs = pcs_q.groupby(["cell_line", "pool"], observed=True)[pc_cols].mean().reset_index()

            combo = pool_level_pcs.merge(props, on=["cell_line", "pool"], how="inner")
            combo["scheme"] = scheme
            combo["repeat"] = f.repeat
            combo["fold"] = f.fold
            combo["split"] = combo["cell_line"].apply(lambda c: "test" if c in held_out else "train")
            all_rows.append(combo)

            elapsed = time.time() - t0
            total_elapsed = time.time() - t_start
            eta = total_elapsed / fold_count * (total_folds - fold_count)
            print(
                f"[{fold_count}/{total_folds}] {scheme} repeat={f.repeat} fold={f.fold} "
                f"took {elapsed:.1f}s (elapsed {total_elapsed / 60:.1f}min, ETA {eta / 60:.1f}min)",
                flush=True,
            )

    result = pd.concat(all_rows, ignore_index=True)
    out_path = OUT_DIR / f"fold_features_D11{out_suffix}.csv"
    result.to_csv(out_path, index=False)
    print(f"\nSaved {len(result)} rows to {out_path}")


if __name__ == "__main__":
    main()
