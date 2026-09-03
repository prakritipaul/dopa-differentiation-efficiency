"""
Extract per-fold D11 PCA + proportion features for the plain and
donor_grouped CV schemes (first-pass prototype scope; LOCO/LODO
deferred -- see modeling/README.md "First-pass prototype scope").

For each fold: PCA is refit on that fold's training-line cells only
(excluding pool11, frozen rule) and projected onto every cell (leakage
safety); proportions don't need per-fold refitting (no fitted parameter,
just tabulating existing celltype labels) so they're computed once and
reused across all folds.

Output: one row per (scheme, repeat, fold, cell_line, pool) combo, with
PC1..PC10 (pool-level means), phat_FPP/phat_NB/phat_P_FPP, and a `split`
column (train/test) -- saved to modeling/fold_features_D11.csv. This is
the leakage-safe per-fold feature table harness.py will consume
(pool-correction, if applied, happens downstream using the `split`
column to compute training-only pool means).
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
    load_lines_with_label,
    persist_folds,
    plain_repeated_kfold,
)

OUT_DIR = Path(__file__).parent
QUALIFYING_COMBOS_CSV = Path(__file__).parent.parent / "metadata_eda" / "qualifying_cell_line_pool_min10_per_timepoint.csv"
N_SPLITS = 5
N_REPEATS = 5
N_PCS = 10
SEED = 0


def main(n_splits: int = N_SPLITS, n_repeats: int = N_REPEATS, out_suffix: str = "") -> None:
    lines = load_lines_with_label()
    folds = {
        "plain": plain_repeated_kfold(lines, n_splits, n_repeats, SEED),
        "donor_grouped": donor_grouped_repeated_kfold(lines, n_splits, n_repeats, SEED),
    }
    persist_folds(folds, OUT_DIR / f"fold_assignments_first_pass{out_suffix}.csv")
    total_folds = sum(len(fs) for fs in folds.values())
    print(f"{total_folds} folds to process (plain + donor_grouped, {n_repeats} repeats each)")

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
