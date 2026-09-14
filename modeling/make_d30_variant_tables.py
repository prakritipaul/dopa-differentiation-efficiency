"""
Derive the D30 "progenitor balance" feature table from the full one.

WHY THIS VARIANT EXISTS
-----------------------
The D30 annotation set contains DA and Sert -- the same two cell types whose
D52 fraction IS the outcome. So a D30 model that keeps them is largely
answering "how much of the D52 phenotype is already established at D30?"
rather than "what forecasts it?". (The raw D30 DA+Sert fraction alone scores
ROC-AUC 0.990 against D52 success; see d30_single_feature_benchmarks.py.)

But D30 is not purely an early readout of the answer. Marker expression over
all 250,923 D30 cells splits the seven annotated types three ways:

  already arrived   DA, Sert                 48% of cells -- this IS the outcome
  still undecided   FPP, P_FPP               31% -- SOX2+/HES1+/VIM-high
                                                   progenitors, not yet
                                                   neuronal; P_FPP actively
                                                   cycling (MKI67, TOP2A)
  off-target        Epen1, U_Neur1, U_Neur2  20% -- Epen1 is ciliated
                                                   ependymal/choroid-plexus-like
                                                   (TTR, FOXJ1, PIFO, RSPH1);
                                                   U_Neur1/2 are post-mitotic
                                                   neurons (STMN2 ~40, no
                                                   proliferation, no DA/Sert
                                                   markers)

The undecided 31% have not committed, and what they become between D30 and D52
is genuine prediction -- the same thing the D11 model does, where all three
annotated types are upstream. This table isolates that question.

WHY RENORMALISE
---------------
Naively deleting phat_DA and phat_Sert does NOT remove them. The five
remaining proportions would then sum to (1 - DA - Sert), so their common scale
still encodes exactly the quantity being removed, and a model can recover it
from any two of them. Renormalising the five to sum to 1 removes the DA+Sert
magnitude and leaves the composition AMONG cells that have not become
target-like -- i.e. of the cells still in play, how many are competent
progenitors (FPP, P_FPP) versus committed off-target fates (Epen1, U_Neur1,
U_Neur2). That balance is the actual quantity of interest, hence the name.

The transcriptomic PCs are carried over untouched. They are NOT independent of
composition (expression reflects which cells are present), so this is a
sensitivity analysis, not a clean decomposition -- see modeling/README.md.

Cross-sectional caveat: marker expression establishes cell IDENTITY, not FATE.
That FPP/P_FPP carry progenitor and cell-cycle programs while U_Neur1/2 carry
neither is strong evidence, not lineage tracing, which this data cannot do.

Usage:
    uv run python -m modeling.make_d30_variant_tables

Reads  modeling/fold_data/fold_features_D30_qualonly.csv
Writes modeling/fold_data/fold_features_D30_progbal_qualonly.csv
"""

from pathlib import Path

import pandas as pd

OUT_DIR = Path(__file__).parent

# The outcome-defining types, removed here. Named for what they are rather
# than listed inline, because the same split drives
# d30_single_feature_benchmarks.py and the two must not drift apart.
TARGET_LIKE = ["phat_DA", "phat_Sert"]
# Still-competent progenitors, among the types that remain after the removal.
PROGENITOR = ["phat_FPP", "phat_P_FPP"]
# Committed off-target fates, the rest of what remains.
OFF_TARGET = ["phat_Epen1", "phat_U_Neur1", "phat_U_Neur2"]

# Kept as a column so the renormalisation is auditable after the fact rather
# than being an invisible transformation: it is the exact mass left in play.
UNDECIDED_MASS_COL = "undecided_mass"


def build_progenitor_balance_table(full: pd.DataFrame) -> pd.DataFrame:
    """Drop the target-like proportions and renormalise the rest to sum to 1.

    Rows where no non-target cells exist cannot be renormalised (0/0). They
    are not silently zero-filled -- that would invent a uniform composition --
    but raise, because a D30 line made up purely of DA/Sert is a real edge
    case the caller must decide about rather than discover as a NaN in a
    results table."""
    missing = [c for c in TARGET_LIKE if c not in full.columns]
    if missing:
        raise ValueError(f"{missing} absent from the D30 table; nothing to remove")

    keep = sorted(c for c in full.columns if c.startswith("phat_") and c not in TARGET_LIKE)
    unexpected = set(keep) - set(PROGENITOR) - set(OFF_TARGET)
    if unexpected:
        raise ValueError(
            f"unclassified D30 cell types {sorted(unexpected)}: every remaining type must be "
            "declared progenitor or off-target, or the variant's interpretation is undefined"
        )

    undecided_mass = full[keep].sum(axis=1)
    degenerate = undecided_mass <= 0
    if degenerate.any():
        raise ValueError(
            f"{int(degenerate.sum())} rows have no non-target-like cells at all, so the "
            "remaining composition is undefined (0/0). Decide explicitly how to treat "
            "them rather than letting them become NaN in a results table."
        )

    out = full.drop(columns=TARGET_LIKE).copy()
    out[keep] = full[keep].div(undecided_mass, axis=0)
    out[UNDECIDED_MASS_COL] = undecided_mass
    return out


def main() -> None:
    src = OUT_DIR / "fold_data/fold_features_D30_qualonly.csv"
    dst = OUT_DIR / "fold_data/fold_features_D30_progbal_qualonly.csv"
    full = pd.read_csv(src)
    out = build_progenitor_balance_table(full)

    kept = sorted(c for c in out.columns if c.startswith("phat_"))
    rowsums = out[kept].sum(axis=1)
    assert rowsums.sub(1.0).abs().max() < 1e-9, f"renormalisation failed: max dev {rowsums.sub(1.0).abs().max()}"
    assert len(out) == len(full), "row count must not change"

    out.to_csv(dst, index=False)
    print(f"read  {src.name}: {len(full)} rows, proportions "
          f"{sorted(c for c in full.columns if c.startswith('phat_'))}")
    print(f"wrote {dst.name}: {len(out)} rows, proportions {kept}")
    print(f"renormalised row sums: min={rowsums.min():.12f} max={rowsums.max():.12f}")
    print(f"undecided mass (1 - DA - Sert): mean={out[UNDECIDED_MASS_COL].mean():.3f} "
          f"min={out[UNDECIDED_MASS_COL].min():.3f} max={out[UNDECIDED_MASS_COL].max():.3f}")


if __name__ == "__main__":
    main()
