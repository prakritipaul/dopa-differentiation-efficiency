"""
What ARE the seven D30 cell types, and which of them can still become DA/Sert?

This exists because the D30 -> D52 model is only interpretable once you know
which D30 types are upstream of the outcome and which are terminal. D11 needed
no such script: all three of its annotated types (FPP, P_FPP, NB) are
progenitor/early, so the whole D11 composition is upstream by construction. At
D30 the annotation set already contains DA and Sert -- the outcome itself --
alongside genuine progenitors and committed off-target fates, and those three
groups have to be told apart from evidence rather than from their names.

Method: mean log-normalised expression (the h5's `X`, which 005 established is
already log-normalised; `raw/X` is counts) of canonical marker panels, averaged
per annotated cell type over all D30 cells. No fitting, no subsetting to the
138 study lines -- this is an annotation-level question, not a modelling one.

Findings are written up in modeling/docs/D30_celltype_interpretation.md; this
script is what regenerates the numbers quoted there.

Usage:
    uv run python -m modeling.d30_celltype_markers
Writes modeling/results/d30_celltype_markers.csv
"""

import importlib.util
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

from modeling.features import TIMEPOINT_FILES

REPO_ROOT = Path(__file__).parent.parent
OUT_DIR = Path(__file__).parent

# Panels chosen a priori from the midbrain-differentiation literature, not
# selected after looking at which ones separated the clusters.
PANELS = {
    "progenitor / floor plate": ["LMX1A", "FOXA2", "CORIN", "SOX2", "NES", "VIM", "HES1", "SOX9"],
    "proliferation": ["MKI67", "TOP2A", "CCNB1"],
    "pan-neuronal": ["STMN2", "DCX", "MAP2", "SYT1", "ELAVL4"],
    "DA neuron": ["TH", "NR4A2", "PITX3", "DDC", "SLC6A3", "KCNJ6"],
    "Sert neuron": ["FEV", "TPH1", "SLC6A4", "GATA2", "GATA3"],
    "ependymal / choroid plexus / glial": ["FOXJ1", "S100B", "PIFO", "RSPH1", "SPAG6", "HDC", "TTR"],
}


def _load_005():
    """005 already streams selected gene columns out of the CSR matrix in
    chunks; reused rather than reimplemented (it is the same operation the
    HVG step performs)."""
    spec = importlib.util.spec_from_file_location("m005", REPO_ROOT / "eda/005_d11_pca_features.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def marker_means(timepoint: str = "D30") -> pd.DataFrame:
    """Mean expression per (cell type, marker). Rows: cell type. Columns:
    MultiIndex (panel, gene)."""
    m005 = _load_005()
    path = TIMEPOINT_FILES[timepoint]

    with h5py.File(path, "r") as f:
        genes = np.array([g.decode() if isinstance(g, bytes) else g for g in f["var/index"][:]])
        celltype = pd.Categorical.from_codes(
            f["obs/celltype"][:], [c.decode() for c in f["obs/__categories/celltype"][:]]
        )

    located, absent = [], []
    for panel, gene_list in PANELS.items():
        for gene in gene_list:
            hit = np.where(genes == gene)[0]
            if len(hit):
                located.append((panel, gene, hit[0]))
            else:
                absent.append(gene)
    if absent:
        # Reported, never silently skipped: a missing marker changes which
        # conclusions the panel can support.
        print(f"WARNING: markers absent from var/index and omitted: {absent}")

    # extract_hvg_matrix requires ascending column indices.
    located.sort(key=lambda t: t[2])
    idx = np.array([t[2] for t in located])
    X = m005.extract_hvg_matrix(path, idx)

    frame = pd.DataFrame(X, columns=[t[1] for t in located])
    frame["celltype"] = np.asarray(celltype)
    means = frame.groupby("celltype", observed=True).mean()
    means.columns = pd.MultiIndex.from_tuples(
        [(panel, gene) for panel, gene, _ in located], names=["panel", "gene"]
    )
    return means


def main() -> None:
    means = marker_means("D30")
    out_path = OUT_DIR / "results/d30_celltype_markers.csv"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    means.to_csv(out_path)

    for panel in PANELS:
        cols = [c for c in means.columns if c[0] == panel]
        if cols:
            print(f"\n=== {panel} ===")
            print(means[cols].droplevel(0, axis=1).round(2).to_string())
    print(f"\nSaved {out_path}")


if __name__ == "__main__":
    main()
