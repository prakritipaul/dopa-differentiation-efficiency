"""
End-to-end smoke test on a TINY synthetic dataset.

Everything else in this suite either uses small in-memory fixtures or reads
the committed result tables. Neither tells a newcomer whether the pipeline
actually *runs*, because the real inputs are three multi-GB h5 files that are
not in the repo. Someone who clones this cannot execute a single stage.

So this builds a miniature dataset in the same on-disk shape as the real one
-- legacy-AnnData h5 with a CSR `X`, `raw/X`, and categorical obs columns --
and drives the real code over it: cohort construction, label build, per-fold
PCA feature extraction, and a CV fit. It is deliberately the REAL functions,
not reimplementations, so it fails if any of them break.

Runs in a couple of seconds. Marked slow only because it writes to a tmp dir
and exercises several stages; it is still cheap enough for every run.
"""

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

h5py = pytest.importorskip("h5py")

N_GENES = 60
N_LINES = 12
N_POOLS = 2
CELLS_PER_COMBO = 40
CELLTYPES = {"D11": ["FPP", "NB", "P_FPP"], "D30": ["DA", "Sert", "FPP", "P_FPP", "Epen1"]}


def _write_h5(path, rng, timepoint, lines, pools, celltypes, treatments=("NONE",)):
    """Write a miniature file in the same layout the real data uses."""
    rows = []
    for line in lines:
        for pool in pools:
            for tx in treatments:
                for _ in range(CELLS_PER_COMBO):
                    rows.append((line, pool, rng.choice(celltypes), tx))
    obs = pd.DataFrame(rows, columns=["donor_id", "pool_id", "celltype", "treatment"])
    n = len(obs)

    # Sparse, non-negative, with a little structure so PCA has something to find.
    dense = rng.random((n, N_GENES)).astype(np.float32)
    dense[dense < 0.7] = 0.0
    X = sp.csr_matrix(dense)

    with h5py.File(path, "w") as f:
        for grp, mat in (("X", X), ("raw/X", X)):
            g = f.create_group(grp)
            g.create_dataset("data", data=mat.data)
            g.create_dataset("indices", data=mat.indices)
            g.create_dataset("indptr", data=mat.indptr)
        f.create_dataset("var/index", data=np.array([f"GENE{i}".encode() for i in range(N_GENES)]))
        f.create_dataset("obs/index", data=np.array([f"cell{i}".encode() for i in range(n)]))
        for col in ["donor_id", "pool_id", "celltype", "treatment"]:
            cats = sorted(obs[col].unique())
            codes = obs[col].map({c: i for i, c in enumerate(cats)}).to_numpy()
            f.create_dataset(f"obs/{col}", data=codes.astype(np.int16))
            f.create_dataset(f"obs/__categories/{col}", data=np.array([c.encode() for c in cats]))
        f.create_dataset("obs/time_point", data=np.zeros(n, dtype=np.int8))
        f.create_dataset("obs/__categories/time_point", data=np.array([timepoint.encode()]))
    return obs


@pytest.fixture(scope="module")
def tiny(tmp_path_factory):
    """A miniature three-timepoint dataset plus its cohort and label."""
    d = tmp_path_factory.mktemp("tiny")
    rng = np.random.default_rng(0)
    lines = [f"HPSI0000i-line_{i}" for i in range(N_LINES)]
    pools = [f"pool{i + 1}" for i in range(N_POOLS)]

    files = {}
    for tp in ["D11", "D30", "D52"]:
        ct = CELLTYPES.get(tp, ["DA", "Sert", "FPP", "Astro"])
        tx = ("NONE", "ROT") if tp == "D52" else ("NONE",)
        files[tp] = d / f"{tp}.h5"
        _write_h5(files[tp], rng, tp, lines, pools, ct, tx)

    cohort = pd.DataFrame([(l, "HPSI0000i", p) for l in lines for p in pools],
                          columns=["cell_line", "donor", "pool"])
    cohort_csv = d / "cohort.csv"
    cohort.to_csv(cohort_csv, index=False)

    eff = rng.random(N_LINES)
    label_csv = d / "label.csv"
    pd.DataFrame({"cell_line": lines, "diff_efficiency": eff}).to_csv(label_csv, index=False)
    return {"dir": d, "files": files, "cohort_csv": cohort_csv, "label_csv": label_csv,
            "lines": lines, "pools": pools}


@pytest.mark.slow
def test_pipeline_runs_end_to_end_on_a_tiny_dataset(tiny, monkeypatch):
    """Cohort -> label -> per-fold PCA features -> CV fit, on the real code."""
    from modeling import features, folds, harness
    from modeling.models import models_for_task

    # Register the tiny dataset as a label variant, so every stage resolves it
    # through the same registry the real run uses rather than a special path.
    monkeypatch.setitem(folds.VARIANTS, "tiny", folds.LabelVariant(
        suffix="_tiny", cohort_csv=tiny["cohort_csv"], label_csv=tiny["label_csv"],
        threshold=0.5, n_lines=N_LINES, description="synthetic smoke-test dataset",
    ))
    monkeypatch.setitem(features.TIMEPOINT_FILES, "D11", str(tiny["files"]["D11"]))
    monkeypatch.setitem(features.DEPTH_OUTLIER_POOLS, "D11", frozenset({"pool2"}))

    lines = folds.load_lines_with_label("tiny")
    assert len(lines) == N_LINES

    fold_list = folds.plain_repeated_kfold(lines, n_splits=3, n_repeats=1, seed=0)
    assert len(fold_list) == 3
    held_out = set(fold_list[0].test_lines)

    meta = features.load_cell_metadata("D11")
    qualifying = pd.read_csv(tiny["cohort_csv"])[["cell_line", "pool"]]

    # The real per-fold extraction: HVGs + PCA fit excluding held-out lines,
    # every cell projected.
    pcs = features.compute_pca_features_for_fold(
        "D11", held_out, meta=meta, n_pcs=3, restrict_to_combos=qualifying)
    assert len(pcs) == len(meta)
    assert not pcs[[f"PC{i}" for i in range(1, 4)]].isna().any().any()

    props = features.compute_proportion_features(
        meta.merge(qualifying, on=["cell_line", "pool"], how="inner"))
    phat = [c for c in props.columns if c.startswith("phat_")]
    assert np.allclose(props[phat].sum(axis=1), 1.0)

    # Assemble a fold table in the shape the harness consumes, then fit.
    pool_pcs = (pcs.merge(qualifying, on=["cell_line", "pool"])
                .groupby(["cell_line", "pool"], observed=True)[pcs.columns[2:].tolist()]
                .mean().reset_index())
    ff = pool_pcs.merge(props, on=["cell_line", "pool"])
    for i in range(4, 11):
        ff[f"PC{i}"] = 0.0
    ff["scheme"], ff["repeat"], ff["fold"] = "plain", 0, 0
    ff["split"] = np.where(ff.cell_line.isin(held_out), "test", "train")

    train, test = harness.build_line_level_for_fold(ff, "plain", 0, 0, False)
    assert set(train.cell_line) | set(test.cell_line) == set(tiny["lines"])
    assert not (set(train.cell_line) & set(test.cell_line))

    y = lines.set_index("cell_line")["diff_efficiency"]
    spec = next(m for m in models_for_task("regression") if m.name == "ridge")
    X_tr = harness.build_feature_matrix(train, 3)
    X_te = harness.build_feature_matrix(test, 3)
    preds, _ = harness.fit_predict(spec, {"alpha": 1.0}, X_tr, y.loc[train.cell_line].to_numpy(), X_te)
    assert len(preds) == len(test) and np.isfinite(preds).all()


@pytest.mark.slow
def test_held_out_lines_do_not_shape_the_basis_on_the_tiny_dataset(tiny, monkeypatch):
    """The central leakage invariant, on data small enough to perturb freely.

    Changing a held-out line's expression must not move a training line's PC
    coordinates. Cheap here, so it runs as an ordinary check rather than
    needing the multi-GB files."""
    from modeling import features

    monkeypatch.setitem(features.TIMEPOINT_FILES, "D11", str(tiny["files"]["D11"]))
    monkeypatch.setitem(features.DEPTH_OUTLIER_POOLS, "D11", frozenset())

    meta = features.load_cell_metadata("D11")
    held_out = {tiny["lines"][0], tiny["lines"][1]}
    train_line = tiny["lines"][5]

    base = features.compute_pca_features_for_fold("D11", held_out, meta=meta, n_pcs=3)

    # Rewrite the held-out lines' expression, then refit.
    perturbed = tiny["dir"] / "D11_perturbed.h5"
    import shutil
    shutil.copy(tiny["files"]["D11"], perturbed)
    with h5py.File(perturbed, "r+") as f:
        codes = f["obs/donor_id"][:]
        cats = [c.decode() for c in f["obs/__categories/donor_id"][:]]
        mask = np.isin(codes, [cats.index(l) for l in held_out])
        indptr = f["X/indptr"][:]
        data = f["X/data"][:]
        for i in np.where(mask)[0]:
            data[indptr[i]:indptr[i + 1]] *= 50.0
        f["X/data"][:] = data

    monkeypatch.setitem(features.TIMEPOINT_FILES, "D11", str(perturbed))
    after = features.compute_pca_features_for_fold("D11", held_out, meta=meta, n_pcs=3)

    pc_cols = [f"PC{i}" for i in range(1, 4)]
    keep = (base.cell_line == train_line).to_numpy()
    # PCA sign/rotation is arbitrary, so compare the subspace via absolute
    # correlation of each component across the two fits on TRAINING cells.
    a = base.loc[keep, pc_cols].to_numpy()
    b = after.loc[keep, pc_cols].to_numpy()
    for j in range(len(pc_cols)):
        r = abs(np.corrcoef(a[:, j], b[:, j])[0, 1])
        assert r > 0.99, f"held-out perturbation moved training {pc_cols[j]} (|r|={r:.3f})"
