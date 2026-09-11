# AGENTS.md

Guidance for Codex working in this repository.

## Project status

Predicting D52 dopaminergic differentiation efficiency from D11 single-cell
features across 138 iPSC lines (Jerber et al. 2021). The pipeline is
implemented and has produced results: ROC-AUC 0.947 under never-seen-donor
cross-validation. Work on a D30 -> D52 model is starting.

- `README.md` — project overview and headline results.
- `FINDINGS.md` — results with literature context.
- `modeling/README.md` — methods, CV design, and a "D30 readiness" section
  listing what still hardcodes D11.
- `modeling/docs/` — investigations (pool correction, feature importance).

`src/pluricon_prototype/` is a leftover `uv init` scaffold with a placeholder
`main()`. Nothing in the analysis imports it.

## Commands

Dependency management is [uv](https://docs.astral.sh/uv/); Python is pinned to
3.11 via `.python-version`.

- Tests: `uv run pytest modeling/ -m "not slow"` (51 fast tests)
- Full suite including the ~3 min integration test: `uv run pytest modeling/`
- Run any module: `uv run python -m modeling.run_regression`
- Add a dependency: `uv add <package>`

`-m "not slow"` is the everyday command, so the slow tier is easy to forget.
It has previously hidden a test that was failing for months
(`modeling/tests/test_features.py:60-65`). Run the full suite before treating
results as final.

## Structure

- `modeling/` — the pipeline: `folds.py`, `features.py`, `models.py`,
  `harness.py`, `feature_importance.py`, and `run_*.py` entry points.
- `modeling/tests/` — pytest suite. All new tests go here.
- `modeling/fold_data/`, `modeling/results/` — per-fold tables and CV grids.
  `_qualonly` suffixes are a variant basis, not duplicates.
- `modeling/archive/` — superseded artifacts kept for provenance. Nothing
  reads these; do not delete them.
- `metadata_eda/` — EDA outputs from the numbered `0NN_*.py` scripts.
- `0NN_*.py` — EDA scripts, run in order, streaming multi-GB h5 files. Not
  covered by tests.

## Your role as independent auditor

You are often asked to audit this pipeline independently — see the
`continuous-independent-audit` skill. When that is the task:

**Write tests from the contract before reading the implementation.** Tests
written after reading the code inherit its assumptions and will confirm a
wrong pipeline.

**The contract here is statistical, not just I/O.** A test asserting a
function returns the right shape is nearly worthless. The invariants that
matter:

1. No held-out line contributes to any fitted quantity — the HVG/PCA basis,
   scaler statistics, or pool means.
2. Pool-level corrections are computed from training rows only.
3. LOCO leaves exactly one cell line out; LODO leaves exactly one donor out,
   including all of that donor's lines.
4. Feature tables and label tables join on cell line without silent row loss
   or duplication.
5. A variant flag (e.g. `--suffix`) must reach every layer it affects, not
   just the first one.

Invariant 5 is not hypothetical: it has already produced mixed-basis
importance tables here, where `--suffix` redirected the input and output paths
but not `load_full_fit_features()`.

**Failures here are usually silent.** Most defects in this repo do not raise —
they produce a plausible, non-crashing, wrong number that lands in a results
table. Prefer oracle tests: construct input whose correct answer is known
independently and check the pipeline recovers it. Comparing against recorded
outputs only proves the code is deterministic.

**Do not weaken a test to make it pass.** If a test fails, decide whether the
code is wrong, the test is wrong, or the requirement is ambiguous — and say
which.
