# AGENTS.md

Guidance for Codex working in this repository.

**This file is deliberately not a copy of `CLAUDE.md`.** That file is Claude's
brief — how to plan, how lazy to be, when to ask for a second opinion. Reading
it would tell you what Claude was told to optimise for, and you would end up
reviewing against Claude's instructions instead of against the contract. That is
the same failure as writing tests after reading the implementation, one level up.
Your contract is the statistical invariants below. `CLAUDE.md` points here for
them rather than restating them, so there is one copy.

## Project status

Predicting D52 dopaminergic differentiation efficiency from D11 and D30
single-cell features across 136 iPSC lines, 20 donors (Jerber et al. 2021).
Both timepoints are implemented and have produced results.

- Outcome: `DA / all D52 cells`, untreated cells only, binarised at 0.2 →
  61 success / 75 failure. This is **not** the published metric; the reason is
  in `FINDINGS.md` §1 and it changes what every number means.
- D11 → D52: ROC-AUC 0.906, R² 0.503, under donor-grouped CV.
- D30 → D52: ROC-AUC 0.945, R² 0.802 — higher, but largely because D30's
  annotation already contains the outcome's numerator.

Where things are written down:

- `README.md` — project overview and results in summary.
- `FINDINGS.md` — results with literature context.
- `modeling/METHODS.md` — features, models, CV design, metrics, feature
  importance.
- `modeling/README.md` — audit trail, decision records, label-variant machinery.
- `modeling/docs/` — investigations (PCA interpretation, D30 cell types, pool
  correction) and `audit_ledger.md`, the status table for audit gates.

## Commands

Dependency management is [uv](https://docs.astral.sh/uv/); Python is pinned to
3.11 via `.python-version`.

- Tests: `uv run pytest modeling/ -m "not slow"` — 120 pass, 3 deselected as
  slow, 1 strict `xfail` encoding a known defect (unseen test pools map to NaN
  during pool correction, `modeling/tests/test_harness.py:164`). Strict means it
  fails the suite if it ever starts passing, so do not "fix" it silently.
- Full suite, 124 tests including the 3 slow ones (~3 min):
  `uv run pytest modeling/`
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
- `metadata_eda/` — EDA outputs from the numbered `eda/0NN_*.py` scripts.
- `eda/0NN_*.py` — ten EDA scripts, run in order, streaming multi-GB h5 files.
  Not covered by tests; their output paths are checked statically in
  `modeling/tests/test_imports.py`.

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
