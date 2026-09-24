# AGENTS.md

Guidance for Codex working in this repository.

**This file is deliberately not a copy of `CLAUDE.md`.** That file is Claude's
brief — how to plan, how lazy to be, when to ask for a second opinion. Two
reasons they stay separate:

1. **One copy of the contract.** `CLAUDE.md` points here for the statistical
   invariants rather than restating them, so the two cannot drift apart. This
   repo has already been bitten by one rule living in two places.
2. **One instruction there is self-defeating if you read it.** Claude is told to
   ask for your opinion *before* revealing its own, so that your answer is not
   agreement coloured by its reasoning. Knowing you are the designated second
   opinion invites the mirror-image bias — disagreeing because dissent is what
   feels useful. Give the answer you would give if you had never been told a
   comparison was coming.

Nothing else in `CLAUDE.md` is withheld from you, and reading it for project
context is fine. Just note that its coding-style rules are Claude's brief, not
your review criteria: a correctness review of this pipeline is about the
invariants below, not about whether a diff was minimal.

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
function returns the right shape is nearly worthless — a shape check passes
identically whether the PCA basis was fit on training lines only or on all 136.

> **What "statistical invariant" means here.** A property that must hold for a
> reported number to estimate the quantity its name claims. Violating one does
> not produce an error or an impossible value; it produces a number of the right
> type and a plausible magnitude that is no longer an estimate of what the label
> says. Break invariant 1 and ROC-AUC 0.906 stops being an out-of-sample
> forecast and becomes an optimistically biased in-sample one — same range, same
> type, different meaning. They are *invariants* because they must hold on every
> fold, every repeat and every variant, not be demonstrated once.
>
> Items 4 and 5 below are data-plumbing properties rather than statistical ones
> in themselves; they are on this list because their consequences are
> statistical — a duplicated join silently reweights lines, and a variant flag
> that stops halfway produces a table whose columns come from different bases.

The invariants that matter:

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
independently and check the pipeline recovers it.

### Where an oracle's "correct answer" comes from

**Independently means: known without executing the code under test.** You are
not looking the answer up anywhere. You derive it, or you force it by
construction. Five sources, none of which is the implementation:

| source | how it yields a known answer | example in this repo |
|---|---|---|
| **Deduction from the invariant** | the invariant's own logic implies an observable equality | invariant 1 says no fitted quantity may depend on a held-out line, so perturbing that line's expression must leave training coordinates *unchanged* |
| **Construction / symmetry** | build two inputs the spec forces to agree, or to differ | `test_progenitor_balance_removes_target_magnitude_oracle` — see below |
| **Counting** | the answer is an integer you can state in advance | LOCO has exactly one fold per line, each holds out exactly one, and the union of held-out lines equals the whole cohort — no cohort size needs to be known |
| **Linear algebra** | a mathematical fact about the design | 7 compositional types must give a full-rank 6-column model matrix, checked with `matrix_rank` |
| **An independent implementation** | a separate computation of the same thing | `proportion_cols` must reproduce the literal pre-refactor D11 constants |

The construction case is the one worth studying, because it shows you can test
a transform without knowing what its output *should be*. `_variant_input()`
builds two lines with the **same relative balance** among the five progenitor
types (1:2:3:4:5) but very different non-target mass (0.9 against 0.2, with
DA+Sert at 0.10 against 0.80). The contract says the progenitor-balance variant
removes the *magnitude* of maturation, not just the DA and Sert columns.
Therefore two lines differing only in that magnitude must come out identical —
and the test asserts exactly that, to 1e-12, without ever stating what the
correct proportions are. Drop the columns without renormalising and line `a`
keeps 0.9 of its mass while `b` keeps 0.2, so the rows differ and the test
fails.

**What is not an oracle: this pipeline's own saved output.** Comparing against
recorded results only proves the code is deterministic, and a leaking pipeline
is perfectly deterministic.

**Do not weaken a test to make it pass.** If a test fails, decide whether the
code is wrong, the test is wrong, or the requirement is ambiguous — and say
which.
