# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project status

Predicting D52 dopaminergic differentiation efficiency from D11 single-cell
features across 138 iPSC lines (Jerber et al. 2021). See `README.md` for the
overview, `FINDINGS.md` for results, `modeling/README.md` for methods.

Pipeline lives in `modeling/`; numbered `0NN_*.py` scripts are the EDA phase
that produced `metadata_eda/`. Test suite: `uv run pytest modeling/ -m "not slow"`.

## Commands

The project uses [uv](https://docs.astral.sh/uv/) for dependency management and packaging (build backend: `uv_build`).

- Run the CLI entry point: `uv run pluricon-prototype`
- Run any script/module in the project's environment: `uv run <command>`
- Add a dependency: `uv add <package>`
- Sync the environment: `uv sync`

Python version is pinned via `.python-version` to 3.11.

## Structure

- `modeling/` — the pipeline (see `modeling/README.md` for the subdirectory layout).
- `metadata_eda/` — EDA outputs, organised into `cohort/ pca/ proportions/ qc/ technical/ plots/`.
- `0NN_*.py` — numbered EDA scripts, run in order. Not covered by tests (they stream multi-GB h5 files); their output paths are checked statically in `modeling/tests/test_imports.py`.

## Planning workflow: get Codex's second opinion

Whenever we're finalizing a non-trivial technical or methodology decision
(e.g. a modeling approach, a CV/train-test scheme, an interpretation of
results, a statistical method) — not for small/obvious implementation
choices — do this before proceeding with implementation:

1. State your own opinion/recommendation first, with reasoning.
2. Get Codex's independent opinion on the same question via the
   `codex:rescue` skill (see `.claude/plugins/openai-codex` if present).
   Pose the question to Codex *without* revealing your own answer first,
   so its take is genuinely independent, not just agreement colored by
   what you said.
3. Present a comparison table (your take vs. Codex's take) plus bullet
   points calling out where you agree, where you disagree, and why —
   before implementing anything based on the discussion.

This applies by default without the user needing to ask each time.

## Track multi-step work with a visible checklist

For any task with more than ~3 steps, or any task spanning background jobs,
open with a checklist and restate it as items complete. Plain markdown in the
response -- there is no todo tool in this project.

```
- [x] regenerate results
- [ ] verify numbers against the README   <- current
- [ ] commit
```

Rules that make it worth doing:
- Write it **before** starting, not retroactively. Its job is to catch the
  step you would otherwise drop, which only works if it exists first.
- One item per verifiable outcome, not per action. "Verify tables reproduce
  byte-identically" beats "check the tables".
- Restate the full list when it changes, so the current state is always in
  the latest message and survives context compaction.
- Add discovered work to the list rather than doing it silently. If a fix
  spawns two more fixes, they become items.
- Do not close an item without the evidence it claims. If it says verified,
  show what verified it.

Why: this project has repeatedly produced work that ran to completion while
silently dropping a step -- permutation importance computed then discarded,
a "fully fixed" defect that was fixed in three of four places, a regenerated
table nobody diffed. Each was a missing checklist item, not a coding error.
