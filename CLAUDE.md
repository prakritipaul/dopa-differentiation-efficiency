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

Deliberately here and not a skill: a skill has to be invoked, and the moment
a checklist is most needed is the moment it is least likely to be reached
for. CLAUDE.md loads every session automatically. There is also no todo tool
in this project, so a skill would have nothing to call -- it would emit the
same markdown either way. Don't re-open this unless a todo tool appears.

## Ponytail overrides for this repo

The `ponytail` skill (lazy/minimal coding mode) earns its keep here — keep the
ladder, the root-cause rule ("grep every caller before you edit"), and
`ponytail:` markers for deliberate shortcuts. These five carve-outs override
the skill where its defaults fight this repo.

1. **Tests go in `modeling/tests/`, never `__main__` self-checks.** The skill's
   "assert-based `demo()`, no frameworks" minimum hand-rolls a worse pytest,
   which is already here — reusing it is rung 2 of ponytail's own ladder.
   Related: `-m "not slow"` hides the slow tier, so run the suite in full
   before declaring results final (`modeling/tests/test_features.py:60-65`
   documents a test that failed unnoticed for exactly this reason).
2. **Methods rationale is the deliverable, not prose debt.** "At most three
   short lines / delete the explanation" applies to defending a *code*
   simplification. Why a statistical choice was made belongs in `FINDINGS.md`
   and `modeling/README.md` at whatever length it takes.
3. **Never delete provenance to shorten a diff.** `archive/`, suffixed variant
   outputs (`_qualonly`), saved fold assignments and seeds are not speculative
   flexibility. Mutating an analysis script in place is the shortest diff and
   silently invalidates every result already published from it.
4. **Timepoint parameterization is not YAGNI.** D30 is the second
   implementation, so parameterizing beats forking `*_d11_*` into `*_d30_*`.
   "No config for a value that never changes" does not apply to a value that
   is about to change.
5. **Ponytail does not cover statistical correctness.** `ponytail-review` and
   `ponytail-audit` scope out correctness by design, and neither names the
   defects that matter here (fold leakage, PCA basis fit on held-out lines,
   label contamination across pools). Those stay with
   `continuous-independent-audit` / `codex:rescue` — an audit pass is not
   correctness coverage.

## Continuous independent audit: on by default

Invoke the `continuous-independent-audit` skill **automatically at the start of
any session that will change `modeling/`** — without being asked. Say in one
line that you have done so. This is distinct from the Codex second-opinion rule
above: that one reviews a *decision* before implementing, this one audits *code*
after writing it.

Why on by default: the defects this repo produces are silent. They do not raise
— they emit a plausible, non-crashing, wrong number into a results table. A
second independent Codex review already caught a case where a previous
`modeling/README.md` entry overstated what had been fixed
(`modeling/README.md:474`). An audit protocol you have to remember to invoke is
one you will skip on exactly the change that needed it.

The skill is tracked at `.claude/skills/continuous-independent-audit/`. That
copy is canonical; a global one may also exist in `~/.claude/skills/`. Edit the
repo copy, and if you touch the global one, sync it back here in the same
change -- two homes drifting apart is the failure this project keeps repeating.

Repo-specific settings, which override the skill's own defaults:

- **Ledger:** `modeling/docs/audit_ledger.md`, not `audit/ledger.md`. It holds
  the **status table only**. Narrative findings stay in `modeling/README.md`;
  ledger rows link to them. Two homes for findings is how "overstated what had
  been fixed" happens again.
- **Cadence:** per milestone, not per change — a CV run is 40+ minutes. The
  three gates are: feature extraction correct → harness/CV correct → results
  table correct.
- **Handoff:** the contract given to Codex must state the *statistical*
  invariants, not just signatures. They are listed in `AGENTS.md` under "Your
  role as independent auditor" — keep that list and this rule in sync.
- **Isolation:** `EnterWorktree`. **Codex:** the `codex:rescue` skill. If
  either is unavailable, the audit is 🚧 blocked — your own checks do not
  substitute.
- **Slow tier:** `-m "not slow"` is the everyday command, so no module reaches
  ✅ on its evidence alone. Run the full suite before a completion gate.
