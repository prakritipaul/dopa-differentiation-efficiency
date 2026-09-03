# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project status

This is a brand-new Python project scaffold (created via `uv init`) with no functional implementation yet. The only code is a placeholder `main()` in `src/pluricon_prototype/__init__.py` that prints a greeting. There is no test suite, linter config, or CI configured. As the project grows, update this file with real architecture notes and commands.

## Commands

The project uses [uv](https://docs.astral.sh/uv/) for dependency management and packaging (build backend: `uv_build`).

- Run the CLI entry point: `uv run pluricon-prototype`
- Run any script/module in the project's environment: `uv run <command>`
- Add a dependency: `uv add <package>`
- Sync the environment: `uv sync`

Python version is pinned via `.python-version` to 3.11.

## Structure

- `src/pluricon_prototype/` — the package. `__init__.py` currently defines `main()`, wired up as the `pluricon-prototype` console script in `pyproject.toml` (`[project.scripts]`).

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
