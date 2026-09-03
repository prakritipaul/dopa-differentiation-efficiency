# AGENTS.md

This file provides guidance to Codex (Codex.ai/code) when working with code in this repository.

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
