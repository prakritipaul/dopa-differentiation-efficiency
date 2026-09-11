---
name: continuous-independent-audit
description: Coordinate independent Codex testing and review during Python development, verify fixes, and maintain a concise per-module audit table. Use when the user requests continuous auditing.
---

# Continuous Independent Audit

Claude writes the code. Codex independently writes tests **before seeing the
implementation**, then reviews it. Neither agent's agreement is evidence --
reproducers are. The user settles disagreements.

Requires a working Codex integration (`codex:rescue`) and an isolated
workspace (`EnterWorktree`). If either is unavailable, mark the audit
🚧 blocked; Claude's own checks cannot substitute for an independent one.

Audit each coherent, runnable change to meaningful behavior, batching small
related edits. Claude may continue unrelated work while Codex audits a fixed
snapshot, but must not declare the audited work complete before verification.

## The loop

1. **Hand off the contract, not the code.** Fresh Codex context: the user's
   request and acceptance criteria, signatures and schemas, and the
   **invariants the code must satisfy** -- derived from requirements, never
   read off the implementation. Withhold implementation bodies, your own
   reasoning, and any suspected fix. Label assumptions; ask the user when an
   ambiguity materially changes correctness, and continue checks that do not
   depend on the answer.
2. **Codex writes tests** covering normal behavior, boundaries, invalid input,
   failure modes, and module interactions. Test observable outcomes; never
   mock away the behavior under test. Prefer an *oracle test* where one
   exists: construct input whose correct answer is known independently, and
   check the pipeline recovers it.
3. **Preserve that blind suite**, then run it against a fixed snapshot.
   Record passed / failed / skipped. Unexecuted tests are pending, not passing.
4. **Only then reveal the implementation.** Codex reviews correctness, missing
   cases, integration behavior, and test quality. Keep implementation-informed
   tests distinct from the blind suite.
5. **Resolve findings.** Accept and fix, or dispute with evidence. Neither
   agent may weaken an assertion, skip a failure, or move an expectation to
   obtain a pass. Codex re-runs its reproducer plus relevant regression checks
   against the updated snapshot before anything is marked verified. After one
   evidence-based exchange, escalate unresolved disagreement to the user with
   both positions and a focused decision.

A failing test means one of three things: wrong code, wrong test, or an
ambiguous requirement. Decide which before fixing.

## Where expected behavior comes from

Step 1 only helps if "expected" has a source. In descending order of strength:

1. **A closed-form or independently computed answer** for input you construct
   -- the oracle route. Strongest, and nearly always available for numerical
   code.
2. **A reference table produced by a separate path**: an EDA output, a
   published figure, a value from the work being replicated. Cite file and row
   in the test.
3. **An explicit requirement** stated in the project's docs or the user's
   request.

Never the recorded output of the code under audit. That establishes
determinism, not correctness. If a number has no source, say so and ask --
do not invent a tolerance that the current output happens to satisfy.

## Ledger

Maintain `audit/ledger.md`, or the project's existing ledger. One row per
meaningful behavior, grouped by module; name both modules for an integration
row.

| Module / file | Behavior | Tests | Review | Status |
|---|---|---|---|---|
| `feature_importance.py` | Ranking order | 12 passed | Reviewed; 1 fixed | ✅ |
| `report.py` + `feature_importance.py` | Renders rankings | pending | pending | 🔄 |

Illustrative rows -- never invent results.

✅ verified · ⏳ incomplete · ❌ confirmed defect · 🔄 stale (code, tests,
requirements or deps moved) · 🚧 blocked

Verification covers only the stated behavior. When code changes, mark
dependent and integration rows 🔄 and recheck; broaden when impact is
unclear. Update at meaningful transitions; explain problems and decisions
rather than narrating routine passes.

## Before declaring completion

Reconcile every changed behavior against the ledger, refresh 🔄 rows, resolve
findings, and run the combined regression + integration suites. Retain
commands, exit results, and snapshot ids.

**A module cannot be ✅ while any required test lacks evidence.** Any
unverified row means verification is incomplete -- say that plainly instead of
reporting that everything passed. Passing tests, coverage numbers, and agent
agreement are not a correctness guarantee.

Deep criteria -- regression sensitivity via controlled defects, integration
boundary validity, per-test checklists, snapshot hygiene -- are in
`reference.md`. Load it when auditing critical logic or before a completion
gate on whole-project scope.
