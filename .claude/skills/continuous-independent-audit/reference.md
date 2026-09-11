# Deep verification criteria

Load when auditing critical logic, or before a completion gate that claims
whole-project scope. Not needed for routine per-change audits.

## Account for every test individually

For each regression and integration test in the agreed scope, record which
criteria apply and the evidence for each. Explain any criterion that does not
apply. If required evidence cannot be obtained, the test stays unverified --
and so does its module.

**Correct expectations.** Assertions trace to an explicit requirement or an
independently justified reference result, not to observed implementation
behavior. A test that asserts what the code currently prints proves only that
the code is deterministic.

**Meaningful execution.** The test actually runs, reaches the behavior it
names, and makes substantive assertions. Skipped, blocked, deselected or
never-run tests are unverified. Check the runner's own selection flags: a
test excluded by a default marker filter can fail silently for months.

**Regression sensitivity.** A regression test must reproduce its target defect
on a known-bad version or a controlled faulty copy, and pass on the corrected
version. A regression test never seen to fail is not known to test anything.

**Integration validity.** The test exercises the component boundary it claims.
Mocks must not stand in for that boundary. Disclose every test double and any
external behavior left untested.

**Independent review.** Codex reviews pre-existing tests; Claude reviews
Codex-authored tests against requirements. Neither agent's approval alone
establishes correctness.

## Controlled defects

To check whether a test is sensitive to the bug it targets, introduce the
defect deliberately and confirm the test fails.

Do this on a disposable copy. Never modify the working implementation for
these experiments. Apply it to critical logic wherever practical, not only to
tests with a named prior defect.

## Oracle tests

The strongest available check for numerical and statistical code: construct
input whose correct output is known by independent means -- a closed-form
answer, a value computed by a separate path, a published reference figure, a
simulation with a planted ground truth -- and assert the pipeline recovers it.

An oracle test can refute a plausible-but-wrong explanation in a way that
comparing against recorded outputs never can, because recorded outputs
inherit whatever the implementation already does.

## Snapshot hygiene

Record the exact code snapshot under audit, including uncommitted files, plus
the requirements and dependency versions in force. A finding without a
snapshot id cannot be re-verified later.

Run against temporary data and test services. Never production credentials.
Import only reviewed test changes; preserve concurrent work.

Codex may edit tests and fixtures only. After the blind suite is preserved,
expose the implementation snapshot for execution and read-only review.

## Faulty tests

Codex may correct a test it wrote, but only with a recorded, contract-based
reason. "The test failed" is not a reason. Silent weakening of assertions is
the failure mode this whole protocol exists to prevent.

Accepted risks stay explicitly unverified in the ledger. They are not ✅.
