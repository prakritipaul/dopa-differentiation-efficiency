# Audit ledger

Status table only. Narrative findings live in `modeling/README.md` and
`modeling/docs/D30_celltype_interpretation.md`; rows here link to them. Two
homes for findings is how "overstated what had been fixed" happens again.

Protocol: `continuous-independent-audit`. Claude writes the code; Codex writes
tests from the contract **before** seeing the implementation, then reviews it.
Neither agent's agreement is evidence — reproducers are.

Legend: ✅ verified · ⏳ incomplete · ❌ confirmed defect · 🔄 stale · 🚧 blocked

Cadence is per milestone, not per change (a CV run is 40+ minutes). The three
gates for the D30 work:

1. feature extraction correct
2. harness / CV correct
3. results table correct

`-m "not slow"` is the everyday command, so **no module reaches ✅ on that
evidence alone** — the full suite must run before a completion gate.

## D30 → D52 milestone

**Gate 1 (feature extraction correct): passed**, 2025-09-11. Codex wrote
`tests/test_d30_contract_blind.py` from the contract without reading the
implementation, then reviewed the code. Full suite (slow tier included):
**77 passed, 1 xfailed** in 5m13s.

| Module / file | Behavior | Tests | Review | Status |
|---|---|---|---|---|
| `features.py` | Per-timepoint depth-outlier exclusion resolves from the timepoint, not a global constant (D30 → pool5, not pool11) | blind: passed | reviewed, no finding | ✅ |
| `features.py` | No held-out line shapes the HVG/PCA basis; all cells still projected | blind: passed — oracle perturbs held-out expression, training coordinates unchanged | reviewed, no finding | ✅ |
| `harness.py` | `proportion_cols` reproduces the pre-refactor D11 constants exactly | oracle vs the literal old constants; passed | reviewed, no finding | ✅ |
| `harness.py` | 7-type D30 composition → 6 free coords, reference dropped, full rank | blind + own: passed | reviewed, no finding | ✅ |
| `harness.py` | Missing reference / no proportions raises rather than fitting rank-deficient | passed | reviewed, no finding | ✅ |
| `folds.py` | Deterministic, timepoint-independent; LOCO/LODO/donor-grouping invariants | blind: passed | reviewed, no finding | ✅ |
| `run_feature_extraction.py` | Fold assignments identical across timepoints; original left intact on mismatch | runtime check fires (log line); blind: passed | reviewed, no finding | ✅ |
| `run_feature_extraction.py` | `fit_population` recorded per row; resume across differing **or unrecorded** populations refused | **❌ found, fixed, re-verified** — see below | Codex finding #1 | ✅ |
| `run_feature_extraction.py` | `append_fold_rows` never changes the column set (only order) | **❌ found, fixed, re-verified** — 3 new tests | Codex finding #1 (second half) | ✅ |
| `make_d30_variant_tables.py` | Renormalisation removes the DA+Sert magnitude, not just the columns | blind: passed — oracle with equal balance, different totals | reviewed, no finding | ✅ |
| `d30_single_feature_benchmarks.py` | Donor-level (not line-level) bootstrap; a-priori metric directions; single-class resamples → NaN, reported | blind: passed | reviewed, no finding | ✅ |
| `d30_celltype_markers.py` | Marker means reproduce the values quoted in the interpretation doc | reruns identically; no unit test | not audited | ⏳ |
| `harness.py` + `run_experiment.py` | D30 CV runs end to end (gate 2) | not started — extraction running | not started | ⏳ |
| results tables | D30 results correct (gate 3) | not started | not started | ⏳ |

Gates 2 and 3 are **not** verified. Only gate 1 is.

### Finding #1 — missing fit-population provenance accepted silently

Two halves, both wrong-code (not wrong-test, not ambiguous spec):

1. `if seen and seen != {fit_population}` — an **empty** `seen` (a table with no
   `fit_population` column) bypassed validation entirely, so a legacy table of
   unknown provenance was treated as compatible and its folds counted as done.
   Unknown provenance is not matching provenance: such a table could have been
   fit on either population. Now refused.
2. Worse, and not spotted when writing the guard: the append path reindexed each
   chunk to the header already on disk. Selecting those columns **silently
   dropped** `fit_population` from the new rows, so a resumed run could produce a
   table carrying no provenance anywhere, with nothing raised. Fixed at root —
   `append_fold_rows` now rejects **any** column-set difference in either
   direction, not just this one column — and given its own reproducers in
   `tests/test_run_feature_extraction.py`.

One test in Codex's informed suite then failed on the error message's wording
after the fix. Behaviour was correct (it raised before any PCA work), so the
message was made to name `fit_population` explicitly rather than the assertion
being relaxed.

## Open items

- **pool5 annotation bias (D30).** Raised by the independent methodology
  review, not yet addressed. pool5's ~10× sequencing-depth deficit may have
  biased D30 *cell-type annotation* itself; excluding it from the HVG/PCA fit
  does nothing to correct biased labels. It holds 18 of 159 qualifying
  (line, pool) combos. See `D30_celltype_interpretation.md` §5.
- **Shared annotation machinery.** Whether D30 and D52 cell types were
  annotated independently or via shared clustering/markers is unverified. Not
  leakage either way, but correlated labelling error would inflate the apparent
  D30→D52 continuity. Would need the upstream paper's methods to settle.
- **Incremental value over `DA+Sert`.** At ROC-AUC 0.990 for a single raw
  column there is no AUC headroom left, so whether the PCs add anything must be
  judged on Brier / log-loss / continuous R², not AUC.
