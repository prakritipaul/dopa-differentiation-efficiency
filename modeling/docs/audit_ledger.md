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
| `harness.py` + `run_experiment.py` | Gate 2: train-only scaling/tuning; per-repeat metrics; nested model-only grouping; pooled LOCO/LODO; compositional rank; provenance; headline selection | blind gate-2 checks passed; full suite 105 passed, 4 failed, 1 xfailed | reviewed, no wrong-code finding | ✅ |
| results tables | Gate 3: four 448-row D11/D30 tables, complete grid, provenance, ranges, unique headline, finite means | blind: passed after the spec was corrected — 19 passed | contract conflict resolved, see below | ✅ |

**Both gates are verified.** Gate 3's four failures were an ambiguous
specification, not a defect, and the resolution tightened the test rather than
relaxing it.

The contract said "metric columns finite". But `loco` and `lodo` are exhaustive
enumerations with exactly one repeat, and the sample SD of a single observation
is undefined — NaN is the correct value. A finite sentinel such as `0.0` would
be worse than useless: it asserts zero across-repeat variability, a precision
claim nobody measured and one that could be quoted as if it had been.

So the spec was corrected to: every metric **mean** finite everywhere; SDs
finite for the repeated schemes; SDs **NaN for the single-repeat schemes**. That
last clause is an assertion, not a tolerance — NaN is now *required* where there
is one repeat. The replacement test is strictly stronger than the original,
which would have passed had a sentinel been written into those cells.

Independent review recommended permitting NaN; the bidirectional form goes
further than either side initially specified.

### Gates 2/3 audit — 2026-09-12

Blind suite: `tests/test_gates_2_3_blind.py`, written before reading
`harness.py`, `run_experiment.py`, or `run_variant.py`. Snapshot:
`e9c0f851bad75d3b78d59be18d7062b8bfad412d` plus the new blind test and this
ledger update. Final blind suite: **15 passed, 4 failed, 0 skipped**. Full
suite, including the slow tier before the result assertions were split for
independent execution: **105 passed, 4 failed, 1 xfailed** in 3m41s. The strict
xfail is the known unseen-pool correction defect and was not changed.

The four failures are one ambiguous specification repeated across the four
committed result tables. Every metric mean is finite; all missing metric cells
are SDs on LOCO/LODO rows, and those schemes have exactly one repeat. Pandas'
sample SD for one observation is correctly NaN, consistent with the explicit
gate-2 single-repeat contract and the pre-existing test for that behavior, but
inconsistent with gate 3's literal requirement that metric columns be finite.
No test was weakened and no result table was rewritten.

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
- **RESOLVED.** The DA-only D30 model vs one raw column was tested properly in
  `modeling/test_incremental_value.py` (both feature sets, identical
  donor-grouped folds, nested inner selection on training lines only, paired
  per repeat). The answer is split and both halves are unanimous across all 10
  repeats: `phat_DA` alone RANKS better (ROC-AUC 0.9563 vs 0.9393; full better
  in 0/10), while the full model is far better CALIBRATED (log loss 0.328 vs
  0.571, Brier 0.093 vs 0.190; full better in 10/10). Original note follows.
- **The DA-only D30 model does not beat one raw column on ranking.** Under the DA-only outcome, D30 `phat_DA` alone
  scores ROC-AUC **0.960** [0.922, 0.986] while the full fitted model —
  6 free proportions + PC1–10, donor-grouped nested CV over 55 configurations —
  scores **0.945**. The two are not strictly like-for-like (the benchmark is
  in-sample with no fitted parameter; the model is out-of-fold), so this is not
  proof the model is worthless, but it is not evidence it adds anything either.
  Settling it needs a paired comparison of out-of-fold predictions against the
  benchmark on Brier / log-loss, which has not been run.
- **A mis-specified benchmark nearly hid that.** The benchmark set was carried
  over from the DA+Sert outcome and compared `phat_DA + phat_Sert` (0.886)
  rather than `phat_DA` (0.960). Against the wrong comparator the model appeared
  to add ~6 AUC points. Fixed in `e9c0f85`; both rows are now reported, each
  labelled with the outcome it matches. **A benchmark that does not match the
  outcome's numerator is not a weaker check, it is a misleading one.**
