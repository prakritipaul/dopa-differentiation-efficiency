# Pool-correction: why we dropped it (for later revisiting)

**Decision**: pool-correction removed from the active pipeline (`run_experiment.py`
`CORRECTIONS = [False]`). Code still supports it (`pool_correction=True` in
`features.py`/`harness.py`) if revisited later.

**The puzzle**: correction reduced regression performance a lot (R² ~0.65 →
~0.55-0.59) but barely touched classification (AUC ~0.94 either way).

## My investigation
- Proportions (`phat_FPP`/`phat_P_FPP`) lose 70-80% of their variance to
  correction — expected, they were the most pool-confounded features (`008`).
- But per-feature correlation with the outcome mostly got *stronger* after
  correction (in-sample check) — contradicts "correction destroys signal."
- The MAE penalty from correction appears already at k=1-2 PCs and plateaus —
  doesn't keep growing with the later, more heavily-corrected PCs (6/8/9).
- Landed on: probably estimation noise from small per-pool training samples
  (some pools had only 4-7 training lines in a fold).

## Codex's independent investigation
- **Ran an oracle test that refutes the estimation-noise theory**: pool means
  computed from all rows (train+test, deliberately leaky) performed
  essentially identically to training-only means (MAE 0.157 vs 0.156). If
  noisy means were the cause, perfect means should have recovered the loss —
  they didn't.
- **Alternative explanation**: the outcome is naturally bimodal — 42 failures
  (mean 0.070, max 0.185) vs. 96 successes (mean 0.619, min 0.221), an actual
  gap between 0.185-0.221. Classification only has to separate two
  well-clustered groups; regression also has to resolve *within-class*
  variation, which pool-aligned PC1-5 apparently helps with even though pool
  barely predicts the binary label (η²≈0.12, matching `010`).
- **Corrected a factual error**: the data has 10 D11 pools, not 12-17 as I'd
  assumed.
- **Two real implementation issues found** (unfixed): (1) inner-CV
  leakage — pool-correction/PCA are computed once per outer fold, not
  refit per inner fold, so inner hyperparameter selection is mildly
  transductive; (2) training rows contribute to their own pool's mean
  (residuals shrink toward zero) while test rows don't — an asymmetry
  between how train vs. test features get treated.
- Suggested follow-ups if revisited: permutation negative-control (shuffle
  pool labels, see how much correction alone costs), decompose each
  feature into between/within-pool components and test each separately,
  sweep correction strength λ instead of all-or-nothing, cross-fit
  training-row pool means to fix the train/test asymmetry.

## Follow-up: checked the within-clump R² (confirms Codex's hypothesis)

The outcome is bimodal: 42 failures (0-0.185) and 96 successes
(0.221-~0.92), a real gap at the threshold. Split the headline nested-CV
predictions (donor_grouped, no correction) by true clump:

| | overall R² | within-success R² | within-failure R² |
|---|---|---|---|
| lasso | 0.667 | 0.222 | **-14.3** |
| ridge | 0.652 | 0.174 | **-14.5** |

The headline R²≈0.667 is mostly "correctly tells success from failure,"
not fine-grained precision. Within the success clump there's modest real
signal (R²≈0.17-0.22). Within the failure clump the model is far worse
than a trivial baseline (large negative R²) — absolute error (~0.135) is
almost as big as that clump's entire true range (0-0.185), so it isn't
resolving "barely failed" from "badly failed" at all. Practical
takeaway: trust this model for success/failure classification; don't
trust the regression output as a precise efficiency estimate, especially
for lines it predicts will fail.
