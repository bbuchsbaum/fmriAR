# Changelog

## fmriAR 0.4.2

### Fixes

- AR order selection (`p = "auto"`) scored BIC on the frame count alone,
  treating many pooled voxels as one series’ worth of evidence. A
  voxel-pooled AR(2) with phi = (0.5, -0.2) could come back AR(1) (as in
  the package vignette), leaving residual autocorrelation behind. The
  BIC sample size is now frames times the effective number of
  independent voxels, as the ARMA search already did. AR(2) over 60
  voxels: residual autocorrelation after whitening 0.025 -\> 0.004;
  ARMA(1,1) noise fitted with AR: 0.053 -\> 0.002. Single-voxel fits,
  white noise, and shared-signal data are unchanged, and
  [`fit_noise()`](https://bbuchsbaum.github.io/fmriAR/reference/fit_noise.md)
  timing is unchanged.

### Documentation

- The vignette described parcel pooling and ARMA estimation as working
  on mean time series; both now pool over voxels and the text says so.
  New sections cover censoring (and dropping censored rows before the
  final fit), residual-bias correction, automatic ARMA order selection,
  AFNI’s `vrt`, and HAC standard errors. Its diagnostics use
  [`acorr_diagnostics()`](https://bbuchsbaum.github.io/fmriAR/reference/acorr_diagnostics.md)
  with run boundaries respected.
- README updated to match: OpenMP is on whenever the toolchain supports
  it, bias correction covers every pooling mode, and `fmriAR.use_cpp_hr`
  affects only the internal single-series helper.

### Tests

- 16 further regression tests close the remaining gaps from the audit:
  the C++ lag-sum and Hannan-Rissanen kernels against plain-R
  references, censor validation, HAC run boundaries, parcel censor
  reporting, the `compat$whiten_with_phi()` default, ARMA failure
  fallback and refinement iterations, AFNI spec validation, and an
  end-to-end check that GLS t-tests keep a nominal false-positive rate
  for AR and ARMA noise.

## fmriAR 0.4.1

Numbers are 0.4.0 -\> 0.4.1 from
`tools/validation/accuracy_benchmark.R`; every metric not listed is
unchanged.

### Fixes

- Residual-bias correction (`design =` / `acvf_correction =`) was
  unstable under censoring. With drift regressors the bias map has a
  near-null direction (roughly a constant offset across lags); without
  censoring the raw autocovariance carries almost no noise along it, but
  censoring breaks the constraint that guarantees this, and the exact
  solve amplified that noise: phi RMSE 0.45 at 10% censoring for a truth
  of 0.4. Near-null directions are now fixed by the short-memory
  assumption the truncated correction already makes (autocovariance ~0
  at the end of the lag budget); well-determined directions are solved
  exactly, so well-conditioned designs are unaffected. Near-null means
  small relative to the median singular value, so designs that leave few
  residual degrees of freedom (where the whole spectrum is small) still
  get the exact solve. RMSE at 10% censoring 0.446 -\> 0.072; without
  censoring 0.035 -\> 0.020.
- Parcel pooling estimated each parcel from its parcel-mean time series,
  whose autocorrelation is dominated by whatever the voxels share.
  Parcels are now fitted from their voxels’ pooled autocovariance, like
  global/run pooling (a single all-voxel parcel reproduces global
  pooling exactly). phi error against the voxels’ own lag-1
  autocorrelation 0.061 -\> 0.006, and with a shared slow component
  0.228 -\> 0.007 (residual autocorrelation after whitening 0.129 -\>
  0.068). Multiscale pooling improves too: RMSE 0.062 -\> 0.038.
- The C++ stationarity/invertibility repair in single-series
  Hannan-Rissanen mis-mapped coefficients when the highest-order
  coefficient was exactly zero (`arma::roots()` drops it): (1.5, 0)
  became (-1, 0.67).
- [`afni_restricted_plan()`](https://bbuchsbaum.github.io/fmriAR/reference/afni_restricted_plan.md)
  rejects negative `a`/`r1`/`r2`, as AFNI does, instead of silently
  clamping them to 0. `vrt` (AFNI’s signal-to-total variance ratio) was
  documented but ignored; with `estimate_ma1 = FALSE` and `vrt < 1` the
  plan is now the exact ARMA(p, p) equivalent of AR(p) plus white noise
  (its autocorrelation is `vrt` times the AR one at every lag, matching
  AFNI). With `estimate_ma1 = TRUE`, a `vrt < 1` warns that it is
  ignored.

### New

- Residual-bias correction is supported for `pooling = "parcel"` in
  [`fit_noise()`](https://bbuchsbaum.github.io/fmriAR/reference/fit_noise.md)
  and
  [`noise_acvf()`](https://bbuchsbaum.github.io/fmriAR/reference/noise_acvf.md)
  (RMSE 0.020 uncensored, 0.062 at 10% censoring, matching global
  pooling).

### Performance

- [`acvf_bias_matrix()`](https://bbuchsbaum.github.io/fmriAR/reference/acvf_bias_matrix.md)
  exploits the low rank of the residual projection: O(L^2 n r) instead
  of O(L^2 n^2); 3.3 s -\> 0.19 s at n = 800, identical to 1e-12.

## fmriAR 0.4.0

Validated with `tools/validation/accuracy_benchmark.R`, which runs the
same seeded scenarios against two installed versions; numbers below are
0.3.3 -\> this version.

### Breaking changes

- `whiten_apply(inplace =)` is deprecated and ignored, with a warning.
  The inputs were never modified (R’s copy-on-modify semantics prevent
  doing so safely), so the argument only made the result invisible.
- `fit_noise(method = "arma", p = "auto")` now chooses the AR order by
  BIC from `0:p_max` instead of always fitting AR order 2 (see below).
- Results change where the fixes below say so: exact start-up for
  AR(p)/ARMA, pooled ARMA estimation, per-voxel
  [`acorr_diagnostics()`](https://bbuchsbaum.github.io/fmriAR/reference/acorr_diagnostics.md),
  the BIC penalty, and `compat$whiten_with_phi(exact_first = TRUE)`.

### Fixes

- An explicit `p` larger than `p_max` is now honoured;
  `fit_noise(p = 8)` returned an AR(6) because the order was capped at
  the default `p_max = 6`.
- Automatic order selection used `2 n log(sigma2) + k log(n)`, half the
  standard BIC penalty. For a single-voxel AR(1) it chose p \> 1 17% of
  the time, now 3% (power to detect a true AR(2) with phi = (0.5, 0.25)
  at n = 200 goes from 97% to 89%, the expected price of the correct
  penalty).
- [`whiten_apply()`](https://bbuchsbaum.github.io/fmriAR/reference/whiten_apply.md)
  accepts a logical `censor` mask. It was passed through
  [`as.integer()`](https://rdrr.io/r/base/integer.html), censoring only
  timepoint 1, which also broke
  [`whiten()`](https://bbuchsbaum.github.io/fmriAR/reference/whiten.md).
- [`whiten_apply()`](https://bbuchsbaum.github.io/fmriAR/reference/whiten_apply.md)
  uses the plan’s `censor` set when none is given and the data have the
  length the plan was fitted on (plans now record `n_time`), mirroring
  the existing `runs` fallback.
- [`whiten_apply()`](https://bbuchsbaum.github.io/fmriAR/reference/whiten_apply.md)
  errors when a per-run plan is applied to data with a different number
  of runs instead of silently pairing runs with the wrong coefficients.
  When censoring is active the result carries `censor`, the rows to drop
  before fitting.
- `enforce_invertible_ma()` moved roots inside the unit circle by
  reflection, which leaves a unit root where it is (`theta = -1` stayed
  `-1`). Roots are now kept at modulus \>= 1/0.99, mirroring the AR PACF
  bound.
- [`sandwich_from_whitened_resid()`](https://bbuchsbaum.github.io/fmriAR/reference/sandwich_from_whitened_resid.md)
  no longer fails on rank-deficient designs; non-estimable coefficients
  get `NA` standard errors, as in
  [`lm()`](https://rdrr.io/r/stats/lm.html).
- Multiscale parcel weights depended on the units of the data (the
  dispersion term was a raw variance). They are now unit-free.
- Multiscale `acvf_pooled` with `p = "auto"` fitted every parcel at
  `p_max`; it now uses the largest order selected at any scale (mean
  fitted order in the benchmark 6 -\> 1.4 for a true AR(1), with
  unchanged accuracy).
- `acorr_diagnostics(aggregate = "mean"/"median")` aggregates per-voxel
  autocorrelations. It used the autocorrelation of the voxel-mean
  series, which measures shared signal and reported 0.86 lag-1
  autocorrelation for voxels whose own was 0.25.
- `compat$update_plan()` keeps the previous plan’s `censor` and
  `exact_first`; `compat$whiten_with_phi()` now defaults to
  `exact_first = TRUE`, matching `plan_from_phi()`.
- The parallel-determinism test passed its input matrices to the
  in-place kernel, so all outputs aliased one buffer and the test
  compared it to itself.

### Algorithms

- Exact stationary start-up for any AR(p)/ARMA(p,q) (Ansley’s banded
  Cholesky). `exact_first = "ar1"` previously rescaled only AR(1)
  starts; for AR(2) (1.2, -0.5) the first two whitened samples of every
  run or censor segment had variance 3.73 and 1.93, now 1.01 and 1.01.
  Whitening now equals dense-Cholesky GLS within each segment to ~1e-14;
  AR(1) output is unchanged. GLS slope variance in short runs: AR(2)
  4x60 runs 50.3 -\> 48.5, ARMA(1,1) 6x40 runs 43.1 -\> 39.5 (x1e-3);
  type-I rates stay nominal.
- ARMA is estimated by a Hannan-Rissanen regression pooled over voxels
  and confined to contiguous segments, instead of on the voxel-mean
  series with censoring gaps spliced together. ARMA(1,1) RMSE
  (phi/theta) 0.079/0.088 -\> 0.014/0.014; with 20% censoring
  0.112/0.165 -\> 0.035/0.040; mean residual autocorrelation after
  whitening drops by 22-56% across scenarios. The censoring warning is
  gone, and the plan’s `sigma2` is now reported for ARMA at voxel scale.
  `afni_restricted_plan(estimate_ma1 = TRUE)` uses the same pooled
  estimator.
- Global AR pooling fits Yule-Walker on the frame-weighted average of
  the runs’ autocorrelations rather than averaging per-run coefficients
  (identical for AR(1); slightly better for AR(p)).
- ARMA orders can be selected automatically: `q = "auto"` (new `q_max`,
  default 2) and `p = "auto"` search the order grid by BIC on the pooled
  Hannan-Rissanen regression. All candidates are column subsets of one
  regression, so selection costs one extra pass over the data. The BIC
  sample size is frames times the effective number of independent voxels
  (Kish design effect from the mean inter-voxel correlation): counting
  frames alone under-fitted shared slow noise, counting every voxel
  over-fitted when voxels share fluctuations. Order recovery: 100% for
  white noise and AR(1), 80% for ARMA(1,1) (30 voxels x 240 frames).
- ARMA(1,1) noise plus a shared slow AR(0.95) component is not an
  ARMA(1,1) process, and the single shared realisation keeps even the
  empirical voxel autocorrelation ~0.14 from the theoretical one, so no
  estimator can recover “the” (1,1) parameters there. What can be fixed
  is how white the voxels end up. Residual autocorrelation after
  whitening: 0.135 in 0.3.3, 0.061 with a fixed (1,1) fit now, 0.055
  with automatic orders, against 0.052 for a correctly specified model.
- Single-series Hannan-Rissanen skips the long-AR burn-in.
- `sandwich_from_whitened_resid(type = "hac")`: Newey-West standard
  errors within runs, for when the noise model leaves autocorrelation
  behind.

### Performance

- OpenMP is now actually compiled (`SHLIB_OPENMP_CXXFLAGS`); thread
  count honours `options(fmriAR.max_threads)`. R API calls were moved
  out of the parallel region, where an interrupt would have terminated
  R.
- Segment-aware lag sums and the Hannan-Rissanen normal equations run in
  C++. At 400 x 20000: `fit_noise(p = "auto")` 1.07 s -\> 0.29 s, HC0
  sandwich 1.07 s -\> 0.30 s (vectorised, identical results). Pooled
  ARMA estimation does strictly more work than the old single-series fit
  (0.47 s -\> 0.93 s).
- Parcel plans whiten the design once per distinct filter.

## fmriAR 0.3.3

### Fixes

- `fit_noise(design =)` and `noise_acvf(design =)` now reject residuals
  that are not numerically orthogonal to the supplied design. The
  correction is derived for OLS residuals from that design’s column
  space; applying it to raw series or residuals from a different
  operator can silently undo a projection that never occurred and
  substantially inflate estimated autocorrelation. Orthogonality is a
  necessary check, not proof of provenance, so callers still must supply
  the residual-forming design that was actually used.

### New

- [`fit_noise()`](https://bbuchsbaum.github.io/fmriAR/reference/fit_noise.md)
  can now correct the residual bias in its autocovariance, via a new
  `design` argument. Autocovariance taken from GLM residuals is biased
  because `E[ehat ehat'] = M Sigma M` for the residual-forming
  projection `M`, which pushes `phi` down by an amount that grows with
  the number of regressors: at T = 300 with true AR(1) rho = 0.4, a
  9-column design returns about 0.36 and a 28-column design about 0.25.
  Passing `design = X` builds the linear map `A` with
  `E[gamma_raw] = A gamma_true` and solves it, which removed essentially
  all of that bias in simulation. The correction is opt-in: it changes
  estimates, and it needs the design that actually produced the
  residuals. Currently supported for `pooling = "global"` and `"run"`
  with `method = "ar"`; other combinations raise an error rather than
  quietly skipping the correction.
  [`acvf_bias_matrix()`](https://bbuchsbaum.github.io/fmriAR/reference/acvf_bias_matrix.md)
  is exported so the matrices can be built once and reused across
  datasets sharing a design, and fed back through `acvf_correction`.
  `correction_max_lag` (default 25) sets the lag budget; the correction
  is exact for noise whose autocovariance dies within it and partial
  under long memory, so with no high-pass filtering the required budget
  exceeds 150 and the approach is not practical. Two guards keep it from
  returning nonsense at the edges. A budget approaching the run length
  makes the system ill-conditioned – at n = 300 the reciprocal condition
  number falls from 0.28 at lag 25 to 5e-9 at lag 250, and solving there
  returned `phi` = 0.96 for a truth of 0.5 – so a correction that
  ill-conditioned is skipped with a warning and that run is left
  uncorrected. And a design leaving fewer residual degrees of freedom
  than the lag budget cannot support it, so the budget is reduced to
  what the data carries, again with a warning. Both failures returned
  finite, positive, plausible-looking numbers, so neither is caught by
  checking for `NaN`.

- New
  [`noise_acvf()`](https://bbuchsbaum.github.io/fmriAR/reference/noise_acvf.md)
  exports the run- and censor-aware autocovariance estimator the package
  already used internally. Lag products never cross a run boundary or a
  censoring gap, and the mean is removed per run rather than per
  fragment. It returns covariances on the scale of the data, along with
  the pair count behind each lag and the segmentation actually used, so
  a consumer can tell an estimate backed by thousands of pairs from one
  backed by three. Previously the only exported route to autocorrelation
  was
  [`acorr_diagnostics()`](https://bbuchsbaum.github.io/fmriAR/reference/acorr_diagnostics.md),
  which returns normalized ACF for inspection rather than covariance.
  For a run-stationary process the covariance of any contrast factors as
  `sum_h gamma_h * L T_h L'`, so `gamma` plus a design gives exact
  contrast variances without refitting.

- `fmriAR_plan` objects now carry the noise scale as well as its shape:
  `gamma` (autocovariance, lags 0..p) and `sigma2` (innovation variance)
  per pooling unit, plus `gamma_by_parcel` and `sigma2_by_parcel` for
  `pooling = "parcel"`. Previously a plan recorded only the correlation
  structure, so two datasets differing 100-fold in variance produced
  identical plans and the magnitude of the noise was unrecoverable
  without refitting. `gamma` is reported at voxel scale for every
  pooling mode, and `sigma2` is derived as
  `gamma_0 - sum_k phi_k gamma_k` from the coefficients stored on the
  plan, so the two are always mutually consistent. Under global pooling,
  runs that reached different lags are truncated to the shortest before
  averaging rather than zero-padded: a zero-padded autocovariance is not
  a valid covariance, and building `Sigma` from one could yield negative
  contrast variances. Truncation can leave `gamma` shorter than
  `length(phi)` when censoring is heavy, and `sigma2` is then `NA`
  rather than a partial sum, which would overstate the innovation
  variance. The condition does not arise below roughly 25% censoring.
  The addition is purely additive; existing fields are unchanged. For
  `method = "arma"`, `gamma` is the voxel-scale noise autocovariance
  pooled the same way as for AR, and `sigma2` is `NA`: Hannan-Rissanen’s
  own innovation variance is that of the run-mean series, smaller than
  the per-voxel value by roughly the number of voxels averaged, so
  reporting it would understate the noise by that factor.

### Fixes

This release fixes several defects that silently produced wrong results
rather than errors. Analyses run with `pooling = "parcel"` or with
`censor` under 0.3.2 should be rerun.

- Global pooling now weights each run by the number of uncensored
  observations that actually contributed. Previously
  [`fit_noise()`](https://bbuchsbaum.github.io/fmriAR/reference/fit_noise.md)
  weighted by the original run length, so a nearly empty run could carry
  the same influence as a complete run and disagree with
  [`noise_acvf()`](https://bbuchsbaum.github.io/fmriAR/reference/noise_acvf.md).

- Run labels are validated and encoded consistently across fitting, ACVF
  estimation, bias correction, whitening,
  [`acorr_diagnostics()`](https://bbuchsbaum.github.io/fmriAR/reference/acorr_diagnostics.md),
  and
  [`afni_restricted_plan()`](https://bbuchsbaum.github.io/fmriAR/reference/afni_restricted_plan.md).
  Wrong-length vectors, missing labels, and labels reused in
  non-contiguous blocks now fail at the boundary instead of recycling or
  silently dropping timepoints; contiguous character labels are
  supported. Previously
  [`acorr_diagnostics()`](https://bbuchsbaum.github.io/fmriAR/reference/acorr_diagnostics.md)
  coerced character labels to `NA` and silently ignored the run split,
  and
  [`afni_restricted_plan()`](https://bbuchsbaum.github.io/fmriAR/reference/afni_restricted_plan.md)
  ordered its run starts by sorted label rather than by time, so runs
  labelled out of time order produced descending starts.

- Under global pooling, a run censored down to one or zero surviving
  frames no longer erases the pooled autocovariance. Such a run carries
  no `gamma` of its own; it previously set the common truncation length
  to zero, so the plan reported `gamma` of length 0 and `sigma2 = NA`
  for the whole fit even though `phi` was pooled correctly and
  [`noise_acvf()`](https://bbuchsbaum.github.io/fmriAR/reference/noise_acvf.md)
  returned the full answer.

- [`fit_noise()`](https://bbuchsbaum.github.io/fmriAR/reference/fit_noise.md)
  rejects residuals containing `NA`, `NaN`, or `Inf`, matching
  [`noise_acvf()`](https://bbuchsbaum.github.io/fmriAR/reference/noise_acvf.md).
  It previously returned an order-0 plan with empty coefficients, which
  [`whiten_apply()`](https://bbuchsbaum.github.io/fmriAR/reference/whiten_apply.md)
  accepted and applied as a no-op.

- Lag budgets are validated the same way in
  [`fit_noise()`](https://bbuchsbaum.github.io/fmriAR/reference/fit_noise.md)
  and
  [`noise_acvf()`](https://bbuchsbaum.github.io/fmriAR/reference/noise_acvf.md):
  both reject a non-finite or non-positive `correction_max_lag`
  ([`fit_noise()`](https://bbuchsbaum.github.io/fmriAR/reference/fit_noise.md)
  used to clamp silently to 1), and budgets past the series length are
  clamped to it rather than overflowing integer coercion into an
  unrelated error.

- [`noise_acvf()`](https://bbuchsbaum.github.io/fmriAR/reference/noise_acvf.md)
  under `pooling = "global"` keys its single pooled unit as `"1"` even
  when only one run survives censoring; it previously leaked the
  surviving run’s label into the name.

- [`noise_acvf()`](https://bbuchsbaum.github.io/fmriAR/reference/noise_acvf.md)
  now has a separate `correction_max_lag` argument (default 25),
  matching
  [`fit_noise()`](https://bbuchsbaum.github.io/fmriAR/reference/fit_noise.md),
  so requesting a short output ACVF no longer weakens residual-bias
  correction. Its `corrected` field now reports whether correction was
  actually applied rather than whether a design was merely supplied.

- Parcel labels containing `NA` or fractional numeric values now fail at
  the boundary instead of dropping voxels or colliding after integer
  coercion.

- Parcel labels that do not survive
  [`as.integer()`](https://rdrr.io/r/base/integer.html) are now refused
  at the boundary with a message naming the offending values. Character
  labels became `NA`, matched no voxel, and surfaced far downstream as
  `invalid K`, which names nothing the caller passed. The check covers
  every exported entry point that accepts parcels:
  [`fit_noise()`](https://bbuchsbaum.github.io/fmriAR/reference/fit_noise.md)
  (including `parcel_sets`),
  [`whiten_apply()`](https://bbuchsbaum.github.io/fmriAR/reference/whiten_apply.md),
  [`afni_restricted_plan()`](https://bbuchsbaum.github.io/fmriAR/reference/afni_restricted_plan.md),
  and `compat$plan_from_phi()`. Integer, numeric, and factor labels are
  unaffected. Character labels are refused rather than mapped to codes
  the way `runs` are, because
  [`fit_noise()`](https://bbuchsbaum.github.io/fmriAR/reference/fit_noise.md)
  and
  [`whiten_apply()`](https://bbuchsbaum.github.io/fmriAR/reference/whiten_apply.md)
  would each have to derive the same mapping from their own copy of the
  vector and nothing on the plan records it; convert once with
  `as.integer(factor(parcels))` and reuse that coding.

- Fixed order selection for `pooling = "parcel"`. The innovation
  sequence used to score BIC was built with a feedback filter rather
  than the intended FIR filter, which inflated the variance at every
  order and made order 0 always win. The effect was that parcel pooling
  returned all-zero AR coefficients:
  [`whiten_apply()`](https://bbuchsbaum.github.io/fmriAR/reference/whiten_apply.md)
  returned its inputs unchanged while the plan reported a non-zero
  order. Order selection now scores the Levinson-Durbin prediction
  error, matching the global and run paths.

- Fixed AR estimation under censoring. Each scrubbing fragment was
  centred on its own mean, which removes the autocorrelation being
  measured – a two-frame fragment yields a lag-1 correlation of exactly
  -1 regardless of the data. The estimate was attenuated toward zero as
  censoring increased (a true phi of 0.6 was recovered as approximately
  0 at 40% censoring) and the pooled autocovariance could lose positive
  semi-definiteness, so
  [`fit_noise()`](https://bbuchsbaum.github.io/fmriAR/reference/fit_noise.md)
  returned non-stationary coefficients and
  [`whiten_apply()`](https://bbuchsbaum.github.io/fmriAR/reference/whiten_apply.md)
  then amplified variance instead of whitening. The mean is now
  estimated per run, while lag products remain confined to contiguous
  valid segments.

- `enforce_stationary_ar()` is now applied on the global and run
  Yule-Walker paths, which previously returned raw coefficients with no
  stationarity check.

- The autocovariance is now positive definite by construction. The
  unbiased pair-count normalization is retained wherever it is valid,
  and otherwise the non-zero lags are shrunk toward white noise only as
  far as needed. Positive definiteness is required with a relative
  margin rather than mere non-negativity: on the boundary a reflection
  coefficient is exactly 1, which collapses the Levinson prediction
  error to its floor, and BIC reads that as a perfect fit and selects
  the maximum order. The resulting filters amplified variance under
  censoring instead of whitening.

- `pooling = "parcel"` now honours `censor`, which was previously
  discarded internally, and no longer estimates across run boundaries.
  Between-run mean offsets previously registered as near-perfect
  autocorrelation.

- Fixed
  [`whiten_apply()`](https://bbuchsbaum.github.io/fmriAR/reference/whiten_apply.md)
  for parcel plans, which passed the caller’s design matrix to an
  in-place routine. The caller’s `X` was overwritten, every `X_by` entry
  aliased a single matrix, and that matrix had been filtered once per
  parcel in sequence rather than once with each parcel’s own
  coefficients.

- Parcel plans now report the AR order actually fitted instead of the
  padded length of the coefficient vector.

- [`acorr_diagnostics()`](https://bbuchsbaum.github.io/fmriAR/reference/acorr_diagnostics.md)
  now uses its `runs` argument, which was documented and accepted but
  never referenced, so results were not run-aware.

- Multiscale pooling sizes the autocovariance to the pooling target
  rather than the selected order, removing a zero-filling step that
  drove Yule-Walker to coefficients pinned at the stationarity boundary.
  The multiscale autocovariance is also estimated with a per-run mean
  rather than a per-segment one, so the `p_target` and `acvf_pooled`
  paths no longer reintroduce the censoring defect described above.

- [`whiten_apply()`](https://bbuchsbaum.github.io/fmriAR/reference/whiten_apply.md)
  returns rows in input order for any run labelling. Results were
  reassembled in sorted run-label order, so runs labelled in any order
  other than ascending-in-time silently produced a row permutation of
  the correct answer for both `X` and `Y`.

- Combining `parcel_sets` with `censor` no longer fails with “incorrect
  length for ‘group’”.

- `enforce_stationary_ar()` now guarantees characteristic roots strictly
  outside the unit circle. Clamping reflection coefficients alone left
  the roots on the circle at high order.

- Order *selection* is bounded by the available sample size, so BIC can
  no longer choose AR(8) from eleven observations. An explicitly
  requested `p` is still honoured as given.

- `p_max` at or near the series length no longer fails with “missing
  value where TRUE/FALSE needed”.

- `method = "arma"` now warns when combined with `censor`. Censored
  frames are excluded, but Hannan-Rissanen then runs on the surviving
  frames spliced together, so its regressions span the gaps and bias
  both the AR and MA coefficients. Prefer `method = "ar"` when censoring
  is present.

- [`whiten_apply()`](https://bbuchsbaum.github.io/fmriAR/reference/whiten_apply.md)
  validates `runs`: a length mismatch was silently recycled by
  [`split()`](https://rdrr.io/r/base/split.html) and an `NA` left that
  row unwritten in both `X` and `Y`. Both now raise an error, and
  character run labels are accepted.

- `compat$plan_from_phi()` works with its documented default
  `theta = NULL` for global and run pooling. It previously reported an
  MA order of `-Inf` and produced a plan that
  [`whiten_apply()`](https://bbuchsbaum.github.io/fmriAR/reference/whiten_apply.md)
  rejected with “subscript out of bounds” – the plan built by the
  example in
  [`?compat`](https://bbuchsbaum.github.io/fmriAR/reference/compat.md).

## fmriAR 0.3.2

CRAN release: 2026-04-15

- [`fit_noise()`](https://bbuchsbaum.github.io/fmriAR/reference/fit_noise.md)
  gains a `censor` parameter for motion scrubbing support. Censored
  timepoints (e.g., frames with high framewise displacement) are
  excluded from

  AR parameter estimation. Accepts either integer indices or a logical
  vector. The time series is segmented at censor points and ACVF is
  pooled across valid segments with proper length-weighting.

- [`whiten()`](https://bbuchsbaum.github.io/fmriAR/reference/whiten.md)
  now passes `censor` to both
  [`fit_noise()`](https://bbuchsbaum.github.io/fmriAR/reference/fit_noise.md)
  and
  [`whiten_apply()`](https://bbuchsbaum.github.io/fmriAR/reference/whiten_apply.md).

- The returned `fmriAR_plan` object now includes the `censor` indices
  for downstream reference.

- Fixed a bug where `fit_noise(method = "ar", p = <integer>)` could
  ignore the requested fixed AR order and still perform order selection,
  sometimes returning a different order.

- Fixed parcel-mode AR fitting so fixed-order requests only trigger
  multiscale pooling when that mode is explicitly requested, preserving
  the expected non-multiscale behavior by default.

- Fixed the parcel `pacf_weighted` multiscale path for `p_target = 1`,
  where a dimension drop could break coefficient averaging.

## fmriAR 0.3.0

- Vignette corrections and clarity improvements.
- Diagnostics: call
  [`acorr_diagnostics()`](https://bbuchsbaum.github.io/fmriAR/reference/acorr_diagnostics.md)
  on innovations (whitened residuals).
- ARMA section: added note about Hannan-Rissanen estimation on run-mean
  series.
- Parcel pooling: clarified `X_by[[pid]]` dimensions.

## fmriAR 0.2.0

CRAN release: 2025-11-03

- Initial CRAN release.
- Core functionality:
  [`fit_noise()`](https://bbuchsbaum.github.io/fmriAR/reference/fit_noise.md),
  [`whiten_apply()`](https://bbuchsbaum.github.io/fmriAR/reference/whiten_apply.md),
  [`whiten()`](https://bbuchsbaum.github.io/fmriAR/reference/whiten.md).
- AR and ARMA(p,q) noise model estimation.
- Run-aware and parcel-aware pooling options.
- Multi-scale parcel pooling with PACF-weighted and ACVF-pooled modes
- AFNI-compatible restricted AR estimation via
  [`afni_restricted_plan()`](https://bbuchsbaum.github.io/fmriAR/reference/afni_restricted_plan.md).
- Autocorrelation diagnostics via
  [`acorr_diagnostics()`](https://bbuchsbaum.github.io/fmriAR/reference/acorr_diagnostics.md).
- Robust standard errors via
  [`sandwich_from_whitened_resid()`](https://bbuchsbaum.github.io/fmriAR/reference/sandwich_from_whitened_resid.md).
