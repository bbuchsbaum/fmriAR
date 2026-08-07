# Changelog

## fmriAR 0.3.3

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
  estimation, bias correction, and whitening. Wrong-length vectors,
  missing labels, and labels reused in non-contiguous blocks now fail at
  the boundary instead of recycling or silently dropping timepoints;
  contiguous character labels are supported.

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
