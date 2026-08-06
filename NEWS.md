# fmriAR 0.3.3

This release fixes several defects that silently produced wrong results rather
than errors. Analyses run with `pooling = "parcel"` or with `censor` under 0.3.2
should be rerun.

* Fixed order selection for `pooling = "parcel"`. The innovation sequence used to
  score BIC was built with a feedback filter rather than the intended FIR filter,
  which inflated the variance at every order and made order 0 always win. The
  effect was that parcel pooling returned all-zero AR coefficients: `whiten_apply()`
  returned its inputs unchanged while the plan reported a non-zero order. Order
  selection now scores the Levinson-Durbin prediction error, matching the
  global and run paths.

* Fixed AR estimation under censoring. Each scrubbing fragment was centred on its
  own mean, which removes the autocorrelation being measured -- a two-frame
  fragment yields a lag-1 correlation of exactly -1 regardless of the data. The
  estimate was attenuated toward zero as censoring increased (a true phi of 0.6
  was recovered as approximately 0 at 40% censoring) and the pooled autocovariance
  could lose positive semi-definiteness, so `fit_noise()` returned non-stationary
  coefficients and `whiten_apply()` then amplified variance instead of whitening.
  The mean is now estimated per run, while lag products remain confined to
  contiguous valid segments.

* `enforce_stationary_ar()` is now applied on the global and run Yule-Walker
  paths, which previously returned raw coefficients with no stationarity check.

* The autocovariance is now positive semi-definite by construction: the unbiased
  pair-count normalization is retained where it is valid and shrunk toward the
  common-divisor form only as far as needed.

* `pooling = "parcel"` now honours `censor`, which was previously discarded
  internally, and no longer estimates across run boundaries. Between-run mean
  offsets previously registered as near-perfect autocorrelation.

* Fixed `whiten_apply()` for parcel plans, which passed the caller's design
  matrix to an in-place routine. The caller's `X` was overwritten, every `X_by`
  entry aliased a single matrix, and that matrix had been filtered once per
  parcel in sequence rather than once with each parcel's own coefficients.

* Parcel plans now report the AR order actually fitted instead of the padded
  length of the coefficient vector.

* `acorr_diagnostics()` now uses its `runs` argument, which was documented and
  accepted but never referenced, so results were not run-aware.

* Multiscale pooling sizes the autocovariance to the pooling target rather than
  the selected order, removing a zero-filling step that drove Yule-Walker to
  coefficients pinned at the stationarity boundary.

# fmriAR 0.3.2

* `fit_noise()` gains a `censor` parameter for motion scrubbing support. Censored
  timepoints (e.g., frames with high framewise displacement) are excluded from

  AR parameter estimation. Accepts either integer indices or a logical vector.
  The time series is segmented at censor points and ACVF is pooled across valid
  segments with proper length-weighting.

* `whiten()` now passes `censor` to both `fit_noise()` and `whiten_apply()`.

* The returned `fmriAR_plan` object now includes the `censor` indices for
  downstream reference.

* Fixed a bug where `fit_noise(method = "ar", p = <integer>)` could ignore the
  requested fixed AR order and still perform order selection, sometimes
  returning a different order.

* Fixed parcel-mode AR fitting so fixed-order requests only trigger multiscale
  pooling when that mode is explicitly requested, preserving the expected
  non-multiscale behavior by default.

* Fixed the parcel `pacf_weighted` multiscale path for `p_target = 1`, where a
  dimension drop could break coefficient averaging.

# fmriAR 0.3.0

* Vignette corrections and clarity improvements.
* Diagnostics: call `acorr_diagnostics()` on innovations (whitened residuals).
* ARMA section: added note about Hannan-Rissanen estimation on run-mean series.
* Parcel pooling: clarified `X_by[[pid]]` dimensions.

# fmriAR 0.2.0

* Initial CRAN release.
* Core functionality: `fit_noise()`, `whiten_apply()`, `whiten()`.
* AR and ARMA(p,q) noise model estimation.
* Run-aware and parcel-aware pooling options.
* Multi-scale parcel pooling with PACF-weighted and ACVF-pooled modes
* AFNI-compatible restricted AR estimation via `afni_restricted_plan()`.
* Autocorrelation diagnostics via `acorr_diagnostics()`.
* Robust standard errors via `sandwich_from_whitened_resid()`.
