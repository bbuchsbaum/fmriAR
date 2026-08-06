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
