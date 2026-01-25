# Changelog

## fmriAR 0.3.1

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
