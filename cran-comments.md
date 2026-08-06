# CRAN Submission Notes (0.3.3)

- This is a bug-fix update to the CRAN version `0.3.2`.

## Changes since 0.3.2

This release corrects defects that produced silently wrong results rather than
errors, so it is submitted as a prompt follow-up to `0.3.2`.

- Fixed AR order selection under `pooling = "parcel"`. The BIC innovation
  sequence was formed with a feedback filter rather than the intended FIR
  filter, so order 0 always won and parcel pooling returned all-zero
  coefficients: whitening was a silent no-op while the plan reported a
  non-zero order.
- Fixed AR estimation under `censor`. Centring each scrubbing fragment on its
  own mean removed the autocorrelation being estimated, attenuating the
  coefficients toward zero as censoring increased and yielding an autocovariance
  that was not positive semi-definite, hence non-stationary coefficients that
  amplified variance during whitening. The mean is now estimated per run while
  lag products stay within contiguous valid segments.
- Stationarity is now enforced on the global and run estimation paths.
- `pooling = "parcel"` now honours `censor` and no longer estimates across run
  boundaries.
- Fixed `whiten_apply()` for parcel plans, which passed the caller's design
  matrix to an in-place routine, overwriting it and returning one repeatedly
  filtered matrix shared by all parcels.
- `acorr_diagnostics()` now uses its documented `runs` argument.
- Added a test file asserting recovery of known AR coefficients for every
  pooling mode, across censoring levels and run structures.

## R CMD check results

- Local `R CMD check --as-cran --no-manual` on macOS Sonoma 14.3, R 4.5.1
  (`aarch64-apple-darwin20`): `0 errors | 1 warning | 1 note`
- All examples, tests, and vignette rebuilds pass.
- The warning and note are unchanged from the `0.3.2` submission and are
  described below; neither originates in package code.

## Notes on warning/note

- The warning is an external compiler message from
  `R_ext/Boolean.h` in the R 4.5.1 headers on this local macOS/Homebrew clang
  setup (`unknown warning group '-Wfixed-enum-extension'`); it is not emitted
  by package source code.
- The note is `unable to verify current time`, which appears to be a
  machine-specific environment issue during the local check.

## Additional notes

- Acronyms such as AR and ARMA are defined in `DESCRIPTION`.
- The package uses optional OpenMP support where available.
