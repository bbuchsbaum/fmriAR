# CRAN Submission Notes (0.3.2)

- This is an update to the CRAN version `0.2.0`.

## Changes since last CRAN release (0.2.0)

- Added motion-scrubbing support via a new `censor` argument in `fit_noise()`,
  with matching propagation through `whiten()` and `whiten_apply()`.
- `fmriAR_plan` objects now retain censor indices for downstream reference.
- Fixed a bug where fixed-order AR fits could still trigger order selection and
  return a different order than requested.
- Fixed parcel-mode and multiscale edge cases affecting fixed-order AR fitting
  and `pacf_weighted` pooling when `p_target = 1`.
- Retained the `0.3.0` vignette and diagnostics clarifications introduced after
  the initial CRAN release.

## R CMD check results

- Local `R CMD check --as-cran --no-manual` on macOS Sonoma 14.3, R 4.5.1
  (`aarch64-apple-darwin20`): `0 errors | 1 warning | 1 note`
- All examples, tests, and vignette rebuilds pass.

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
