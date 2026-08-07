# CRAN Submission Notes (fmriAR 0.3.3)

Update to the CRAN version `0.3.2`, containing new functionality and a set of
correctness fixes.

## Changes since 0.3.2

New functionality (all additions are backward compatible):

- `fit_noise()` gains an opt-in correction for the residual bias in the
  autocovariance (new `design` argument); `acvf_bias_matrix()` is exported so
  the correction matrices can be reused across datasets sharing a design.
- New `noise_acvf()` exports the run- and censor-aware autocovariance
  estimator the package already used internally.
- `fmriAR_plan` objects now record the noise scale (`gamma`, `sigma2`) in
  addition to its shape.

Bug fixes. This release corrects defects that produced silently wrong results
rather than errors, including:

- AR order selection under `pooling = "parcel"` always chose order 0, making
  whitening a silent no-op.
- AR estimation under `censor` centred each scrubbing fragment on its own
  mean, attenuating coefficients toward zero and allowing non-stationary
  fits that amplified variance instead of whitening.
- `whiten_apply()` reassembled rows in sorted run-label order, silently
  permuting the output for non-ascending run labellings.
- Parcel-plan whitening overwrote the caller's design matrix via an in-place
  routine, returning one repeatedly filtered matrix shared by all parcels.
- Run and parcel labels are now validated at the boundary instead of being
  recycled, dropped, or coerced to `NA` downstream.

The complete list is in `NEWS.md`.

## R CMD check results

- Local `R CMD check --as-cran --no-manual` on macOS (aarch64-apple-darwin20,
  R 4.5.1): `0 errors | 1 warning | 0 notes`.
- The warning is unchanged from the accepted `0.3.2` submission and does not
  originate in package code: it is an external compiler message from
  `R_ext/Boolean.h` in the R 4.5.1 headers on this local macOS/Homebrew clang
  toolchain (`unknown warning group '-Wfixed-enum-extension'`).

## Reverse dependencies

There are no reverse dependencies on CRAN.
