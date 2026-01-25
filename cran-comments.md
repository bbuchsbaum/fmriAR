# CRAN Submission Notes (0.3.1)

- This is a resubmission (previous CRAN version: 0.2.0).

## Changes since last CRAN release (0.2.0)

### New in 0.3.1
- `fit_noise()` gains a `censor` parameter for motion scrubbing support.
  Censored timepoints are excluded from AR parameter estimation, improving
  robustness when motion-corrupted frames are present.

### From 0.3.0
- Vignette corrections and clarity improvements; no API changes.

## R CMD check results
- 0 errors | 0 warnings | 0 notes

## Test environments
- Local: macOS (darwin), R 4.x
- GitHub Actions: ubuntu-latest, windows-latest, macOS-latest

## Additional notes
- Acronyms (AR, ARMA) are defined in DESCRIPTION.
- Examples avoid the `:::` operator; no use of non-exported functions in examples.
