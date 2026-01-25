# fmriAR: Fast AR and ARMA Noise Whitening for Functional MRI (fMRI) Design and Data

Lightweight utilities to estimate autoregressive (AR) and autoregressive
moving average (ARMA) noise models from residuals and apply matched
generalized least squares to whiten functional magnetic resonance
imaging (fMRI) design and data matrices. The ARMA estimator follows a
classic 1982 approach
[doi:10.1093/biomet/69.1.81](https://doi.org/10.1093/biomet/69.1.81) ,
and a restricted AR family mirrors workflows described by Cox (2012)
[doi:10.1016/j.neuroimage.2011.08.056](https://doi.org/10.1016/j.neuroimage.2011.08.056)
.

Estimate AR/ARMA noise models from residuals and apply matched GLS
prewhitening to fMRI design and data matrices. Run-aware and
censor-aware.

## Details

The fmriAR package provides efficient implementations for:

- AR and ARMA model estimation from fMRI residuals

- Run-aware and censor-aware whitening transformations

- Parcel-based parameter pooling

- Sandwich standard error computation

## See also

Useful links:

- <https://bbuchsbaum.github.io/fmriAR/>

Useful links:

- [`fit_noise`](https://bbuchsbaum.github.io/fmriAR/reference/fit_noise.md)
  for noise model estimation

- [`whiten_apply`](https://bbuchsbaum.github.io/fmriAR/reference/whiten_apply.md)
  for applying whitening transformations

- [`whiten`](https://bbuchsbaum.github.io/fmriAR/reference/whiten.md)
  for one-step whitening

## Author

**Maintainer**: Bradley Buchsbaum <brad.buchsbaum@gmail.com>
