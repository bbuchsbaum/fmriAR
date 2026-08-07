# Run- and censor-aware noise autocovariance

Estimates the autocovariance of fMRI noise directly, without fitting an
AR model first. Lag products never cross a run boundary or a censoring
gap, and the mean is removed per run rather than per fragment, so
scrubbed data does not have its autocorrelation destroyed by the
centering.

## Usage

``` r
noise_acvf(
  resid,
  runs = NULL,
  censor = NULL,
  max_lag = 20L,
  pooling = c("global", "run", "parcel"),
  parcels = NULL,
  design = NULL,
  correction_max_lag = 25L
)
```

## Arguments

- resid:

  Numeric matrix of residuals (timepoints x voxels).

- runs:

  Optional run labels, length `nrow(resid)`. Each label must occupy one
  contiguous block and may not be missing.

- censor:

  Optional indices of censored timepoints (1-based).

- max_lag:

  Highest lag to return.

- pooling:

  `"global"` for one autocovariance over all runs, `"run"` for one per
  run, `"parcel"` for one per parcel.

- parcels:

  Parcel labels (length `ncol(resid)`), required when
  `pooling = "parcel"`.

- design:

  Optional design matrix whose projection produced `resid`. When
  supplied the residual bias is corrected; see
  [`acvf_bias_matrix()`](https://bbuchsbaum.github.io/fmriAR/reference/acvf_bias_matrix.md).

- correction_max_lag:

  Lag budget used to estimate and undo residual projection bias when
  `design` is supplied. This is independent of `max_lag`: correction
  needs enough tail information even when only a few output lags are
  requested.

## Value

An object of class `fmriAR_acvf`: a list with

- `acvf`: named list of autocovariance vectors, lags 0 upward.

- `pairs`: matching list of the pair count behind each lag, so a
  consumer can tell an estimate backed by thousands of pairs from one
  backed by three.

- `n_segments`, `segment_lengths`: the segmentation actually used.

- `max_lag`, `pooling`, `corrected`.

Lags with no surviving pairs are dropped rather than reported as zero.

## Details

This is the estimator
[`fit_noise()`](https://bbuchsbaum.github.io/fmriAR/reference/fit_noise.md)
uses internally, exposed so consumers can get the covariance itself. For
a run-stationary process with autocovariance `gamma`, the covariance of
any linear contrast `L` factors as `sum_h gamma_h * L T_h L'` with `T_h`
the lag-`h` indicator, so `gamma` plus a design is enough to compute
exact contrast variances without refitting.

Unlike
[`acorr_diagnostics()`](https://bbuchsbaum.github.io/fmriAR/reference/acorr_diagnostics.md),
which returns normalized autocorrelation for inspection, this returns
covariances on the scale of the data.

## See also

[`fit_noise()`](https://bbuchsbaum.github.io/fmriAR/reference/fit_noise.md),
[`acvf_bias_matrix()`](https://bbuchsbaum.github.io/fmriAR/reference/acvf_bias_matrix.md)

## Examples

``` r
set.seed(1)
r <- matrix(rnorm(400 * 5), 400, 5)
a <- noise_acvf(r, runs = rep(1:2, each = 200), max_lag = 5, pooling = "run")
a$acvf[["1"]]
#> [1]  1.059682131 -0.023317336 -0.032081492 -0.036966387  0.007025246
#> [6]  0.103570754
```
