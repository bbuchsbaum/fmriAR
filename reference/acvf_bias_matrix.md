# Bias matrix relating raw residual autocovariance to the truth

Builds the linear map `A` satisfying `E[gamma_raw] = A gamma_true`
induced by projecting a design out of the data. Autocovariance estimated
from GLM residuals is biased low, increasingly so for richer designs,
because `E[ehat ehat'] = M Sigma M` for the residual-forming projection
`M`. Solving `A gamma = gamma_raw` removes that bias.

## Usage

``` r
acvf_bias_matrix(design, runs = NULL, censor = NULL, max_lag = 20L)
```

## Arguments

- design:

  Numeric design matrix (timepoints x regressors), the one whose
  projection produced the residuals.

- runs:

  Optional run labels, length `nrow(design)`. Lag products never cross a
  run boundary.

- censor:

  Optional indices of censored timepoints (1-based), excluded exactly as
  in
  [`fit_noise()`](https://bbuchsbaum.github.io/fmriAR/reference/fit_noise.md).

- max_lag:

  Highest lag to correct.

## Value

A named list of `(max_lag + 1)` square matrices, one per run, rows
indexed by the lag of the raw estimate and columns by the lag of the
truth.

## Details

Ordinarily you do not need this: pass `design` to
[`fit_noise()`](https://bbuchsbaum.github.io/fmriAR/reference/fit_noise.md)
and it builds and applies the matrices itself. Use this directly to
inspect the bias, or to cache the matrices when many datasets share one
design and segmentation and feed them back as `acvf_correction`.

One matrix is returned per run, because
[`fit_noise()`](https://bbuchsbaum.github.io/fmriAR/reference/fit_noise.md)
estimates each run separately and the residual operator restricted to a
run still involves the whole design.

The correction is exact for noise whose autocovariance dies within
`max_lag`. Under long memory it is partial, since truncating the system
aliases in structure past the budget. With high-pass filtering (DCT, 128
s) a budget near 20-25 lags is typically enough; without high-pass it
exceeds 150 and the approach stops being practical.

Cost is `O(max_lag^2 * n_valid * n)` per run, so it is worth caching.

## See also

[`fit_noise()`](https://bbuchsbaum.github.io/fmriAR/reference/fit_noise.md),
[`noise_acvf()`](https://bbuchsbaum.github.io/fmriAR/reference/noise_acvf.md)

## Examples

``` r
n <- 120
X <- cbind(1, poly(seq_len(n), 3))
A <- acvf_bias_matrix(X, max_lag = 5)
round(A[[1]][1:3, 1:3], 3)
#>        [,1]   [,2]   [,3]
#> [1,]  0.967 -0.064 -0.062
#> [2,] -0.032  0.935 -0.063
#> [3,] -0.032 -0.063  0.937
```
