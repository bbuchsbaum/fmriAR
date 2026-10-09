# GLS standard errors from whitened residuals

GLS standard errors from whitened residuals

## Usage

``` r
sandwich_from_whitened_resid(
  Xw,
  Yw,
  beta = NULL,
  type = c("iid", "hc0", "hac"),
  df_mode = c("rankX", "n-p"),
  runs = NULL,
  hac_lag = NULL
)
```

## Arguments

- Xw:

  Whitened design matrix.

- Yw:

  Whitened data matrix (time x voxels).

- beta:

  Optional coefficients (p x v); estimated if `NULL`.

- type:

  `"iid"` (default) assumes the whitened errors are white and
  homoskedastic. `"hc0"` is a heteroskedasticity-robust sandwich.
  `"hac"` is a Newey-West (Bartlett kernel) sandwich computed within
  runs, which stays valid when the noise model leaves some
  autocorrelation behind.

- df_mode:

  Degrees-of-freedom mode: "rankX" (default) or "n-p".

- runs:

  Optional run labels (length `nrow(Xw)`). Used by `type = "hac"` so
  that no lag product spans a run boundary.

- hac_lag:

  Bandwidth for `type = "hac"`. Defaults to
  `floor(4 * (n_run / 100)^(2/9))` using the median run length.

## Value

List containing standard errors, innovation variances, and XtX inverse.
Coefficients that are not estimable because `Xw` is rank deficient get
`NA` standard errors (and `NA` rows/columns in `XtX_inv`), following
[`stats::lm()`](https://rdrr.io/r/stats/lm.html).

## Examples

``` r
# Generate example whitened data
n_time <- 200
n_pred <- 3
n_voxels <- 50
Xw <- matrix(rnorm(n_time * n_pred), n_time, n_pred)
Yw <- matrix(rnorm(n_time * n_voxels), n_time, n_voxels)

# Compute standard errors
se_result <- sandwich_from_whitened_resid(Xw, Yw, type = "iid")

# Extract standard errors for first voxel
se_voxel1 <- se_result$se[, 1]
```
