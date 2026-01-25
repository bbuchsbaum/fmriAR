# GLS standard errors from whitened residuals

GLS standard errors from whitened residuals

## Usage

``` r
sandwich_from_whitened_resid(
  Xw,
  Yw,
  beta = NULL,
  type = c("iid", "hc0"),
  df_mode = c("rankX", "n-p"),
  runs = NULL
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

  Either "iid" (default) or "hc0" for a robust sandwich.

- df_mode:

  Degrees-of-freedom mode: "rankX" (default) or "n-p".

- runs:

  Optional run labels (reserved for future per-run scaling).

## Value

List containing standard errors, innovation variances, and XtX inverse.

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
