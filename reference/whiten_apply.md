# Apply a whitening plan to design and data matrices

Apply a whitening plan to design and data matrices

## Usage

``` r
whiten_apply(
  plan,
  X,
  Y,
  runs = NULL,
  run_starts = NULL,
  censor = NULL,
  parcels = NULL,
  inplace = FALSE,
  parallel = TRUE
)
```

## Arguments

- plan:

  Whitening plan from
  [`fit_noise()`](https://bbuchsbaum.github.io/fmriAR/reference/fit_noise.md).

- X:

  Numeric matrix of predictors (time x regressors).

- Y:

  Numeric matrix of data (time x voxels).

- runs:

  Optional run labels, one per row. Each label must occupy one
  contiguous block and may not be missing.

- run_starts:

  Optional 0-based run start indices (alternative to `runs`).

- censor:

  Optional indices of censored TRs (1-based); filter resets after gaps.

- parcels:

  Optional parcel labels (length = ncol(Y)) when using parcel plans.

- inplace:

  Modify inputs in place (logical).

- parallel:

  Use OpenMP parallelism if available.

## Value

List with whitened data. Parcel plans return `X_by` per parcel; others
return a single `X` matrix.

## Examples

``` r
# Create example design matrix and data
n_time <- 200
n_pred <- 3
n_voxels <- 50
X <- matrix(rnorm(n_time * n_pred), n_time, n_pred)
Y <- X %*% matrix(rnorm(n_pred * n_voxels), n_pred, n_voxels) +
     matrix(rnorm(n_time * n_voxels), n_time, n_voxels)

# Fit noise model from residuals
residuals <- Y - X %*% solve(crossprod(X), crossprod(X, Y))
plan <- fit_noise(residuals, method = "ar", p = 2)

# Apply whitening
whitened <- whiten_apply(plan, X, Y)
Xw <- whitened$X
Yw <- whitened$Y
```
