# Fit and apply whitening in one call

Fit and apply whitening in one call

## Usage

``` r
whiten(X, Y, runs = NULL, censor = NULL, ...)
```

## Arguments

- X:

  Design matrix (time x regressors).

- Y:

  Data matrix (time x voxels).

- runs:

  Optional run labels.

- censor:

  Optional censor indices.

- ...:

  Additional parameters passed to
  [`fit_noise()`](https://bbuchsbaum.github.io/fmriAR/reference/fit_noise.md).

## Value

List with whitened `X` and `Y` matrices.

## Examples

``` r
# Create example data
n_time <- 200
n_pred <- 3
n_voxels <- 50
X <- matrix(rnorm(n_time * n_pred), n_time, n_pred)
Y <- X %*% matrix(rnorm(n_pred * n_voxels), n_pred, n_voxels) +
     matrix(rnorm(n_time * n_voxels, sd = 2), n_time, n_voxels)

# One-step whitening
whitened <- whiten(X, Y, method = "ar", p = 2)
```
