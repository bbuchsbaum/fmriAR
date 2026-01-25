# Autocorrelation diagnostics for residuals

Autocorrelation diagnostics for residuals

## Usage

``` r
acorr_diagnostics(
  resid,
  runs = NULL,
  max_lag = 20L,
  aggregate = c("mean", "median", "none")
)
```

## Arguments

- resid:

  Numeric matrix (time x voxels), typically whitened residuals.

- runs:

  Optional run labels.

- max_lag:

  Maximum lag to evaluate.

- aggregate:

  Aggregation across voxels: "mean", "median", or "none".

## Value

List of autocorrelation values and nominal confidence interval.

## Examples

``` r
# Generate example residuals with some autocorrelation
n_time <- 200
n_voxels <- 50
resid <- matrix(rnorm(n_time * n_voxels), n_time, n_voxels)

# Add some AR(1) structure
for (v in 1:n_voxels) {
  resid[, v] <- filter(resid[, v], filter = 0.3, method = "recursive")
}

# Check autocorrelation
acorr_check <- acorr_diagnostics(resid, max_lag = 10, aggregate = "mean")

# Examine lag-1 autocorrelation
lag1_acorr <- acorr_check$acf[2]  # First element is lag-0 (always 1)
```
