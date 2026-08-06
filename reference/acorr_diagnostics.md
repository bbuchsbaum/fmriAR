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

  Optional run labels, length `nrow(resid)`. When supplied, each run is
  centred separately and no lag product spans a run boundary.

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

# Examine lag-1 autocorrelation (lag 0 is not returned, so acf[1] is lag 1)
lag1_acorr <- acorr_check$acf[1]
```
