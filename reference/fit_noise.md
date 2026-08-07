# Fit an AR/ARMA noise model (run-aware) and return a whitening plan

Fit an AR/ARMA noise model (run-aware) and return a whitening plan

## Usage

``` r
fit_noise(
  resid = NULL,
  Y = NULL,
  X = NULL,
  runs = NULL,
  censor = NULL,
  method = c("ar", "arma"),
  p = "auto",
  q = 0L,
  p_max = 6L,
  exact_first = c("ar1", "none"),
  pooling = c("global", "run", "parcel"),
  parcels = NULL,
  parcel_sets = NULL,
  multiscale = c("pacf_weighted", "acvf_pooled"),
  ms_mode = NULL,
  p_target = NULL,
  beta = 0.5,
  hr_iter = 0L,
  step1 = c("burg", "yw"),
  parallel = FALSE
)
```

## Arguments

- resid:

  Numeric matrix (time x voxels) of residuals from an initial OLS fit.

- Y:

  Optional data matrix used to compute residuals when `resid` is
  omitted.

- X:

  Optional design matrix used with `Y` to compute residuals.

- runs:

  Optional integer vector of run identifiers.

- censor:

  Optional integer vector of 1-based timepoint indices to exclude from
  AR parameter estimation, or a logical vector of length `nrow(resid)`
  where `TRUE`

  indicates censored timepoints. Censored frames (e.g.,
  motion-corrupted) are excluded when computing autocorrelations. Each
  run's estimation uses only its own valid (non-censored) segments.

- method:

  Either "ar" or "arma".

- p:

  AR order (integer or "auto" if method == "ar").

- q:

  MA order (integer).

- p_max:

  Maximum AR order when `p = "auto"`.

- exact_first:

  Apply exact AR(1) scaling at segment starts ("ar1" or "none").

- pooling:

  Combine parameters across runs or parcels ("global", "run", "parcel").

- parcels:

  Integer vector (length = ncol(resid)) giving fine parcel memberships
  when `pooling = "parcel"`.

- parcel_sets:

  Optional named list with entries `coarse`, `medium`, `fine` of equal
  length specifying nested parcel labels for multi-scale pooling.

- multiscale:

  Multi-scale pooling mode when `parcel_sets` is supplied
  ("pacf_weighted" or "acvf_pooled"), or `TRUE/FALSE` to toggle pooling.

- ms_mode:

  Explicit multiscale mode when `multiscale` is logical.

- p_target:

  Target AR order for multi-scale pooling (defaults to `p_max`).

- beta:

  Size exponent for multi-scale weights (default 0.5).

- hr_iter:

  Number of Hannan–Rissanen refinement iterations for ARMA.

- step1:

  Preliminary high-order AR fit method for HR ("burg" or "yw").

- parallel:

  Reserved for future parallel estimation (logical).

## Value

An object of class `fmriAR_plan` used by
[`whiten_apply()`](https://bbuchsbaum.github.io/fmriAR/reference/whiten_apply.md).
Besides the AR/MA coefficients the plan carries the noise scale and
shape it was fitted from, so consumers can reconstruct the covariance it
implies rather than only its correlation structure:

- `gamma`: list of autocovariance vectors, one per pooling unit – a
  single entry for `pooling = "global"`, one per run for
  `pooling = "run"`. Lags run 0 to the highest the data supported, which
  is governed by `p_max` and the run length rather than by `p`, so
  `fit_noise(p = 1, p_max = 6)` returns seven values, not two. Under
  global pooling every run is truncated to the shortest available length
  before averaging, since a zero-padded autocovariance is not a valid
  covariance.

- `sigma2`: list of innovation variances, matching `gamma`, derived as
  `gamma_0 - sum_k phi_k gamma_k` from the coefficients stored on the
  plan so the two are always mutually consistent. `NA` for
  `method = "arma"`, where no comparably cheap voxel-scale innovation
  variance is available, and `NA` whenever `gamma` does not reach lag
  `length(phi)` – heavy censoring can truncate it that far, and a
  partial sum would overstate the innovation variance rather than report
  that it is unavailable.

- `gamma_by_parcel`, `sigma2_by_parcel`: the same quantities per parcel
  when `pooling = "parcel"`, keyed like `phi_by_parcel`.

For a run-stationary noise process with autocovariance `gamma`, the
covariance of the data within a run is the Toeplitz matrix built from
it, which is what makes design-specific variance calculations possible
downstream without refitting.

## Examples

``` r
# Generate example data with AR(1) structure
n_time <- 200
n_voxels <- 50
phi_true <- 0.5

# Simulate residuals with AR(1) structure
resid <- matrix(0, n_time, n_voxels)
for (v in 1:n_voxels) {
  e <- rnorm(n_time)
  resid[1, v] <- e[1]
  for (t in 2:n_time) {
    resid[t, v] <- phi_true * resid[t-1, v] + e[t]
  }
}

# Fit AR model
plan <- fit_noise(resid, method = "ar", p = 1)

# With multiple runs
runs <- rep(1:2, each = 100)
plan_runs <- fit_noise(resid, runs = runs, method = "ar", pooling = "run")
```
