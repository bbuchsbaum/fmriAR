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
  q_max = 2L,
  exact_first = c("ar1", "none"),
  pooling = c("global", "run", "parcel"),
  parcels = NULL,
  parcel_sets = NULL,
  multiscale = c("pacf_weighted", "acvf_pooled"),
  ms_mode = NULL,
  p_target = NULL,
  beta = 0.5,
  ms_weights = c("adaptive", "size"),
  hr_iter = 0L,
  step1 = c("burg", "yw"),
  design = NULL,
  acvf_correction = NULL,
  correction_max_lag = "auto",
  data_norm = NULL,
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

  Optional run labels, one per timepoint. Each label must occupy one
  contiguous block and may not be missing.

- censor:

  Optional integer vector of 1-based timepoint indices to exclude from
  AR parameter estimation, or a logical vector of length `nrow(resid)`
  where `TRUE`

  indicates censored timepoints. Censored frames (e.g.,
  motion-corrupted) are excluded when computing autocorrelations. Each
  run's estimation uses only its own valid (non-censored) segments.

- method:

  Either "ar" or "arma". ARMA models are estimated by a Hannan–Rissanen
  regression pooled over voxels (every voxel is treated as a replicate
  series sharing the coefficients), with all lags kept inside contiguous
  run/censor segments. With `design`, the estimate is then refitted to
  the design-corrected autocovariance. The automatic ARMA order search
  (`p`/`q = "auto"`) tends to choose the largest orders offered when
  many voxels are pooled, and warns when it does; set the orders when
  theory suggests them (AR(1) noise plus white noise is ARMA(1,1)).

- p:

  AR order (integer or "auto" if method == "ar"). An explicit order is
  honoured even when it exceeds `p_max`.

- q:

  MA order (integer), or `"auto"` (ARMA only) to choose it by BIC from
  `0:q_max`.

- p_max:

  Maximum AR order when `p = "auto"`. The order is chosen by BIC whose
  sample size is the number of frames times the effective number of
  independent voxels (a design effect from the mean inter-voxel
  correlation), so pooling many independent voxels can support a higher
  order than one series would. For `method = "arma"`, `p = "auto"`
  searches `0:p_max` jointly with the MA order by BIC on the pooled
  Hannan–Rissanen regression.

- q_max:

  Maximum MA order when `q = "auto"`.

- exact_first:

  How
  [`whiten_apply()`](https://bbuchsbaum.github.io/fmriAR/reference/whiten_apply.md)
  treats the start of each run or post-censoring segment. `"ar1"`
  (default; the name is historical) applies the exact stationary
  initialisation of the fitted AR/ARMA model, i.e. exact GLS whitening
  within each segment; for AR(1) this is the familiar `sqrt(1 - phi^2)`
  scaling of the first sample. `"none"` uses the conditional filter with
  a zero pre-sample.

- pooling:

  Combine parameters across runs or parcels ("global", "run", "parcel").

- parcels:

  Integer vector (length = ncol(resid)) giving fine parcel memberships
  when `pooling = "parcel"`.

- parcel_sets:

  Optional named list with entries `coarse`, `medium`, `fine` of equal
  length specifying nested parcel labels for multi-scale pooling.

- multiscale:

  Multi-scale pooling under `pooling = "parcel"`: `"pacf_weighted"`
  (pool partial autocorrelations, which always yields a stationary
  filter) or `"acvf_pooled"` (pool autocorrelations and refit), or
  `TRUE`/`FALSE`. By default, pooling is on (`"pacf_weighted"`) when
  `parcel_sets` is supplied and off otherwise, whatever `p` is. With
  `multiscale = TRUE` and no `parcel_sets`, each parcel is pooled toward
  the estimate over all voxels.

- ms_mode:

  Explicit multiscale mode when `multiscale` is logical.

- p_target:

  Target AR order for multi-scale pooling (defaults to `p_max`).

- beta:

  Size exponent for `ms_weights = "size"` (default 0.5).

- ms_weights:

  How finer scales are pooled toward coarser ones. `"adaptive"`
  (default) uses empirical-Bayes weights: each estimate is shrunk toward
  its parent's by `tau2 / (tau2 + s2)`, where `s2` is its sampling
  variance (from the number of frames and the effective number of
  independent voxels) and `tau2` the spread of estimates around their
  parents, estimated from the data. Parcels that genuinely differ keep
  their estimates; small, noisy parcels are pulled toward their parents.
  `"size"` uses the weights of earlier versions, proportional to
  `size^beta / (1 + dispersion)`, which do not depend on whether the
  scales agree.

- hr_iter:

  Number of Hannan–Rissanen refinement iterations for ARMA.

- step1:

  Preliminary high-order AR fit method for HR ("burg" or "yw").

- design:

  Design matrix (timepoints x regressors) whose projection produced
  `resid`. Supplying it corrects the downward bias that projecting a
  design out of the data puts into the autocovariance, and hence into
  the coefficients; drift terms cause most of it. Pass it whenever
  `resid` came from a regression: without it, pooling many voxels can
  also lead BIC to choose too high an order. When `resid` is omitted and
  `Y` and `X` are given, `X` is used automatically. Without a design,
  `fit_noise()` still corrects the part of the bias its own per-run
  centring introduces, which needs no design; `design = FALSE` turns
  every correction off and returns the plain estimator. A numerical
  orthogonality check rejects detectable mismatches, including typical
  raw, GLS, robust, or unrelated residuals. Directions an
  ill-conditioned design barely determines (singular values below 1e-6
  of the largest) are left out of the check, since solvers disagree on
  whether they belong to the fit; for a fit at a lower rank than
  `design`, passing a basis of the fitted column space makes the
  correction exact. Orthogonality cannot prove the exact
  residual-forming provenance, so the caller remains responsible for
  supplying the matching design. Supported for every pooling mode with
  `method = "ar"`, where under `pooling = "parcel"` the per-run bias
  maps are combined in proportion to the lag pairs each run contributes,
  and for `method = "arma"`, where the Hannan–Rissanen estimate is
  refitted by minimum distance to the corrected autocovariance.

- acvf_correction:

  Precomputed bias matrices from
  [`acvf_bias_matrix()`](https://bbuchsbaum.github.io/fmriAR/reference/acvf_bias_matrix.md),
  as an alternative to `design` when many datasets share one design. A
  single matrix is applied to every run; a list is matched against the
  runs in order. Mutually exclusive with `design`.

- correction_max_lag:

  Lag budget for the bias correction: `"auto"` (default) or a number.
  The correction solves a system truncated at this lag, assuming the
  noise autocovariance has died out by then, so too short a budget
  leaves bias behind; too long a one makes the system ill-conditioned,
  faster the more drift terms the design has. `"auto"` starts at 25 and
  lowers the budget as far as needed for a well-conditioned system
  (reciprocal condition number at least 1e-3). A number is used as given
  unless it exceeds a third of a run's valid volumes or makes the system
  ill-conditioned; it is then lowered the same way, with a warning. A
  warning also reports when the fitted noise model is still clearly
  autocorrelated at the end of the budget, where the correction becomes
  partial.

- data_norm:

  Optional norm of the data `resid` was computed from: one number for
  the whole matrix (`sqrt(sum(Y^2))`) or one per column
  (`sqrt(colSums(Y^2))`). It gives the orthogonality check behind
  `design` a scale for residuals the design explains almost exactly.
  Those residuals are floating-point error, and unless they come from
  the same QR factorisation the check uses, they are not orthogonal to
  the design relative to their own size. With `data_norm`, a projection
  below `sqrt(.Machine$double.eps) * data_norm` is accepted as rounding
  error. Computed from `Y` when the residuals are formed from `Y` and
  `X`, and by
  [`whiten()`](https://bbuchsbaum.github.io/fmriAR/reference/whiten.md).
  Ignored unless `design` is supplied.

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
  covariance. Global AR coefficients are fitted by Yule–Walker on the
  frame-weighted average of the runs' autocorrelations.

- `sigma2`: list of innovation variances, matching `gamma`. For
  `method = "ar"` it is derived as `gamma_0 - sum_k phi_k gamma_k` from
  the coefficients stored on the plan so the two are always mutually
  consistent, and is `NA` whenever `gamma` does not reach lag
  `length(phi)` – heavy censoring can truncate it that far, and a
  partial sum would overstate the innovation variance rather than report
  that it is unavailable. For `method = "arma"` it is the mean squared
  exact-whitened residual over all voxels and valid frames.

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
