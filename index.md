# fmriAR

Fast AR and ARMA prewhitening for fMRI GLM workflows. Estimate a noise
model from residuals, apply the same filter to the design matrix and the
data, then fit a standard linear model on the whitened series.

The C++ core (RcppArmadillo) keeps large-scale whitening fast.
Estimation is run-aware and censor-aware: lag products never cross a run
boundary or a scrubbed frame.

## Installation

``` r

install.packages("fmriAR")
```

Development version:

``` r

# install.packages("remotes")
remotes::install_github("bbuchsbaum/fmriAR")
```

Documentation: <https://bbuchsbaum.github.io/fmriAR/>

## Typical workflow

1.  Fit an initial OLS model and take residuals.
2.  Estimate an AR/ARMA plan with
    [`fit_noise()`](https://bbuchsbaum.github.io/fmriAR/reference/fit_noise.md).
3.  Whiten both `X` and `Y` with
    [`whiten_apply()`](https://bbuchsbaum.github.io/fmriAR/reference/whiten_apply.md)
    (or do both steps with
    [`whiten()`](https://bbuchsbaum.github.io/fmriAR/reference/whiten.md)).
4.  Fit the linear model on the whitened matrices.
5.  Optionally compute sandwich standard errors and residual
    autocorrelation.

``` r

library(fmriAR)

# X: design (n x p), Y: data (n x voxels), runs: one label per timepoint
resid <- Y - X %*% qr.solve(X, Y)

plan <- fit_noise(
  resid,
  runs = runs,
  method = "ar",
  p = "auto",
  pooling = "global",
  design = X            # optional: undo residual-projection bias
)

xyw <- whiten_apply(plan, X, Y, runs = runs)
fit <- lm.fit(xyw$X, xyw$Y)
se  <- sandwich_from_whitened_resid(xyw$X, xyw$Y, beta = fit$coefficients)
ac  <- acorr_diagnostics(xyw$Y - xyw$X %*% fit$coefficients, runs = runs)
```

With censoring, the scrubbed rows stay in the whitened output and are
listed in `xyw$censor`; drop them before fitting, e.g.
`lm.fit(xyw$X[-xyw$censor, ], xyw$Y[-xyw$censor, ])`.

One-step shortcut (fits the plan from `Y` and `X` internally):

``` r

xyw <- whiten(X, Y, runs = runs, method = "ar", p = "auto")
```

## What the package does

- **AR and ARMA plans.** `method = "ar"` fits Yule–Walker on
  autocovariances pooled over voxels and selects the order by BIC
  (`p = "auto"`), counting pooled voxels at their effective number.
  `method = "arma"` uses Hannan–Rissanen (1982) pooled over voxels, with
  `p`/`q = "auto"` for BIC order selection.
- **Exact whitening at run starts.** Every run and post-censoring
  segment is whitened with the exact stationary start of the fitted
  AR/ARMA model, so the result equals exact GLS within each segment
  (`exact_first = "ar1"`, the default).
- **Pooling.** `"global"` (one filter), `"run"` (one per run), or
  `"parcel"` (one per parcel, estimated from its voxels, with optional
  multiscale shrinkage via `parcel_sets`).
- **Censoring.** Pass motion-scrubbed frames as indices or a logical
  mask; they are dropped from estimation and treated as segment breaks
  when whitening.
  [`whiten_apply()`](https://bbuchsbaum.github.io/fmriAR/reference/whiten_apply.md)
  reuses the plan’s censor set for data of the same length.
- **Residual-bias correction.** Autocovariance from GLM residuals is
  biased low (`E[ehat ehat'] = M Sigma M`). Pass `design = X` to
  [`fit_noise()`](https://bbuchsbaum.github.io/fmriAR/reference/fit_noise.md)
  or
  [`noise_acvf()`](https://bbuchsbaum.github.io/fmriAR/reference/noise_acvf.md)
  to undo that bias (`method = "ar"`, any pooling, with or without
  censoring). Cache the map with
  [`acvf_bias_matrix()`](https://bbuchsbaum.github.io/fmriAR/reference/acvf_bias_matrix.md)
  when many datasets share a design.
- **Standard errors.**
  [`sandwich_from_whitened_resid()`](https://bbuchsbaum.github.io/fmriAR/reference/sandwich_from_whitened_resid.md)
  gives iid, HC0, or Newey–West HAC (`type = "hac"`, within runs)
  standard errors, and handles rank-deficient designs.
- **Noise scale on the plan.** An `fmriAR_plan` now stores `gamma`
  (voxel-scale autocovariance) and `sigma2` (innovation variance) per
  pooling unit, not just the AR/MA coefficients.
- **Autocovariance without a model.**
  [`noise_acvf()`](https://bbuchsbaum.github.io/fmriAR/reference/noise_acvf.md)
  returns the same run- and censor-aware covariances
  [`fit_noise()`](https://bbuchsbaum.github.io/fmriAR/reference/fit_noise.md)
  uses internally, plus pair counts so you can see how much data backs
  each lag.
- **AFNI-style restricted AR.**
  [`afni_restricted_plan()`](https://bbuchsbaum.github.io/fmriAR/reference/afni_restricted_plan.md)
  builds a plan from AFNI root parameters (Cox, 2012) for pipeline
  comparison, including AFNI’s additive white-noise ratio `vrt`.

See
[`vignette("fmriAR-introduction")`](https://bbuchsbaum.github.io/fmriAR/articles/fmriAR-introduction.md)
and
[`?fit_noise`](https://bbuchsbaum.github.io/fmriAR/reference/fit_noise.md)
for parcel pooling, ARMA, and AFNI examples.

## Options

- `options(fmriAR.max_threads = n)` — cap the OpenMP threads used for
  whitening. OpenMP is used automatically when the compiler supports it;
  results are identical for any thread count.
- `options(fmriAR.use_cpp_hr = TRUE)` — use the C++ single-series
  Hannan–Rissanen estimator (default) in the internal `hr_arma()`
  helper. Set `FALSE` to fall back to the R implementation.
  [`fit_noise()`](https://bbuchsbaum.github.io/fmriAR/reference/fit_noise.md)
  uses the pooled estimator regardless.

## Validation

`tools/validation/accuracy_benchmark.R` (in the source repository) runs
seeded accuracy, calibration, and timing scenarios through the exported
API, so any two installed versions can be compared directly:

``` sh
Rscript tools/validation/accuracy_benchmark.R <lib_path_or_-> results.rds
```

## References

- Hannan, E. J., & Rissanen, J. (1982). Recursive estimation of mixed
  autoregressive-moving average order. *Biometrika*, 69(1), 81–94.
  <https://doi.org/10.1093/biomet/69.1.81>
- Cox, R. W. (2012). AFNI: What, where, how? *NeuroImage*, 62(2),
  743–747. <https://doi.org/10.1016/j.neuroimage.2011.08.056>
