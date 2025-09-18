# fmriAR

fmriAR provides fast AR/ARMA-based prewhitening for fMRI GLM workflows. It estimates voxel-wise or parcel-based noise models, applies segment-aware whitening, and exposes diagnostics that make it easy to confirm residual independence.

## Key capabilities
- Automatic AR/ARMA order selection via Hannan–Rissanen initialization and iterative refinement
- Segment-aware whitening that respects run boundaries and optional multiscale pooling across parcels
- Convenience helpers to whiten design matrices, refit GLMs, and inspect autocorrelation diagnostics

## Installation

```r
# install.packages("remotes")  # only needed once
remotes::install_github("bbuchsbaum/fmriAR")
```

For local development you can load the package with devtools:

```r
devtools::document()
devtools::load_all()
```

## Quick start

```r
# X: design matrix (n x p), Y: voxel data (n x v), runs: factor or integer run labels
res   <- Y - X %*% qr.solve(X, Y)                      # pre-fit residuals
plan  <- fit_noise(res, runs = runs, method = "ar",    # estimate AR model
                   p = "auto", pooling = "global")
xyw   <- whiten_apply(plan, X, Y, runs = runs)         # whiten design and data
fit   <- lm.fit(xyw$X, xyw$Y)
se    <- sandwich_from_whitened_resid(xyw$X, xyw$Y, beta = fit$coefficients)
ac    <- acorr_diagnostics(xyw$Y - xyw$X %*% fit$coefficients)
```

See `vignettes/` and `?fit_noise` for more detailed workflows, including multiscale pooling and ARMA whitening.
