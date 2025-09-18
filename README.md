# fmriAR

Fast AR/ARMA prewhitening for fMRI design and data.

## Install (dev)

```r
Rcpp::compileAttributes()
devtools::document()
devtools::load_all()
```

## Example

```r
# X: (n x p), Y: (n x v), runs: length n
res   <- Y - X %*% qr.solve(X, Y)
plan  <- fit_noise(res, runs = runs, method = "ar", p = "auto", pooling = "global")
xyw   <- whiten_apply(plan, X, Y, runs = runs)
fit   <- lm.fit(xyw$X, xyw$Y)
se    <- sandwich_from_whitened_resid(xyw$X, xyw$Y, beta = fit$coefficients)
ac    <- acorr_diagnostics(xyw$Y - xyw$X %*% fit$coefficients)
```
