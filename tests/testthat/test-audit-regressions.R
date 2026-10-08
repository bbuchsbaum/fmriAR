# Regression and property tests for the correctness/performance audit.
# Each block names the defect it guards against.

sim_arma_mat <- function(n, v, phi = numeric(0), theta = numeric(0)) {
  vapply(seq_len(v), function(j) {
    as.numeric(stats::arima.sim(list(ar = phi, ma = theta), n, n.start = 200L))
  }, numeric(n))
}

# Exact GLS oracle: forward-solve with the Cholesky factor of the stationary
# covariance (unit innovation variance).
gls_oracle <- function(y, phi, theta) {
  n <- length(y)
  g0 <- sum(c(1, stats::ARMAtoMA(phi, theta, 5000L))^2)
  S <- stats::toeplitz(stats::ARMAacf(ar = phi, ma = theta, lag.max = n - 1L) * g0)
  forwardsolve(t(chol(S)), y)
}

# --- order selection ----------------------------------------------------------

test_that("an explicit p above p_max is honoured, not silently capped", {
  set.seed(9101)
  R <- sim_arma_mat(400, 10, c(0.3, 0.1, 0, 0, 0, 0, 0, 0.1))
  plan <- fit_noise(R, p = 8L)
  expect_equal(plan$order[["p"]], 8L)
  expect_length(plan$phi[[1]], 8L)
  plan_par <- fit_noise(R, p = 8L, pooling = "parcel", parcels = rep(1:2, 5))
  expect_equal(plan_par$order[["p"]], 8L)
  expect_error(fit_noise(R, p = 1.5), "whole number")
})

test_that("BIC uses the standard penalty and does not over-select", {
  # With the fit term doubled, BIC chose p > 1 for a true AR(1) 17% of the time.
  set.seed(9102)
  sel1 <- replicate(150, length(fit_noise(sim_arma_mat(200, 1, 0.5), p = "auto")$phi[[1]]))
  expect_lt(mean(sel1 > 1), 0.08)
  expect_gt(mean(sel1 == 1), 0.9)
  # ...while still detecting a genuine second lag.
  sel2 <- replicate(100, length(fit_noise(sim_arma_mat(200, 1, c(0.5, 0.25)), p = "auto")$phi[[1]]))
  expect_gt(mean(sel2 >= 2), 0.8)
})

# --- censoring ----------------------------------------------------------------

test_that("whiten_apply treats a logical censor mask like the equivalent indices", {
  set.seed(9103)
  R <- sim_arma_mat(200, 4, 0.5)
  X <- cbind(1, rnorm(200))
  mask <- rep(FALSE, 200); mask[c(50, 120, 121)] <- TRUE
  plan <- fit_noise(R, p = 1, censor = mask)
  a <- whiten_apply(plan, X, R, censor = mask)
  b <- whiten_apply(plan, X, R, censor = which(mask))
  expect_identical(a$Y, b$Y)
  expect_identical(a$X, b$X)
  expect_identical(a$censor, which(mask))
  expect_error(whiten_apply(plan, X, R, censor = mask[-1]), "one entry per timepoint")
  # whiten() forwards the mask to both steps.
  w1 <- whiten(X, R, censor = mask, p = 1)
  w2 <- whiten(X, R, censor = which(mask), p = 1)
  expect_identical(w1$Y, w2$Y)
})

test_that("whiten_apply falls back to the plan's censor set for matching data", {
  set.seed(9104)
  R <- sim_arma_mat(200, 4, 0.5)
  X <- cbind(1, rnorm(200))
  cens <- c(40L, 41L, 150L)
  plan <- fit_noise(R, p = 1, censor = cens)
  expect_identical(whiten_apply(plan, X, R)$Y,
                   whiten_apply(plan, X, R, censor = cens)$Y)
  # A different dataset (different length) does not inherit the training set.
  out <- whiten_apply(plan, X[1:150, ], R[1:150, ])
  expect_null(out$censor)
})

test_that("censored rows are reported so they can be dropped", {
  set.seed(9105)
  R <- sim_arma_mat(100, 2, 0.4)
  plan <- fit_noise(R, p = 1)
  out <- whiten_apply(plan, matrix(1, 100, 1), R, censor = c(10, 60))
  expect_identical(out$censor, c(10L, 60L))
  expect_null(whiten_apply(plan, matrix(1, 100, 1), R)$censor)
})

# --- plan / data consistency ------------------------------------------------------

test_that("a per-run plan refuses data with a different number of runs", {
  set.seed(9106)
  R <- sim_arma_mat(300, 3, 0.5)
  plan <- fit_noise(R, runs = rep(1:3, each = 100), pooling = "run", p = 1)
  expect_error(
    whiten_apply(plan, matrix(1, 200, 1), R[1:200, ], runs = rep(1:2, each = 100)),
    "3 'phi' sets"
  )
  # A global plan broadcasts to any number of runs.
  gplan <- fit_noise(R, runs = rep(1:3, each = 100), p = 1)
  expect_silent(whiten_apply(gplan, matrix(1, 200, 1), R[1:200, ], runs = rep(1:2, each = 100)))
})

# --- MA invertibility -------------------------------------------------------------

test_that("enforce_invertible_ma moves unit roots off the unit circle", {
  min_mod <- function(th) min(Mod(polyroot(c(1, th))))
  expect_gt(min_mod(enforce_invertible_ma(-1)), 1)
  expect_gt(min_mod(enforce_invertible_ma(1)), 1)
  expect_gt(min_mod(enforce_invertible_ma(c(0, -1))), 1)
  set.seed(9107)
  for (i in 1:50) {
    th <- enforce_invertible_ma(rnorm(3, sd = 1.5))
    expect_gte(min_mod(th), 1 / 0.99 - 1e-8)
  }
  # Comfortably invertible coefficients are returned untouched.
  expect_identical(enforce_invertible_ma(c(0.4, 0.1)), c(0.4, 0.1))
})

# --- exact stationary initialisation -------------------------------------------------

test_that("arma_acvf_cpp matches stats::ARMAacf", {
  cases <- list(list(0.5, numeric(0)), list(c(1.2, -0.5), numeric(0)),
                list(0.5, 0.4), list(c(0.6, 0.2), c(0.3, -0.2)),
                list(numeric(0), c(0.5, 0.3)))
  for (cs in cases) {
    ref <- stats::ARMAacf(cs[[1]], cs[[2]], 8) *
      sum(c(1, stats::ARMAtoMA(cs[[1]], cs[[2]], 5000L))^2)
    expect_equal(arma_acvf_cpp(cs[[1]], cs[[2]], 8L), unname(ref), tolerance = 1e-10)
  }
  expect_true(all(is.na(arma_acvf_cpp(1.1, numeric(0), 3L))))
})

test_that("exact whitening equals dense-Cholesky GLS for AR, MA and ARMA", {
  set.seed(9108)
  cases <- list(list(0.5, numeric(0)), list(c(1.2, -0.5), numeric(0)),
                list(c(0.3, 0.2, 0.1), numeric(0)), list(0.5, 0.4),
                list(c(0.6, 0.2), c(0.3, -0.2)), list(numeric(0), 0.7),
                list(c(0.3, 0.2, 0.1), 0.8))
  for (cs in cases) {
    y <- rnorm(70)
    starts <- c(0L, 30L, 33L)  # includes a 3-frame segment
    out <- arma_whiten_inplace(matrix(y, ncol = 1), matrix(0, 70, 0), cs[[1]], cs[[2]],
                               starts, exact_first = TRUE, parallel = FALSE)$Y
    ref <- c(gls_oracle(y[1:30], cs[[1]], cs[[2]]),
             gls_oracle(y[31:33], cs[[1]], cs[[2]]),
             gls_oracle(y[34:70], cs[[1]], cs[[2]]))
    expect_equal(drop(out), ref, tolerance = 1e-10)
  }
})

test_that("exact AR(1) start is the Prais-Winsten scaling, as before", {
  set.seed(9109)
  y <- rnorm(50); phi <- 0.7
  out <- drop(arma_whiten_inplace(matrix(y, ncol = 1), matrix(0, 50, 0), phi, numeric(0),
                                  0L, exact_first = TRUE, parallel = FALSE)$Y)
  expect_equal(out[1], y[1] * sqrt(1 - phi^2))
  expect_equal(out[-1], y[-1] - phi * y[-50])
})

test_that("whitened noise has unit variance at every timepoint, run starts included", {
  # The old filter whitened only AR(1) starts exactly: for AR(2) (1.2, -0.5)
  # the first two whitened samples of each run had variance 3.7 and 1.9.
  set.seed(9110)
  for (model in list(list(c(1.2, -0.5), numeric(0)), list(0.6, 0.5))) {
    Y <- sim_arma_mat(30, 3000, model[[1]], model[[2]])
    plan <- compat$plan_from_phi(model[[1]], theta = model[[2]], exact_first = TRUE,
                                 method = if (length(model[[2]])) "arma" else "ar")
    W <- whiten_apply(plan, matrix(1, 30, 1), Y)$Y
    v <- apply(W, 1, var)
    expect_true(all(abs(v - 1) < 0.1), info = paste(round(v[1:4], 2), collapse = " "))
    # and successive whitened samples are uncorrelated
    expect_lt(abs(cor(W[1, ], W[2, ])), 0.06)
  }
})

test_that("a non-stationary filter falls back to the conditional recursion", {
  y <- c(1, 2, 3, 4, 5)
  out <- drop(arma_whiten_inplace(matrix(y, ncol = 1), matrix(0, 5, 0), 1.0, numeric(0),
                                  0L, exact_first = TRUE, parallel = FALSE)$Y)
  expect_equal(out, c(1, 1, 1, 1, 1))
})

# --- OpenMP ----------------------------------------------------------------------

test_that("whitening is identical across thread counts and never aliases inputs", {
  set.seed(9111)
  Y <- matrix(rnorm(300 * 97), 300, 97)
  X <- matrix(rnorm(300 * 3), 300, 3)
  Y0 <- Y + 0; X0 <- X + 0
  ref <- arma_whiten_inplace(Y + 0, X + 0, c(0.5, 0.2), 0.3, c(0L, 120L),
                             exact_first = TRUE, parallel = FALSE)
  for (nt in c(1L, 2L, 4L)) {
    out <- arma_whiten_inplace(Y + 0, X + 0, c(0.5, 0.2), 0.3, c(0L, 120L),
                               exact_first = TRUE, parallel = TRUE, n_threads = nt)
    expect_identical(out$Y, ref$Y)
    expect_identical(out$X, ref$X)
  }
  # whiten_apply never modifies its arguments.
  plan <- compat$plan_from_phi(c(0.5, 0.2))
  withr::with_options(list(fmriAR.max_threads = 2L), {
    w <- whiten_apply(plan, X, Y)
  })
  expect_identical(Y, Y0)
  expect_identical(X, X0)
  expect_false(isTRUE(all.equal(w$Y, Y0)))
})

# --- pooled ARMA ---------------------------------------------------------------------

test_that("ARMA is estimated from voxel-level structure, not the voxel mean", {
  # Near-white voxels plus a weak shared slow signal: the voxel-mean series is
  # dominated by the shared component (phi ~ 0.86 under the old estimator).
  set.seed(9112)
  n <- 400; v <- 200
  g <- as.numeric(stats::filter(rnorm(n), 0.9, method = "recursive"))
  R <- matrix(rnorm(n * v), n, v) + 0.3 * g
  plan <- fit_noise(R, method = "arma", p = 1, q = 1)
  lag1 <- mean(apply(R, 2, function(y) stats::acf(y, 1, plot = FALSE)$acf[2]))
  # ARMA(1,1) lag-1 autocorrelation implied by the estimate
  rho1 <- stats::ARMAacf(plan$phi[[1]], plan$theta[[1]], 1)[2]
  expect_equal(unname(rho1), lag1, tolerance = 0.05)
  # whitened voxels are white
  w <- whiten_apply(plan, matrix(1, n, 1), R)$Y
  expect_lt(abs(acorr_diagnostics(w, max_lag = 1)$acf), 0.03)
})

test_that("pooled ARMA recovers parameters and agrees across pooling modes", {
  set.seed(9113)
  R <- sim_arma_mat(400, 30, 0.5, 0.4)
  runs <- rep(1:2, each = 200)
  pg <- fit_noise(R, runs = runs, method = "arma", p = 1, q = 1)
  pr <- fit_noise(R, runs = runs, method = "arma", p = 1, q = 1, pooling = "run")
  expect_equal(pg$phi[[1]], 0.5, tolerance = 0.1)
  expect_equal(pg$theta[[1]], 0.4, tolerance = 0.1)
  for (i in 1:2) {
    expect_equal(pr$phi[[i]], pg$phi[[1]], tolerance = 0.1)
    expect_equal(pr$theta[[i]], pg$theta[[1]], tolerance = 0.1)
  }
  expect_equal(pg$sigma2[[1]], 1, tolerance = 0.1)
})

test_that("single-series HR skips the long-AR burn-in", {
  set.seed(9114)
  est <- replicate(200, {
    y <- as.numeric(stats::arima.sim(list(ar = 0.5, ma = 0.4), 150))
    unlist(hr_arma(y, 1, 1)[c("phi", "theta")])
  })
  expect_lt(abs(mean(est[1, ]) - 0.5), 0.04)
  expect_lt(abs(mean(est[2, ]) - 0.4), 0.04)
  # R fallback agrees with the C++ path.
  y <- as.numeric(stats::arima.sim(list(ar = 0.5, ma = 0.4), 1000))
  cpp <- hr_arma(y, 1, 1)
  r <- withr::with_options(list(fmriAR.use_cpp_hr = FALSE), hr_arma(y, 1, 1, step1 = "yw"))
  expect_lt(abs(cpp$phi - r$phi), 0.05)
  expect_lt(abs(cpp$theta - r$theta), 0.05)
})

test_that("AFNI MA(1) is estimated from the voxels, not their mean", {
  set.seed(9115)
  n <- 300; v <- 100
  # A weak shared component: negligible in each voxel (lag-1 autocorrelation
  # ~0.02) but dominant in the 100-voxel mean.
  g <- as.numeric(stats::filter(rnorm(n), 0.95, method = "recursive"))
  R <- matrix(rnorm(n * v), n, v) + 0.05 * g
  plan <- afni_restricted_plan(R, p = 3L, roots = list(a = 0, r1 = 0, t1 = 0))
  # AR part is the identity here, so the voxel noise is ~white: theta ~ 0.
  expect_lt(abs(plan$theta[[1]]), 0.06)
})

# --- global pooling ------------------------------------------------------------------

test_that("global AR pooling is invariant to per-run scale", {
  set.seed(9116)
  R <- sim_arma_mat(300, 10, c(0.5, 0.2))
  runs <- rep(1:2, each = 150)
  R2 <- R; R2[runs == 2, ] <- R2[runs == 2, ] * 10
  p1 <- fit_noise(R, runs = runs, p = 2)
  p2 <- fit_noise(R2, runs = runs, p = 2)
  expect_equal(p1$phi[[1]], p2$phi[[1]], tolerance = 0.02)
})

test_that("global AR(p) pooling is at least as accurate as per-run averaging", {
  set.seed(9117)
  truth <- c(0.5, 0.2)
  err <- t(replicate(40, {
    lens <- c(120, 200, 80)
    runs <- rep(seq_along(lens), lens)
    R <- do.call(rbind, lapply(lens, function(L) sim_arma_mat(L, 10, truth)))
    g <- fit_noise(R, runs = runs, p = 2)$phi[[1]]
    pr <- fit_noise(R, runs = runs, p = 2, pooling = "run")$phi
    avg <- colSums(do.call(rbind, pr) * lens) / sum(lens)
    c(sum((g - truth)^2), sum((avg - truth)^2))
  }))
  expect_lte(mean(err[, 1]), mean(err[, 2]) * 1.05)
})

# --- multiscale ------------------------------------------------------------------------

test_that("multiscale pooling is invariant to the units of the data", {
  set.seed(9118)
  n <- 200; v <- 64
  R <- sim_arma_mat(n, v, 0.4) * rep(runif(v, 0.5, 2), each = n)
  ps <- list(coarse = rep(1:2, each = 32), medium = rep(1:4, each = 16),
             fine = rep(1:8, each = 8))
  for (mode in c("pacf_weighted", "acvf_pooled")) {
    a <- fit_noise(R, pooling = "parcel", parcels = ps$fine, parcel_sets = ps,
                   p = 1, multiscale = mode)
    b <- fit_noise(R * 100, pooling = "parcel", parcels = ps$fine, parcel_sets = ps,
                   p = 1, multiscale = mode)
    expect_equal(unlist(a$phi_by_parcel), unlist(b$phi_by_parcel), tolerance = 1e-8)
  }
})

test_that("acvf_pooled under p = 'auto' fits the selected order, not p_max", {
  set.seed(9119)
  R <- sim_arma_mat(300, 32, 0.5)
  ps <- list(coarse = rep(1:2, each = 16), medium = rep(1:4, each = 8),
             fine = rep(1:8, each = 4))
  plan <- fit_noise(R, pooling = "parcel", parcels = ps$fine, parcel_sets = ps,
                    p = "auto", p_max = 6, multiscale = "acvf_pooled")
  expect_lt(plan$order[["p"]], 6L)
})

# --- parcel design cache ----------------------------------------------------------------

test_that("parcels sharing a filter share one whitened design with identical values", {
  set.seed(9120)
  R <- sim_arma_mat(120, 6, 0.4)
  X <- cbind(1, rnorm(120))
  plan <- compat$plan_from_phi(list(`1` = 0.4, `2` = 0.4, `3` = 0.2),
                               parcels = c(1, 1, 2, 2, 3, 3), pooling = "parcel")
  out <- whiten_apply(plan, X, R)
  expect_identical(out$X_by[["1"]], out$X_by[["2"]])
  ref3 <- whiten_apply(compat$plan_from_phi(0.2), X, R)$X
  expect_equal(out$X_by[["3"]], ref3)
})

# --- diagnostics --------------------------------------------------------------------

test_that("acorr_diagnostics aggregates per-voxel autocorrelations", {
  set.seed(9121)
  n <- 300; v <- 50
  g <- as.numeric(stats::filter(rnorm(n), 0.9, method = "recursive"))
  R <- matrix(rnorm(n * v), n, v) + 0.5 * g
  per_vox <- vapply(seq_len(v), function(j) stats::acf(R[, j], 3, plot = FALSE)$acf[-1], numeric(3))
  expect_equal(acorr_diagnostics(R, max_lag = 3, aggregate = "none")$acf, per_vox,
               tolerance = 1e-12)
  expect_equal(acorr_diagnostics(R, max_lag = 3)$acf, rowMeans(per_vox), tolerance = 1e-12)
  expect_equal(acorr_diagnostics(R, max_lag = 3, aggregate = "median")$acf,
               apply(per_vox, 1, median), tolerance = 1e-12)
  # With runs, no product spans the boundary: two independent constant-offset
  # runs of white noise show no lag-1 structure.
  runs <- rep(1:2, each = 150)
  W <- matrix(rnorm(n * v), n, v) + ifelse(runs == 1, -5, 5)
  expect_lt(abs(acorr_diagnostics(W, runs = runs, max_lag = 1)$acf), 0.02)
})

# --- standard errors ------------------------------------------------------------------

test_that("sandwich handles rank-deficient designs like lm", {
  set.seed(9122)
  x <- rnorm(100)
  X <- cbind(1, x, 2 * x)
  y <- 1 + x + rnorm(100)
  out <- sandwich_from_whitened_resid(X, matrix(y, ncol = 1))
  fit <- summary(stats::lm(y ~ x))
  expect_true(is.na(out$se[3, 1]))
  expect_equal(out$se[1:2, 1], unname(fit$coefficients[, 2]), tolerance = 1e-10)
  expect_equal(out$df, 98)
})

test_that("vectorised HC0 equals the explicit sandwich", {
  set.seed(9123)
  X <- cbind(1, matrix(rnorm(200 * 3), 200))
  Y <- matrix(rnorm(200 * 7), 200) * (1 + abs(X[, 2]))
  out <- sandwich_from_whitened_resid(X, Y, type = "hc0")
  XtXi <- solve(crossprod(X))
  E <- Y - X %*% XtXi %*% crossprod(X, Y)
  ref <- vapply(1:7, function(j) {
    sqrt(diag(XtXi %*% crossprod(X * E[, j]) %*% XtXi))
  }, numeric(4))
  expect_equal(out$se, ref, tolerance = 1e-10)
  # HAC with zero bandwidth is HC0.
  expect_equal(sandwich_from_whitened_resid(X, Y, type = "hac", hac_lag = 0)$se, out$se)
})

test_that("HAC standard errors account for residual autocorrelation", {
  # Under-whitened noise (AR(1) left in) with a slow regressor: iid SEs are too
  # small, HAC SEs are not.
  set.seed(9124)
  n <- 400
  x <- rep(rep(c(0, 1), each = 20), length.out = n)
  X <- cbind(1, x)
  tt <- replicate(300, {
    e <- as.numeric(stats::arima.sim(list(ar = 0.5), n))
    s_iid <- sandwich_from_whitened_resid(X, matrix(e, ncol = 1))$se[2, 1]
    s_hac <- sandwich_from_whitened_resid(X, matrix(e, ncol = 1), type = "hac",
                                          hac_lag = 10)$se[2, 1]
    b <- qr.solve(X, e)[2]
    c(b / s_iid, b / s_hac)
  })
  expect_gt(mean(abs(tt[1, ]) > 1.96), 0.12)
  expect_lt(mean(abs(tt[2, ]) > 1.96), 0.09)
})

# --- compat ---------------------------------------------------------------------------

test_that("compat$update_plan keeps the previous plan's censor and exact_first", {
  set.seed(9125)
  R <- sim_arma_mat(200, 5, 0.5)
  cens <- c(20L, 21L, 100L)
  prev <- fit_noise(R, p = 1, censor = cens, exact_first = "none")
  upd <- compat$update_plan(prev, R)
  expect_identical(upd$censor, cens)
  expect_false(upd$exact_first)
})
