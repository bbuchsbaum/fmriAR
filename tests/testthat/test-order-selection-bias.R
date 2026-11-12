# Order selection and the bias that fitting puts into residual autocorrelation.
#
# With thousands of voxels pooled, BIC has enough evidence to treat the small,
# systematic negative bias that a fitted model leaves at every lag as real
# structure, and chooses too high an order. Two sources: the per-run centring
# the estimator applies itself (known exactly, corrected by default), and the
# user's design, above all its drift terms (corrected by `design = X`, which
# whiten() and the Y/X path now supply on their own).

sim_runs <- function(n_vox, ar, nr = 150L, seed = 1) {
  set.seed(seed)
  replicate(n_vox, c(stats::arima.sim(list(ar = ar), nr),
                     stats::arima.sim(list(ar = ar), nr)))
}
drift_design <- function(k, nr = 150L) {
  runs <- rep(1:2, each = nr)
  cosines <- sapply(seq_len(k), function(j) cos(pi * j * (seq_len(nr) - 0.5) / nr))
  z <- matrix(0, nr, k)
  cbind(runs == 1, runs == 2, rbind(cosines, z), rbind(z, cosines))
}

test_that("the per-run centring is corrected by default", {
  skip_on_cran()
  runs <- rep(1:2, each = 150L)
  E <- sim_runs(200L, 0.5, seed = 3)
  R <- E - rep(colMeans(E), each = nrow(E))
  dflt <- fit_noise(R, runs = runs, p = 1L)
  plain <- fit_noise(R, runs = runs, p = 1L, design = FALSE)
  # Centring alone pulls lag-1 autocorrelation down by roughly 1/150.
  expect_lt(abs(dflt$phi[[1]] - 0.5), abs(plain$phi[[1]] - 0.5))
  expect_lt(abs(dflt$phi[[1]] - 0.5), 0.01)
})

test_that("without drift terms the default no longer over-selects the order", {
  skip_on_cran()
  runs <- rep(1:2, each = 150L)
  E <- sim_runs(2000L, c(0.4, 0.2), seed = 11)
  R <- qr.resid(qr(cbind(runs == 1, runs == 2)), E)
  expect_identical(fit_noise(R, runs = runs)$order[["p"]], 2L)
})

test_that("drift residuals without design warn at the order ceiling", {
  skip_on_cran()
  runs <- rep(1:2, each = 150L)
  X <- drift_design(3L)
  E <- sim_runs(2000L, c(0.4, 0.2), seed = 11)
  R <- qr.resid(qr(X), E)
  expect_warning(plan <- fit_noise(R, runs = runs), class = "fmriAR_order_at_ceiling")
  expect_identical(plan$order[["p"]], 6L)
  expect_silent(plan_d <- fit_noise(R, runs = runs, design = X))
  expect_identical(plan_d$order[["p"]], 2L)
})

test_that("whiten() and the Y/X path correct for their own design", {
  skip_on_cran()
  runs <- rep(1:2, each = 150L)
  X <- drift_design(3L)
  E <- sim_runs(2000L, c(0.4, 0.2), seed = 11)
  expect_silent(yx <- fit_noise(Y = E, X = X, runs = runs))
  expect_identical(yx$order[["p"]], 2L)
  ref <- fit_noise(qr.resid(qr(X), E), runs = runs, design = X)
  expect_equal(yx$phi, ref$phi)

  w <- whiten(X, E, runs = runs)
  w_ref <- whiten_apply(ref, X, E, runs = runs)
  expect_equal(w$Y, w_ref$Y)
  # design = FALSE still turns every correction off.
  w_off <- whiten(X, E, runs = runs, design = FALSE, p = 2L)
  plan_off <- fit_noise(qr.resid(qr(X), E), runs = runs, design = FALSE, p = 2L)
  expect_equal(w_off$Y, whiten_apply(plan_off, X, E, runs = runs)$Y)
})

test_that("design = FALSE reproduces the plain estimator exactly", {
  set.seed(5)
  R <- matrix(rnorm(200 * 6), 200, 6)
  a <- noise_acvf(R, max_lag = 4L, design = FALSE)
  expect_identical(a$correction, "none")
  # Plain estimator: centred, lag products over pairs.
  Rc <- sweep(R, 2, colMeans(R))
  manual <- vapply(0:4, function(k) {
    mean(colSums(Rc[(1 + k):200, , drop = FALSE] * Rc[1:(200 - k), , drop = FALSE]) / (200 - k))
  }, 0)
  expect_equal(a$acvf[[1]], manual, tolerance = 1e-10)
})
