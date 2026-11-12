# Parcel pooling defaults and adaptive (empirical-Bayes) weights.

ar1_runs <- function(phi, nr = 150L) {
  c(stats::arima.sim(list(ar = phi), nr), stats::arima.sim(list(ar = phi), nr))
}
first <- function(l) unname(vapply(l, function(p) if (length(p)) p[1] else 0, 0))

test_that("large regions that really differ keep their own estimates", {
  skip_on_cran()
  set.seed(101)
  runs <- rep(1:2, each = 150L)
  region_phi <- c(0.2, 0.45, 0.7)
  region <- rep(1:3, each = 600L)
  Y <- sapply(region_phi[region], ar1_runs)
  X <- cbind(runs == 1, runs == 2)
  R <- qr.resid(qr(X), Y)
  # Default: no pooling across parcels without parcel_sets, for any p.
  for (pp in list("auto", 1L)) {
    plan <- fit_noise(R, runs = runs, design = X, p = pp, pooling = "parcel",
                      parcels = region)
    expect_equal(first(plan$phi_by_parcel), region_phi, tolerance = 0.03)
  }
  # Asking for pooling without a hierarchy pools toward all voxels, but only
  # as far as the evidence allows: these regions stay put.
  pooled <- fit_noise(R, runs = runs, design = X, pooling = "parcel",
                      parcels = region, multiscale = TRUE)
  expect_equal(first(pooled$phi_by_parcel), region_phi, tolerance = 0.03)
})

test_that("small noisy parcels are steadied by adaptive pooling", {
  skip_on_cran()
  set.seed(103)
  runs <- rep(1:2, each = 150L)
  fine <- rep(1:24, each = 4L)
  medium <- (fine - 1L) %/% 4L + 1L
  coarse <- (fine - 1L) %/% 8L + 1L
  truth <- c(0.2, 0.45, 0.7)[(1:24 - 1L) %/% 8L + 1L] + rnorm(24, sd = 0.03)
  Y <- sapply(truth[fine], ar1_runs)
  X <- cbind(runs == 1, runs == 2)
  R <- qr.resid(qr(X), Y)
  sets <- list(fine = fine, medium = medium, coarse = coarse)
  own <- first(fit_noise(R, runs = runs, design = X, p = 1L, pooling = "parcel",
                         parcels = fine)$phi_by_parcel)
  # With parcel_sets, pooling is on by default, also for an explicit p.
  ada <- first(fit_noise(R, runs = runs, design = X, p = 1L, pooling = "parcel",
                         parcels = fine, parcel_sets = sets)$phi_by_parcel)
  rms <- function(x) sqrt(mean((x - truth)^2))
  expect_lt(rms(ada), rms(own))
  # The legacy size-based weights remain available.
  size <- first(fit_noise(R, runs = runs, design = X, p = 1L, pooling = "parcel",
                          parcels = fine, parcel_sets = sets,
                          ms_weights = "size")$phi_by_parcel)
  expect_length(size, 24L)
  expect_false(isTRUE(all.equal(size, ada)))
})

test_that("empirical-Bayes shrinkage follows the evidence", {
  K <- matrix(c(0.2, 0.7), 1L)
  parent <- matrix(c(0.45, 0.45), 1L)
  # Precise, different children: kept.
  out <- fmriAR:::.ms_eb_shrink(K, matrix(1e-6, 1L, 2L), parent)
  expect_equal(out, K, tolerance = 1e-3)
  # Imprecise children that differ only by noise: pulled to the parent.
  out2 <- fmriAR:::.ms_eb_shrink(K, matrix(1, 1L, 2L), parent)
  expect_equal(out2, parent, tolerance = 1e-12)
})
