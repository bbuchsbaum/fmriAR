# ARMA estimates corrected for the design that produced the residuals.
#
# Hannan-Rissanen runs on the residual series and so inherits the bias that
# projecting out slow regressors puts into residual autocorrelation. With a
# design, fit_noise() refits the ARMA coefficients to the design-corrected
# autocovariance by minimum distance.

arma_case <- function(n_vox = 1000L, seed = 3) {
  set.seed(seed)
  nr <- 150L
  runs <- rep(1:2, each = nr)
  cosines <- sapply(1:3, function(j) cos(pi * j * (seq_len(nr) - 0.5) / nr))
  z <- matrix(0, nr, 3L)
  X <- cbind(runs == 1, runs == 2, rbind(cosines, z), rbind(z, cosines))
  # AR(1) phi = 0.8 plus white noise: ARMA(1,1) with theta about -0.343
  one <- function() as.numeric(stats::arima.sim(list(ar = 0.8), nr)) + rnorm(nr)
  E <- replicate(n_vox, c(one(), one()))
  list(runs = runs, X = X, E = E, R = qr.resid(qr(X), E))
}

test_that("the design correction removes ARMA projection bias", {
  skip_on_cran()
  d <- arma_case()
  hr <- fit_noise(d$R, runs = d$runs, method = "arma", p = 1L, q = 1L)
  dc <- fit_noise(d$R, runs = d$runs, method = "arma", p = 1L, q = 1L, design = d$X)
  expect_lt(abs(dc$phi[[1]] - 0.8), abs(hr$phi[[1]] - 0.8))
  expect_lt(abs(dc$phi[[1]] - 0.8), 0.02)
  expect_lt(abs(dc$theta[[1]] + 0.343), 0.04)
  # Invertible and stationary by construction.
  expect_true(all(Mod(polyroot(c(1, -dc$phi[[1]]))) > 1))
  expect_true(all(Mod(polyroot(c(1, dc$theta[[1]]))) > 1))
  # sigma2 on the voxel scale: AR(1) innovations plus white noise.
  expect_gt(dc$sigma2[[1]], 1)
})

test_that("per-run ARMA pooling is corrected too", {
  skip_on_cran()
  d <- arma_case()
  dc <- fit_noise(d$R, runs = d$runs, method = "arma", p = 1L, q = 1L,
                  pooling = "run", design = d$X)
  for (k in seq_along(dc$phi)) expect_lt(abs(dc$phi[[k]] - 0.8), 0.03)
})

test_that("whiten() and the Y/X path correct ARMA for their own design", {
  skip_on_cran()
  d <- arma_case(n_vox = 400L)
  ref <- fit_noise(d$R, runs = d$runs, method = "arma", p = 1L, q = 1L, design = d$X)
  yx <- fit_noise(Y = d$E, X = d$X, runs = d$runs, method = "arma", p = 1L, q = 1L)
  expect_equal(yx$phi, ref$phi)
  w <- whiten(d$X, d$E, runs = d$runs, method = "arma", p = 1L, q = 1L)
  expect_equal(w$Y, whiten_apply(ref, d$X, d$E, runs = d$runs)$Y)
})

test_that("the refinement recovers coefficients from an exact autocovariance", {
  g <- fmriAR:::arma_acvf_cpp(0.8, -0.343, 20L)
  ref <- fmriAR:::.arma_refine_acvf(2 * g, 0.6, -0.1)
  expect_equal(ref$phi, 0.8, tolerance = 1e-4)
  expect_equal(ref$theta, -0.343, tolerance = 1e-4)
  expect_equal(ref$sigma2, 2, tolerance = 1e-4)
  g2 <- fmriAR:::arma_acvf_cpp(c(0.4, 0.2), numeric(0), 20L)
  ref2 <- fmriAR:::.arma_refine_acvf(g2, c(0.3, 0.1), numeric(0))
  expect_equal(ref2$phi, c(0.4, 0.2), tolerance = 1e-4)
  expect_null(fmriAR:::.arma_refine_acvf(c(1, 0.5), 0.5, numeric(0)))
})

test_that("an ARMA order search that hits its ceiling warns", {
  skip_on_cran()
  d <- arma_case(n_vox = 800L)
  expect_warning(fit_noise(d$R, runs = d$runs, method = "arma", p = "auto",
                           q = "auto", design = d$X),
                 class = "fmriAR_order_at_ceiling")
})
