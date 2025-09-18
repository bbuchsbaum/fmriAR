
test_that("HR ARMA(1,1) recovers coefficients and whitening yields near-white residuals", {
  skip_if_not_installed("fmriAR")
  # Simulate ARMA(1,1)
  gen <- simulate_arma11(n = 8000, phi = 0.6, theta = 0.4, sigma = 1.0, burnin = 300, seed = 123)
  y <- gen$y
  res <- fmriAR:::hr_arma_fit_cpp(y, p = 1L, q = 1L, p_big = 0L, iter = 2L)
  expect_true(isTRUE(res$ok))
  expect_equal(as.numeric(res$phi), 0.6, tolerance = 0.06)
  expect_equal(as.numeric(res$theta), 0.4, tolerance = 0.07)
  expect_equal(as.numeric(res$sigma2), 1.0, tolerance = 0.15)

  # Whiten with estimated parameters; residuals should be close to white
  yw <- do_whiten_Y(y, phi = res$phi, theta = res$theta, run_starts0 = 0L, exact_first_ar1 = FALSE, parallel = FALSE)
  ac <- as.numeric(stats::acf(yw, plot = FALSE, lag.max = 5L)$acf)[-1L]
  expect_lt(max(abs(ac)), 0.05)
})
