
test_that("Estimated AR/MA polynomials are (likely) stationary/invertible", {
  skip_if_not_installed("fmriAR")
  # This test will pass after enforcement patch; before that it may fail for some seeds.
  set.seed(42)
  n <- 600
  phi_true <- c(1.60, -0.64)  # stable AR(2), roots at 1.25 (double)
  theta_true <- 0.5           # invertible MA(1)
  # simulate ARMA(2,1) by building from innovations
  e <- rnorm(n + 200L, 0, 1)
  y <- numeric(n + 200L)
  for (t in seq_len(n + 200L)) {
    y[t] <- (if (t > 1) (phi_true[1] * y[t - 1] + phi_true[2] * (if (t > 2) y[t - 2] else 0)) else 0) +
             e[t] + (if (t > 1) theta_true * e[t - 1] else 0)
  }
  y <- y[201:(200 + n)] - mean(y[201:(200 + n)])

  fit <- fmriAR:::hr_arma_fit_cpp(y, p = 2L, q = 1L, p_big = 0L, iter = 1L)
  expect_true(isTRUE(fit$ok))

  # Check root moduli (should exceed 1 by a margin)
  ar_mod <- ar_root_moduli(as.numeric(fit$phi))
  ma_mod <- ma_root_moduli(as.numeric(fit$theta))

  expect_gt(min(ar_mod), 1.0 + 1e-6)
  expect_gt(min(ma_mod), 1.0 + 1e-6)
})
