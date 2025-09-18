test_that("hr_arma_fit_cpp recovers ARMA(1,1) structure", {
  skip_if_not_installed("stats")
  set.seed(123)
  n <- 4000L
  phi <- 0.6
  theta <- 0.4
  e <- rnorm(n)
  y <- numeric(n)
  y[1] <- e[1]
  for (t in 2:n) {
    y[t] <- phi * y[t - 1] + e[t] + theta * e[t - 1]
  }

  est <- fmriAR:::hr_arma_fit_cpp(y, p = 1L, q = 1L, p_big = 12L, iter = 1L)
  phi_hat <- fmriAR:::enforce_stationary_ar(est$phi, bound = 0.99)
  theta_hat <- fmriAR:::enforce_invertible_ma(est$theta)

  expect_lt(abs(phi_hat[1] - phi), 0.10)
  expect_lt(abs(theta_hat[1] - theta), 0.15)

  Y <- cbind(y)
  X <- cbind(1)
  out <- arma_whiten_inplace(Y, X,
                             phi = phi_hat,
                             theta = theta_hat,
                             run_starts = 0L,
                             exact_first_ar1 = FALSE,
                             parallel = FALSE)
  innovations <- drop(out$Y)
  ac_vals <- stats::acf(innovations, plot = FALSE, lag.max = 12, demean = TRUE)$acf[-1L]
  ci <- 1.96 / sqrt(length(innovations))
  expect_lt(mean(abs(ac_vals[1:5])), 2 * ci)
})
