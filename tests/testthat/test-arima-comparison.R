
test_that("HR ARMA(1,1) aligns with stats::arima estimates", {
  skip_if_not_installed("fmriAR")
  set.seed(202401)
  sim <- simulate_arma11(n = 4000, phi = 0.6, theta = 0.4, sigma = 1.0, burnin = 500)
  y <- sim$y

  hr_fit <- fmriAR:::hr_arma(y, p = 1L, q = 1L, iter = 2L, step1 = "yw", enforce = TRUE)
  arima_fit <- stats::arima(y, order = c(1L, 0L, 1L), include.mean = FALSE, method = "ML")

  expect_equal(as.numeric(hr_fit$phi), unname(arima_fit$coef["ar1"]), tolerance = 0.05)
  expect_equal(as.numeric(hr_fit$theta), unname(arima_fit$coef["ma1"]), tolerance = 0.05)
})

test_that("Pure AR fits agree with stats::arima", {
  skip_if_not_installed("fmriAR")
  set.seed(202402)
  phi_true <- c(0.55, -0.25)
  y <- as.numeric(stats::arima.sim(list(ar = phi_true), n = 5000, sd = 1))

  hr_fit <- fmriAR:::hr_arma(y, p = 2L, q = 0L, iter = 1L, step1 = "yw", enforce = TRUE)
  arima_fit <- stats::arima(y, order = c(2L, 0L, 0L), include.mean = FALSE, method = "ML")

  expect_equal(as.numeric(hr_fit$phi), unname(arima_fit$coef[c("ar1", "ar2")]), tolerance = 0.04)
  expect_length(hr_fit$theta, 0L)
})
