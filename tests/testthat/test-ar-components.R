test_that("fit_noise recovers AR coefficients for single series", {
  set.seed(42)
  n <- 100

  # AR(1)
  true_phi1 <- 0.7
  x <- numeric(n)
  x[1] <- rnorm(1)
  for (i in 2:n) x[i] <- true_phi1 * x[i - 1] + rnorm(1)

  plan1 <- fit_noise(resid = matrix(x, ncol = 1), method = "ar", p = 1L)
  expect_equal(plan1$order[["p"]], 1L)
  expect_equal(length(plan1$phi[[1]]), 1L)
  expect_equal(plan1$phi[[1]][1], true_phi1, tolerance = 0.1)

  # AR(2)
  true_phi2 <- c(0.5, 0.3)
  y <- numeric(n)
  y[1:2] <- rnorm(2)
  for (i in 3:n) y[i] <- true_phi2[1] * y[i - 1] + true_phi2[2] * y[i - 2] + rnorm(1)

  plan2 <- fit_noise(resid = matrix(y, ncol = 1), method = "ar", p = 2L)
  expect_equal(plan2$order[["p"]], 2L)
  expect_equal(length(plan2$phi[[1]]), 2L)
  expect_equal(plan2$phi[[1]], true_phi2, tolerance = 0.3)
})

test_that("arma_whiten_inplace applies AR(1) transform correctly", {
  set.seed(43)
  n <- 50
  p <- 3
  phi <- 0.6

  X_orig <- cbind(1, matrix(rnorm(n * (p - 1)), n, p - 1))
  Y_orig <- matrix(rnorm(n * 2), n, 2)

  X_exact <- X_orig + 0
  Y_exact <- Y_orig + 0
  res_exact <- arma_whiten_inplace(Y_exact, X_exact, phi = c(phi), theta = numeric(),
                                   run_starts = 0L, exact_first_ar1 = TRUE,
                                   parallel = FALSE)

  expect_equal(dim(res_exact$X), dim(X_orig))
  expect_equal(dim(res_exact$Y), dim(Y_orig))
  expect_false(isTRUE(all.equal(res_exact$X, X_orig)))
  expect_false(isTRUE(all.equal(res_exact$Y, Y_orig)))

  scale <- sqrt(1 - phi^2)
  expect_equal(as.numeric(res_exact$X[1, ]), as.numeric(X_orig[1, ] * scale), tolerance = 1e-10)
  expect_equal(as.numeric(res_exact$Y[1, ]), as.numeric(Y_orig[1, ] * scale), tolerance = 1e-10)
  expect_equal(as.numeric(res_exact$X[2, ]), as.numeric(X_orig[2, ] - phi * X_orig[1, ]), tolerance = 1e-10)

  X_plain <- X_orig + 0
  Y_plain <- Y_orig + 0
  res_plain <- arma_whiten_inplace(Y_plain, X_plain, phi = c(phi), theta = numeric(),
                                    run_starts = 0L, exact_first_ar1 = FALSE,
                                    parallel = FALSE)
  expect_equal(dim(res_plain$X), dim(X_orig))
  expect_equal(dim(res_plain$Y), dim(Y_orig))
  expect_equal(as.numeric(res_plain$X[1, ]), as.numeric(X_orig[1, ]), tolerance = 1e-10)
  expect_equal(as.numeric(res_plain$Y[1, ]), as.numeric(Y_orig[1, ]), tolerance = 1e-10)
})

test_that("arma_whiten_inplace handles AR(2) filters", {
  set.seed(44)
  n <- 100
  p <- 2
  phi <- c(0.5, 0.2)

  X_orig <- matrix(rnorm(n * p), n, p)
  Y_orig <- matrix(rnorm(n * 3), n, 3)

  X_run <- X_orig + 0
  Y_run <- Y_orig + 0
  res <- arma_whiten_inplace(Y_run, X_run, phi = phi, theta = numeric(),
                             run_starts = 0L, exact_first_ar1 = TRUE,
                             parallel = FALSE)

  expect_equal(dim(res$X), dim(X_orig))
  expect_equal(dim(res$Y), dim(Y_orig))
  expect_false(isTRUE(all.equal(res$X, X_orig)))
  expect_false(isTRUE(all.equal(res$Y, Y_orig)))

  expected_row3 <- X_orig[3, ] - phi[1] * X_orig[2, ] - phi[2] * X_orig[1, ]
  expect_equal(as.numeric(res$X[3, ]), as.numeric(expected_row3), tolerance = 1e-10)
})

test_that("whitening reduces residual autocorrelation", {
  set.seed(45)
  n <- 200
  p <- 5
  phi_true <- 0.8

  X <- cbind(1, matrix(rnorm(n * (p - 1)), n, p - 1))
  beta_true <- rnorm(p)

  errors <- numeric(n)
  errors[1] <- rnorm(1)
  for (i in 2:n) errors[i] <- phi_true * errors[i - 1] + rnorm(1)

  Y_vec <- as.numeric(X %*% beta_true + errors)

  beta_ols <- qr.solve(X, Y_vec)
  resid_ols <- Y_vec - X %*% beta_ols

  plan_hat <- fit_noise(resid = matrix(resid_ols, ncol = 1), method = "ar", p = 1L)
  phi_hat <- plan_hat$phi[[1]][1]
  expect_gt(phi_hat, 0.5)

  X_in <- X
  Y_in <- matrix(Y_vec, ncol = 1)
  whitened <- arma_whiten_inplace(Y = Y_in, X = X_in,
                                  phi = c(phi_hat), theta = numeric(),
                                  run_starts = 0L, exact_first_ar1 = TRUE,
                                  parallel = FALSE)

  beta_gls <- qr.solve(whitened$X, drop(whitened$Y))
  resid_gls <- drop(whitened$Y) - whitened$X %*% beta_gls

  plan_resid <- fit_noise(resid = matrix(resid_gls, ncol = 1), method = "ar", p = 1L)
  phi_resid <- if (plan_resid$order[["p"]] > 0L) plan_resid$phi[[1]][1] else 0
  expect_lt(abs(phi_resid), abs(phi_hat))
  expect_lt(abs(phi_resid), 0.2)
})

test_that("fit_noise and arma_whiten_inplace handle edge cases", {
  set.seed(46)
  expect_error(fit_noise(resid = matrix(numeric(0), nrow = 0, ncol = 1),
                         method = "ar", p = 1L))

  const_series <- matrix(rep(5, 100), ncol = 1)
  plan_const <- fit_noise(resid = const_series, method = "ar", p = "auto", p_max = 5L)
  expect_equal(plan_const$order[["p"]], 0L)

  X <- matrix(1:12, 4, 3)
  Y <- matrix(1:8, 4, 2)
  phi <- 0.5

  res <- arma_whiten_inplace(Y, X, phi = c(phi), theta = numeric(),
                              run_starts = 0L, exact_first_ar1 = TRUE,
                              parallel = FALSE)
  expect_equal(dim(res$X), dim(X))
  expect_equal(dim(res$Y), dim(Y))
})
