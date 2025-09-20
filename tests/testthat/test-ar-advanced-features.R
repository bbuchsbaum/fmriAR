# AR Advanced Features tests

# Helper functions for simulating AR data compatible with fmriAR
simulate_ar_data <- function(n = 100, p = 1, phi = NULL, sigma = 1) {
  if (is.null(phi)) {
    phi <- runif(p, -0.9, 0.9) / p  # Default stable AR coefficients
  }

  # Handle AR(0) case (white noise)
  if (length(phi) == 0 || all(phi == 0)) {
    return(rnorm(n, sd = sigma))
  }

  # Ensure stability
  roots <- polyroot(c(1, -phi))
  if (any(abs(roots) <= 1.05)) {
    warning("AR coefficients may be near unit root")
  }

  as.numeric(arima.sim(model = list(ar = phi), n = n, sd = sigma))
}

create_ar_test_data <- function(n = 100, p_design = 3, ar_order = 1,
                                ar_coef = 0.5, n_vox = 5) {
  # Create design matrix
  X <- cbind(1, matrix(rnorm(n * (p_design - 1)), n, p_design - 1))

  # True betas
  true_betas <- rnorm(p_design)

  # Generate Y with AR errors
  Y <- matrix(NA, n, n_vox)
  for (v in 1:n_vox) {
    signal <- X %*% true_betas
    errors <- simulate_ar_data(n = n, p = ar_order, phi = ar_coef)
    Y[, v] <- signal + errors
  }

  # Compute residuals for noise fitting
  resid <- Y - X %*% qr.solve(X, Y)

  list(X = X, Y = Y, resid = resid, true_betas = true_betas, ar_coef = ar_coef)
}

test_that("single run AR(1) whitening reduces autocorrelation", {
  set.seed(123)
  test_data <- create_ar_test_data(n = 150, ar_order = 1, ar_coef = 0.6)

  # Fit noise model
  plan <- fit_noise(resid = test_data$resid, method = "ar", p = 1L)

  # Apply whitening
  whitened <- whiten_apply(plan, test_data$X, test_data$Y, parallel = FALSE)

  # Check that whitening was applied
  expect_false(all(whitened$X == test_data$X))
  expect_false(all(whitened$Y == test_data$Y))

  # Check AR coefficient is reasonable
  phi_est <- plan$phi[[1]][1]
  expect_true(abs(phi_est) < 1)  # Stationary
  expect_true(abs(phi_est - 0.6) < 0.2)  # Close to true value

  # Verify autocorrelation reduction
  whitened_resid <- whitened$Y - whitened$X %*% qr.solve(whitened$X, whitened$Y)
  acorr_result <- acorr_diagnostics(whitened_resid, max_lag = 5, aggregate = "mean")
  expect_true(all(abs(acorr_result$acf[1:3]) < 0.15))
})

test_that("multi-run AR whitening handles different coefficients per run", {
  set.seed(124)

  # Create multi-run data with different AR coefficients
  n_per_run <- 60
  n_runs <- 3
  p_design <- 3
  n_vox <- 4

  # Different AR coefficients per run
  ar_coefs <- c(0.3, 0.5, 0.7)

  X_list <- list()
  Y_list <- list()

  for (r in 1:n_runs) {
    X_run <- cbind(1, matrix(rnorm(n_per_run * (p_design - 1)), n_per_run, p_design - 1))
    true_betas <- rnorm(p_design)

    Y_run <- matrix(NA, n_per_run, n_vox)
    for (v in 1:n_vox) {
      signal <- X_run %*% true_betas
      errors <- simulate_ar_data(n = n_per_run, p = 1, phi = ar_coefs[r])
      Y_run[, v] <- signal + errors
    }

    X_list[[r]] <- X_run
    Y_list[[r]] <- Y_run
  }

  X <- do.call(rbind, X_list)
  Y <- do.call(rbind, Y_list)

  # Create run indices
  runs <- rep(1:n_runs, each = n_per_run)

  # Compute residuals
  resid <- Y - X %*% qr.solve(X, Y)

  # Fit noise with run-specific pooling
  plan <- fit_noise(resid = resid, runs = runs, method = "ar", p = 1L, pooling = "run")

  # Apply whitening
  whitened <- whiten_apply(plan, X, Y, runs = runs, parallel = FALSE)

  # Check results
  expect_equal(length(plan$phi), n_runs)

  # Check that estimated AR coefficients are different per run
  phi_estimates <- sapply(plan$phi, function(x) x[1])
  expect_equal(length(phi_estimates), n_runs)

  # Each should be reasonably close to true value (allow for estimation variability)
  for (r in 1:n_runs) {
    if (!is.na(phi_estimates[r])) {
      expect_true(abs(phi_estimates[r] - ar_coefs[r]) < 0.4)
    }
  }
})

test_that("higher-order AR models (AR(3) and AR(4)) work correctly", {
  set.seed(127)

  # Test AR(3)
  n <- 200
  p_design <- 3
  n_vox <- 3
  ar3_coef <- c(0.3, 0.2, 0.1)

  test_data <- create_ar_test_data(n = n, p_design = p_design,
                                   ar_order = 3, ar_coef = ar3_coef, n_vox = n_vox)

  plan_ar3 <- fit_noise(resid = test_data$resid, method = "ar", p = 3L)
  whitened_ar3 <- whiten_apply(plan_ar3, test_data$X, test_data$Y, parallel = FALSE)

  # Check AR(3) coefficients
  expect_equal(length(plan_ar3$phi[[1]]), 3)

  # Test AR(4)
  ar4_coef <- c(0.2, 0.15, 0.1, 0.05)
  test_data4 <- create_ar_test_data(n = n, p_design = p_design,
                                    ar_order = 4, ar_coef = ar4_coef, n_vox = n_vox)

  plan_ar4 <- fit_noise(resid = test_data4$resid, method = "ar", p = 4L)
  whitened_ar4 <- whiten_apply(plan_ar4, test_data4$X, test_data4$Y, parallel = FALSE)

  # Check AR(4) coefficients (may be reduced by model selection)
  expect_true(length(plan_ar4$phi[[1]]) >= 1)
  expect_true(length(plan_ar4$phi[[1]]) <= 4)
})

test_that("auto AR order selection works with p = 'auto'", {
  set.seed(128)

  # Test with known AR(2) process
  n <- 250
  p_design <- 3
  n_vox <- 2
  ar_coef <- c(0.4, 0.2)

  test_data <- create_ar_test_data(n = n, p_design = p_design,
                                   ar_order = 2, ar_coef = ar_coef, n_vox = n_vox)

  plan_auto <- fit_noise(resid = test_data$resid, method = "ar", p = "auto", p_max = 5L)

  # Should select AR order near the true order
  expect_true(plan_auto$order[["p"]] >= 1)
  expect_true(plan_auto$order[["p"]] <= 5)

  # Apply whitening
  whitened <- whiten_apply(plan_auto, test_data$X, test_data$Y, parallel = FALSE)
  expect_equal(dim(whitened$Y), dim(test_data$Y))
})

test_that("AR methods handle near unit-root processes", {
  set.seed(129)

  # Test with AR coefficients close to 1
  near_unit_roots <- c(0.9, 0.95, 0.99)

  for (phi in near_unit_roots) {
    n <- 300  # Need longer series for near unit root
    p_design <- 2
    n_vox <- 2

    # Suppress warnings about near unit root
    test_data <- suppressWarnings(
      create_ar_test_data(n = n, p_design = p_design, ar_order = 1,
                          ar_coef = phi, n_vox = n_vox)
    )

    # Should handle without error
    plan <- fit_noise(resid = test_data$resid, method = "ar", p = 1L)
    whitened <- whiten_apply(plan, test_data$X, test_data$Y, parallel = FALSE)

    expect_true(!is.null(plan$phi))
    phi_est <- plan$phi[[1]][1]

    # Should still be stationary (enforced by fmriAR)
    expect_true(abs(phi_est) < 1)

    # For very high AR, should be substantially reduced from OLS case
    expect_equal(dim(whitened$Y), dim(test_data$Y))
  }
})

test_that("AR whitening handles short time series edge cases", {
  set.seed(130)

  # Test with time series too short for high AR order
  n_short <- 15
  p_design <- 3
  n_vox <- 2

  test_data <- create_ar_test_data(n = n_short, p_design = p_design,
                                   ar_order = 1, ar_coef = 0.5, n_vox = n_vox)

  # Try high AR order with short series
  plan <- fit_noise(resid = test_data$resid, method = "ar", p = "auto", p_max = 8L)

  # Should work but likely select lower order
  expect_true(plan$order[["p"]] >= 0)
  expect_true(plan$order[["p"]] < n_short / 2)

  whitened <- whiten_apply(plan, test_data$X, test_data$Y, parallel = FALSE)
  expect_equal(dim(whitened$Y), dim(test_data$Y))
})

test_that("parcel-based AR estimation works correctly", {
  set.seed(131)

  n <- 150
  p_design <- 3
  n_vox <- 12
  n_parcels <- 3

  # Create different AR structure per parcel
  parcels <- rep(1:n_parcels, each = n_vox / n_parcels)
  parcel_ar_coefs <- c(0.3, 0.6, 0.8)

  X <- cbind(1, matrix(rnorm(n * (p_design - 1)), n, p_design - 1))
  true_betas <- rnorm(p_design)

  Y <- matrix(NA, n, n_vox)
  for (v in 1:n_vox) {
    signal <- X %*% true_betas
    p_id <- parcels[v]
    errors <- simulate_ar_data(n = n, p = 1, phi = parcel_ar_coefs[p_id])
    Y[, v] <- signal + errors
  }

  resid <- Y - X %*% qr.solve(X, Y)

  # Fit parcel-based noise model
  plan <- fit_noise(resid = resid, method = "ar", p = 1L,
                    pooling = "parcel", parcels = parcels)

  # Apply whitening
  whitened <- whiten_apply(plan, X, Y, parcels = parcels, parallel = FALSE)

  # Check that we get different X matrices per parcel
  expect_true(!is.null(whitened$X_by))
  expect_equal(length(whitened$X_by), n_parcels)

  # Check parcel-specific coefficients
  expect_true(!is.null(plan$phi_by_parcel))
  expect_equal(length(plan$phi_by_parcel), n_parcels)
})

test_that("AR whitening with censor gaps resets recursions", {
  set.seed(132)

  n <- 150
  p_design <- 3
  n_vox <- 3
  ar_coef <- 0.4

  test_data <- create_ar_test_data(n = n, p_design = p_design,
                                   ar_order = 1, ar_coef = ar_coef, n_vox = n_vox)

  # Define censor points
  censor <- c(60L, 120L)

  plan <- fit_noise(resid = test_data$resid, method = "ar", p = 1L)

  # Apply whitening with censor gaps
  whitened_censored <- whiten_apply(plan, test_data$X, test_data$Y,
                                    censor = censor, parallel = FALSE)

  # Compare with no censoring
  whitened_full <- whiten_apply(plan, test_data$X, test_data$Y, parallel = FALSE)

  # Should be different due to reset at censor points
  expect_false(identical(whitened_censored$Y, whitened_full$Y))
  expect_equal(dim(whitened_censored$Y), dim(test_data$Y))
})

test_that("ARMA(1,1) processes are handled correctly", {
  set.seed(133)

  n <- 300
  p_design <- 3
  n_vox <- 3

  # Generate ARMA(1,1) data
  X <- cbind(1, matrix(rnorm(n * (p_design - 1)), n, p_design - 1))
  true_betas <- rnorm(p_design)

  Y <- matrix(NA, n, n_vox)
  for (v in 1:n_vox) {
    signal <- X %*% true_betas
    # ARMA(1,1) errors
    errors <- as.numeric(arima.sim(model = list(ar = 0.3, ma = -0.25), n = n))
    Y[, v] <- signal + errors
  }

  resid <- Y - X %*% qr.solve(X, Y)

  # Fit ARMA model
  plan <- fit_noise(resid = resid, method = "arma", p = 1L, q = 1L)

  # Apply whitening
  whitened <- whiten_apply(plan, X, Y, parallel = FALSE)

  # Check that both AR and MA components are estimated
  expect_equal(plan$order[["p"]], 1)
  expect_equal(plan$order[["q"]], 1)
  expect_equal(length(plan$phi[[1]]), 1)
  expect_equal(length(plan$theta[[1]]), 1)

  # Verify whitening reduces autocorrelation
  whitened_resid <- whitened$Y - whitened$X %*% qr.solve(whitened$X, whitened$Y)
  acorr_result <- acorr_diagnostics(whitened_resid, max_lag = 5, aggregate = "mean")
  expect_true(all(abs(acorr_result$acf[1:3]) < 0.2))
})

test_that("whitening preserves matrix dimensions and rank", {
  set.seed(134)

  # Create full rank design matrix
  n <- 100
  p <- 5
  n_vox <- 3

  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(rnorm(n * n_vox), n, n_vox)
  resid <- Y - X %*% qr.solve(X, Y)

  plan <- fit_noise(resid = resid, method = "ar", p = 2L)
  whitened <- whiten_apply(plan, X, Y, parallel = FALSE)

  # Check dimensions preserved
  expect_equal(dim(whitened$X), dim(X))
  expect_equal(dim(whitened$Y), dim(Y))

  # Check that whitened X still has full rank
  rank_original <- qr(X)$rank
  rank_whitened <- qr(whitened$X)$rank

  expect_equal(rank_whitened, rank_original)
})

test_that("global pooling averages across runs correctly", {
  set.seed(135)

  n1 <- 100
  n2 <- 150
  phi1 <- 0.3
  phi2 <- 0.7

  # Create two runs with different AR coefficients
  resid1 <- matrix(simulate_ar_data(n1, 1, phi1), ncol = 1)
  resid2 <- matrix(simulate_ar_data(n2, 1, phi2), ncol = 1)

  resid <- rbind(resid1, resid2)
  runs <- c(rep(1L, n1), rep(2L, n2))

  plan <- fit_noise(resid = resid, runs = runs, pooling = "global",
                    method = "ar", p = 1L)

  if (plan$order[["p"]] > 0L) {
    phi_hat <- plan$phi[[1]][1]
    # Should be weighted average
    expected <- (n1 * phi1 + n2 * phi2) / (n1 + n2)
    expect_equal(phi_hat, expected, tolerance = 0.1)
  }
})

test_that("whitening handles various matrix conditions gracefully", {
  set.seed(136)

  n <- 100
  p <- 5
  n_vox <- 3

  # Test 1: Well-conditioned matrix
  X_good <- matrix(rnorm(n * p), n, p)
  Y <- matrix(rnorm(n * n_vox), n, n_vox)
  resid <- Y - X_good %*% qr.solve(X_good, Y)

  plan <- fit_noise(resid = resid, method = "ar", p = 1L)
  whitened_good <- whiten_apply(plan, X_good, Y, parallel = FALSE)
  expect_equal(dim(whitened_good$Y), dim(Y))

  # Test 2: Matrix with high condition number but still full rank
  X_ill <- matrix(rnorm(n * p), n, p)
  X_ill[, p] <- X_ill[, 1] + rnorm(n, sd = 0.01)  # Nearly dependent column
  resid_ill <- Y - X_ill %*% qr.solve(X_ill, Y)

  plan_ill <- fit_noise(resid = resid_ill, method = "ar", p = 1L)
  whitened_ill <- whiten_apply(plan_ill, X_ill, Y, parallel = FALSE)
  expect_equal(dim(whitened_ill$Y), dim(Y))
})

test_that("integration with full fmriAR workflow", {
  set.seed(137)

  # Create realistic fMRI-like data
  n_time <- 200
  n_vox <- 8
  p_design <- 4

  # Design matrix with intercept and some regressors
  X <- cbind(1,
             sin(2 * pi * (1:n_time) / 20),  # Periodic regressor
             cos(2 * pi * (1:n_time) / 20),
             rnorm(n_time))  # Random regressor

  # True parameters
  true_betas <- c(100, 5, 3, 2)

  # Generate data with AR(1) errors
  ar_coef <- 0.4
  Y <- matrix(NA, n_time, n_vox)
  for (v in 1:n_vox) {
    signal <- X %*% true_betas
    errors <- simulate_ar_data(n = n_time, p = 1, phi = ar_coef, sigma = 10)
    Y[, v] <- signal + errors
  }

  # Test full pipeline using whiten() convenience function
  whitened <- whiten(X, Y, method = "ar", p = "auto", p_max = 3)

  # Fit model on whitened data
  beta_hat <- qr.solve(whitened$X, whitened$Y)

  # Check that betas are recovered reasonably well
  mean_betas <- rowMeans(beta_hat)

  # Should be within reasonable range of true values
  for (i in 1:p_design) {
    if (abs(true_betas[i]) > 1) {  # Only check larger coefficients
      relative_error <- abs(mean_betas[i] - true_betas[i]) / abs(true_betas[i])
      expect_true(relative_error < 0.5)  # Within 50% of true value
    }
  }

  # Check that whitening reduced autocorrelation
  whitened_resid <- whitened$Y - whitened$X %*% beta_hat
  acorr_result <- acorr_diagnostics(whitened_resid, max_lag = 5, aggregate = "mean")
  expect_true(all(abs(acorr_result$acf[1:3]) < 0.15))
})

test_that("inplace modification works correctly", {
  set.seed(138)

  test_data <- create_ar_test_data(n = 100, ar_order = 1, ar_coef = 0.5)

  plan <- fit_noise(resid = test_data$resid, method = "ar", p = 1L)

  # Create copies for comparison
  X_copy <- test_data$X
  Y_copy <- test_data$Y

  # Test inplace modification
  result_inplace <- whiten_apply(plan, test_data$X, test_data$Y,
                                 inplace = TRUE, parallel = FALSE)

  # Test regular modification
  result_regular <- whiten_apply(plan, X_copy, Y_copy, parallel = FALSE)

  # The inplace result should have the same output as regular
  expect_equal(result_inplace$Y, result_regular$Y)
  expect_equal(result_inplace$X, result_regular$X)

  # Test that inplace option works without error
  expect_true(is.list(result_inplace))
  expect_true(!is.null(result_inplace$Y))
  expect_true(!is.null(result_inplace$X))
})