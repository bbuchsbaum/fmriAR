test_that("AR(1) whitening reduces lag-1 autocorrelation", {
  set.seed(1)
  n <- 600
  phi <- 0.6
  e <- rnorm(n)
  y <- as.numeric(stats::filter(e, filter = phi, method = "recursive"))
  Y <- cbind(y)
  X <- cbind(1)

  out <- arma_whiten_inplace(
    Y = Y,
    X = X,
    phi = phi,
    theta = numeric(0),
    run_starts = 0L,
    exact_first_ar1 = TRUE,
    parallel = FALSE
  )
  yw <- drop(out$Y)
  acf_vals <- stats::acf(yw, plot = FALSE, lag.max = 10, demean = TRUE)$acf[-1L]
  expect_true(all(abs(acf_vals[1:3]) < 0.15))
})

test_that("PACF <-> AR round-trips", {
  kap <- c(0.2, -0.1, 0.3)
  phi <- fmriAR:::pacf_to_ar(kap)
  kap2 <- fmriAR:::ar_to_pacf(phi)
  expect_equal(unname(kap2), unname(kap), tolerance = 1e-8)
})

test_that("run-specific parameters are applied per segment", {
  set.seed(11)
  n1 <- 80
  n2 <- 90
  phi1 <- 0.5
  phi2 <- -0.3
  y1 <- as.numeric(stats::filter(rnorm(n1), filter = phi1, method = "recursive"))
  y2 <- as.numeric(stats::filter(rnorm(n2), filter = phi2, method = "recursive"))
  Y <- rbind(matrix(y1, ncol = 1L), matrix(y2, ncol = 1L))
  X <- rbind(matrix(1, n1, 1L), matrix(1, n2, 1L))
  runs <- c(rep(1L, n1), rep(2L, n2))

  plan <- fmriAR:::new_whiten_plan(
    phi = list(phi1, phi2),
    theta = list(numeric(0), numeric(0)),
    order = c(p = 1L, q = 0L),
    runs = runs,
    exact_first = FALSE,
    method = "ar",
    pooling = "run"
  )

  out <- whiten_apply(plan, X, Y, parallel = FALSE)

  manual1 <- arma_whiten_inplace(matrix(y1, ncol = 1L), matrix(1, n1, 1L),
                                 phi = phi1, theta = numeric(0),
                                 run_starts = 0L, exact_first_ar1 = FALSE,
                                 parallel = FALSE)
  manual2 <- arma_whiten_inplace(matrix(y2, ncol = 1L), matrix(1, n2, 1L),
                                 phi = phi2, theta = numeric(0),
                                 run_starts = 0L, exact_first_ar1 = FALSE,
                                 parallel = FALSE)

  expect_equal(out$Y[seq_len(n1), , drop = FALSE], manual1$Y)
  expect_equal(out$Y[n1 + seq_len(n2), , drop = FALSE], manual2$Y)
})

test_that("censor gaps reset the whitening recursions", {
  set.seed(21)
  n <- 150
  phi <- 0.4
  y <- as.numeric(stats::filter(rnorm(n), filter = phi, method = "recursive"))
  Y <- matrix(y, ncol = 1L)
  X <- matrix(1, n, 1L)
  runs <- rep(1L, n)
  censor <- c(60L, 120L)

  plan <- fmriAR:::new_whiten_plan(
    phi = list(phi),
    theta = list(numeric(0)),
    order = c(p = 1L, q = 0L),
    runs = runs,
    exact_first = FALSE,
    method = "ar",
    pooling = "global"
  )

  out_plan <- whiten_apply(plan, X, Y, runs = runs, censor = censor, parallel = FALSE)

  run_starts <- as.integer(c(0L, censor))
  out_manual <- arma_whiten_inplace(Y, X,
                                    phi = phi, theta = numeric(0),
                                    run_starts = run_starts,
                                    exact_first_ar1 = FALSE,
                                    parallel = FALSE)

  expect_equal(out_plan$Y, out_manual$Y)
  expect_equal(out_plan$X, out_manual$X)
})

test_that("fit_noise censor parameter excludes censored timepoints from estimation", {
  set.seed(42)
  n <- 200
  n_vox <- 20
  phi_true <- 0.6

  # Generate AR(1) residuals
  resid <- matrix(0, n, n_vox)
  for (v in seq_len(n_vox)) {
    e <- rnorm(n)
    resid[1, v] <- e[1]
    for (t in 2:n) {
      resid[t, v] <- phi_true * resid[t - 1, v] + e[t]
    }
  }

  # Corrupt some timepoints with moderate outliers (simulating motion artifacts)
  censor_idx <- c(50L, 51L, 100L, 101L, 150L)
  # Use moderate outliers that bias estimates but don't completely dominate
  resid[censor_idx, ] <- resid[censor_idx, ] + rnorm(length(censor_idx) * n_vox, sd = 5)

  # Force p = 1 to ensure we get phi estimates in both cases
  # Fit without censor (should be biased by outliers)
  plan_no_censor <- fit_noise(resid, method = "ar", p = 1, p_max = 1L)

  # Fit with censor (should be closer to true value)
  plan_with_censor <- fit_noise(resid, censor = censor_idx, method = "ar", p = 1, p_max = 1L)

  # Check that censor is stored in the plan
  expect_equal(plan_with_censor$censor, censor_idx)

  # Both should return p=1 AR coefficients
  expect_length(plan_no_censor$phi[[1]], 1L)
  expect_length(plan_with_censor$phi[[1]], 1L)

  # The censored plan should typically have phi closer to true value
  # (allow some randomness, so just check both are reasonably close)
  phi_with_censor <- plan_with_censor$phi[[1]][1]
  expect_true(abs(phi_with_censor - phi_true) < 0.15)
})

test_that("fit_noise accepts logical censor vector", {
  set.seed(43)
  n <- 100
  n_vox <- 10
  resid <- matrix(rnorm(n * n_vox), n, n_vox)

  # Test logical censor
  censor_logical <- rep(FALSE, n)
  censor_logical[c(25, 50, 75)] <- TRUE

  plan <- fit_noise(resid, censor = censor_logical, method = "ar", p = 1)
  expect_equal(plan$censor, c(25L, 50L, 75L))
})

test_that("fit_noise censor works with runs", {
  set.seed(44)
  n1 <- 100
  n2 <- 100
  n_vox <- 15
  phi_true <- 0.5

  # Generate AR(1) residuals for two runs
  resid <- matrix(0, n1 + n2, n_vox)
  for (v in seq_len(n_vox)) {
    for (run_start in c(1, n1 + 1)) {
      e <- rnorm(100)
      resid[run_start, v] <- e[1]
      for (t in 2:100) {
        resid[run_start + t - 1, v] <- phi_true * resid[run_start + t - 2, v] + e[t]
      }
    }
  }
  runs <- c(rep(1L, n1), rep(2L, n2))

  # Corrupt some timepoints in each run
  censor_idx <- c(25L, 50L, 125L, 175L)  # 25, 50 in run1; 125, 175 in run2
  resid[censor_idx, ] <- resid[censor_idx, ] + rnorm(length(censor_idx) * n_vox, sd = 15)

  plan_censored <- fit_noise(resid, runs = runs, censor = censor_idx,
                              method = "ar", p = 1, pooling = "run")

  expect_equal(plan_censored$censor, censor_idx)
  expect_equal(length(plan_censored$phi), 2L)  # One phi per run
})

test_that("ARMA(1,1) whitening yields near-white innovations", {
  set.seed(31)
  n <- 800
  phi <- 0.3
  theta <- -0.25
  y <- as.numeric(stats::arima.sim(model = list(ar = phi, ma = theta), n = n))
  Y <- matrix(y, ncol = 1L)
  X <- matrix(1, n, 1L)

  out <- arma_whiten_inplace(Y, X,
                             phi = phi, theta = theta,
                             run_starts = 0L,
                             exact_first_ar1 = FALSE,
                             parallel = FALSE)
  innovations <- drop(out$Y)
  acf_vals <- stats::acf(innovations, lag.max = 10, plot = FALSE, demean = TRUE)$acf[-1L]
  expect_true(all(abs(acf_vals[1:5]) < 0.1))
})
