test_that("fit_noise with Y/X matches precomputed residuals", {
  set.seed(42)
  n <- 160
  runs <- rep(1:2, each = n / 2)
  X <- cbind(1, rnorm(n))
  beta <- c(0.5, -0.2)
  eps <- as.numeric(stats::arima.sim(list(ar = 0.4), n = n))
  Y <- matrix(X %*% beta + eps, ncol = 1L)
  resid <- Y - X %*% qr.solve(X, Y)

  plan_resid <- fit_noise(resid = resid, runs = runs, method = "ar", p = "auto", p_max = 4)
  plan_yx    <- fit_noise(Y = Y, X = X, runs = runs, method = "ar", p = "auto", p_max = 4)

  expect_equal(plan_resid$order, plan_yx$order)
  expect_equal(plan_resid$phi, plan_yx$phi)
  expect_equal(plan_resid$theta, plan_yx$theta)
})

test_that("fit_noise ignores ms_mode when multiscale = FALSE", {
  set.seed(9)
  n <- 120
  v <- 20
  X <- cbind(1, rnorm(n))
  Y <- matrix(rnorm(n * v), n, v)
  resid <- Y - X %*% qr.solve(X, Y)
  parcels <- rep(1:5, length.out = v)

  plan_plain <- fit_noise(resid, pooling = "parcel", parcels = parcels,
                          method = "ar", p = 2L)
  plan_off   <- fit_noise(resid, pooling = "parcel", parcels = parcels,
                          method = "ar", p = 2L,
                          multiscale = FALSE, ms_mode = "acvf_pooled")

  out_plain <- whiten_apply(plan_plain, X, Y, parcels = parcels, parallel = FALSE)
  out_off   <- whiten_apply(plan_off, X, Y, parcels = parcels, parallel = FALSE)

  expect_equal(out_plain$Y, out_off$Y)
})

test_that("whiten_apply honors run_starts input", {
  set.seed(7)
  n1 <- 60
  n2 <- 80
  phi1 <- 0.5
  phi2 <- -0.2
  y1 <- as.numeric(stats::filter(rnorm(n1), filter = phi1, method = "recursive"))
  y2 <- as.numeric(stats::filter(rnorm(n2), filter = phi2, method = "recursive"))
  Y <- rbind(matrix(y1, ncol = 1L), matrix(y2, ncol = 1L))
  X <- matrix(1, nrow(Y), 1L)

  plan <- fmriAR:::new_whiten_plan(
    phi = list(phi1, phi2),
    theta = list(numeric(0), numeric(0)),
    order = c(p = 1L, q = 0L),
    runs = NULL,
    exact_first = FALSE,
    method = "ar",
    pooling = "run"
  )

  out_runs <- whiten_apply(plan, X, Y, runs = c(rep(1L, n1), rep(2L, n2)), parallel = FALSE)
  out_starts <- whiten_apply(plan, X, Y, run_starts = c(0L, n1), parallel = FALSE)

  expect_equal(out_runs$Y, out_starts$Y)
  expect_equal(out_runs$X, out_starts$X)
})

test_that("whiten_apply inplace modifies matrices", {
  set.seed(12)
  n <- 90
  phi <- 0.4
  y <- as.numeric(stats::filter(rnorm(n), filter = phi, method = "recursive"))
  Y <- matrix(y, ncol = 1L)
  X <- matrix(1, n, 1L)
  plan <- fmriAR:::new_whiten_plan(
    phi = list(phi),
    theta = list(numeric(0)),
    order = c(p = 1L, q = 0L),
    runs = NULL,
    exact_first = FALSE,
    method = "ar",
    pooling = "global"
  )
  Y_ref <- Y
  X_ref <- X
  res_inplace <- whiten_apply(plan, X, Y, inplace = TRUE, parallel = FALSE)
  res_regular <- whiten_apply(plan, X_ref, Y_ref, parallel = FALSE)

  expect_equal(res_inplace$Y, res_regular$Y)
  expect_equal(res_inplace$X, res_regular$X)
})

test_that("global pooling averages run-level coefficients", {
  set.seed(23)
  n1 <- 400
  n2 <- 600
  phi1 <- 0.3
  phi2 <- 0.7
  y1 <- as.numeric(stats::filter(rnorm(n1), filter = phi1, method = "recursive"))
  y2 <- as.numeric(stats::filter(rnorm(n2), filter = phi2, method = "recursive"))
  resid <- rbind(cbind(y1, y1), cbind(y2, y2))
  runs <- c(rep(1L, n1), rep(2L, n2))

  plan <- fit_noise(resid = resid, runs = runs, pooling = "global",
                    method = "ar", p = 1L)
  if (plan$order[["p"]] == 0L) testthat::skip("Estimated AR order was zero")
  phi_hat <- plan$phi[[1]][1]
  expected <- (n1 * phi1 + n2 * phi2) / (n1 + n2)
  expect_equal(phi_hat, expected, tolerance = 0.05)
})

test_that("multiscale acvf_pooled combines parcel scales", {
  set.seed(19)
  n <- 300
  v <- 96
  parcels_coarse <- rep(1:4, each = v / 4)
  parcels_medium <- rep(1:8, each = v / 8)
  parcels_fine <- rep(1:16, each = v / 16)

  phi_coarse <- runif(4, 0.3, 0.7)
  phi_vox <- phi_coarse[parcels_coarse] + rnorm(v, 0, 0.04)

  noise <- matrix(rnorm(n * v), n, v)
  Y <- noise
  for (j in seq_len(v)) {
    for (t in 2:n) {
      Y[t, j] <- phi_vox[j] * Y[t - 1, j] + noise[t, j]
    }
  }
  X <- cbind(1, rnorm(n))
  resid <- Y - X %*% qr.solve(X, Y)

  plan_fine <- fit_noise(resid, pooling = "parcel", parcels = parcels_fine,
                         method = "ar", p = "auto", p_max = 4)
  plan_ms <- fit_noise(resid, pooling = "parcel", parcels = parcels_fine,
                       parcel_sets = list(coarse = parcels_coarse,
                                          medium = parcels_medium,
                                          fine = parcels_fine),
                       method = "ar", p = "auto", p_max = 4,
                       multiscale = TRUE, ms_mode = "acvf_pooled")

  expect_equal(plan_ms$order[["p"]], plan_fine$order[["p"]])
  expect_true(any(abs(unlist(plan_ms$phi_by_parcel) - unlist(plan_fine$phi_by_parcel)) > 1e-6))

  out_ms <- whiten_apply(plan_ms, X, Y, parcels = parcels_fine, parallel = FALSE)
  expect_equal(dim(out_ms$Y), dim(Y))
})

test_that("parcel means handles NA removal", {
  set.seed(33)
  resid <- matrix(rnorm(30), nrow = 5)
  resid[2, 3] <- NA_real_
  parcels <- c(1, 1, 2, 2, 3, 3)

  fast <- fmriAR:::`.parcel_means`(resid, parcels, na.rm = TRUE)
  manual <- matrix(0, nrow = 5, ncol = 3)
  manual[, 1] <- rowMeans(resid[, 1:2], na.rm = TRUE)
  manual[, 2] <- rowMeans(resid[, 3:4], na.rm = TRUE)
  manual[, 3] <- rowMeans(resid[, 5:6], na.rm = TRUE)
  colnames(manual) <- colnames(fast)

  expect_equal(fast, manual)
})

test_that("segmented_acvf supports unbiased and no-centering", {
  y <- c(1, 3, 5, 7)
  rs0 <- 0L
  g_biased <- fmriAR:::`.segment_acvf`(y, rs0, lag_max = 2L, unbiased = FALSE, center = FALSE)
  g_unbiased <- fmriAR:::`.segment_acvf`(y, rs0, lag_max = 2L, unbiased = TRUE, center = FALSE)

  expect_equal(g_biased[1], mean(y^2))
  expect_equal(g_unbiased[2], sum(y[-1] * y[-length(y)]) / (length(y) - 1))
})

test_that("whiten() matches fit_noise + whiten_apply", {
  set.seed(27)
  n <- 180
  runs <- rep(1:3, each = 60)
  X <- cbind(1, rnorm(n))
  beta <- c(0.2, 0.5)
  eps <- as.numeric(stats::arima.sim(list(ar = 0.3), n = n))
  Y <- matrix(X %*% beta + eps, ncol = 1L)

  plan <- fit_noise(resid = Y - X %*% qr.solve(X, Y), runs = runs,
                    method = "ar", p = "auto", p_max = 4)
  manual <- whiten_apply(plan, X, Y, runs = runs, parallel = FALSE)
  shortcut <- whiten(X, Y, runs = runs, method = "ar", p = "auto", p_max = 4)

  expect_equal(manual$X, shortcut$X)
  expect_equal(manual$Y, shortcut$Y)
})

test_that("compat plan_info exposes parcel coefficients", {
  set.seed(101)
  n <- 160
  v <- 40
  parcels <- rep(1:8, each = v / 8)
  X <- cbind(1, rnorm(n))
  Y <- matrix(rnorm(n * v), n, v)
  resid <- Y - X %*% qr.solve(X, Y)

  plan <- fit_noise(resid, pooling = "parcel", parcels = parcels,
                    method = "ar", p = 2L, p_max = 4)
  info <- fmriAR:::compat$plan_info(plan)

  expect_true(!is.null(info$phi_by_parcel))
  expect_equal(dim(info$phi_by_parcel), c(plan$order[["p"]], length(unique(parcels))))
})

test_that("compat plan_info omits parcel coefficients for global plans", {
  set.seed(55)
  resid <- matrix(rnorm(200), nrow = 100, ncol = 2)
  plan <- fit_noise(resid, method = "ar", p = 1L, pooling = "global")
  info <- fmriAR:::compat$plan_info(plan)
  expect_true(is.null(info$phi_by_parcel))
})
