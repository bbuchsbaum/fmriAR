test_that("parcel whitening returns expected structure", {
  set.seed(3)
  n <- 200
  v <- 60
  parcels <- rep(1:12, each = v / 12)
  X <- cbind(1, rnorm(n))
  Y <- matrix(rnorm(n * v), n, v)
  res <- Y - X %*% qr.solve(X, Y)

  plan <- fit_noise(res, runs = NULL, method = "ar", p = "auto", p_max = 4,
                    pooling = "parcel", parcels = parcels)
  out <- whiten_apply(plan, X, Y, parcels = parcels)

  expect_null(out$X)
  expect_true(is.list(out$X_by))
  expect_equal(length(out$X_by), length(unique(parcels)))
  for (pid in names(out$X_by)) {
    expect_equal(dim(out$X_by[[pid]]), dim(X))
  }
  expect_equal(dim(out$Y), dim(Y))
})

test_that("multiscale pooling improves whiteness", {
  set.seed(11)
  n <- 400
  v <- 120
  parcels_coarse <- rep(1:6, each = v / 6)
  parcels_medium <- rep(1:12, each = v / 12)
  parcels_fine <- rep(1:24, each = v / 24)

  phi_coarse <- runif(6, 0.3, 0.8)
  phi_vox <- phi_coarse[parcels_coarse] + rnorm(v, 0, 0.03)

  noise <- matrix(rnorm(n * v), n, v)
  Y <- noise
  for (j in seq_len(v)) {
    for (t in 2:n) {
      Y[t, j] <- phi_vox[j] * Y[t - 1, j] + noise[t, j]
    }
  }
  X <- cbind(1, rnorm(n))
  res <- Y - X %*% qr.solve(X, Y)

  plan_f <- fit_noise(res, runs = NULL, method = "ar", p = "auto", p_max = 4,
                      pooling = "parcel", parcels = parcels_fine)
  out_f <- whiten_apply(plan_f, X, Y, parcels = parcels_fine)

  plan_ms <- fit_noise(res, runs = NULL, method = "ar", p = "auto", p_max = 4,
                       pooling = "parcel", parcels = parcels_fine,
                       parcel_sets = list(coarse = parcels_coarse,
                                          medium = parcels_medium,
                                          fine = parcels_fine),
                       multiscale = "pacf_weighted", beta = 0.5)
  out_ms <- whiten_apply(plan_ms, X, Y, parcels = parcels_fine)

  whiteness <- function(Yw) {
    apply(Yw, 2, function(y) {
      ac <- stats::acf(y, lag.max = 3, plot = FALSE, demean = TRUE)$acf[-1L]
      mean(abs(ac))
    }) |> mean()
  }

  m_f <- whiteness(out_f$Y)
  m_ms <- whiteness(out_ms$Y)

  expect_lte(m_ms, m_f)
})

test_that("parcel means helper matches brute-force computation", {
  set.seed(51)
  n <- 30
  v <- 25
  resid <- matrix(rnorm(n * v), nrow = n, ncol = v)
  parcels <- sample(1:8, size = v, replace = TRUE)

  fast <- fmriAR:::`.parcel_means`(resid, parcels)
  ref <- sapply(sort(unique(parcels)), function(pid) {
    idx <- which(parcels == pid)
    rowMeans(resid[, idx, drop = FALSE])
  })
  colnames(ref) <- as.character(sort(unique(parcels)))

  expect_equal(fast, ref, tolerance = 1e-12)
})

test_that("segment-aware ACVF pools within runs", {
  set.seed(99)
  n <- 150
  runs <- rep(1:3, c(40, 60, 50))
  y <- as.numeric(stats::arima.sim(list(ar = 0.6), n = n))
  run_starts0 <- fmriAR:::`.full_run_starts`(runs, censor = NULL, n = n)
  lag_max <- 4L
  acvf_seg <- fmriAR:::`.segment_acvf`(y, run_starts0, lag_max)

  segments <- split(y, runs)
  manual <- numeric(lag_max + 1L)
  total_n <- 0L
  for (seg in segments) {
    L <- length(seg)
    mu <- mean(seg)
    total_n <- total_n + L
    for (lag in 0:lag_max) {
      if (lag < L) {
        seg1 <- seg[(lag + 1L):L] - mu
        seg2 <- seg[1:(L - lag)] - mu
        manual[lag + 1L] <- manual[lag + 1L] + sum(seg1 * seg2)
      }
    }
  }
  manual <- manual / total_n

  expect_equal(acvf_seg, manual, tolerance = 1e-10)
})
