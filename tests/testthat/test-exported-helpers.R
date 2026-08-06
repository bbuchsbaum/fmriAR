test_that("acorr_diagnostics aggregate = 'none' matches per-voxel acf", {
  resid <- cbind(
    c(0, 1, 0, -1, 0, 1),
    c(1, 0, -1, 0, 1, 0)
  )

  out <- acorr_diagnostics(resid, max_lag = 3L, aggregate = "none")
  ref <- vapply(
    seq_len(ncol(resid)),
    function(j) stats::acf(resid[, j], lag.max = 3L, plot = FALSE, demean = TRUE)$acf[-1L],
    numeric(3L)
  )

  expect_equal(out$acf, ref, tolerance = 1e-12)
  expect_equal(out$lags, 1:3)
})

test_that("acorr_diagnostics aggregate = 'median' matches median trajectory acf", {
  resid <- cbind(
    c(0, 2, 0, -2, 0, 2),
    c(1, 1, -1, -1, 1, 1),
    c(-1, 0, 1, 0, -1, 0)
  )

  out <- acorr_diagnostics(resid, max_lag = 2L, aggregate = "median")
  ref_series <- apply(resid, 1L, stats::median)
  ref <- stats::acf(ref_series, lag.max = 2L, plot = FALSE, demean = TRUE)$acf[-1L]

  expect_equal(out$acf, ref, tolerance = 1e-12)
})

test_that("sandwich_from_whitened_resid iid matches closed-form GLS standard errors", {
  Xw <- cbind(1, c(-1, 0, 1, 2, 3))
  beta <- rbind(c(0.5, -0.25), c(1.0, 0.75))
  E <- matrix(c(
    0.10, -0.20,
   -0.05,  0.10,
    0.08, -0.12,
   -0.02,  0.05,
    0.04,  0.15
  ), nrow = 5, byrow = TRUE)
  Yw <- Xw %*% beta + E

  out <- sandwich_from_whitened_resid(Xw, Yw, beta = beta, type = "iid")

  XtX_inv <- solve(crossprod(Xw))
  sigma2 <- colSums(E^2) / (nrow(Xw) - qr(Xw)$rank)
  ref_se <- sqrt(outer(diag(XtX_inv), sigma2))

  expect_equal(out$XtX_inv, XtX_inv, tolerance = 1e-12)
  expect_equal(out$sigma2, sigma2, tolerance = 1e-12)
  expect_equal(out$se, ref_se, tolerance = 1e-12)
})

test_that("sandwich_from_whitened_resid hc0 matches manual sandwich construction", {
  Xw <- cbind(1, c(0, 1, 2, 3, 4))
  beta <- rbind(c(0.2, -0.1), c(0.4, 0.6))
  E <- matrix(c(
    0.20, -0.10,
   -0.10,  0.15,
    0.05, -0.20,
    0.10,  0.05,
   -0.15,  0.10
  ), nrow = 5, byrow = TRUE)
  Yw <- Xw %*% beta + E

  out <- sandwich_from_whitened_resid(Xw, Yw, beta = beta, type = "hc0")
  XtX_inv <- solve(crossprod(Xw))

  ref_se <- matrix(NA_real_, nrow = ncol(Xw), ncol = ncol(Yw))
  for (j in seq_len(ncol(Yw))) {
    Xe <- Xw * E[, j]
    meat <- crossprod(Xe, Xe)
    ref_se[, j] <- sqrt(diag(XtX_inv %*% meat %*% XtX_inv))
  }

  expect_equal(out$XtX_inv, XtX_inv, tolerance = 1e-12)
  expect_equal(out$se, ref_se, tolerance = 1e-12)
  expect_equal(out$type, "hc0")
})
