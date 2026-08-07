library(testthat)
library(fmriAR)

# Exact expectation of the raw estimator, computed by brute force from the
# residual operator and an explicitly built Sigma. This is the definitive check
# on the bias matrix: no simulation, no tolerance for Monte Carlo error.
#
# Sigma is given finite support within the lag budget, so Sigma = sum_k gamma_k
# S_k holds exactly and A gamma must reproduce the expectation to machine
# precision. A gamma with tails past the budget would only match approximately,
# which would test nothing sharply.
exact_expectation <- function(X, gam, runs = NULL, censor = NULL) {
  n <- nrow(X)
  H <- length(gam) - 1L
  run_vec <- if (is.null(runs)) rep(1L, n) else as.integer(match(runs, unique(runs)))
  M <- diag(n) - X %*% solve(crossprod(X), t(X))

  Sig <- matrix(0, n, n)
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      k <- abs(i - j)
      if (k <= H && run_vec[i] == run_vec[j]) Sig[i, j] <- gam[k + 1L]
    }
  }

  Rsets <- if (is.null(runs)) list(seq_len(n)) else split(seq_len(n), run_vec)
  cens <- if (is.null(censor)) integer(0) else as.integer(censor)

  lapply(Rsets, function(ridx) {
    keep <- setdiff(ridx, cens)
    nv <- length(keep)
    seg_id <- cumsum(c(1L, as.integer(diff(keep) > 1L)))
    R <- M[keep, , drop = FALSE]
    R <- R - rep(colMeans(R), each = nv)
    RSR <- R %*% Sig %*% t(R)
    out <- numeric(H + 1L)
    out[1L] <- mean(diag(RSR))
    for (h in seq_len(H)) {
      hi <- seq.int(h + 1L, nv); lo <- seq.int(1L, nv - h)
      ok <- seg_id[hi] == seg_id[lo]
      out[h + 1L] <- if (any(ok)) mean(RSR[cbind(hi[ok], lo[ok])]) else 0
    }
    out
  })
}

test_that("the bias matrix reproduces the exact expectation of the raw estimator", {
  n <- 60L; H <- 6L
  gam <- c(1, 0.5, 0.2, 0, 0, 0, 0)   # finite support inside the budget
  X <- cbind(1, poly(seq_len(n), 3))

  A <- acvf_bias_matrix(X, max_lag = H)
  ex <- exact_expectation(X, gam)
  expect_equal(as.numeric(A[[1]] %*% gam), ex[[1]], tolerance = 1e-10)

  # The whole point: the raw expectation is NOT the truth, so this test would be
  # vacuous if A were the identity.
  expect_false(isTRUE(all.equal(as.numeric(A[[1]] %*% gam), gam, tolerance = 1e-3)))
})

test_that("the bias matrix is exact with runs and censoring too", {
  n <- 80L; H <- 5L
  gam <- c(1, 0.45, 0.15, 0.05, 0, 0)
  runs <- rep(1:2, each = 40L)
  censor <- c(5:9, 50:52, 70L)
  X <- cbind(1, poly(seq_len(n), 2), matrix(sin(seq_len(n) / 7), n, 1))

  A <- acvf_bias_matrix(X, runs = runs, censor = censor, max_lag = H)
  ex <- exact_expectation(X, gam, runs = runs, censor = censor)
  expect_length(A, 2L)
  for (i in seq_along(A)) {
    expect_equal(as.numeric(A[[i]] %*% gam), ex[[i]], tolerance = 1e-10,
                 info = paste("run", i))
  }
})

test_that("lag products never cross a run boundary in the bias matrix", {
  # Two runs whose concatenation would create spurious lag-1 pairs at the seam.
  # Building A run by run must give the same answer as building it for each run
  # analysed entirely on its own would, for the seam pairs specifically: the
  # count of pairs behind each lag must equal the within-run count.
  n <- 40L; H <- 3L
  runs <- rep(1:2, each = 20L)
  X <- cbind(1, seq_len(n) / n)
  A <- acvf_bias_matrix(X, runs = runs, max_lag = H)
  expect_length(A, 2L)
  expect_equal(dim(A[[1]]), c(H + 1L, H + 1L))
  # A run of 20 frames supports 20 - h pairs at lag h; a seam-crossing estimator
  # would find 40 - h across the concatenation.
  gam <- c(1, rep(0, H))
  for (i in seq_along(A)) expect_equal(as.numeric(A[[i]] %*% gam)[1], mean(diag(
    { M <- diag(n) - X %*% solve(crossprod(X), t(X))
      keep <- which(runs == i); nv <- length(keep)
      R <- M[keep, , drop = FALSE]; R <- R - rep(colMeans(R), each = nv)
      R %*% t(R) })), tolerance = 1e-10)
})

test_that("correction reduces AR bias and does not invent autocorrelation", {
  skip_on_cran()
  set.seed(4242)
  n <- 300L; nvox <- 60L
  X <- cbind(1, poly(seq_len(n), 3), matrix(rnorm(n * 20), n, 20))
  M <- diag(n) - X %*% solve(crossprod(X), t(X))

  sim <- function(phi) {
    R <- matrix(0, n, nvox)
    for (v in seq_len(nvox)) {
      e <- if (phi == 0) stats::rnorm(n) else as.numeric(stats::arima.sim(list(ar = phi), n))
      R[, v] <- M %*% e
    }
    R
  }

  R <- sim(0.4)
  raw <- fit_noise(R, method = "ar", p = 1L, pooling = "global")
  cor <- fit_noise(R, method = "ar", p = 1L, pooling = "global", design = X)
  expect_lt(abs(cor$phi[[1]][1] - 0.4), abs(raw$phi[[1]][1] - 0.4))
  expect_lt(abs(cor$phi[[1]][1] - 0.4), 0.06)

  # White noise must stay white: the correction must not manufacture structure
  # where there is none. This is the failure mode that would make it dangerous.
  Rw <- sim(0)
  cw <- fit_noise(Rw, method = "ar", p = 1L, pooling = "global", design = X)
  expect_lt(abs(cw$phi[[1]][1]), 0.08)
})

test_that("corrected gamma remains a valid covariance", {
  skip_on_cran()
  set.seed(99)
  n <- 240L
  X <- cbind(1, poly(seq_len(n), 3), matrix(rnorm(n * 10), n, 10))
  M <- diag(n) - X %*% solve(crossprod(X), t(X))
  R <- matrix(0, n, 30L)
  for (v in 1:30) R[, v] <- M %*% as.numeric(stats::arima.sim(list(ar = 0.6), n))

  for (runs in list(NULL, rep(1:2, each = 120L))) {
    for (cens in list(NULL, c(10:20, 150:160))) {
      pl <- fit_noise(R, runs = runs, censor = cens, method = "ar", p = "auto",
                      pooling = "global", design = X)
      g <- pl$gamma[[1]]
      expect_true(length(g) >= 1L)
      expect_true(all(is.finite(g)))
      expect_gt(g[1], 0)
      if (length(g) > 1L) {
        ev <- eigen(stats::toeplitz(g), symmetric = TRUE, only.values = TRUE)$values
        expect_gt(min(ev), 0)
      }
      expect_true(all(abs(fmriAR:::ar_to_pacf(pl$phi[[1]])) < 1))
    }
  }
})

test_that("acvf_correction reproduces what design computes", {
  skip_on_cran()
  set.seed(808)
  n <- 200L
  X <- cbind(1, poly(seq_len(n), 3), matrix(rnorm(n * 6), n, 6))
  M <- diag(n) - X %*% solve(crossprod(X), t(X))
  R <- matrix(0, n, 20L)
  for (v in 1:20) R[, v] <- M %*% as.numeric(stats::arima.sim(list(ar = 0.5), n))

  A <- acvf_bias_matrix(X, max_lag = 25L)
  by_design <- fit_noise(R, method = "ar", p = 2L, pooling = "global", design = X)
  by_matrix <- fit_noise(R, method = "ar", p = 2L, pooling = "global",
                         acvf_correction = A)
  expect_equal(by_design$phi[[1]], by_matrix$phi[[1]], tolerance = 1e-12)
  expect_equal(by_design$gamma[[1]], by_matrix$gamma[[1]], tolerance = 1e-12)

  # A bare matrix is broadcast to every run.
  by_bare <- fit_noise(R, method = "ar", p = 2L, pooling = "global",
                       acvf_correction = A[[1]])
  expect_equal(by_design$phi[[1]], by_bare$phi[[1]], tolerance = 1e-12)
})

test_that("unsupported correction combinations are refused, not silently ignored", {
  set.seed(2)
  n <- 200L
  R <- matrix(rnorm(n * 8), n, 8)
  X <- cbind(1, poly(seq_len(n), 3))

  expect_error(fit_noise(R, design = X, acvf_correction = acvf_bias_matrix(X, max_lag = 3L)),
               "not both")
  expect_error(fit_noise(R, pooling = "parcel", parcels = rep(1:4, each = 2), design = X),
               "not yet supported")
  expect_error(fit_noise(R, method = "arma", p = 1L, q = 1L, design = X),
               "method = 'ar' only")
  expect_error(fit_noise(R, design = X[1:50, , drop = FALSE]),
               "one row per timepoint")
  expect_error(fit_noise(R, acvf_correction = list(diag(4), diag(4), diag(4))),
               "3 matrices but there are 1 runs")
})

test_that("correction restores the noise scale, not just its shape, for AR(2)", {
  skip_on_cran()
  # gamma_0 matters on its own: it is the noise variance a consumer builds
  # standard errors from. A 28-column design understates it by about 24%
  # (1.72 against a truth of 2.24), which understates every standard error
  # derived from it. Both AR(2) coefficients must improve as well -- a
  # correction that fixed lag 1 by distorting lag 2 would be no use.
  set.seed(5150)
  n <- 300L; nvox <- 60L
  phi <- c(0.5, 0.3)
  g0_true <- (1 - phi[2]) / ((1 + phi[2]) * ((1 - phi[2])^2 - phi[1]^2))

  X <- cbind(1, poly(seq_len(n), 3), matrix(rnorm(n * 24), n, 24))
  M <- diag(n) - X %*% solve(crossprod(X), t(X))
  R <- matrix(0, n, nvox)
  for (v in seq_len(nvox)) R[, v] <- M %*% as.numeric(stats::arima.sim(list(ar = phi), n))

  raw <- fit_noise(R, method = "ar", p = 2L, pooling = "global")
  cor <- fit_noise(R, method = "ar", p = 2L, pooling = "global", design = X)

  expect_lt(abs(cor$gamma[[1]][1] - g0_true), abs(raw$gamma[[1]][1] - g0_true))
  expect_lt(abs(cor$gamma[[1]][1] - g0_true) / g0_true, 0.12)
  for (j in 1:2) {
    expect_lt(abs(cor$phi[[1]][j] - phi[j]), abs(raw$phi[[1]][j] - phi[j]) + 1e-8,
              label = paste("coefficient", j))
  }
})

test_that("correction still helps across runs and censoring", {
  skip_on_cran()
  set.seed(6161)
  n <- 300L; nvox <- 40L; phi <- 0.5
  runs <- rep(1:4, each = 75L)
  X <- cbind(1, poly(seq_len(n), 3), matrix(rnorm(n * 12), n, 12))
  M <- diag(n) - X %*% solve(crossprod(X), t(X))
  R <- matrix(0, n, nvox)
  for (v in seq_len(nvox)) {
    e <- numeric(n)
    for (r in 1:4) {
      idx <- which(runs == r)
      e[idx] <- as.numeric(stats::arima.sim(list(ar = phi), length(idx)))
    }
    R[, v] <- M %*% e
  }
  for (cf in list(NULL, sort(sample(seq_len(n), 60L)), sort(sample(seq_len(n), 120L)))) {
    raw <- fit_noise(R, runs = runs, censor = cf, method = "ar", p = 1L,
                     pooling = "global")
    cor <- fit_noise(R, runs = runs, censor = cf, method = "ar", p = 1L,
                     pooling = "global", design = X)
    lbl <- paste("censored", length(cf))
    expect_lt(abs(cor$phi[[1]][1] - phi), abs(raw$phi[[1]][1] - phi) + 1e-8, label = lbl)
    expect_lt(abs(cor$phi[[1]][1] - phi), 0.1, label = lbl)
    expect_true(all(is.finite(cor$gamma[[1]])), label = lbl)
  }
})

test_that("an ill-conditioned correction is refused rather than trusted", {
  skip_on_cran()
  set.seed(77)
  n <- 300L
  X <- cbind(1, poly(seq_len(n), 3), matrix(rnorm(n * 4), n, 4))
  M <- diag(n) - X %*% solve(crossprod(X), t(X))
  R <- matrix(0, n, 20L)
  for (v in 1:20) R[, v] <- M %*% as.numeric(stats::arima.sim(list(ar = 0.5), n))

  # A lag budget near the run length drives rcond(A) to ~5e-9, and solving
  # there returned phi = 0.96 for a truth of 0.5 -- finite and positive, so a
  # NaN check would have passed it straight through.
  expect_warning(bad <- fit_noise(R, method = "ar", p = 1L, pooling = "global",
                                  design = X, correction_max_lag = 250L),
                 "ill-conditioned")
  plain <- fit_noise(R, method = "ar", p = 1L, pooling = "global")
  # Refused means it fell back to the uncorrected estimate exactly.
  expect_equal(bad$phi[[1]], plain$phi[[1]], tolerance = 1e-12)
  expect_lt(abs(bad$phi[[1]][1] - 0.5), 0.15)

  # The sane budget is untouched and still corrects.
  good <- fit_noise(R, method = "ar", p = 1L, pooling = "global", design = X)
  expect_lt(abs(good$phi[[1]][1] - 0.5), abs(plain$phi[[1]][1] - 0.5))
})

test_that("too few residual degrees of freedom caps the lag budget", {
  skip_on_cran()
  set.seed(1234)
  n <- 300L
  X <- cbind(1, matrix(rnorm(n * 280), n, 280))   # 281 columns, 19 residual df
  M <- diag(n) - X %*% solve(crossprod(X), t(X))
  R <- matrix(0, n, 20L)
  for (v in 1:20) R[, v] <- M %*% as.numeric(stats::arima.sim(list(ar = 0.5), n))

  expect_warning(out <- fit_noise(R, method = "ar", p = 1L, pooling = "global",
                                  design = X),
                 "residual degrees of freedom")
  plain <- fit_noise(R, method = "ar", p = 1L, pooling = "global")

  # Uncorrected, a design this rich flattens the noise almost completely: the
  # correction has to restore both the shape and the scale.
  expect_lt(abs(plain$phi[[1]][1]), 0.1)
  expect_gt(out$phi[[1]][1], 0.15)
  expect_gt(out$gamma[[1]][1], 5 * plain$gamma[[1]][1])
  # And must not pin at the stationarity clamp, which is what it did uncapped.
  expect_lt(out$phi[[1]][1], 0.95)
  expect_true(all(is.finite(out$gamma[[1]])))
})

test_that("omitting design leaves estimates bit-identical", {
  # The correction is opt-in. Nothing about adding the arguments may perturb the
  # default path.
  set.seed(31337)
  R <- matrix(rnorm(300 * 12), 300, 12)
  runs <- rep(1:2, each = 150L)
  a <- fit_noise(R, runs = runs, method = "ar", p = "auto", pooling = "global")
  b <- fit_noise(R, runs = runs, method = "ar", p = "auto", pooling = "global",
                 design = NULL, acvf_correction = NULL, correction_max_lag = 25L)
  expect_identical(a$phi, b$phi)
  expect_identical(a$gamma, b$gamma)
  expect_identical(a$sigma2, b$sigma2)
})
