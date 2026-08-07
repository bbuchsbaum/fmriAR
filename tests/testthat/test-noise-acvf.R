library(testthat)
library(fmriAR)

ar1 <- function(n, phi, nvox = 1L, seed = NULL, sd = 1) {
  if (!is.null(seed)) set.seed(seed)
  burn <- 200L
  out <- matrix(0, n, nvox)
  for (j in seq_len(nvox)) {
    e <- stats::rnorm(n + burn, sd = sd)
    y <- numeric(n + burn)
    for (t in 2:(n + burn)) y[t] <- phi * y[t - 1L] + e[t]
    out[, j] <- y[(burn + 1L):(burn + n)]
  }
  out
}

test_that("noise_acvf returns covariances on the scale of the data", {
  A <- ar1(400L, 0.5, nvox = 30L, seed = 11, sd = 3)
  out <- noise_acvf(A, max_lag = 6L)
  expect_s3_class(out, "fmriAR_acvf")
  expect_length(out$acvf, 1L)
  g <- out$acvf[[1L]]
  expect_length(g, 7L)
  # Covariance, not correlation: gamma_0 tracks the empirical variance.
  expect_equal(g[1], mean(apply(A, 2, var)), tolerance = 0.05)
  expect_equal(g[2] / g[1], 0.5, tolerance = 0.08)
  # A valid covariance.
  expect_gt(min(eigen(stats::toeplitz(g), symmetric = TRUE, only.values = TRUE)$values), 0)
})

test_that("lag products never cross a run boundary", {
  # The pair counts are the tell. Two runs of 100 support (100 - h) pairs each;
  # an estimator that ran across the seam would report (200 - h).
  A <- ar1(200L, 0.4, nvox = 8L, seed = 21)
  runs <- rep(1:2, each = 100L)

  per_run <- noise_acvf(A, runs = runs, max_lag = 4L, pooling = "run")
  for (u in per_run$pairs) expect_equal(u, c(100, 99, 98, 97, 96))

  glob <- noise_acvf(A, runs = runs, max_lag = 4L, pooling = "global")
  expect_equal(glob$pairs[[1L]], c(200, 198, 196, 194, 192))

  # Without run labels the same data is treated as one series, which is exactly
  # the seam-crossing behaviour, and the counts differ.
  flat <- noise_acvf(A, max_lag = 4L)
  expect_equal(flat$pairs[[1L]], c(200, 199, 198, 197, 196))
  expect_false(isTRUE(all.equal(flat$acvf[[1L]], glob$acvf[[1L]])))
})

test_that("run labels are validated and encoded consistently", {
  A <- ar1(200L, 0.4, nvox = 6L, seed = 22)
  numeric_runs <- rep(c(20L, 10L), each = 100L)
  character_runs <- rep(c("run-b", "run-a"), each = 100L)

  num <- noise_acvf(A, runs = numeric_runs, max_lag = 3L, pooling = "run")
  chr <- noise_acvf(A, runs = character_runs, max_lag = 3L, pooling = "run")
  expect_equal(unname(chr$acvf), unname(num$acvf), tolerance = 1e-12)
  expect_named(chr$acvf, c("run-b", "run-a"))

  expect_error(noise_acvf(A, runs = rep(1:2, each = 20L)), "length")
  bad_na <- numeric_runs
  bad_na[50L] <- NA
  expect_error(noise_acvf(A, runs = bad_na), "NA")
  expect_error(noise_acvf(A, runs = rep(c("a", "b", "a", "b"), each = 50L)),
               "contiguous")
})

test_that("censoring segments the series and is reported", {
  A <- ar1(300L, 0.5, nvox = 10L, seed = 31)
  cens <- c(50:59, 200:204)
  out <- noise_acvf(A, censor = cens, max_lag = 5L)
  mask <- rep(FALSE, nrow(A))
  mask[cens] <- TRUE
  out_logical <- noise_acvf(A, censor = mask, max_lag = 5L)

  expect_equal(out$n_segments[[1L]], 3L)
  expect_equal(sum(out$segment_lengths[[1L]]), 300L - length(cens))
  # Lag-0 pairs are exactly the surviving frames; higher lags lose one pair per
  # segment per lag because products cannot span a gap.
  expect_equal(out$pairs[[1L]][1], 300 - length(cens))
  expect_lt(out$pairs[[1L]][2], out$pairs[[1L]][1])
  expect_equal(out$acvf[[1L]][2] / out$acvf[[1L]][1], 0.5, tolerance = 0.1)
  expect_equal(out_logical$acvf, out$acvf, tolerance = 1e-12)
  expect_equal(out_logical$pairs, out$pairs)
})

test_that("global pooling weights runs by surviving observations", {
  # Run 2 has 100x the variance but only two surviving frames. It must carry
  # 2/22 of the pooled estimate, not half merely because both original runs had
  # 20 frames.
  A <- matrix(c(rep(c(-1, 1), 10L), rep(c(-10, 10), 10L)), ncol = 1L)
  runs <- rep(1:2, each = 20L)
  censor <- 23:40

  per <- noise_acvf(A, runs = runs, censor = censor, max_lag = 1L,
                    pooling = "run")
  glob <- noise_acvf(A, runs = runs, censor = censor, max_lag = 1L,
                     pooling = "global")
  expected_g0 <- weighted.mean(vapply(per$acvf, `[`, 0, 1L), c(20, 2))
  expect_equal(glob$acvf[[1L]][1], expected_g0, tolerance = 1e-12)

  plan <- fit_noise(A, runs = runs, censor = censor, method = "ar", p = 1L,
                    pooling = "global")
  expect_equal(plan$gamma[[1L]], glob$acvf[[1L]], tolerance = 1e-12)
})

test_that("a run with no surviving frames does not erase the pooled gamma", {
  # A run censored down to <= 1 frame carries no autocovariance at all. It must
  # simply drop out of the pooled gamma; letting its empty vector set the
  # common truncation length collapsed gamma to numeric(0) and sigma2 to NA
  # for the whole global fit, while phi was reported correctly.
  A <- ar1(200L, 0.5, nvox = 8L, seed = 44)
  runs <- rep(1:2, each = 100L)
  for (keep in 0:1) {
    censor <- 101:(200L - keep)
    plan <- fit_noise(A, runs = runs, censor = censor, method = "ar", p = 2L,
                      pooling = "global")
    ac <- noise_acvf(A, runs = runs, censor = censor, max_lag = 2L,
                     pooling = "global")
    expect_gt(length(plan$gamma[[1L]]), 0)
    expect_true(is.finite(plan$sigma2[[1L]]))
    expect_equal(plan$gamma[[1L]][seq_along(ac$acvf[[1L]])], ac$acvf[[1L]],
                 tolerance = 1e-12, info = paste("surviving frames:", keep))
  }
})

test_that("noise_acvf agrees with the autocovariance fit_noise uses", {
  A <- ar1(300L, 0.55, nvox = 20L, seed = 41)
  pl <- fit_noise(A, method = "ar", p = "auto", p_max = 6L, pooling = "global")
  ac <- noise_acvf(A, max_lag = 6L)
  expect_equal(ac$acvf[[1L]], pl$gamma[[1L]], tolerance = 1e-12)

  # And with runs and censoring in play.
  runs <- rep(1:2, each = 150L)
  cens <- c(20:24, 200:203)
  pl2 <- fit_noise(A, runs = runs, censor = cens, method = "ar", p = "auto",
                   p_max = 6L, pooling = "run")
  ac2 <- noise_acvf(A, runs = runs, censor = cens, max_lag = 6L, pooling = "run")
  for (i in seq_along(pl2$gamma)) {
    expect_equal(ac2$acvf[[i]], pl2$gamma[[i]], tolerance = 1e-12,
                 info = paste("run", i))
  }
})

test_that("noise_acvf and fit_noise use the same correction lag budget", {
  set.seed(43)
  n <- 120L
  X <- cbind(1, poly(seq_len(n), 3), matrix(rnorm(n * 6L), n, 6L))
  M <- diag(n) - X %*% solve(crossprod(X), t(X))
  R <- M %*% ar1(n, 0.5, nvox = 8L, seed = 44)

  plan <- fit_noise(R, method = "ar", p = 2L, p_max = 6L, pooling = "global",
                    design = X, correction_max_lag = 25L)
  ac <- noise_acvf(R, max_lag = 6L, design = X, correction_max_lag = 25L)
  expect_equal(ac$acvf[[1L]], plan$gamma[[1L]], tolerance = 1e-12)
  expect_true(ac$corrected)
})

test_that("corrected ACVF recovers downstream contrast variance", {
  # End-to-end oracle for the documented use case: gamma plus a design should
  # recover the sampling variance of a linear-model contrast, not merely the AR
  # coefficient or the shape of the ACVF.
  n <- 100L
  phi <- 0.5
  X <- cbind(1, sin(seq_len(n) * 2 * pi / n))
  contrast <- c(0, 1)
  XtX_inv <- solve(crossprod(X))
  weights <- as.numeric(X %*% XtX_inv %*% contrast)

  E <- ar1(n, phi, nvox = 500L, seed = 45)
  beta_hat <- as.numeric(crossprod(weights, E))
  M <- diag(n) - X %*% XtX_inv %*% t(X)
  resid <- M %*% E

  ac <- noise_acvf(resid, max_lag = 25L, design = X,
                   correction_max_lag = 25L)
  g_hat <- c(ac$acvf[[1L]], rep(0, n - length(ac$acvf[[1L]])))
  predicted <- as.numeric(crossprod(weights, stats::toeplitz(g_hat) %*% weights))

  g_true <- phi^(0:(n - 1L)) / (1 - phi^2)
  analytic <- as.numeric(crossprod(weights, stats::toeplitz(g_true) %*% weights))
  empirical <- stats::var(beta_hat)

  expect_true(ac$corrected)
  expect_equal(predicted / analytic, 1, tolerance = 0.08)
  expect_equal(empirical / analytic, 1, tolerance = 0.12)
})

test_that("parcel pooling is keyed by parcel and refuses bad labels", {
  A <- ar1(300L, 0.5, nvox = 16L, seed = 51)
  pf <- rep(1:4, length.out = 16L)
  out <- noise_acvf(A, max_lag = 4L, pooling = "parcel", parcels = pf)
  expect_named(out$acvf, c("1", "2", "3", "4"))
  expect_length(out$pairs, 4L)
  for (g in out$acvf) expect_equal(g[2] / g[1], 0.5, tolerance = 0.12)

  expect_error(noise_acvf(A, max_lag = 4L, pooling = "parcel"),
               "'parcels' is required")
  expect_error(noise_acvf(A, max_lag = 4L, pooling = "parcel",
                          parcels = rep(c("a", "b"), length.out = 16L)),
               "must be integer, numeric, or factor labels")
})

test_that("noise_acvf can correct the residual bias", {
  skip_on_cran()
  set.seed(61)
  n <- 300L
  X <- cbind(1, poly(seq_len(n), 3), matrix(rnorm(n * 18), n, 18))
  M <- diag(n) - X %*% solve(crossprod(X), t(X))
  R <- matrix(0, n, 40L)
  for (v in 1:40) R[, v] <- M %*% as.numeric(stats::arima.sim(list(ar = 0.5), n))

  raw <- noise_acvf(R, max_lag = 25L)
  cor <- noise_acvf(R, max_lag = 25L, design = X)
  expect_false(raw$corrected)
  expect_true(cor$corrected)
  expect_lt(abs(cor$acvf[[1L]][2] / cor$acvf[[1L]][1] - 0.5),
            abs(raw$acvf[[1L]][2] / raw$acvf[[1L]][1] - 0.5))

  expect_error(noise_acvf(R, pooling = "parcel", parcels = rep(1:4, length.out = 40L),
                          design = X),
               "not yet supported")
})

test_that("noise_acvf reports whether correction was actually applied", {
  set.seed(62)
  n <- 100L
  X <- cbind(1, poly(seq_len(n), 3))
  R <- matrix(rnorm(n * 5L), n, 5L)

  expect_warning(
    out <- noise_acvf(R, max_lag = 5L, design = X, correction_max_lag = 90L),
    "ill-conditioned"
  )
  raw <- noise_acvf(R, max_lag = 5L)
  expect_false(out$corrected)
  expect_equal(out$acvf, raw$acvf, tolerance = 1e-12)
})

test_that("noise_acvf rejects non-finite data and malformed lag budgets", {
  A <- ar1(20L, 0.4, nvox = 2L, seed = 63)
  zero <- noise_acvf(A, max_lag = 0L)
  expect_length(zero$acvf[[1L]], 1L)
  expect_length(zero$pairs[[1L]], 1L)

  bad <- A
  bad[3L, 1L] <- Inf
  expect_error(noise_acvf(bad), "Inf")
  expect_error(fit_noise(bad), "Inf")
  expect_error(noise_acvf(A, max_lag = NA_real_), "max_lag")
  expect_error(noise_acvf(A, max_lag = -1L), "max_lag")
  expect_error(noise_acvf(A, design = cbind(1, seq_len(20L)),
                          correction_max_lag = 0L), "correction_max_lag")
  expect_error(fit_noise(A, design = cbind(1, seq_len(20L)),
                         correction_max_lag = 0L), "correction_max_lag")

  # Budgets past .Machine$integer.max used to overflow as.integer() to NA and
  # die inside vector allocation; they are clamped to what the data supports.
  big <- noise_acvf(A, max_lag = 1e10)
  expect_identical(big$acvf, noise_acvf(A, max_lag = 19L)$acvf)
})

test_that("print method summarises without dumping the object", {
  A <- ar1(200L, 0.4, nvox = 6L, seed = 71)
  out <- noise_acvf(A, runs = rep(1:2, each = 100L), max_lag = 5L, pooling = "run")
  txt <- paste(capture.output(print(out)), collapse = "\n")
  expect_match(txt, "fmriAR noise autocovariance")
  expect_match(txt, "pooling:\\s+run")
  expect_match(txt, "gamma_0")
  capture.output(expect_invisible(print(out)))
})
