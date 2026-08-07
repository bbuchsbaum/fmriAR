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

test_that("censoring segments the series and is reported", {
  A <- ar1(300L, 0.5, nvox = 10L, seed = 31)
  cens <- c(50:59, 200:204)
  out <- noise_acvf(A, censor = cens, max_lag = 5L)

  expect_equal(out$n_segments[[1L]], 3L)
  expect_equal(sum(out$segment_lengths[[1L]]), 300L - length(cens))
  # Lag-0 pairs are exactly the surviving frames; higher lags lose one pair per
  # segment per lag because products cannot span a gap.
  expect_equal(out$pairs[[1L]][1], 300 - length(cens))
  expect_lt(out$pairs[[1L]][2], out$pairs[[1L]][1])
  expect_equal(out$acvf[[1L]][2] / out$acvf[[1L]][1], 0.5, tolerance = 0.1)
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

test_that("print method summarises without dumping the object", {
  A <- ar1(200L, 0.4, nvox = 6L, seed = 71)
  out <- noise_acvf(A, runs = rep(1:2, each = 100L), max_lag = 5L, pooling = "run")
  txt <- paste(capture.output(print(out)), collapse = "\n")
  expect_match(txt, "fmriAR noise autocovariance")
  expect_match(txt, "pooling:\\s+run")
  expect_match(txt, "gamma_0")
  capture.output(expect_invisible(print(out)))
})
