test_that("whitening reduces low-lag autocorrelation when voxel means cancel", {
  set.seed(827)

  make_run <- function(n, phi) {
    burn <- 100L
    e <- rnorm(n + burn)
    y <- numeric(n + burn)
    for (t in 2:(n + burn)) {
      y[t] <- phi * y[t - 1L] + e[t]
    }
    core <- y[(burn + 1L):(burn + n)]
    z <- rnorm(n, sd = 0.05)
    cbind(core + z, -(core + z))
  }

  n_per_run <- 100L
  resid_run1 <- make_run(n_per_run, phi = 0.6)
  resid_run2 <- make_run(n_per_run, phi = 0.55)
  resid <- rbind(resid_run1, resid_run2)
  runs <- rep(1:2, each = n_per_run)

  plan <- fit_noise(resid,
                    runs = runs,
                    method = "ar",
                    p = "auto",
                    p_max = 4L,
                    pooling = "run")

  expect_gt(plan$order[["p"]], 0L)
  phi_lengths <- vapply(plan$phi, length, 0L)
  expect_true(all(phi_lengths > 0L))
  expect_gt(mean(abs(unlist(plan$phi))), 0.2)

  lag_summary <- function(mat) {
    apply(mat, 2, function(y) {
      ac <- stats::acf(y, plot = FALSE, lag.max = 5)$acf[-1L]
      mean(abs(ac))
    })
  }

  X_zero <- matrix(0, nrow(resid), 1L)
  whitened <- whiten_apply(plan, X_zero, resid, runs = runs, parallel = FALSE)
  innov <- whitened$Y

  raw_stats <- lag_summary(resid)
  innov_stats <- lag_summary(innov)

  expect_lt(mean(innov_stats), 0.5 * mean(raw_stats))
})
