
test_that("Whitening is deterministic across thread counts and parallel flag", {
  skip_if_not_installed("fmriAR")
  set.seed(1)
  n  <- 512L
  v  <- 64L
  Y  <- matrix(rnorm(n * v), nrow = n, ncol = v)
  X  <- matrix(rnorm(n * 3L), nrow = n, ncol = 3L)
  phi <- c(0.7, -0.2)   # stable AR(2)
  theta <- 0.5          # MA(1)
  rs <- 0L

  # Serial
  out1 <- fmriAR:::arma_whiten_inplace(Y, X, phi = phi, theta = theta, run_starts = rs,
                                       exact_first_ar1 = FALSE, parallel = FALSE)

  # Parallel (threads may be ignored on systems without OpenMP; results should still be identical)
  old <- Sys.getenv("OMP_NUM_THREADS", unset = NA_character_)
  on.exit(if (!is.na(old)) Sys.setenv(OMP_NUM_THREADS = old), add = TRUE)

  Sys.setenv(OMP_NUM_THREADS = "1")
  out2 <- fmriAR:::arma_whiten_inplace(Y, X, phi = phi, theta = theta, run_starts = rs,
                                       exact_first_ar1 = FALSE, parallel = TRUE)
  Sys.setenv(OMP_NUM_THREADS = "2")
  out3 <- fmriAR:::arma_whiten_inplace(Y, X, phi = phi, theta = theta, run_starts = rs,
                                       exact_first_ar1 = FALSE, parallel = TRUE)
  Sys.setenv(OMP_NUM_THREADS = "8")
  out4 <- fmriAR:::arma_whiten_inplace(Y, X, phi = phi, theta = theta, run_starts = rs,
                                       exact_first_ar1 = FALSE, parallel = TRUE)

  expect_equal(out1$Y, out2$Y, tolerance = 0)  # exact equality
  expect_equal(out1$X, out2$X, tolerance = 0)
  expect_equal(out1$Y, out3$Y, tolerance = 0)
  expect_equal(out1$X, out3$X, tolerance = 0)
  expect_equal(out1$Y, out4$Y, tolerance = 0)
  expect_equal(out1$X, out4$X, tolerance = 0)
})
