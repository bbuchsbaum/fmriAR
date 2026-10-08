test_that("Whitening is deterministic across thread counts and parallel flag", {
  set.seed(1)
  n  <- 512L
  v  <- 64L
  Y  <- matrix(rnorm(n * v), nrow = n, ncol = v)
  X  <- matrix(rnorm(n * 3L), nrow = n, ncol = 3L)
  phi <- c(0.7, -0.2)   # stable AR(2)
  theta <- 0.5          # MA(1)
  rs <- 0L

  # arma_whiten_inplace() writes through its arguments, so every call gets
  # fresh copies. Passing Y itself made all outputs alias one buffer that was
  # whitened repeatedly, and the equality checks below compared it to itself.
  run <- function(parallel, nt) {
    fmriAR:::arma_whiten_inplace(Y + 0, X + 0, phi = phi, theta = theta, run_starts = rs,
                                 exact_first = FALSE, parallel = parallel,
                                 n_threads = nt)
  }
  out1 <- run(FALSE, 0L)
  expect_false(isTRUE(all.equal(out1$Y, Y)))
  for (nt in c(1L, 2L, 8L)) {
    outp <- run(TRUE, nt)
    expect_identical(out1$Y, outp$Y)
    expect_identical(out1$X, outp$X)
  }
})
