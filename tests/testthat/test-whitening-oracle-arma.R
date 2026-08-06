manual_arma_whiten <- function(M, phi, theta, run_starts) {
  M <- as.matrix(M)
  out <- matrix(0, nrow = nrow(M), ncol = ncol(M))
  p <- length(phi)
  q <- length(theta)
  seg_beg <- as.integer(run_starts) + 1L
  seg_end <- c(seg_beg[-1L] - 1L, nrow(M))

  for (col in seq_len(ncol(M))) {
    for (s in seq_along(seg_beg)) {
      prev_in <- rep(0, max(1L, p))
      prev_out <- rep(0, max(1L, q))

      for (t in seg_beg[s]:seg_end[s]) {
        orig <- M[t, col]
        st <- orig
        if (p > 0L) {
          for (k in seq_len(p)) st <- st - phi[k] * prev_in[k]
        }
        zt <- st
        if (q > 0L) {
          for (j in seq_len(q)) zt <- zt - theta[j] * prev_out[j]
        }

        out[t, col] <- zt

        if (p > 0L) {
          if (p > 1L) prev_in[2:p] <- prev_in[1:(p - 1L)]
          prev_in[1] <- orig
        }
        if (q > 0L) {
          if (q > 1L) prev_out[2:q] <- prev_out[1:(q - 1L)]
          prev_out[1] <- zt
        }
      }
    }
  }

  out
}

test_that("arma_whiten_inplace matches manual ARMA recursion across segments", {
  phi <- c(0.4, -0.15)
  theta <- 0.25
  run_starts <- c(0L, 3L)

  Y <- matrix(c(
    1,  2,
    3,  5,
    8, 13,
    2,  4,
    6,  8,
    1,  3
  ), nrow = 6, byrow = TRUE)
  X <- matrix(c(
    1, 0, 2,
    1, 1, 3,
    1, 2, 5,
    1, 3, 1,
    1, 4, 0,
    1, 5, 2
  ), nrow = 6, byrow = TRUE)
  Y_ref <- Y
  X_ref <- X
  manual_Y <- manual_arma_whiten(Y_ref, phi, theta, run_starts)
  manual_X <- manual_arma_whiten(X_ref, phi, theta, run_starts)

  out <- arma_whiten_inplace(
    Y = Y,
    X = X,
    phi = phi,
    theta = theta,
    run_starts = run_starts,
    exact_first_ar1 = FALSE,
    parallel = FALSE
  )

  expect_equal(out$Y, manual_Y, tolerance = 1e-12)
  expect_equal(out$X, manual_X, tolerance = 1e-12)
})

test_that("whiten_apply matches manual ARMA whitening with run and censor resets", {
  phi <- c(0.35, -0.1)
  theta <- 0.2
  runs <- c(rep(1L, 4L), rep(2L, 4L))
  censor <- c(2L, 7L)

  Y <- matrix(c(
    1,  0,
    2,  1,
    4,  1,
    7,  2,
    1,  1,
    3,  2,
    5,  3,
    8,  5
  ), nrow = 8, byrow = TRUE)
  X <- cbind(1, seq_len(8), c(0, 1, 0, 1, 1, 0, 1, 0))

  plan <- fmriAR:::new_whiten_plan(
    phi = list(phi, phi),
    theta = list(theta, theta),
    order = c(p = 2L, q = 1L),
    runs = runs,
    exact_first = FALSE,
    method = "arma",
    pooling = "run"
  )

  idx1 <- which(runs == 1L)
  idx2 <- which(runs == 2L)
  manual_X <- rbind(
    manual_arma_whiten(X[idx1, , drop = FALSE], phi, theta, fmriAR:::`.sub_run_starts`(length(idx1), 2L)),
    manual_arma_whiten(X[idx2, , drop = FALSE], phi, theta, fmriAR:::`.sub_run_starts`(length(idx2), 3L))
  )
  manual_Y <- rbind(
    manual_arma_whiten(Y[idx1, , drop = FALSE], phi, theta, fmriAR:::`.sub_run_starts`(length(idx1), 2L)),
    manual_arma_whiten(Y[idx2, , drop = FALSE], phi, theta, fmriAR:::`.sub_run_starts`(length(idx2), 3L))
  )
  out <- whiten_apply(plan, X, Y, runs = runs, censor = censor, parallel = FALSE)

  expect_equal(out$X, manual_X, tolerance = 1e-12)
  expect_equal(out$Y, manual_Y, tolerance = 1e-12)
})
