
# Helper utilities for fmriAR tests

simulate_arma11 <- function(n, phi = 0.6, theta = 0.4, sigma = 1.0, burnin = 200, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  e <- rnorm(n + burnin, 0, sigma)
  y <- numeric(n + burnin)
  for (t in seq_len(n + burnin)) {
    y[t] <- (if (t > 1) phi * y[t - 1] else 0) + e[t] + (if (t > 1) theta * e[t - 1] else 0)
  }
  list(y = y[(burnin + 1):(burnin + n)], e = e[(burnin + 1):(burnin + n)])
}

simulate_ar1_runs <- function(S, L, phi = 0.8, sigma = 1.0, seed = 1L) {
  set.seed(seed)
  n <- S * L
  y <- numeric(n)
  starts <- seq.int(1L, n, by = L)
  idx <- 1L
  for (s in seq_len(S)) {
    # Stationary initial state for an AR(1):
    y_prev <- rnorm(1L, mean = 0, sd = sigma / sqrt(1 - phi^2))
    for (t in seq_len(L)) {
      e_t <- rnorm(1L, 0, sigma)
      y[idx] <- phi * y_prev + e_t
      y_prev <- y[idx]
      idx <- idx + 1L
    }
  }
  list(y = y, run_starts0 = as.integer(starts - 1L))
}

# Roots (moduli) for AR and MA polynomials used by fmriAR:
# AR: A(z) = 1 - sum_k phi_k z^k
# MA: M(z) = 1 + sum_j theta_j z^j
ar_root_moduli <- function(phi) {
  if (length(phi) == 0L) return(numeric())
  Mod(polyroot(c(1, -as.numeric(phi))))
}
ma_root_moduli <- function(theta) {
  if (length(theta) == 0L) return(numeric())
  Mod(polyroot(c(1, as.numeric(theta))))
}

do_whiten_Y <- function(y, phi, theta = numeric(), run_starts0 = 0L, exact_first_ar1 = FALSE, parallel = TRUE, X_cols = 2L) {
  stopifnot(is.numeric(y), is.numeric(phi), is.numeric(theta))
  Y <- matrix(y, ncol = 1L)
  X <- matrix(rnorm(length(y) * X_cols), ncol = X_cols)
  # Call the Rcpp whitening (exported by the package)
  out <- fmriAR:::arma_whiten_inplace(Y, X, phi = phi, theta = theta, run_starts = as.integer(run_starts0),
                                      exact_first_ar1 = exact_first_ar1, parallel = parallel)
  as.numeric(out$Y[, 1L])
}
