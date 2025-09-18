
# helper-ms-diagnostics.R
# Utilities for multiscale diagnostics tests

# Ensure fmriAR is available
skip_if_no_fmriAR <- function() {
  if (!requireNamespace("fmriAR", quietly = TRUE)) {
    testthat::skip("fmriAR not installed")
  }
  ns <- asNamespace("fmriAR")
  needed <- c("fit_noise", "whiten_apply")
  missing <- !vapply(needed, function(fn) exists(fn, envir = ns, inherits = FALSE), logical(1))
  if (any(missing)) {
    testthat::skip(paste("Missing functions in fmriAR:", paste(needed[missing], collapse = ", ")))
  }
}

# Enforce AR stationarity (reflect roots outside the unit circle)
enforce_ar_stationary <- function(phi) {
  p <- length(phi)
  if (p == 0L) return(phi)
  r <- polyroot(c(1, -phi))
  changed <- FALSE
  for (i in seq_along(r)) {
    if (Mod(r[i]) <= 1) {
      r[i] <- (1/conj(r[i])) * 1.05
      changed <- TRUE
    }
  }
  if (changed) {
    cf <- poly(r) # leading 1
    phi <- -Re(cf[-1L])
  }
  as.numeric(phi)
}

# Build a tidy multiscale hierarchy (coarse -> medium -> fine), and voxel labels
make_hierarchy <- function(n_coarse = 4L, medium_per_coarse = 3L, fine_per_medium = 3L, vox_per_fine = 4L) {
  n_medium <- n_coarse * medium_per_coarse
  n_fine   <- n_medium * fine_per_medium
  # parent maps (1-based)
  parent_medium <- rep(seq_len(n_medium), each = fine_per_medium)
  parent_coarse <- rep(rep(seq_len(n_coarse), each = medium_per_coarse), each = fine_per_medium)
  # voxels per fine
  v <- n_fine * vox_per_fine
  parcels_fine <- rep(seq_len(n_fine), each = vox_per_fine)
  list(
    n_coarse = n_coarse,
    n_medium = n_medium,
    n_fine   = n_fine,
    parent_medium = parent_medium,
    parent_coarse = parent_coarse,
    parcels_fine = parcels_fine,
    vox_per_fine = vox_per_fine
  )
}

# Simulate hierarchical AR(2) noise with small parcel-level deviations around parents
# Returns Train/Test Y matrices and design matrices (X), plus run_starts (0-based)
simulate_hier_ar2 <- function(h, n_train_per_run = 150L, n_test = 150L, runs_train = 2L, seed = 1L) {
  set.seed(seed)
  # True AR(2) at coarse level, safe/stable
  phi_coarse <- replicate(h$n_coarse, enforce_ar_stationary(c(0.55, -0.18)))
  # Medium inherit + small jitter
  phi_medium <- matrix(NA_real_, 2L, h$n_medium)
  for (m in seq_len(h$n_medium)) {
    cix <- ceiling(m / (h$n_medium / h$n_coarse)) # parent coarse index
    base <- phi_coarse[, cix]
    phi_medium[, m] <- enforce_ar_stationary(base + rnorm(2L, 0, 0.04))
  }
  # Fine inherit + a bit more jitter
  phi_fine <- matrix(NA_real_, 2L, h$n_fine)
  for (f in seq_len(h$n_fine)) {
    mix <- h$parent_medium[f]
    base <- phi_medium[, mix]
    phi_fine[, f] <- enforce_ar_stationary(base + rnorm(2L, 0, 0.05))
  }

  # Simulate runs for TRAIN (runs_train) and one TEST run
  v <- length(h$parcels_fine)
  # Map fine parcel -> voxel phi
  phi_vox <- phi_fine[, h$parcels_fine]

  sim_run <- function(n) {
    Y <- matrix(0.0, n, v)
    for (j in seq_len(v)) {
      phi <- phi_vox[, j]
      e <- rnorm(n + 5L) # small burn-in margin
      y <- numeric(n + 5L)
      for (t in seq_len(n + 5L)) {
        y[t] <- (if (t > 1) phi[1] * y[t - 1] else 0) + (if (t > 2) phi[2] * y[t - 2] else 0) + e[t]
      }
      Y[, j] <- y[(5L + 1):(5L + n)]
    }
    Y
  }

  # Train concatenated over runs
  Y_train <- do.call(rbind, replicate(runs_train, sim_run(n_train_per_run), simplify = FALSE))
  Y_test  <- sim_run(n_test)

  # Simple design (intercept) to satisfy interfaces
  X_train <- matrix(1.0, nrow(Y_train), 1L)
  X_test  <- matrix(1.0, nrow(Y_test), 1L)

  # 0-based run_starts vectors
  rs_train0 <- as.integer(seq.int(1L, nrow(Y_train), by = n_train_per_run) - 1L)
  rs_test0  <- as.integer(0L)

  list(
    Y_train = Y_train, X_train = X_train, run_starts_train0 = rs_train0,
    Y_test = Y_test,   X_test  = X_test,  run_starts_test0  = rs_test0,
    parcels_fine = h$parcels_fine
  )
}

# Ljung-Box p-values and KS distance to Uniform(0,1)
lb_pvals <- function(E, lag = 10L) {
  apply(E, 2L, function(e) stats::Box.test(e, lag = lag, type = "Ljung-Box")$p.value)
}
ks_to_uniform <- function(p) {
  p <- sort(p[is.finite(p) & p >= 0 & p <= 1])
  if (!length(p)) return(NA_real_)
  n <- length(p)
  max(abs(p - (seq_len(n)/n)))
}

# Negative log-likelihood per series under Gaussian with estimated sigma^2
series_nll <- function(E) {
  # MLE sigma^2 per series
  s2 <- colMeans(E^2)
  n <- nrow(E)
  0.5 * n * (log(2 * pi) + log(s2) + 1)
}

# Extract phi_by_parcel (p x K_fine) from a plan (robust to internal changes)
plan_phi_by_parcel <- function(plan) {
  # Try compat$plan_info
  pi <- try(fmriAR:::compat$plan_info(plan), silent = TRUE)
  if (!inherits(pi, "try-error") && !is.null(pi$phi_by_parcel)) return(pi$phi_by_parcel)
  # Try plan_info
  pi2 <- try(fmriAR:::plan_info(plan), silent = TRUE)
  if (!inherits(pi2, "try-error") && !is.null(pi2$phi_by_parcel)) return(pi2$phi_by_parcel)
  # Try direct field
  if (is.list(plan) && !is.null(plan$phi_by_parcel)) return(plan$phi_by_parcel)
  stop("Cannot extract phi_by_parcel from plan")
}
