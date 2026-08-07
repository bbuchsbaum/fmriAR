# Internal helpers -------------------------------------------------------------

.sub_run_starts <- function(n_run, censor_idx_rel = integer()) {
  starts <- 1L
  if (length(censor_idx_rel)) {
    add <- censor_idx_rel + 1L
    add <- add[add <= n_run]
    starts <- sort(unique(c(starts, add)))
  }
  as.integer(starts - 1L)
}

.split_runs <- function(runs) {
  if (is.null(runs)) return(list(seq_along(integer(0))))
  runs <- as.integer(runs)
  split(seq_along(runs), runs)
}

new_whiten_plan <- function(phi, theta, order, runs, exact_first, method, pooling,
                            parcels = NULL, parcel_ids = NULL,
                            phi_by_parcel = NULL, theta_by_parcel = NULL,
                            censor = NULL,
                            gamma = NULL, sigma2 = NULL,
                            gamma_by_parcel = NULL, sigma2_by_parcel = NULL) {
  structure(
    list(
      phi = phi,
      theta = theta,
      order = order,
      runs = runs,
      exact_first = exact_first,
      method = method,
      pooling = pooling,
      parcels = parcels,
      parcel_ids = parcel_ids,
      phi_by_parcel = phi_by_parcel,
      theta_by_parcel = theta_by_parcel,
      censor = censor,
      # Noise scale and shape. Without these a plan describes the correlation
      # structure but not its magnitude, so a consumer cannot recover the
      # covariance the plan implies -- two datasets differing 100-fold in
      # variance produced identical plans.
      gamma = gamma,
      sigma2 = sigma2,
      gamma_by_parcel = gamma_by_parcel,
      sigma2_by_parcel = sigma2_by_parcel
    ),
    class = "fmriAR_plan"
  )
}

.arma_innovations <- function(y, phi, theta) {
  Y <- matrix(as.numeric(y), ncol = 1L)
  X <- matrix(0, nrow = length(y), ncol = 1L)
  out <- arma_whiten_inplace(Y, X, phi, theta, run_starts = 0L,
                             exact_first_ar1 = FALSE, parallel = FALSE)
  drop(out$Y)
}

.run_avg_acvf <- function(mat, max_lag) {
  if (!is.matrix(mat)) mat <- as.matrix(mat)
  run_avg_acvf_cpp(mat, as.integer(max_lag))
}

# run_avg_acvf_cpp normalizes lag h by the pair count (n - h), the unbiased
# estimator, which is not positive semi-definite. Rescaling each lag by
# (n - h) / n converts it to the biased estimator, whose Toeplitz matrix is PSD
# by construction. Yule-Walker on a non-PSD acvf can return explosive
# coefficients, so every path that feeds yw_from_acvf_fast() must use this form.
.run_avg_acvf_psd <- function(mat, max_lag) {
  if (!is.matrix(mat)) mat <- as.matrix(mat)
  max_lag <- as.integer(max_lag)
  gamma <- run_avg_acvf_cpp(mat, max_lag)
  n <- nrow(mat)
  if (n <= 0L) return(gamma)
  lags <- seq_along(gamma) - 1L
  gamma * pmax(n - lags, 0) / n
}

# Pooled autocovariance over valid, possibly fragmented data.
#
# Two groupings, deliberately distinct:
#
#   center_id -- the unit the mean is estimated over. The mean is a per-run
#     nuisance parameter (offsets and drift differ between runs), so it must be
#     estimated per run, NOT per segment and NOT across runs. Centering each
#     scrubbing fragment on its own mean removes the very signal being measured:
#     a 2-frame fragment so centered has lag-1 product -(x1-x2)^2/4, negative by
#     construction, so rho1 = -1 exactly. Pooling those drives phi toward and
#     past zero as censoring rises, and is what makes the pooled acvf non-PSD.
#     Centering across runs instead leaks between-run offsets into every lag and
#     pushes rho toward 1.
#
#   seg_id -- the unit lag products are confined to, so no product spans a run
#     boundary or a scrubbed frame.
#
# Returns per-lag product sums (averaged over columns) and their pair counts.
.pooled_acvf_segments <- function(mat, seg_id, max_lag, center_id = NULL) {
  if (!is.matrix(mat)) mat <- as.matrix(mat)
  max_lag <- max(0L, as.integer(max_lag))
  nv <- nrow(mat)
  nc <- ncol(mat)
  num <- numeric(max_lag + 1L)
  pairs <- numeric(max_lag + 1L)
  if (nv == 0L || nc == 0L) return(list(num = num, pairs = pairs))

  if (is.null(center_id)) {
    mat <- mat - rep(colMeans(mat), each = nv)
  } else {
    center_id <- as.integer(center_id)
    grp <- match(center_id, sort(unique(center_id)))
    mu <- rowsum(mat, grp, reorder = TRUE) / as.numeric(table(grp))
    mat <- mat - mu[grp, , drop = FALSE]
  }
  num[1L] <- sum(mat * mat) / nc
  pairs[1L] <- nv
  for (lg in seq_len(max_lag)) {
    if (nv <= lg) break
    hi <- seq.int(lg + 1L, nv)
    lo <- seq.int(1L, nv - lg)
    ok <- seg_id[hi] == seg_id[lo]
    if (!any(ok)) next
    num[lg + 1L] <- sum(mat[hi[ok], , drop = FALSE] * mat[lo[ok], , drop = FALSE]) / nc
    pairs[lg + 1L] <- sum(ok)
  }
  list(num = num, pairs = pairs)
}

# Highest lag with any contributing pairs; -1 when the pooled acvf is empty.
.acvf_max_lag <- function(pooled) {
  usable <- which(pooled$pairs > 0)
  if (!length(usable)) return(-1L)
  max(usable) - 1L
}

# Positive definite with a strict relative margin, not merely non-negative.
# Sitting on the cone boundary makes a reflection coefficient exactly +/-1, which
# drives the Levinson prediction error to the 1e-12 floor; BIC then reads that as
# a perfect fit and always selects p_max.
.acvf_is_psd <- function(gamma, tol = 1e-6) {
  if (!length(gamma) || !all(is.finite(gamma)) || gamma[1] <= 0) return(FALSE)
  if (length(gamma) < 2L) return(TRUE)
  ev <- tryCatch(
    eigen(stats::toeplitz(gamma), symmetric = TRUE, only.values = TRUE)$values,
    error = function(e) NA_real_
  )
  if (anyNA(ev)) return(FALSE)
  min(ev) >= tol * gamma[1]
}

# Prefer the pair-count (unbiased) normalization, which is not attenuated when
# censoring thins the pairs available at each lag. Where that is not positive
# definite, shrink the non-zero lags toward white noise only as far as needed,
# since Yule-Walker on a non-PD acvf returns explosive coefficients.
.acvf_from_pooled <- function(pooled, order = NULL, tol = 1e-6) {
  num <- pooled$num
  pairs <- pooled$pairs
  usable <- which(pairs > 0)
  if (!length(usable)) return(numeric(0))
  avail <- max(usable)
  # Correct only the Toeplitz matrix actually being inverted. Fragmented data
  # leaves few pairs at long lags, so those entries are noisy and frequently
  # break PSD; shrinking the whole vector to repair them attenuates the short
  # lags that carry the signal (phi = 0.95 came back as 0.85 at 40% censoring).
  keep <- if (is.null(order)) seq_len(max(1L, avail))
          else seq_len(max(1L, min(as.integer(order) + 1L, avail)))
  num <- num[keep]
  pairs <- pairs[keep]
  if (!length(num) || pairs[1L] <= 0) return(numeric(0))

  g_unb <- ifelse(pairs > 0, num / pairs, 0)
  if (!is.finite(g_unb[1]) || g_unb[1] <= 0) return(numeric(0))
  if (.acvf_is_psd(g_unb, tol)) return(g_unb)

  # Shrink the non-zero lags toward white noise until the Toeplitz matrix is
  # positive definite. Shrinking toward the common-divisor form instead is not
  # guaranteed to terminate: that form can itself be non-PSD, leaving both ends
  # of the search invalid. At lambda = 0 the matrix is gamma[1] * I, so a valid
  # point always exists and the bisection cannot fail.
  lo <- 0
  hi <- 1
  for (i in seq_len(50L)) {
    mid <- (lo + hi) / 2
    if (.acvf_is_psd(c(g_unb[1L], g_unb[-1L] * mid), tol)) lo <- mid else hi <- mid
  }
  c(g_unb[1L], g_unb[-1L] * lo)
}

# Innovation variance implied by an autocovariance and an AR coefficient vector:
# sigma2 = gamma_0 - sum_k phi_k gamma_k. Deriving it here rather than carrying
# the value from estimation keeps sigma2 consistent with the phi actually stored
# on the plan, which differs after stationarity clamping or multiscale pooling.
.sigma2_from_gamma_phi <- function(gamma, phi) {
  if (!length(gamma) || !is.finite(gamma[1]) || gamma[1] <= 0) return(NA_real_)
  if (!length(phi)) return(gamma[1])
  # Dropping the non-finite entries would slide every later coefficient onto the
  # wrong lag, so a phi that is not wholly finite yields no answer at all.
  if (!all(is.finite(phi))) return(NA_real_)
  # The sum needs gamma out to lag p. Global pooling truncates gamma to the
  # shortest run, so a single heavily censored run can leave it shorter than
  # that, and the terms lost are exactly the ones that pull sigma2 below
  # gamma_0. A partial sum therefore overstates the innovation variance -- at
  # the limit, a length-1 gamma returns gamma_0 and declares white noise while
  # phi on the same plan says otherwise. Report the absence instead.
  if (length(gamma) - 1L < length(phi)) return(NA_real_)
  s2 <- gamma[1] - sum(phi * gamma[seq_along(phi) + 1L])
  # No stationary process has innovation variance above its total variance.
  # phi that was stationarity-clamped or multiscale-pooled need not match this
  # unit's own gamma, which let the ratio reach 1.06.
  min(max(s2, 1e-12), gamma[1])
}

# Indices of non-censored timepoints plus the 0-based starts of the contiguous
# segments they form. A segment breaks at a run boundary or wherever censoring
# removed a frame, so no lag product ever spans a discontinuity.
.valid_segments <- function(n, runs = NULL, censor = NULL) {
  n <- as.integer(n)
  valid <- rep(TRUE, n)
  if (!is.null(censor) && length(censor)) {
    censor <- as.integer(censor)
    censor <- censor[censor >= 1L & censor <= n]
    valid[censor] <- FALSE
  }
  idx <- which(valid)
  if (!length(idx)) {
    return(list(idx = integer(0), starts0 = integer(0), run_id = integer(0)))
  }
  r <- if (is.null(runs)) rep(1L, n) else as.integer(runs)
  brk <- c(TRUE, diff(idx) != 1L | r[idx[-1L]] != r[idx[-length(idx)]])
  list(idx = idx, starts0 = as.integer(which(brk) - 1L), run_id = r[idx])
}

# Order selection and fitting for a single series, aware of the segment
# structure it was drawn from. The autocovariance comes from segmented_acvf_cpp
# (per-segment centering, PSD-safe normalization) and BIC scores the
# Levinson-Durbin prediction error, matching the global/run path exactly.
.estimate_ar_series <- function(y, p_max, p = "auto", starts0 = 0L, center_id = NULL) {
  y <- as.numeric(y)
  n <- length(y)
  starts0 <- as.integer(starts0)
  if (!length(starts0)) starts0 <- 0L

  seg_len <- diff(c(starts0, n))
  p_cap <- min(as.integer(p_max), max(seg_len) - 1L, n - 1L)
  empty <- list(phi = numeric(0), order = c(p = 0L, q = 0L),
                gamma = numeric(0), sigma2 = NA_real_)
  if (p_cap < 1L) return(empty)

  seg_id <- cumsum(seq_len(n) %in% (starts0 + 1L))
  pooled <- .pooled_acvf_segments(matrix(y, ncol = 1L), seg_id, p_cap,
                                  center_id = center_id)
  gamma0 <- .acvf_from_pooled(pooled, order = 0L)
  if (!length(gamma0) || !is.finite(gamma0[1]) || gamma0[1] <= 0) return(empty)
  p_cap <- min(p_cap, .acvf_max_lag(pooled))
  if (p_cap < 1L) {
    empty$gamma <- gamma0
    empty$sigma2 <- gamma0[1]
    return(empty)
  }
  gamma_full <- .acvf_from_pooled(pooled, order = p_cap)

  # PSD-correct at the order being fitted, not at p_max.
  fit_order <- function(pp) {
    g <- .acvf_from_pooled(pooled, order = pp)
    yw <- yw_from_acvf_fast(g[seq_len(pp + 1L)], pp)
    list(phi = enforce_stationary_ar(yw$phi, 0.99),
         sigma2 = pmax(yw$sigma2, 1e-12))
  }

  if (!identical(p, "auto")) {
    pp <- min(as.integer(p), p_cap)
    if (pp <= 0L) {
      empty$gamma <- gamma_full
      empty$sigma2 <- gamma_full[1]
      return(empty)
    }
    f <- fit_order(pp)
    return(list(phi = f$phi, order = c(p = pp, q = 0L),
                gamma = gamma_full, sigma2 = f$sigma2))
  }

  # Cap the order BIC may select by the data available. Selection on a handful
  # of points will happily choose AR(8) from 11 observations, and the resulting
  # filter inflates variance instead of whitening. An explicitly requested order
  # is honoured as given; this bounds only the search.
  p_sel <- min(p_cap, floor(n / 5))
  n_log <- log(n)
  best <- list(bic = 2 * n * log(pmax(gamma0[1], 1e-12)) + n_log,
               phi = numeric(0), p = 0L, sigma2 = gamma0[1])
  for (pp in seq_len(max(0L, p_sel))) {
    f <- fit_order(pp)
    # enforce_stationary_ar() returns length 0 when it cannot produce a
    # stationary filter, which must not be recorded as an order-pp fit.
    if (!is.finite(f$sigma2) || length(f$phi) != pp || !all(is.finite(f$phi))) next
    bic <- 2 * n * log(f$sigma2) + (pp + 1L) * n_log
    if (is.finite(bic) && bic < best$bic) {
      best <- list(bic = bic, phi = f$phi, p = pp, sigma2 = f$sigma2)
    }
  }
  list(phi = best$phi, order = c(p = best$p, q = 0L),
       gamma = gamma_full, sigma2 = best$sigma2)
}

.full_run_starts <- function(runs, censor, n) {
  starts <- 1L
  if (!is.null(runs)) {
    runs <- as.integer(runs)
    starts <- sort(unique(c(starts, which(diff(runs) != 0L) + 1L)))
  }
  if (!is.null(censor) && length(censor)) {
    censor <- sort(unique(as.integer(censor)))
    extra <- censor + 1L
    extra <- extra[extra <= n]
    starts <- sort(unique(c(starts, extra)))
  }
  as.integer(starts - 1L)
}

.runs_from_starts0 <- function(run_starts0, n) {
  rs <- sort(unique(as.integer(run_starts0)))
  if (!length(rs) || rs[1] != 0L) stop("run_starts must include 0")
  if (rs[length(rs)] == n) rs <- rs[-length(rs)]
  if (!length(rs) || rs[length(rs)] >= n) stop("run_starts out of bounds")
  rs1 <- rs + 1L
  bounds <- c(rs1, n + 1L)
  runs <- integer(n)
  for (i in seq_along(rs1)) {
    runs[seq(rs1[i], bounds[i + 1L] - 1L)] <- i
  }
  runs
}

# Exported API -----------------------------------------------------------------

#' Fit an AR/ARMA noise model (run-aware) and return a whitening plan
#'
#' @param resid Numeric matrix (time x voxels) of residuals from an initial OLS fit.
#' @param Y Optional data matrix used to compute residuals when `resid` is omitted.
#' @param X Optional design matrix used with `Y` to compute residuals.
#' @param runs Optional integer vector of run identifiers.
#' @param censor Optional integer vector of 1-based timepoint indices to exclude from
#'   AR parameter estimation, or a logical vector of length `nrow(resid)` where `TRUE`
#'
#'   indicates censored timepoints. Censored frames (e.g., motion-corrupted) are excluded
#'   when computing autocorrelations. Each run's estimation uses only its own valid
#'   (non-censored) segments.
#' @param method Either "ar" or "arma".
#' @param p AR order (integer or "auto" if method == "ar").
#' @param q MA order (integer).
#' @param p_max Maximum AR order when `p = "auto"`.
#' @param exact_first Apply exact AR(1) scaling at segment starts ("ar1" or "none").
#' @param pooling Combine parameters across runs or parcels ("global", "run", "parcel").
#' @param parcels Integer vector (length = ncol(resid)) giving fine parcel memberships when `pooling = "parcel"`.
#' @param parcel_sets Optional named list with entries `coarse`, `medium`, `fine` of equal length specifying nested parcel labels for multi-scale pooling.
#' @param multiscale Multi-scale pooling mode when `parcel_sets` is supplied ("pacf_weighted" or "acvf_pooled"), or `TRUE/FALSE` to toggle pooling.
#' @param ms_mode Explicit multiscale mode when `multiscale` is logical.
#' @param p_target Target AR order for multi-scale pooling (defaults to `p_max`).
#' @param beta Size exponent for multi-scale weights (default 0.5).
#' @param hr_iter Number of Hannan--Rissanen refinement iterations for ARMA.
#' @param step1 Preliminary high-order AR fit method for HR ("burg" or "yw").
#' @param parallel Reserved for future parallel estimation (logical).
#' @return An object of class `fmriAR_plan` used by [whiten_apply()]. Besides the
#'   AR/MA coefficients the plan carries the noise scale and shape it was fitted
#'   from, so consumers can reconstruct the covariance it implies rather than
#'   only its correlation structure:
#'   \itemize{
#'     \item `gamma`: list of autocovariance vectors, one per pooling unit -- a
#'       single entry for `pooling = "global"`, one per run for
#'       `pooling = "run"`. Lags run 0 to the highest the data supported, which
#'       is governed by `p_max` and the run length rather than by `p`, so
#'       `fit_noise(p = 1, p_max = 6)` returns seven values, not two. Under
#'       global pooling every run is truncated to the shortest available length
#'       before averaging, since a zero-padded autocovariance is not a valid
#'       covariance.
#'     \item `sigma2`: list of innovation variances, matching `gamma`, derived
#'       as `gamma_0 - sum_k phi_k gamma_k` from the coefficients stored on the
#'       plan so the two are always mutually consistent. `NA` for
#'       `method = "arma"`, where no comparably cheap voxel-scale innovation
#'       variance is available, and `NA` whenever `gamma` does not reach lag
#'       `length(phi)` -- heavy censoring can truncate it that far, and a
#'       partial sum would overstate the innovation variance rather than
#'       report that it is unavailable.
#'     \item `gamma_by_parcel`, `sigma2_by_parcel`: the same quantities per
#'       parcel when `pooling = "parcel"`, keyed like `phi_by_parcel`.
#'   }
#'   For a run-stationary noise process with autocovariance `gamma`, the
#'   covariance of the data within a run is the Toeplitz matrix built from it,
#'   which is what makes design-specific variance calculations possible
#'   downstream without refitting.
#' @examples
#' # Generate example data with AR(1) structure
#' n_time <- 200
#' n_voxels <- 50
#' phi_true <- 0.5
#'
#' # Simulate residuals with AR(1) structure
#' resid <- matrix(0, n_time, n_voxels)
#' for (v in 1:n_voxels) {
#'   e <- rnorm(n_time)
#'   resid[1, v] <- e[1]
#'   for (t in 2:n_time) {
#'     resid[t, v] <- phi_true * resid[t-1, v] + e[t]
#'   }
#' }
#'
#' # Fit AR model
#' plan <- fit_noise(resid, method = "ar", p = 1)
#'
#' # With multiple runs
#' runs <- rep(1:2, each = 100)
#' plan_runs <- fit_noise(resid, runs = runs, method = "ar", pooling = "run")
#' @export
fit_noise <- function(resid = NULL,
                      Y = NULL,
                      X = NULL,
                      runs = NULL,
                      censor = NULL,
                      method = c("ar", "arma"),
                      p = "auto",
                      q = 0L,
                      p_max = 6L,
                      exact_first = c("ar1", "none"),
                      pooling = c("global", "run", "parcel"),
                      parcels = NULL,
                      parcel_sets = NULL,
                      multiscale = c("pacf_weighted", "acvf_pooled"),
                      ms_mode = NULL,
                      p_target = NULL,
                      beta = 0.5,
                      hr_iter = 0L,
                      step1 = c("burg", "yw"),
                      parallel = FALSE) {

  if (is.null(resid)) {
    if (!is.null(Y) && !is.null(X)) {
      if (!is.matrix(Y)) Y <- as.matrix(Y)
      if (!is.matrix(X)) X <- as.matrix(X)
      stopifnot(nrow(Y) == nrow(X))
      coef <- qr.solve(X, Y)
      resid <- Y - X %*% coef
    } else {
      stop("fit_noise: supply 'resid' or both 'Y' and 'X'")
    }
  }

  stopifnot(is.matrix(resid))
  method <- match.arg(method)
  exact_first <- match.arg(exact_first)
  pooling <- match.arg(pooling)
  step1 <- match.arg(step1)

  ms_modes <- c("pacf_weighted", "acvf_pooled")
  multiscale_mode <- NULL
  multiscale_explicit <- isTRUE(multiscale) || (!is.logical(multiscale) && length(multiscale) == 1L)
  if (is.logical(multiscale)) {
    if (isTRUE(multiscale)) {
      multiscale_mode <- if (is.null(ms_mode)) "pacf_weighted" else match.arg(ms_mode, ms_modes)
    }
  } else {
    multiscale_mode <- match.arg(multiscale, ms_modes)
  }
  if (!is.null(ms_mode) && (!is.logical(multiscale) || isTRUE(multiscale))) {
    multiscale_mode <- match.arg(ms_mode, ms_modes)
  }
  if (!multiscale_explicit && identical(pooling, "parcel") && !identical(p, "auto")) {
    multiscale_mode <- NULL
  }

  n <- nrow(resid)
  if (n < 10) stop("series too short")


  # Normalize censor input: convert logical to integer indices
  if (!is.null(censor)) {
    if (is.logical(censor)) {
      stopifnot(length(censor) == n)
      censor <- which(censor)
    }
    censor <- sort(unique(as.integer(censor)))
    censor <- censor[censor >= 1L & censor <= n]
    if (!length(censor)) censor <- NULL
  }

  Rsets <- if (is.null(runs)) list(seq_len(n)) else split(seq_len(n), as.integer(runs))

  # Split censor indices by run (relative to run start)
  censor_by_run <- lapply(Rsets, function(idx) integer(0L))
  if (!is.null(censor)) {
    for (ri in seq_along(Rsets)) {
      idx <- Rsets[[ri]]
      c_in <- intersect(censor, idx)
      if (length(c_in)) {
        censor_by_run[[ri]] <- as.integer(c_in - min(idx) + 1L)
      }
    }
  }

  run_mats <- lapply(Rsets, function(idx) resid[idx, , drop = FALSE])

  est_run <- function(mat, censor_rel = integer(0L)) {
    # Get valid (non-censored) segments for estimation
    # censor_rel contains 1-based indices within this run
    n_run <- nrow(mat)
    if (length(censor_rel)) {
      # Create mask of valid timepoints
      valid <- rep(TRUE, n_run)
      valid[censor_rel] <- FALSE
      valid_idx <- which(valid)
    } else {
      valid_idx <- seq_len(n_run)
    }

    if (method == "ar") {
      n_eff <- length(valid_idx)
      if (n_eff <= 1L) {
        return(list(phi = numeric(0), theta = numeric(0), order = c(p = 0L, q = 0L)))
      }
      p_cap <- min(as.integer(p_max), n_eff - 1L)

      # Pool the autocovariance over the valid segments, centering once across
      # all of this run's surviving frames.
      seg_id <- cumsum(c(1L, as.integer(diff(valid_idx) > 1L)))
      pooled <- .pooled_acvf_segments(mat[valid_idx, , drop = FALSE], seg_id, p_cap)
      gamma0 <- .acvf_from_pooled(pooled, order = 0L)
      null_fit <- list(phi = numeric(0), theta = numeric(0),
                       order = c(p = 0L, q = 0L),
                       gamma = gamma0, sigma2 = if (length(gamma0)) gamma0[1] else NA_real_)
      if (!length(gamma0) || !is.finite(gamma0[1]) || gamma0[1] <= 0) return(null_fit)
      p_cap <- min(p_cap, .acvf_max_lag(pooled))
      if (p_cap < 1L) return(null_fit)
      gamma_full <- .acvf_from_pooled(pooled, order = p_cap)
      null_fit$gamma <- gamma_full
      null_fit$sigma2 <- gamma_full[1]

      if (!identical(p, "auto")) {
        pp <- min(as.integer(p), p_cap)
        if (pp <= 0L) return(null_fit)
        gamma_pp <- .acvf_from_pooled(pooled, order = pp)
        yw <- yw_from_acvf_fast(gamma_pp[seq_len(pp + 1L)], pp)
        return(list(phi = enforce_stationary_ar(yw$phi, 0.99),
                    theta = numeric(0), order = c(p = pp, q = 0L),
                    gamma = gamma_full, sigma2 = pmax(yw$sigma2, 1e-12)))
      }

      best_phi <- numeric(0)
      best_order <- c(p = 0L, q = 0L)
      n_eff_log <- log(n_eff)
      sigma0 <- pmax(gamma0[1], 1e-12)
      best_bic <- 2 * n_eff * log(sigma0) + n_eff_log
      best_sigma2 <- sigma0
      # Bound the order BIC may select by available data; an explicit p is honoured.
      p_sel <- min(p_cap, floor(n_eff / 5))
      if (p_cap >= 1L) {
        for (pp in seq_len(max(0L, p_sel))) {
          gamma_pp <- .acvf_from_pooled(pooled, order = pp)
          yw <- yw_from_acvf_fast(gamma_pp[seq_len(pp + 1L)], pp)
          sigma2 <- pmax(yw$sigma2, 1e-12)
          if (!is.finite(sigma2)) next
          bic <- 2 * n_eff * log(sigma2) + (pp + 1L) * n_eff_log
          if (!is.finite(bic) || bic >= best_bic) next
          phi_pp <- enforce_stationary_ar(yw$phi, 0.99)
          if (length(phi_pp) != pp || !all(is.finite(phi_pp))) next
          best_bic <- bic
          best_phi <- phi_pp
          best_order <- c(p = pp, q = 0L)
          best_sigma2 <- sigma2
        }
      }
      list(phi = best_phi, theta = numeric(0), order = best_order,
           gamma = gamma_full, sigma2 = best_sigma2)
    } else {
      # For ARMA: use valid timepoints only.
      #
      # Censored frames are dropped, but Hannan-Rissanen then runs on the
      # surviving frames spliced together, so its regressions do span the gaps.
      # That biases the estimate (an ARMA(1,1) with theta = 0.4 comes back near
      # 0.30 at 25% censoring). Segmenting it properly is a larger change than
      # this release carries, so warn rather than fail silently.
      if (length(censor_rel) && any(diff(valid_idx) > 1L)) {
        warning("fit_noise: method = 'arma' estimates across censoring gaps; ",
                "AR and MA coefficients will be biased. Prefer method = 'ar' ",
                "when censoring is present.", call. = FALSE)
      }
      mat_valid <- mat[valid_idx, , drop = FALSE]
      y_mean <- rowMeans(mat_valid)
      pp <- if (identical(p, "auto")) min(2L, p_max) else as.integer(p)
      qq <- as.integer(q)
      fit <- hr_arma(y_mean, p = pp, q = qq, iter = as.integer(hr_iter),
                     step1 = step1)

      # Report the noise autocovariance at voxel scale, pooled the same way the
      # AR path does. hr_arma's own sigma2 is the innovation variance of the
      # run-MEAN series, which is smaller than the per-voxel value by roughly
      # the number of voxels averaged -- reporting it as the plan's sigma2 would
      # understate the noise by that factor. There is no comparably cheap
      # voxel-scale innovation variance for ARMA, so it is left NA rather than
      # filled with a number on the wrong scale.
      seg_id_ma <- cumsum(c(1L, as.integer(diff(valid_idx) > 1L)))
      lag_ma <- max(1L, pp + qq)
      pooled_ma <- .pooled_acvf_segments(mat_valid, seg_id_ma, lag_ma)
      fit$gamma <- .acvf_from_pooled(pooled_ma, order = lag_ma)
      fit$sigma2 <- NA_real_
      fit
    }
  }

  if (pooling == "parcel") {
    if (!identical(method, "ar")) stop("Parcel pooling currently supports method = 'ar' only")
    stopifnot(!is.null(parcels))
    parcels <- as.integer(parcels)
    stopifnot(length(parcels) == ncol(resid))

    # Drop censored frames and segment at both run boundaries and censoring
    # gaps, so parcel estimation sees the same valid segment structure the
    # global/run path uses. Passing censor = NULL here (as earlier versions
    # did) silently estimated across scrubbed frames and run boundaries alike.
    seg <- .valid_segments(n, runs = runs, censor = censor)
    if (length(seg$idx) < 2L) stop("no valid timepoints remain after censoring")
    run_starts0 <- seg$starts0

    estimator <- function(y, starts0 = run_starts0) {
      .estimate_ar_series(y, p_max, p = p, starts0 = starts0,
                          center_id = seg$run_id)
    }
    M_fine <- .parcel_means(resid, parcels)[seg$idx, , drop = FALSE]

    target <- if (is.null(p_target)) {
      if (identical(p, "auto")) {
        p_max
      } else if (multiscale_explicit) {
        min(as.integer(p), p_max)
      } else {
        p_max
      }
    } else {
      min(as.integer(p_target), p_max)
    }

    if (is.null(parcel_sets)) {
      est_f <- .ms_estimate_scale(M_fine, estimator, run_starts0, lag_max = target, center_id = seg$run_id)
      if (is.null(multiscale_mode) || target == 0L) {
        phi_parcel <- est_f$phi
      } else if (identical(multiscale_mode, "pacf_weighted")) {
        shrink <- 0.6
        kap_mat <- do.call(cbind, lapply(est_f$phi, function(phi) .ms_pad(ar_to_pacf(phi), target)))
        if (!is.matrix(kap_mat)) kap_mat <- matrix(kap_mat, nrow = target)
        avg_kap <- if (target > 0L) rowMeans(kap_mat, na.rm = TRUE) else numeric(0)
        avg_kap <- pmin(pmax(avg_kap, -0.99), 0.99)
        phi_parcel <- lapply(est_f$phi, function(phi) {
          kap_f <- .ms_pad(ar_to_pacf(phi), target)
          kap_mix <- (1 - shrink) * kap_f + shrink * avg_kap
          pacf_to_ar(pmin(pmax(kap_mix, -0.99), 0.99))
        })
      } else {
        shrink <- 0.6
        acvf_mat <- vapply(est_f$acvf, function(g) .ms_pad(g, target + 1L), numeric(target + 1L))
        avg_g <- rowMeans(acvf_mat, na.rm = TRUE)
        phi_parcel <- lapply(est_f$acvf, function(g) {
          g_pad <- .ms_pad(g, target + 1L)
          g_mix <- (1 - shrink) * g_pad + shrink * avg_g
          yw <- yw_from_acvf_fast(g_mix, target)
          enforce_stationary_ar(yw$phi)
        })
      }
    } else {
      required_keys <- c("coarse", "medium", "fine")
      stopifnot(all(required_keys %in% names(parcel_sets)))
      parcels_coarse <- as.integer(parcel_sets$coarse)
      parcels_medium <- as.integer(parcel_sets$medium)
      parcels_fine <- as.integer(parcel_sets$fine)
      stopifnot(length(parcels_coarse) == ncol(resid))
      stopifnot(length(parcels_medium) == ncol(resid))
      stopifnot(all(parcels_fine == parcels))

      # Subset to the surviving frames exactly as M_fine is, or the centering
      # grouping and the data disagree in length.
      M_coarse <- .parcel_means(resid, parcels_coarse)[seg$idx, , drop = FALSE]
      M_medium <- .parcel_means(resid, parcels_medium)[seg$idx, , drop = FALSE]
      est_c <- .ms_estimate_scale(M_coarse, estimator, run_starts0, lag_max = target, center_id = seg$run_id)
      est_m <- .ms_estimate_scale(M_medium, estimator, run_starts0, lag_max = target, center_id = seg$run_id)
      est_f <- .ms_estimate_scale(M_fine, estimator, run_starts0, lag_max = target, center_id = seg$run_id)

      parents <- .ms_parent_maps(parcels_fine, parcels_medium, parcels_coarse)
      sizes <- list(
        n_t = nrow(resid),
        n_runs = if (is.null(runs)) 1L else length(unique(as.integer(runs))),
        beta = beta,
        coarse = as.list(table(parcels_coarse)),
        medium = as.list(table(parcels_medium)),
        fine = as.list(table(parcels_fine))
      )
      disp_list <- list(
        coarse = .ms_dispersion(resid, parcels_coarse),
        medium = .ms_dispersion(resid, parcels_medium),
        fine = .ms_dispersion(resid, parcels_fine)
      )
      if (is.null(multiscale_mode)) {
        phi_parcel <- est_f$phi
      } else {
        phi_parcel <- .ms_combine_to_fine(
          phi_by_coarse = est_c$phi,
          phi_by_medium = est_m$phi,
          phi_by_fine   = est_f$phi,
          acvf_by_coarse = if (identical(multiscale_mode, "acvf_pooled")) est_c$acvf else NULL,
          acvf_by_medium = if (identical(multiscale_mode, "acvf_pooled")) est_m$acvf else NULL,
          acvf_by_fine   = if (identical(multiscale_mode, "acvf_pooled")) est_f$acvf else NULL,
          parents = parents,
          sizes = sizes,
          disp_list = disp_list,
          p_target = target,
          mode = multiscale_mode
        )
      }
    }

    if (is.null(multiscale_mode) && !is.null(p_target) && target > 0L) {
      phi_parcel <- mapply(function(phi, g) {
        g_pad <- .ms_pad(g, target + 1L)
        yw <- yw_from_acvf_fast(g_pad, target)
        enforce_stationary_ar(yw$phi)
      }, phi_parcel, est_f$acvf, SIMPLIFY = FALSE)
    }

    # .ms_pad() zero-fills each parcel out to the pooling target, so the stored
    # length reports padding rather than the order actually fitted -- an AR(2)
    # fit with p_max = 6 would claim order 6. Trim the trailing zeros that every
    # parcel shares and report the order that was really estimated.
    eff <- vapply(phi_parcel, function(ph) {
      nz <- which(ph != 0)
      if (!length(nz)) 0L else max(nz)
    }, 0L)
    keep <- if (length(eff)) max(eff) else 0L
    phi_parcel <- lapply(phi_parcel, function(ph) {
      if (length(ph) >= keep) ph[seq_len(keep)] else c(ph, rep(0, keep - length(ph)))
    })

    theta_parcel <- setNames(vector("list", length(phi_parcel)), names(phi_parcel))
    order_vec <- c(p = keep, q = 0L)

    # Per-parcel noise scale and shape, keyed like phi_by_parcel.
    #
    # gamma is computed from the VOXELS in each parcel, not from the parcel-mean
    # series that phi was estimated on. The parcel mean is the right basis for
    # phi -- correlation structure is scale-invariant -- but its variance is
    # smaller than a voxel's by the number of voxels averaged, so reporting it
    # as the noise scale understated the noise by exactly that factor (16-fold
    # at 16 voxels per parcel) and made the plan's units depend on the pooling
    # mode. sigma2 is then derived from the phi actually stored, since under
    # multiscale pooling the fine-scale innovation variance no longer matches
    # the pooled coefficients.
    n_valid <- length(seg$idx)
    seg_id_p <- cumsum(seq_len(n_valid) %in% (seg$starts0 + 1L))
    resid_valid <- resid[seg$idx, , drop = FALSE]
    gamma_parcel <- setNames(lapply(names(phi_parcel), function(k) {
      cols <- which(parcels == as.integer(k))
      if (!length(cols)) return(numeric(0))
      lag_k <- max(as.integer(target), length(phi_parcel[[k]]), 1L)
      .acvf_from_pooled(
        .pooled_acvf_segments(resid_valid[, cols, drop = FALSE], seg_id_p,
                              lag_k, center_id = seg$run_id),
        order = lag_k)
    }), names(phi_parcel))
    sigma2_parcel <- setNames(
      lapply(names(phi_parcel), function(k)
        .sigma2_from_gamma_phi(gamma_parcel[[k]], phi_parcel[[k]])),
      names(phi_parcel))

    parcel_ids <- names(phi_parcel)
    if (is.null(parcel_ids)) parcel_ids <- as.character(sort(unique(parcels)))
    return(new_whiten_plan(
      phi = NULL,
      theta = NULL,
      order = order_vec,
      runs = runs,
      exact_first = (exact_first == "ar1"),
      method = method,
      pooling = "parcel",
      parcels = parcels,
      parcel_ids = parcel_ids,
      phi_by_parcel = phi_parcel,
      theta_by_parcel = theta_parcel,
      censor = censor,
      gamma_by_parcel = gamma_parcel,
      sigma2_by_parcel = sigma2_parcel
    ))
  }

  estimates <- mapply(est_run, run_mats, censor_by_run, SIMPLIFY = FALSE)

  if (pooling == "global") {
    lens <- vapply(Rsets, length, 0L)
    pmax_len <- max(vapply(estimates, function(e) length(e$phi), 0L))
    qmax_len <- max(vapply(estimates, function(e) length(e$theta), 0L))
    Phi <- matrix(0, length(estimates), pmax_len)
    Th <- matrix(0, length(estimates), qmax_len)
    for (i in seq_along(estimates)) {
      if (length(estimates[[i]]$phi))   Phi[i, seq_along(estimates[[i]]$phi)]   <- estimates[[i]]$phi
      if (length(estimates[[i]]$theta)) Th[i, seq_along(estimates[[i]]$theta)] <- estimates[[i]]$theta
    }
    w <- lens / sum(lens)
    # Averaging per-run coefficient vectors does not preserve stationarity (or
    # invertibility), and runs that selected a lower order are zero-padded to the
    # longest before averaging, so the pooled vector need not be a valid AR of
    # any run. Re-impose the constraints on the pooled result.
    phi_pooled <- as.numeric(drop(crossprod(w, Phi)))
    theta_pooled <- as.numeric(drop(crossprod(w, Th)))
    if (length(phi_pooled)) phi_pooled <- enforce_stationary_ar(phi_pooled, 0.99)
    if (length(theta_pooled)) theta_pooled <- enforce_invertible_ma(theta_pooled)
    phi_list <- list(phi_pooled)
    theta_list <- list(theta_pooled)
    # Pool the autocovariance length-weighted, but TRUNCATE every run to the
    # shortest available length rather than zero-padding to the longest. Runs
    # legitimately reach different lags when their lengths or censoring differ,
    # and a zero-padded acvf is not a covariance at all: Toeplitz(1, 0.9, 0, 0)
    # has a negative eigenvalue. Averaging equal-length PD Toeplitz matrices is
    # PD by convexity; averaging zero-padded ones is not, and consumers building
    # Sigma from gamma got negative contrast variances.
    glens <- vapply(estimates, function(e) length(e$gamma %||% numeric(0)), 0L)
    gmin <- if (length(glens)) min(glens) else 0L
    gamma_list <- if (gmin > 0L) {
      G <- do.call(rbind, lapply(estimates, function(e) e$gamma[seq_len(gmin)]))
      list(as.numeric(drop(crossprod(w, G))))
    } else {
      list(numeric(0))
    }
    sigma2_list <- list(if (identical(method, "arma")) NA_real_ else
      .sigma2_from_gamma_phi(gamma_list[[1]], phi_pooled))
  } else {
    phi_list <- lapply(estimates, `[[`, "phi")
    theta_list <- lapply(estimates, `[[`, "theta")
    gamma_list <- lapply(estimates, function(e) e$gamma %||% numeric(0))
    sigma2_list <- if (identical(method, "arma")) {
      rep(list(NA_real_), length(estimates))
    } else {
      mapply(.sigma2_from_gamma_phi, gamma_list, phi_list, SIMPLIFY = FALSE)
    }
  }

  order_vec <- if (pooling == "global") {
    c(p = length(phi_list[[1]]), q = length(theta_list[[1]]))
  } else {
    c(p = max(vapply(phi_list, length, 0L)), q = max(vapply(theta_list, length, 0L)))
  }

  new_whiten_plan(
    phi = phi_list,
    theta = theta_list,
    order = order_vec,
    runs = runs,
    exact_first = (exact_first == "ar1"),
    method = method,
    pooling = pooling,
    censor = censor,
    gamma = gamma_list,
    sigma2 = sigma2_list
  )
}

#' Apply a whitening plan to design and data matrices
#'
#' @param plan Whitening plan from [fit_noise()].
#' @param X Numeric matrix of predictors (time x regressors).
#' @param Y Numeric matrix of data (time x voxels).
#' @param runs Optional run labels.
#' @param run_starts Optional 0-based run start indices (alternative to `runs`).
#' @param censor Optional indices of censored TRs (1-based); filter resets after gaps.
#' @param parcels Optional parcel labels (length = ncol(Y)) when using parcel plans.
#' @param inplace Modify inputs in place (logical).
#' @param parallel Use OpenMP parallelism if available.
#' @return List with whitened data. Parcel plans return `X_by` per parcel; others return a single `X` matrix.
#' @examples
#' # Create example design matrix and data
#' n_time <- 200
#' n_pred <- 3
#' n_voxels <- 50
#' X <- matrix(rnorm(n_time * n_pred), n_time, n_pred)
#' Y <- X %*% matrix(rnorm(n_pred * n_voxels), n_pred, n_voxels) +
#'      matrix(rnorm(n_time * n_voxels), n_time, n_voxels)
#'
#' # Fit noise model from residuals
#' residuals <- Y - X %*% solve(crossprod(X), crossprod(X, Y))
#' plan <- fit_noise(residuals, method = "ar", p = 2)
#'
#' # Apply whitening
#' whitened <- whiten_apply(plan, X, Y)
#' Xw <- whitened$X
#' Yw <- whitened$Y
#' @export
whiten_apply <- function(plan, X, Y, runs = NULL, run_starts = NULL, censor = NULL, parcels = NULL,
                         inplace = FALSE, parallel = TRUE) {
  stopifnot(inherits(plan, "fmriAR_plan"))
  if (!is.matrix(X)) X <- as.matrix(X)
  if (!is.matrix(Y)) Y <- as.matrix(Y)
  if (anyNA(X) || anyNA(Y)) {
    stop("whiten_apply: NA values detected in X or Y")
  }
  n <- nrow(X)
  stopifnot(nrow(Y) == n)

  if (!is.null(run_starts)) run_starts <- as.integer(run_starts)
  if (is.null(runs) && !is.null(run_starts)) {
    runs <- .runs_from_starts0(run_starts, n)
  }
  if (is.null(runs) && !is.null(plan$runs) && length(plan$runs) == n) runs <- plan$runs
  if (is.null(runs)) runs <- rep_len(1L, n)

  # Validate rather than let split() recycle or drop. A short `runs` was
  # silently recycled, and an NA left that row unwritten in both X and Y.
  if (length(runs) != n) {
    stop(sprintf("whiten_apply: 'runs' has length %d but X and Y have %d rows",
                 length(runs), n))
  }
  if (anyNA(runs)) stop("whiten_apply: 'runs' contains NA")
  runs <- if (is.numeric(runs)) as.integer(runs) else match(runs, unique(runs))

  if (identical(plan$pooling, "parcel")) {
    parcels_vec <- plan$parcels
    if (!is.null(parcels)) {
      stopifnot(length(parcels) == ncol(Y))
      parcels_vec <- as.integer(parcels)
    }
    stopifnot(!is.null(parcels_vec))
    stopifnot(length(parcels_vec) == ncol(Y))

    run_starts_vec <- .full_run_starts(runs, censor, n)
    parcel_ids <- if (!is.null(plan$parcel_ids)) plan$parcel_ids else sort(unique(parcels_vec))
    phi_by <- plan$phi_by_parcel
    theta_by <- plan$theta_by_parcel

    Yw <- matrix(NA_real_, n, ncol(Y))
    X_by <- setNames(vector("list", length(parcel_ids)), as.character(parcel_ids))
    X_base <- X

    for (pid in parcel_ids) {
      cols <- which(parcels_vec == pid)
      if (!length(cols)) next
      key <- as.character(pid)
      phi <- phi_by[[key]]
      if (is.null(phi)) phi <- numeric(0)
      theta <- theta_by[[key]]
      if (is.null(theta)) theta <- numeric(0)
      Y_sub <- Y[, cols, drop = FALSE]
      # arma_whiten_inplace() writes through its arguments. Assigning X_base
      # without subsetting shares storage with the caller's X, so every parcel
      # filtered the previous parcel's output, all X_by entries aliased one
      # matrix, and the caller's design was silently overwritten.
      X_sub <- X_base[seq_len(n), , drop = FALSE]
      out <- arma_whiten_inplace(
        Y = Y_sub,
        X = X_sub,
        phi = phi,
        theta = theta,
        run_starts = run_starts_vec,
        exact_first_ar1 = isTRUE(plan$exact_first),
        parallel = parallel
      )
      Yw[, cols] <- out$Y
      X_by[[key]] <- out$X
    }

    if (inplace) {
      Y[,] <- Yw
      return(invisible(list(X = NULL, X_by = X_by, Y = Y)))
    }
    return(list(X = NULL, X_by = X_by, Y = Yw))
  }

  rsplits <- split(seq_len(n), as.integer(runs))

  censor_by_run <- lapply(rsplits, function(idx) integer(0L))
  if (!is.null(censor)) {
    censor <- as.integer(censor)
    for (ri in seq_along(rsplits)) {
      idx <- rsplits[[ri]]
      c_in <- intersect(censor, idx)
      if (length(c_in)) censor_by_run[[ri]] <- as.integer(c_in - min(idx) + 1L)
    }
  }

  phi_list <- if (length(plan$phi) == 1L) rep(plan$phi, length(rsplits)) else plan$phi
  theta_list <- if (length(plan$theta) == 1L) rep(plan$theta, length(rsplits)) else plan$theta

  # Scatter each run's result back to the rows it came from. rbind()-ing the
  # split pieces reassembles them in sorted run-label order, so any labelling
  # that is not ascending in time (e.g. run 2 acquired first) silently returned
  # a row permutation of the right answer for both X and Y.
  Xw <- matrix(NA_real_, n, ncol(X), dimnames = dimnames(X))
  Yw <- matrix(NA_real_, n, ncol(Y), dimnames = dimnames(Y))

  for (ri in seq_along(rsplits)) {
    idx <- rsplits[[ri]]
    Xr <- X[idx, , drop = FALSE]
    Yr <- Y[idx, , drop = FALSE]
    rs <- .sub_run_starts(n_run = nrow(Xr), censor_idx_rel = censor_by_run[[ri]])
    out <- arma_whiten_inplace(
      Yr,
      Xr,
      phi = phi_list[[ri]],
      theta = theta_list[[ri]],
      run_starts = rs,
      exact_first_ar1 = isTRUE(plan$exact_first),
      parallel = parallel
    )
    Xw[idx, ] <- out$X
    Yw[idx, ] <- out$Y
  }

  if (inplace) {
    X[,] <- Xw
    Y[,] <- Yw
    invisible(list(X = X, Y = Y))
  } else {
    list(X = Xw, Y = Yw)
  }
}

#' Fit and apply whitening in one call
#'
#' @param X Design matrix (time x regressors).
#' @param Y Data matrix (time x voxels).
#' @param runs Optional run labels.
#' @param censor Optional censor indices.
#' @param ... Additional parameters passed to [fit_noise()].
#' @return List with whitened `X` and `Y` matrices.
#' @examples
#' # Create example data
#' n_time <- 200
#' n_pred <- 3
#' n_voxels <- 50
#' X <- matrix(rnorm(n_time * n_pred), n_time, n_pred)
#' Y <- X %*% matrix(rnorm(n_pred * n_voxels), n_pred, n_voxels) +
#'      matrix(rnorm(n_time * n_voxels, sd = 2), n_time, n_voxels)
#'
#' # One-step whitening
#' whitened <- whiten(X, Y, method = "ar", p = 2)
#' @export
whiten <- function(X, Y, runs = NULL, censor = NULL, ...) {
  res <- Y - X %*% qr.solve(X, Y)
  plan <- fit_noise(resid = res, runs = runs, censor = censor, ...)
  whiten_apply(plan, X, Y, runs = runs, censor = censor)
}
