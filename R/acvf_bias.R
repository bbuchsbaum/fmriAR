# Residual-bias correction for autocovariance estimated from GLM residuals.
#
# Autocovariance estimated from residuals is biased. With residuals
# ehat = M y and M = I - X(X'X)^-1 X' the residual-forming projection,
#
#   E[ehat ehat'] = M Sigma M
#
# so every lag product mixes in structure from the whole design. The bias runs
# downward on the correlations and grows with the number of regressors: for
# AR(1) rho = 0.4 at T = 300, a 9-column design returns about 0.36 and a
# 28-column design about 0.25.
#
# Sigma is Toeplitz within a run, Sigma = sum_k gamma_k S_k, where S_k is the
# lag-k indicator that does not cross a run boundary. Expectation is linear in
# gamma, so what the estimator reports is a fixed linear map of the truth,
#
#   E[gamma_raw] = A gamma_true,
#   A[h,k] = (1/pairs_h) * sum_{(a,b) in P_h} (R S_k R')[a,b]
#
# where R is the operator carrying the data to the residuals the estimator
# actually accumulates, and P_h is the set of lag-h pairs it actually uses.
# A depends only on the design and the segmentation, so it is built once and
# reused; gamma_true is recovered by solving A gamma = gamma_raw.
#
# Two things make this exact rather than approximate. R folds in the per-run
# centering the estimator applies after projection, since Cov(C M y) is
# (CM) Sigma (CM)' and CM differs from M whenever the design carries no per-run
# intercept. And P_h and pairs_h are read off the same segmentation the
# estimator uses, so censoring and run boundaries are handled by construction
# rather than assumed away.
#
# The correction is exact for noise whose autocovariance dies within the lag
# budget. Under long memory it is partial, because truncating the system
# aliases in structure beyond the budget.

# Run membership per original timepoint, as integer codes.
.run_vector <- function(n, runs = NULL) {
  n <- as.integer(n)
  if (is.null(runs)) return(rep(1L, n))
  if (length(runs) != n) stop("'runs' must have one entry per timepoint")
  if (anyNA(runs)) stop("'runs' contains NA")
  as.integer(match(runs, unique(runs)))
}

# S_k applied along the time axis to every row of B (each row an n-vector):
#   (S_k v)[i] = v[i-k] + v[i+k], keeping only neighbours inside the same run.
# S_0 is the identity.
.apply_lag_operator <- function(B, k, run_vec) {
  if (k == 0L) return(B)
  n <- ncol(B)
  out <- matrix(0, nrow(B), n)
  if (k >= n) return(out)
  hi <- seq.int(k + 1L, n)
  lo <- seq.int(1L, n - k)
  same <- run_vec[hi] == run_vec[lo]
  if (!any(same)) return(out)
  hi <- hi[same]; lo <- lo[same]
  out[, hi] <- out[, hi] + B[, lo, drop = FALSE]
  out[, lo] <- out[, lo] + B[, hi, drop = FALSE]
  out
}

# Bias matrix for one pooling unit. `Q` is an orthonormal basis for the design
# column space, `idx` the unit's valid timepoints in original time, `seg_id`
# their segment labels. Centering is a single group over `idx`, matching
# .pooled_acvf_segments() called without center_id, which is how fit_noise
# estimates each run.
.acvf_bias_core <- function(Q, run_vec, idx, seg_id, max_lag) {
  n <- length(run_vec)
  nv <- length(idx)
  max_lag <- max(0L, as.integer(max_lag))
  if (nv < 2L) return(diag(max_lag + 1L))

  R <- -tcrossprod(Q[idx, , drop = FALSE], Q)
  R[cbind(seq_len(nv), idx)] <- R[cbind(seq_len(nv), idx)] + 1
  R <- R - rep(colMeans(R), each = nv)

  a_list <- vector("list", max_lag + 1L)
  b_list <- vector("list", max_lag + 1L)
  npair <- numeric(max_lag + 1L)
  a_list[[1L]] <- seq_len(nv); b_list[[1L]] <- seq_len(nv); npair[1L] <- nv
  for (lg in seq_len(max_lag)) {
    if (nv <= lg) break
    hi <- seq.int(lg + 1L, nv)
    lo <- seq.int(1L, nv - lg)
    ok <- seg_id[hi] == seg_id[lo]
    a_list[[lg + 1L]] <- hi[ok]
    b_list[[lg + 1L]] <- lo[ok]
    npair[lg + 1L] <- sum(ok)
  }

  A <- diag(max_lag + 1L)
  for (k in 0:max_lag) {
    SkR <- .apply_lag_operator(R, k, run_vec)
    for (h in 0:max_lag) {
      if (npair[h + 1L] <= 0) next
      A[h + 1L, k + 1L] <-
        sum(R[a_list[[h + 1L]], , drop = FALSE] *
            SkR[b_list[[h + 1L]], , drop = FALSE]) / npair[h + 1L]
    }
  }
  A
}

# One bias matrix per run, keyed like the run split fit_noise uses. Each run is
# estimated separately, and the residual operator restricted to a run still
# involves the whole design, so the matrices genuinely differ between runs.
.acvf_bias_by_run <- function(design, n, runs = NULL, censor = NULL, max_lag = 20L) {
  if (!is.matrix(design)) design <- as.matrix(design)
  storage.mode(design) <- "double"
  if (nrow(design) != n) {
    stop("'design' must have one row per timepoint (", n, "), not ", nrow(design))
  }
  if (anyNA(design)) stop("'design' contains NA")
  run_vec <- .run_vector(n, runs)
  qrX <- qr(design)
  Q <- qr.Q(qrX)[, seq_len(qrX$rank), drop = FALSE]

  # Recovering max_lag + 1 autocovariance lags needs at least that many residual
  # dimensions to recover them from. Below the boundary the correction degrades
  # and individual fits start pinning phi at the stationarity clamp: at n = 300
  # with 8 residual degrees of freedom against a 25-lag budget the mean estimate
  # went to -0.18 for a truth of 0.5, worse than no correction at all.
  df <- n - qrX$rank
  if (df < as.integer(max_lag) + 1L) {
    capped <- max(1L, df - 1L)
    warning("residual-bias correction: the design leaves ", df,
            " residual degrees of freedom, fewer than the ",
            as.integer(max_lag) + 1L, " autocovariance lags requested. ",
            "Reducing the correction lag budget to ", capped,
            ". Correcting at a budget the data cannot support returns ",
            "coefficients pinned at the stationarity clamp.", call. = FALSE)
    max_lag <- capped
  }

  Rsets <- if (is.null(runs)) list(seq_len(n)) else split(seq_len(n), as.integer(runs))
  cens <- if (is.null(censor)) integer(0) else as.integer(censor)

  setNames(lapply(Rsets, function(ridx) {
    keep <- setdiff(ridx, cens)
    if (length(keep) < 2L) return(diag(as.integer(max_lag) + 1L))
    seg_id <- cumsum(c(1L, as.integer(diff(keep) > 1L)))
    .acvf_bias_core(Q, run_vec, keep, seg_id, max_lag)
  }), names(Rsets))
}

#' Bias matrix relating raw residual autocovariance to the truth
#'
#' Builds the linear map `A` satisfying `E[gamma_raw] = A gamma_true` induced by
#' projecting a design out of the data. Autocovariance estimated from GLM
#' residuals is biased low, increasingly so for richer designs, because
#' `E[ehat ehat'] = M Sigma M` for the residual-forming projection `M`. Solving
#' `A gamma = gamma_raw` removes that bias.
#'
#' Ordinarily you do not need this: pass `design` to [fit_noise()] and it builds
#' and applies the matrices itself. Use this directly to inspect the bias, or to
#' cache the matrices when many datasets share one design and segmentation and
#' feed them back as `acvf_correction`.
#'
#' One matrix is returned per run, because [fit_noise()] estimates each run
#' separately and the residual operator restricted to a run still involves the
#' whole design.
#'
#' The correction is exact for noise whose autocovariance dies within `max_lag`.
#' Under long memory it is partial, since truncating the system aliases in
#' structure past the budget. With high-pass filtering (DCT, 128 s) a budget
#' near 20-25 lags is typically enough; without high-pass it exceeds 150 and the
#' approach stops being practical.
#'
#' Cost is `O(max_lag^2 * n_valid * n)` per run, so it is worth caching.
#'
#' @param design Numeric design matrix (timepoints x regressors), the one whose
#'   projection produced the residuals.
#' @param runs Optional run labels, length `nrow(design)`. Lag products never
#'   cross a run boundary.
#' @param censor Optional indices of censored timepoints (1-based), excluded
#'   exactly as in [fit_noise()].
#' @param max_lag Highest lag to correct.
#' @return A named list of `(max_lag + 1)` square matrices, one per run, rows
#'   indexed by the lag of the raw estimate and columns by the lag of the truth.
#' @seealso [fit_noise()], [noise_acvf()]
#' @examples
#' n <- 120
#' X <- cbind(1, poly(seq_len(n), 3))
#' A <- acvf_bias_matrix(X, max_lag = 5)
#' round(A[[1]][1:3, 1:3], 3)
#' @export
acvf_bias_matrix <- function(design, runs = NULL, censor = NULL, max_lag = 20L) {
  if (!is.matrix(design)) design <- as.matrix(design)
  .acvf_bias_by_run(design, nrow(design), runs = runs, censor = censor,
                    max_lag = max_lag)
}

# Normalize whatever the caller supplied into one matrix per run. A single
# matrix is broadcast; a list is matched positionally against the run split.
.normalize_correction <- function(correction, n_runs) {
  if (is.null(correction)) return(NULL)
  if (is.matrix(correction)) return(rep(list(correction), n_runs))
  if (is.list(correction)) {
    if (length(correction) == 1L) return(rep(correction, n_runs))
    if (length(correction) != n_runs) {
      stop("'acvf_correction' has ", length(correction),
           " matrices but there are ", n_runs, " runs")
    }
    return(correction)
  }
  stop("'acvf_correction' must be a matrix or a list of matrices")
}
