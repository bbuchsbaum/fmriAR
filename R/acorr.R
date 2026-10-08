#' Autocorrelation diagnostics for residuals
#'
#' @param resid Numeric matrix (time x voxels), typically whitened residuals.
#' @param runs Optional run labels, length `nrow(resid)`. When supplied, each run
#'   is centred separately and no lag product spans a run boundary.
#' @param max_lag Maximum lag to evaluate.
#' @param aggregate How per-voxel autocorrelations are combined: "mean",
#'   "median" (across voxels, lag by lag), or "none" (lags x voxels matrix).
#' @return List of autocorrelation values and nominal confidence interval.
#' @examples
#' # Generate example residuals with some autocorrelation
#' n_time <- 200
#' n_voxels <- 50
#' resid <- matrix(rnorm(n_time * n_voxels), n_time, n_voxels)
#'
#' # Add some AR(1) structure
#' for (v in 1:n_voxels) {
#'   resid[, v] <- filter(resid[, v], filter = 0.3, method = "recursive")
#' }
#'
#' # Check autocorrelation
#' acorr_check <- acorr_diagnostics(resid, max_lag = 10, aggregate = "mean")
#'
#' # Examine lag-1 autocorrelation (lag 0 is not returned, so acf[1] is lag 1)
#' lag1_acorr <- acorr_check$acf[1]
#' @export
acorr_diagnostics <- function(resid, runs = NULL, max_lag = 20L,
                              aggregate = c("mean", "median", "none")) {
  stopifnot(is.matrix(resid))
  aggregate <- match.arg(aggregate)
  n <- nrow(resid)
  ci <- 1.96 / sqrt(n)

  # `runs` was previously accepted and documented but never used, so the result
  # looked run-aware while lag products still spanned run boundaries.
  seg_id <- NULL
  if (!is.null(runs)) {
    seg_id <- .run_codes(runs, n)
  }

  # Per-voxel autocorrelation, lags x voxels: sum_t y_t y_{t-k} / sum_t y_t^2
  # with each voxel centred per run (or overall) and no lag product spanning a
  # run boundary. Without runs this is exactly stats::acf().
  if (is.null(seg_id)) seg_id <- rep(1L, n)
  Rc <- resid
  for (r in unique(seg_id)) {
    rows <- which(seg_id == r)
    Rc[rows, ] <- Rc[rows, , drop = FALSE] -
      rep(colMeans(Rc[rows, , drop = FALSE]), each = length(rows))
  }
  den <- colSums(Rc * Rc)
  A <- matrix(NA_real_, max_lag, ncol(resid))
  for (k in seq_len(max_lag)) {
    if (k >= n) break
    hi <- seq.int(k + 1L, n)
    lo <- seq.int(1L, n - k)
    ok <- seg_id[hi] == seg_id[lo]
    A[k, ] <- if (any(ok)) {
      colSums(Rc[hi[ok], , drop = FALSE] * Rc[lo[ok], , drop = FALSE]) / den
    } else 0
  }
  A[, !(den > 0)] <- NA_real_

  if (aggregate == "none") {
    return(list(lags = seq_len(max_lag), acf = A, ci = ci))
  }

  # Aggregate the per-voxel autocorrelations. Taking the autocorrelation of the
  # voxel-mean series instead (as earlier versions did) measures the
  # cross-voxel covariance, i.e. whatever signal is shared, and reported strong
  # residual autocorrelation for voxels that were individually white.
  a <- switch(aggregate,
              mean = rowMeans(A, na.rm = TRUE),
              median = apply(A, 1L, stats::median, na.rm = TRUE))
  list(lags = seq_len(max_lag), acf = a, ci = ci)
}
