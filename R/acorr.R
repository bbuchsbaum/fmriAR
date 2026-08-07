#' Autocorrelation diagnostics for residuals
#'
#' @param resid Numeric matrix (time x voxels), typically whitened residuals.
#' @param runs Optional run labels, length `nrow(resid)`. When supplied, each run
#'   is centred separately and no lag product spans a run boundary.
#' @param max_lag Maximum lag to evaluate.
#' @param aggregate Aggregation across voxels: "mean", "median", or "none".
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

  acf_one <- function(y) {
    if (is.null(seg_id)) {
      return(stats::acf(y, lag.max = max_lag, plot = FALSE, demean = TRUE)$acf[-1L])
    }
    g <- .acvf_from_pooled(
      .pooled_acvf_segments(matrix(as.numeric(y), ncol = 1L), seg_id, max_lag,
                            center_id = seg_id)
    )
    g <- c(g, rep(0, max_lag + 1L - length(g)))[seq_len(max_lag + 1L)]
    if (!is.finite(g[1L]) || g[1L] <= 0) return(rep(NA_real_, max_lag))
    g[-1L] / g[1L]
  }

  if (aggregate == "none") {
    A <- vapply(seq_len(ncol(resid)), function(j) acf_one(resid[, j]), numeric(max_lag))
    return(list(lags = seq_len(max_lag), acf = A, ci = ci))
  }

  ybar <- switch(aggregate,
                 mean = rowMeans(resid),
                 median = apply(resid, 1L, stats::median))
  a <- acf_one(ybar)
  list(lags = seq_len(max_lag), acf = a, ci = ci)
}
