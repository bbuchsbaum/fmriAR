#' Autocorrelation diagnostics for residuals
#'
#' @param resid Numeric matrix (time x voxels), typically whitened residuals.
#' @param runs Optional run labels.
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
#' # Examine lag-1 autocorrelation
#' lag1_acorr <- acorr_check$acf[2]  # First element is lag-0 (always 1)
#' @export
acorr_diagnostics <- function(resid, runs = NULL, max_lag = 20L,
                              aggregate = c("mean", "median", "none")) {
  stopifnot(is.matrix(resid))
  aggregate <- match.arg(aggregate)
  n <- nrow(resid)
  ci <- 1.96 / sqrt(n)

  acf_one <- function(y) stats::acf(y, lag.max = max_lag, plot = FALSE, demean = TRUE)$acf[-1L]

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
