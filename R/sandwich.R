#' GLS standard errors from whitened residuals
#'
#' @param Xw Whitened design matrix.
#' @param Yw Whitened data matrix (time x voxels).
#' @param beta Optional coefficients (p x v); estimated if `NULL`.
#' @param type `"iid"` (default) assumes the whitened errors are white and
#'   homoskedastic. `"hc0"` is a heteroskedasticity-robust sandwich. `"hac"` is
#'   a Newey-West (Bartlett kernel) sandwich computed within runs, which stays
#'   valid when the noise model leaves some autocorrelation behind.
#' @param df_mode Degrees-of-freedom mode: "rankX" (default) or "n-p".
#' @param runs Optional run labels (length `nrow(Xw)`). Used by `type = "hac"`
#'   so that no lag product spans a run boundary.
#' @param hac_lag Bandwidth for `type = "hac"`. Defaults to
#'   `floor(4 * (n_run / 100)^(2/9))` using the median run length.
#' @return List containing standard errors, innovation variances, and XtX
#'   inverse. Coefficients that are not estimable because `Xw` is rank
#'   deficient get `NA` standard errors (and `NA` rows/columns in `XtX_inv`),
#'   following [stats::lm()].
#' @examples
#' # Generate example whitened data
#' n_time <- 200
#' n_pred <- 3
#' n_voxels <- 50
#' Xw <- matrix(rnorm(n_time * n_pred), n_time, n_pred)
#' Yw <- matrix(rnorm(n_time * n_voxels), n_time, n_voxels)
#'
#' # Compute standard errors
#' se_result <- sandwich_from_whitened_resid(Xw, Yw, type = "iid")
#'
#' # Extract standard errors for first voxel
#' se_voxel1 <- se_result$se[, 1]
#' @export
sandwich_from_whitened_resid <- function(Xw, Yw, beta = NULL,
                                         type = c("iid", "hc0", "hac"),
                                         df_mode = c("rankX", "n-p"),
                                         runs = NULL,
                                         hac_lag = NULL) {
  stopifnot(is.matrix(Xw), is.matrix(Yw), nrow(Xw) == nrow(Yw))
  type <- match.arg(type)
  df_mode <- match.arg(df_mode)

  n <- nrow(Xw)
  p <- ncol(Xw)
  v <- ncol(Yw)

  # Pivoted QR identifies the estimable columns. chol(crossprod(Xw)) failed
  # outright on any rank-deficient design, which made df_mode = "rankX" moot.
  qrX <- qr(Xw)
  rankX <- qrX$rank
  est <- sort(qrX$pivot[seq_len(rankX)])
  Xe <- Xw[, est, drop = FALSE]
  XtX_e <- chol2inv(chol(crossprod(Xe)))
  XtX_inv <- matrix(NA_real_, p, p)
  if (!is.null(colnames(Xw))) dimnames(XtX_inv) <- list(colnames(Xw), colnames(Xw))
  XtX_inv[est, est] <- XtX_e

  if (is.null(beta)) {
    beta <- matrix(NA_real_, p, v)
    beta[est, ] <- XtX_e %*% crossprod(Xe, Yw)
    E <- Yw - Xe %*% beta[est, , drop = FALSE]
  } else {
    beta <- as.matrix(beta)
    stopifnot(nrow(beta) == p, ncol(beta) == v)
    E <- Yw - Xe %*% beta[est, , drop = FALSE]
  }
  df <- if (df_mode == "rankX") n - rankX else n - p
  sigma2 <- colSums(E^2) / df

  se <- matrix(NA_real_, p, v)
  if (type == "iid") {
    se[est, ] <- sqrt(outer(diag(XtX_e), sigma2))
    return(list(se = se, sigma2 = sigma2, XtX_inv = XtX_inv, df = df, type = "iid"))
  }

  # H = (X'X)^-1 X'. Var(beta_k) = sum_t sum_s H_kt H_ks Omega_ts, so only the
  # rows of H and the residuals are needed -- no per-voxel p x p products.
  H <- XtX_e %*% t(Xe)
  if (type == "hc0") {
    se[est, ] <- sqrt((H * H) %*% (E * E))
    return(list(se = se, sigma2 = sigma2, XtX_inv = XtX_inv, df = df, type = "hc0"))
  }

  run_codes <- .run_codes(runs, n)
  if (is.null(hac_lag)) {
    run_len <- stats::median(tabulate(run_codes))
    hac_lag <- floor(4 * (run_len / 100)^(2 / 9))
  }
  hac_lag <- as.integer(max(0, min(hac_lag, n - 1L)))
  pairs <- lapply(seq_len(hac_lag), function(l) {
    hi <- seq.int(l + 1L, n)
    lo <- seq.int(1L, n - l)
    ok <- run_codes[hi] == run_codes[lo]
    list(hi = hi[ok], lo = lo[ok], w = 1 - l / (hac_lag + 1))
  })
  for (k in seq_along(est)) {
    U <- H[k, ] * E
    vk <- colSums(U * U)
    for (pr in pairs) {
      if (!length(pr$hi)) next
      vk <- vk + 2 * pr$w * colSums(U[pr$hi, , drop = FALSE] * U[pr$lo, , drop = FALSE])
    }
    se[est[k], ] <- sqrt(pmax(vk, 0))
  }
  list(se = se, sigma2 = sigma2, XtX_inv = XtX_inv, df = df, type = "hac",
       hac_lag = hac_lag)
}
