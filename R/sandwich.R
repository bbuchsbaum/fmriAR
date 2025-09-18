#' GLS standard errors from whitened residuals
#'
#' @param Xw Whitened design matrix.
#' @param Yw Whitened data matrix (time x voxels).
#' @param beta Optional coefficients (p x v); estimated if `NULL`.
#' @param type Either "iid" (default) or "hc0" for a robust sandwich.
#' @param df_mode Degrees-of-freedom mode: "rankX" (default) or "n-p".
#' @param runs Optional run labels (reserved for future per-run scaling).
#' @return List containing standard errors, innovation variances, and XtX inverse.
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
                                         type = c("iid", "hc0"),
                                         df_mode = c("rankX", "n-p"),
                                         runs = NULL) {
  stopifnot(is.matrix(Xw), is.matrix(Yw), nrow(Xw) == nrow(Yw))
  type <- match.arg(type)
  df_mode <- match.arg(df_mode)

  n <- nrow(Xw)
  p <- ncol(Xw)
  v <- ncol(Yw)

  XtX <- crossprod(Xw)
  Rchol <- chol(XtX)
  XtX_inv <- chol2inv(Rchol)

  if (is.null(beta)) beta <- XtX_inv %*% crossprod(Xw, Yw)

  E <- Yw - Xw %*% beta
  rankX <- qr(Xw)$rank
  df <- if (df_mode == "rankX") n - rankX else n - p

  if (type == "iid") {
    sigma2 <- colSums(E^2) / df
    se <- sqrt(outer(diag(XtX_inv), sigma2))
    return(list(se = se, sigma2 = sigma2, XtX_inv = XtX_inv, df = df, type = "iid"))
  }

  se <- matrix(NA_real_, p, v)
  for (j in 1L:v) {
    e <- E[, j]
    Xe <- Xw * as.numeric(e)
    meat <- crossprod(Xe, Xe)
    V <- XtX_inv %*% meat %*% XtX_inv
    se[, j] <- sqrt(diag(V))
  }
  list(se = se, sigma2 = colSums(E^2) / df, XtX_inv = XtX_inv, df = df, type = "hc0")
}
