#' @keywords internal
.lag_matrix <- function(x, k) {
  n <- length(x)
  if (k == 0L) return(matrix(numeric(0), n, 0L))
  M <- matrix(NA_real_, n, k)
  for (i in 1L:k) M[(i + 1L):n, i] <- x[1L:(n - i)]
  M
}

#' @keywords internal
hr_arma <- function(y, p, q,
                    p_big = NULL,
                    step1 = c("burg", "yw"),
                    iter = 0L,
                    ar_pacf_bound = 0.99,
                    enforce = TRUE) {
  stopifnot(is.numeric(y), length(p) == 1L, length(q) == 1L, p >= 0L, q >= 0L)
  use_cpp <- isTRUE(getOption("fmriAR.use_cpp_hr", TRUE))

  res <- NULL
  cpp_ok <- FALSE
  if (use_cpp) {
    res <- try(hr_arma_fit_cpp(as.numeric(y), as.integer(p), as.integer(q),
                               if (is.null(p_big)) 0L else as.integer(p_big),
                               as.integer(iter)),
               silent = TRUE)
    if (!inherits(res, "try-error")) {
      cpp_ok <- isTRUE(res$ok)
    }
  }

  if (!use_cpp || inherits(res, "try-error") || !cpp_ok) {
    res <- hr_arma_R(y, p, q,
                     p_big = p_big,
                     step1 = step1,
                     iter = iter,
                     ar_pacf_bound = ar_pacf_bound,
                     enforce = FALSE)
    res$ok <- TRUE
  }

  phi   <- res$phi
  theta <- res$theta
  if (enforce) {
    if (length(phi))   phi   <- enforce_stationary_ar(phi, bound = ar_pacf_bound)
    if (length(theta)) theta <- enforce_invertible_ma(theta)
  }

  list(phi = as.numeric(phi),
       theta = as.numeric(theta),
       sigma2 = as.numeric(res$sigma2),
       order = c(p = as.integer(p), q = as.integer(q)),
       method = "hr",
       p_big = as.integer(res$p_big),
       iterations = as.integer(res$iter))
}

#' @keywords internal
hr_arma_R <- function(y, p, q,
                      p_big = NULL,
                      step1 = c("burg", "yw"),
                      iter = 0L,
                      ar_pacf_bound = 0.99,
                      enforce = TRUE) {
  step1 <- match.arg(step1)
  y <- as.numeric(y) - mean(y)
  n <- length(y)
  if (n < 10L) stop("Series too short for HR estimation")

  if (is.null(p_big)) {
    p_big <- max(8L, p + q + 5L, ceiling(10 * log10(n)))
    p_big <- min(p_big, max(2L, n - 2L), 40L)
  }

  arfit <- stats::ar(y, aic = FALSE, order.max = p_big, method = step1)
  ehat <- as.numeric(arfit$resid)
  ehat[is.na(ehat)] <- 0

  phi <- if (p > 0L) numeric(p) else numeric(0)
  theta <- if (q > 0L) numeric(q) else numeric(0)

  for (ii in 0L:iter) {
    Ylags <- if (p > 0L) .lag_matrix(y, p) else matrix(numeric(0), n, 0L)
    Elags <- if (q > 0L) .lag_matrix(ehat, q) else matrix(numeric(0), n, 0L)
    m <- max(p, q)
    idx <- seq.int(m + 1L, n)
    Z <- cbind(Ylags[idx, , drop = FALSE], Elags[idx, , drop = FALSE])
    z_y <- y[idx]
    if (nrow(Z) < (ncol(Z) + 1L)) stop("Not enough data for the requested (p,q)")
    coef <- tryCatch(qr.solve(Z, z_y), error = function(e) NA)
    if (anyNA(coef)) stop("HR regression failed")
    if (p > 0L)  phi   <- coef[seq_len(p)]
    if (q > 0L)  theta <- coef[p + seq_len(q)]
    ehat <- .arma_innovations(y, phi, theta)
  }

  if (enforce) {
    if (length(phi))   phi   <- enforce_stationary_ar(phi, bound = ar_pacf_bound)
    if (length(theta)) theta <- enforce_invertible_ma(theta)
  }

  list(phi = phi,
       theta = theta,
       sigma2 = mean(ehat^2),
       p_big = p_big,
       iter = iter)
}
