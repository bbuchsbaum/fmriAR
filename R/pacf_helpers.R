#' @keywords internal
pacf_to_ar <- function(kappa) {
  p <- length(kappa)
  if (p == 0L) return(numeric(0))
  Phi <- vector("list", p)
  Phi[[1L]] <- kappa[1L]
  if (p >= 2L) {
    for (m in 2L:p) {
      km <- kappa[m]
      prev <- Phi[[m-1L]]
      pm1 <- m - 1L
      cur <- numeric(m)
      for (j in 1L:pm1) cur[j] <- prev[j] - km * prev[pm1 - j + 1L]
      cur[m] <- km
      Phi[[m]] <- cur
    }
  }
  as.numeric(Phi[[p]])
}

#' @keywords internal
ar_to_pacf <- function(phi, eps = 1e-12) {
  p <- length(phi)
  if (p == 0L) return(numeric(0))
  a <- as.numeric(phi)
  if (!all(is.finite(a))) return(rep(NA_real_, p))
  kappa <- numeric(p)
  for (m in p:1L) {
    km <- a[m]
    kappa[m] <- km
    if (m == 1L) break
    # A near-singular Levinson step can make km non-finite; comparing NaN here
    # raised "missing value where TRUE/FALSE needed" rather than a usable error.
    den <- 1 - km * km
    if (!is.finite(den) || den < eps) den <- eps
    anew <- numeric(m - 1L)
    for (j in 1L:(m-1L)) anew[j] <- (a[j] + km * a[m - j]) / den
    a <- anew
  }
  kappa
}

#' @keywords internal
enforce_stationary_ar <- function(phi, bound = 0.99, margin = 1e-6) {
  if (length(phi) == 0L) return(phi)
  if (!all(is.finite(phi))) return(numeric(0))
  kap <- ar_to_pacf(phi)
  if (!all(is.finite(kap))) return(numeric(0))
  kap <- pmin(pmax(kap, -bound), bound)
  out <- pacf_to_ar(kap)
  if (!all(is.finite(out))) return(numeric(0))

  # Clamping the reflection coefficients bounds each |kappa| < 1, but at high
  # order the implied roots sit arbitrarily close to the unit circle: a PACF of
  # 0.99 repeated ten times yields coefficients above 45 with min|root| equal to
  # 1 to machine precision, which whitens by amplifying. Shrink geometrically
  # until the characteristic roots are strictly outside by a real margin.
  for (i in seq_len(60L)) {
    roots <- tryCatch(Mod(polyroot(c(1, -out))), error = function(e) NA_real_)
    if (!anyNA(roots) && min(roots) > 1 + margin) break
    out <- out * 0.9
    if (max(abs(out)) < 1e-12) return(numeric(0))
  }
  out
}
