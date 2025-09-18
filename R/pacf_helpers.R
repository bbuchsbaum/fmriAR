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
  kappa <- numeric(p)
  for (m in p:1L) {
    km <- a[m]
    kappa[m] <- km
    if (m == 1L) break
    den <- 1 - km * km
    if (den < eps) den <- eps
    anew <- numeric(m - 1L)
    for (j in 1L:(m-1L)) anew[j] <- (a[j] + km * a[m - j]) / den
    a <- anew
  }
  kappa
}

#' @keywords internal
enforce_stationary_ar <- function(phi, bound = 0.99) {
  if (length(phi) == 0L) return(phi)
  kap <- ar_to_pacf(phi)
  kap <- pmin(pmax(kap, -bound), bound)
  pacf_to_ar(kap)
}
