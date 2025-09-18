#' @keywords internal
.coeff_from_roots <- function(roots) {
  coeffs <- 1 + 0i
  for (r in roots) coeffs <- c(coeffs, 0 + 0i) + (-1 / r) * c(0 + 0i, coeffs)
  Re(coeffs)
}

#' @keywords internal
enforce_invertible_ma <- function(theta, tol = 1e-8) {
  q <- length(theta)
  if (q == 0L) return(theta)
  r <- polyroot(c(1, theta))
  if (!length(r)) return(theta)
  r_new <- r
  for (i in seq_along(r)) {
    if (Mod(r[i]) <= 1 + tol) r_new[i] <- 1 / Conj(r[i])
  }
  as.numeric(.coeff_from_roots(r_new)[-1L])
}
