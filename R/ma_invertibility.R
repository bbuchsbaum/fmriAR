#' @keywords internal
.coeff_from_roots <- function(roots) {
  coeffs <- 1 + 0i
  for (r in roots) coeffs <- c(coeffs, 0 + 0i) + (-1 / r) * c(0 + 0i, coeffs)
  Re(coeffs)
}

#' @keywords internal
enforce_invertible_ma <- function(theta, tol = 1e-8, bound = 0.99) {
  q <- length(theta)
  if (q == 0L) return(theta)
  if (!all(is.finite(theta))) return(numeric(q))
  r <- polyroot(c(1, theta))
  if (!length(r)) return(theta)
  r_new <- r
  # Reflect roots inside the unit circle. Reflection maps a root ON the circle
  # to itself, so it cannot repair a unit root (theta = -1 came back as -1, and
  # the whitening recursion then integrates like a random walk). Push any root
  # still within 1/bound of the origin out to that radius, mirroring the PACF
  # bound applied to the AR side.
  r_min <- 1 / bound
  changed <- FALSE
  for (i in seq_along(r)) {
    if (Mod(r[i]) <= 1 + tol) {
      r_new[i] <- 1 / Conj(r[i])
      changed <- TRUE
    }
    if (Mod(r_new[i]) < r_min) {
      r_new[i] <- r_new[i] * (r_min / Mod(r_new[i]))
      changed <- TRUE
    }
  }
  if (!changed) return(theta)
  as.numeric(.coeff_from_roots(r_new)[-1L])
}
