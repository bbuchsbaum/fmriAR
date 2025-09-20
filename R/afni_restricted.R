
# R/afni_restricted.R — internal AFNI-style restricted AR helpers
# We expose a plan with plan$method = "afni" to make the provenance clear.
# This produces standard fmriAR_plan objects that work with whiten_apply().

#' @keywords internal
.afni_phi_ar3 <- function(a, r1, t1) {
  stopifnot(is.numeric(a), is.numeric(r1), is.numeric(t1))
  a  <- pmin(pmax(a,  0), 0.95)
  r1 <- pmin(pmax(r1, 0), 0.95)
  t1 <- pmin(pmax(t1, 0), pi)
  c1 <- cos(t1)
  p1 <- a + 2 * r1 * c1
  p2 <- -2 * a * r1 * c1 - r1^2
  p3 <- a * r1^2
  c(p1, p2, p3)
}

#' @keywords internal
.afni_phi_ar5 <- function(a, r1, t1, r2, t2) {
  stopifnot(is.numeric(a), is.numeric(r1), is.numeric(t1), is.numeric(r2), is.numeric(t2))
  a  <- pmin(pmax(a,  0), 0.95)
  r1 <- pmin(pmax(r1, 0), 0.95)
  r2 <- pmin(pmax(r2, 0), 0.95)
  t1 <- pmin(pmax(t1, 0), pi)
  t2 <- pmin(pmax(t2, 0), pi)
  c1 <- cos(t1); c2 <- cos(t2)
  # from AFNI code (see user paste)
  p1 <-  2*r1*c1 + 2*r2*c2 + a
  p2 <- -4*r1*r2*c1*c2 - 2*a*(r1*c1 + r2*c2) - r1^2 - r2^2
  p3 <- a*(r1^2 + r2^2 + 4*r1*r2*c1*c2) + 2*r1*r2*(r2*c1 + r1*c2)
  p4 <- -2*a*r1*r2*(r2*c1 + r1*c2) - (r1^2)*(r2^2)
  p5 <- a*(r1^2)*(r2^2)
  c(p1, p2, p3, p4, p5)
}

# local helper: 0-based run starts from integer run labels
.afni_run_starts0 <- function(runs, n) {
  if (is.null(runs)) return(0L)
  lab <- as.integer(runs)
  idx <- split(seq_len(n), lab)
  starts <- vapply(idx, function(ix) min(ix)-1L, 0L)  # 0-based
  as.integer(starts)
}

#' Build an AFNI-style restricted AR plan from root parameters
#'
#' @param resid (n x v) residual matrix (used only if estimate_ma1=TRUE)
#' @param runs integer vector length n (optional)
#' @param parcels integer vector length v (optional; if provided, plan pooling='parcel')
#' @param p either 3 or 5
#' @param roots either a single list with elements named as needed
#'        - for p=3: list(a, r1, t1, vrt = 1.0)
#'        - for p=5: list(a, r1, t1, r2, t2, vrt = 1.0)
#'      or a named list of such lists keyed by parcel id (character) for per-parcel specs.
#' @param estimate_ma1 logical, if TRUE estimate MA(1) on AR residuals to mimic AFNI's additive white
#' @param exact_first apply exact AR(1) scaling at segment starts (harmless here; default TRUE)
#' @keywords internal
afni_restricted_plan <- function(resid, runs = NULL, parcels = NULL,
                                 p = 3L, roots,
                                 estimate_ma1 = TRUE,
                                 exact_first = TRUE) {
  stopifnot(p %in% c(3L,5L))
  n <- nrow(resid); v <- ncol(resid)

  as_phi <- function(spec) {
    if (p == 3L) .afni_phi_ar3(spec$a, spec$r1, spec$t1) else
      .afni_phi_ar5(spec$a, spec$r1, spec$t1, spec$r2, spec$t2)
  }

  # construct phi lists
  if (is.null(parcels)) {
    # global/run plan from a single spec
    stopifnot(is.list(roots), !is.null(roots$a))
    phi <- as_phi(roots)
    phi_list <- list(phi)
    theta_list <- list(numeric(0))  # ensure theta placeholder exists
    order_vec <- c(p = length(phi), q = 0L)
    plan <- new_whiten_plan(phi = phi_list, theta = theta_list, order = order_vec,
                            runs = runs, exact_first = isTRUE(exact_first),
                            method = "afni", pooling = if (is.null(runs)) "global" else "run")
    # optional MA(1) estimation on AR residuals (global)
    if (isTRUE(estimate_ma1)) {
      rs0 <- .afni_run_starts0(runs, n)
      ymean <- rowMeans(resid)  # global mean series
      # AR residuals via C++ whitener with q=0
      out <- arma_whiten_inplace(matrix(ymean, ncol = 1L),
                                 matrix(0, nrow = n, ncol = 1L),
                                 phi = phi, theta = numeric(0),
                                 run_starts = rs0,
                                 exact_first_ar1 = FALSE, parallel = FALSE)
      s <- drop(out$Y)
      est <- hr_arma_fit_cpp(s, p = 0L, q = 1L, p_big = 0L, iter = 0L)
      if (isTRUE(est$ok) || is.null(est$ok)) {
        plan$theta <- list(as.numeric(est$theta))
        plan$order["q"] <- 1L
      } else {
        plan$theta <- list(numeric(0))
      }
    } else {
      plan$theta <- list(numeric(0))
    }
    return(plan)
  }

  # parcel mode
  parcels <- as.integer(parcels); stopifnot(length(parcels) == v)
  Pids <- sort(unique(parcels))
  phi_by <- setNames(vector("list", length(Pids)), as.character(Pids))
  th_by  <- setNames(vector("list", length(Pids)), as.character(Pids))

  # normalize roots input to per-parcel list
  roots_by <- if (is.list(roots) && !is.null(roots$a)) {
    # single spec for all parcels
    setNames(rep(list(roots), length(Pids)), as.character(Pids))
  } else {
    roots  # expect named list keyed by parcel id
  }

  rs0 <- .afni_run_starts0(runs, n)

  for (pid in Pids) {
    key <- as.character(pid)
    spec <- roots_by[[key]]
    if (is.null(spec)) spec <- roots_by[[1]] %||% roots_by[[names(roots_by)[1]]]
    phi <- as_phi(spec)
    phi_by[[key]] <- phi

    if (isTRUE(estimate_ma1)) {
      cols <- which(parcels == pid)
      ymean <- if (length(cols) == 1L) resid[, cols] else rowMeans(resid[, cols, drop = FALSE])
      out <- arma_whiten_inplace(matrix(ymean, ncol = 1L),
                                 matrix(0, nrow = n, ncol = 1L),
                                 phi = phi, theta = numeric(0),
                                 run_starts = rs0,
                                 exact_first_ar1 = FALSE, parallel = FALSE)
      s <- drop(out$Y)
      est <- hr_arma_fit_cpp(s, p = 0L, q = 1L, p_big = 0L, iter = 0L)
      if (isTRUE(est$ok) || is.null(est$ok)) th_by[[key]] <- as.numeric(est$theta)
    }
  }

  order_vec <- c(p = p, q = if (isTRUE(estimate_ma1)) 1L else 0L)
  new_whiten_plan(
    phi = NULL, theta = NULL, order = order_vec, runs = runs,
    exact_first = isTRUE(exact_first), method = "afni", pooling = "parcel",
    parcels = parcels, parcel_ids = Pids,
    phi_by_parcel = phi_by, theta_by_parcel = th_by
  )
}
