
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

# AFNI (armacor.c) rejects negative pole parameters rather than clamping them;
# clamping silently turned a sign error into a different noise model. Values
# above 0.95 and angles outside [0, pi] are clamped, as AFNI does.
.afni_check_spec <- function(spec, p) {
  need <- if (p == 3L) c("a", "r1", "t1") else c("a", "r1", "t1", "r2", "t2")
  miss <- setdiff(need, names(spec))
  if (length(miss)) {
    stop("AFNI root spec is missing: ", paste(miss, collapse = ", "), call. = FALSE)
  }
  vals <- unlist(spec[need])
  if (!is.numeric(vals) || any(!is.finite(vals))) {
    stop("AFNI root parameters must be finite numbers", call. = FALSE)
  }
  mods <- intersect(c("a", "r1", "r2"), need)
  if (any(unlist(spec[mods]) < 0)) {
    stop("AFNI root parameters a, r1, r2 must be >= 0 (AFNI rejects negative ",
         "values)", call. = FALSE)
  }
  vrt <- spec$vrt %||% 1
  if (length(vrt) != 1L || !is.finite(vrt) || vrt <= 0) {
    stop("'vrt' must be a single number in (0, 1]", call. = FALSE)
  }
  invisible(TRUE)
}

# MA part implied by AFNI's additive white noise. AFNI models the noise as an
# AR(p) signal plus white noise, vrt = s^2 / (s^2 + w^2). Then
#   phi(B) y = e + phi(B) w,
# whose right-hand side is an MA(p) with autocovariance
#   c_k = delta_k sigma_e^2 + sigma_w^2 sum_j phi~_j phi~_{j+k},
# phi~ = (1, -phi). Spectral factorisation (roots outside the unit circle)
# gives the invertible theta, so ARMA(p, p) reproduces AFNI's correlations
# exactly: rho_k = vrt * rho_AR(k) for k >= 1.
.afni_theta_from_vrt <- function(phi, vrt) {
  p <- length(phi)
  if (!p || vrt >= 1) return(numeric(0))
  g <- arma_acvf_cpp(phi, numeric(0), 0L)       # AR variance, sigma_e = 1
  if (!is.finite(g[1])) return(numeric(0))
  s2w <- g[1] * (1 - vrt) / vrt
  pt <- c(1, -phi)
  ck <- vapply(0:p, function(k) sum(pt[seq_len(p + 1L - k)] * pt[seq_len(p + 1L - k) + k]), 0)
  ck <- s2w * ck
  ck[1] <- ck[1] + 1
  coefs <- c(rev(ck[-1]), ck)                   # z^0 .. z^(2p), symmetric
  r <- polyroot(coefs)
  r_out <- r[order(-Mod(r))][seq_len(p)]        # the p roots outside |z| = 1
  as.numeric(.coeff_from_roots(r_out)[-1L])
}

# local helper: 0-based run starts from run labels
.afni_run_starts0 <- function(runs, n) {
  if (is.null(runs)) return(0L)
  # .run_codes renumbers labels 1..K in block order, so starts come back
  # ascending-in-time even when labels are not ascending; the old
  # split(, as.integer(runs)) returned them in label-sorted order and
  # collapsed character labels to a single NA group.
  codes <- .run_codes(runs, n)
  as.integer(c(1L, which(diff(codes) != 0L) + 1L) - 1L)
}

# MA(1) for the additive-white-noise part of AFNI's restricted model, estimated
# on the AR-filtered VOXELS, pooled. Fitting it on the filtered voxel-mean
# series (as earlier versions did) measures whatever is shared across voxels
# rather than the voxel-level noise being whitened. The first p frames of each
# run are dropped because the AR filter has no pre-sample there.
.afni_ma1_pooled <- function(resid_cols, phi, rs0) {
  n <- nrow(resid_cols)
  filt <- arma_whiten_inplace(resid_cols + 0, matrix(0, n, 0L),
                              phi = phi, theta = numeric(0),
                              run_starts = rs0, exact_first = FALSE,
                              parallel = TRUE, n_threads = .n_threads())$Y
  ends <- c(rs0[-1L], n)
  units <- lapply(seq_along(rs0), function(i) {
    first <- rs0[i] + 1L + length(phi)
    rows <- if (first <= ends[i]) seq.int(first, ends[i]) else integer(0)
    list(mat = filt[rows, , drop = FALSE], starts0 = 0L)
  })
  units <- units[vapply(units, function(u) nrow(u$mat) >= 10L, TRUE)]
  if (!length(units)) return(NULL)
  est <- suppressWarnings(.hr_arma_pooled(units, p = 0L, q = 1L))
  if (!isTRUE(est$ok)) return(NULL)
  as.numeric(est$theta)
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
#'      As in AFNI, `a`, `r1`, `r2` must be non-negative (values above 0.95
#'      are clamped to 0.95) and angles are clamped to `[0, pi]`. `vrt` is
#'      AFNI's signal-to-total variance ratio `s^2 / (s^2 + w^2)` for additive
#'      white noise `w`; it is used when `estimate_ma1 = FALSE`.
#' @param estimate_ma1 logical. If `TRUE`, estimate an MA(1) term from the
#'   AR-filtered voxels to stand in for AFNI's additive white noise (`vrt` is
#'   then ignored). If `FALSE`, the white-noise share is taken from `vrt`: for
#'   `vrt < 1` the plan is the exact ARMA(p, p) equivalent of AR(p) plus white
#'   noise, so its autocorrelation is `vrt` times the AR autocorrelation at
#'   every non-zero lag; `vrt = 1` (default) gives the pure AR(p).
#' @param exact_first apply exact AR(1) scaling at segment starts (harmless here; default TRUE)
#' @return An `fmriAR_plan` with `method = "afni"` that can be supplied to
#'   [whiten_apply()].
#' @examples NULL
#' @export
afni_restricted_plan <- function(resid, runs = NULL, parcels = NULL,
                                 p = 3L, roots,
                                 estimate_ma1 = TRUE,
                                 exact_first = TRUE) {
  stopifnot(p %in% c(3L,5L))
  n <- nrow(resid); v <- ncol(resid)

  as_phi <- function(spec) {
    .afni_check_spec(spec, p)
    if (p == 3L) .afni_phi_ar3(spec$a, spec$r1, spec$t1) else
      .afni_phi_ar5(spec$a, spec$r1, spec$t1, spec$r2, spec$t2)
  }
  warned <- FALSE
  # MA part from vrt when it is not estimated. AFNI treats vrt <= 0.01 as
  # white noise.
  model_from_spec <- function(spec) {
    phi <- as_phi(spec)
    vrt <- min(spec$vrt %||% 1, 1)
    if (isTRUE(estimate_ma1) && vrt < 1 && !warned) {
      warning("afni_restricted_plan: 'vrt' is ignored when estimate_ma1 = TRUE; ",
              "set estimate_ma1 = FALSE to use it", call. = FALSE)
      warned <<- TRUE
    }
    if (vrt <= 0.01) return(list(phi = numeric(0), theta = numeric(0)))
    list(phi = phi,
         theta = if (isTRUE(estimate_ma1)) numeric(0) else .afni_theta_from_vrt(phi, vrt))
  }

  # construct phi lists
  if (is.null(parcels)) {
    # global/run plan from a single spec
    stopifnot(is.list(roots), !is.null(roots$a))
    mod <- model_from_spec(roots)
    phi <- mod$phi
    phi_list <- list(phi)
    theta_list <- list(mod$theta)
    order_vec <- c(p = length(phi), q = length(mod$theta))
    plan <- new_whiten_plan(phi = phi_list, theta = theta_list, order = order_vec,
                            runs = runs, exact_first = isTRUE(exact_first),
                            method = "afni", pooling = if (is.null(runs)) "global" else "run")
    # optional MA(1) estimation on AR residuals (global)
    if (isTRUE(estimate_ma1)) {
      rs0 <- .afni_run_starts0(runs, n)
      th <- .afni_ma1_pooled(as.matrix(resid), phi, rs0)
      if (!is.null(th)) {
        plan$theta <- list(th)
        plan$order["q"] <- 1L
      } else {
        plan$theta <- list(numeric(0))
      }
    }
    return(plan)
  }

  # parcel mode
  parcels <- .parcel_codes(parcels); stopifnot(length(parcels) == v)
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
    mod <- model_from_spec(spec)
    phi <- mod$phi
    phi_by[[key]] <- phi
    th_by[key] <- list(mod$theta)

    if (isTRUE(estimate_ma1)) {
      cols <- which(parcels == pid)
      th <- .afni_ma1_pooled(as.matrix(resid)[, cols, drop = FALSE], phi, rs0)
      if (!is.null(th)) th_by[[key]] <- th
    }
  }

  order_vec <- c(p = max(vapply(phi_by, length, 0L)),
                 q = max(0L, vapply(th_by, length, 0L)))
  new_whiten_plan(
    phi = NULL, theta = NULL, order = order_vec, runs = runs,
    exact_first = isTRUE(exact_first), method = "afni", pooling = "parcel",
    parcels = parcels, parcel_ids = Pids,
    phi_by_parcel = phi_by, theta_by_parcel = th_by
  )
}
