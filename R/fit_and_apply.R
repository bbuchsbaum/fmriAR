# Internal helpers -------------------------------------------------------------

.sub_run_starts <- function(n_run, censor_idx_rel = integer()) {
  starts <- 1L
  if (length(censor_idx_rel)) {
    add <- censor_idx_rel + 1L
    add <- add[add <= n_run]
    starts <- sort(unique(c(starts, add)))
  }
  as.integer(starts - 1L)
}

.split_runs <- function(runs) {
  if (is.null(runs)) return(list(seq_along(integer(0))))
  runs <- as.integer(runs)
  split(seq_along(runs), runs)
}

new_whiten_plan <- function(phi, theta, order, runs, exact_first, method, pooling,
                            parcels = NULL, parcel_ids = NULL,
                            phi_by_parcel = NULL, theta_by_parcel = NULL) {
  structure(
    list(
      phi = phi,
      theta = theta,
      order = order,
      runs = runs,
      exact_first = exact_first,
      method = method,
      pooling = pooling,
      parcels = parcels,
      parcel_ids = parcel_ids,
      phi_by_parcel = phi_by_parcel,
      theta_by_parcel = theta_by_parcel
    ),
    class = "fmriAR_plan"
  )
}

.arma_innovations <- function(y, phi, theta) {
  Y <- matrix(as.numeric(y), ncol = 1L)
  X <- matrix(0, nrow = length(y), ncol = 1L)
  out <- arma_whiten_inplace(Y, X, phi, theta, run_starts = 0L,
                             exact_first_ar1 = FALSE, parallel = FALSE)
  drop(out$Y)
}

.run_avg_acvf <- function(mat, max_lag) {
  if (!is.matrix(mat)) mat <- as.matrix(mat)
  n <- nrow(mat)
  v <- ncol(mat)
  if (n == 0L || v == 0L) return(rep(0, max(0L, max_lag) + 1L))
  max_lag <- max(0L, min(as.integer(max_lag), n - 1L))
  centered <- sweep(mat, 2L, colMeans(mat))
  out <- numeric(max_lag + 1L)
  out[1] <- mean(centered^2)
  if (max_lag >= 1L) {
    for (lag in seq_len(max_lag)) {
      rows <- seq_len(n - lag)
      out[lag + 1L] <- mean(centered[rows + lag, , drop = FALSE] *
                             centered[rows, , drop = FALSE])
    }
  }
  out
}

.estimate_ar_series <- function(y, p_max) {
  y <- as.numeric(y)
  y_center <- y - mean(y)
  best <- list(score = Inf, phi = numeric(0), p = 0L)
  for (pp in 0L:p_max) {
    if (pp == 0L) {
      e <- y_center
      ll <- -length(e) * log(stats::var(e))
      k <- 1L
      phi <- numeric(0)
    } else {
      acvf <- stats::acf(y_center, lag.max = pp, plot = FALSE, type = "covariance")$acf
      gamma <- as.numeric(acvf)
      R <- stats::toeplitz(gamma[1:pp])
      r <- gamma[2:(pp + 1L)]
      phi_try <- tryCatch(drop(solve(R, r)), error = function(e) rep(0, pp))
      phi_try <- enforce_stationary_ar(phi_try, 0.99)
      e <- stats::filter(y_center, c(1, -phi_try), method = "recursive")
      e <- e[!is.na(e)]
      ll <- -length(e) * log(mean(e^2))
      k <- pp + 1L
      phi <- phi_try
    }
    n0 <- length(y_center)
    bic <- -2 * ll + k * log(n0)
    if (bic < best$score) best <- list(score = bic, phi = if (pp > 0) phi else numeric(0), p = pp)
  }
  list(phi = best$phi, order = c(p = best$p, q = 0L))
}

.full_run_starts <- function(runs, censor, n) {
  starts <- 1L
  if (!is.null(runs)) {
    runs <- as.integer(runs)
    starts <- sort(unique(c(starts, which(diff(runs) != 0L) + 1L)))
  }
  if (!is.null(censor) && length(censor)) {
    censor <- sort(unique(as.integer(censor)))
    extra <- censor + 1L
    extra <- extra[extra <= n]
    starts <- sort(unique(c(starts, extra)))
  }
  as.integer(starts - 1L)
}

.runs_from_starts0 <- function(run_starts0, n) {
  rs <- sort(unique(as.integer(run_starts0)))
  if (!length(rs) || rs[1] != 0L) stop("run_starts must include 0")
  if (rs[length(rs)] == n) rs <- rs[-length(rs)]
  if (!length(rs) || rs[length(rs)] >= n) stop("run_starts out of bounds")
  rs1 <- rs + 1L
  bounds <- c(rs1, n + 1L)
  runs <- integer(n)
  for (i in seq_along(rs1)) {
    runs[seq(rs1[i], bounds[i + 1L] - 1L)] <- i
  }
  runs
}

# Exported API -----------------------------------------------------------------

#' Fit an AR/ARMA noise model (run-aware) and return a whitening plan
#'
#' @param resid Numeric matrix (time x voxels) of residuals from an initial OLS fit.
#' @param Y Optional data matrix used to compute residuals when `resid` is omitted.
#' @param X Optional design matrix used with `Y` to compute residuals.
#' @param runs Optional integer vector of run identifiers.
#' @param method Either "ar" or "arma".
#' @param p AR order (integer or "auto" if method == "ar").
#' @param q MA order (integer).
#' @param p_max Maximum AR order when `p = "auto"`.
#' @param exact_first Apply exact AR(1) scaling at segment starts ("ar1" or "none").
#' @param pooling Combine parameters across runs or parcels ("global", "run", "parcel").
#' @param parcels Integer vector (length = ncol(resid)) giving fine parcel memberships when `pooling = "parcel"`.
#' @param parcel_sets Optional named list with entries `coarse`, `medium`, `fine` of equal length specifying nested parcel labels for multi-scale pooling.
#' @param multiscale Multi-scale pooling mode when `parcel_sets` is supplied ("pacf_weighted" or "acvf_pooled"), or `TRUE/FALSE` to toggle pooling.
#' @param ms_mode Explicit multiscale mode when `multiscale` is logical.
#' @param p_target Target AR order for multi-scale pooling (defaults to `p_max`).
#' @param beta Size exponent for multi-scale weights (default 0.5).
#' @param hr_iter Number of Hannan--Rissanen refinement iterations for ARMA.
#' @param step1 Preliminary high-order AR fit method for HR ("burg" or "yw").
#' @param parallel Reserved for future parallel estimation (logical).
#' @return An object of class `fmriAR_plan` used by [whiten_apply()].
#' @examples
#' # Generate example data with AR(1) structure
#' n_time <- 200
#' n_voxels <- 50
#' phi_true <- 0.5
#'
#' # Simulate residuals with AR(1) structure
#' resid <- matrix(0, n_time, n_voxels)
#' for (v in 1:n_voxels) {
#'   e <- rnorm(n_time)
#'   resid[1, v] <- e[1]
#'   for (t in 2:n_time) {
#'     resid[t, v] <- phi_true * resid[t-1, v] + e[t]
#'   }
#' }
#'
#' # Fit AR model
#' plan <- fit_noise(resid, method = "ar", p = 1)
#'
#' # With multiple runs
#' runs <- rep(1:2, each = 100)
#' plan_runs <- fit_noise(resid, runs = runs, method = "ar", pooling = "run")
#' @export
fit_noise <- function(resid = NULL,
                      Y = NULL,
                      X = NULL,
                      runs = NULL,
                      method = c("ar", "arma"),
                      p = "auto",
                      q = 0L,
                      p_max = 6L,
                      exact_first = c("ar1", "none"),
                      pooling = c("global", "run", "parcel"),
                      parcels = NULL,
                      parcel_sets = NULL,
                      multiscale = c("pacf_weighted", "acvf_pooled"),
                      ms_mode = NULL,
                      p_target = NULL,
                      beta = 0.5,
                      hr_iter = 0L,
                      step1 = c("burg", "yw"),
                      parallel = FALSE) {

  if (is.null(resid)) {
    if (!is.null(Y) && !is.null(X)) {
      if (!is.matrix(Y)) Y <- as.matrix(Y)
      if (!is.matrix(X)) X <- as.matrix(X)
      stopifnot(nrow(Y) == nrow(X))
      coef <- qr.solve(X, Y)
      resid <- Y - X %*% coef
    } else {
      stop("fit_noise: supply 'resid' or both 'Y' and 'X'")
    }
  }

  stopifnot(is.matrix(resid))
  method <- match.arg(method)
  exact_first <- match.arg(exact_first)
  pooling <- match.arg(pooling)
  step1 <- match.arg(step1)

  ms_modes <- c("pacf_weighted", "acvf_pooled")
  multiscale_mode <- NULL
  if (is.logical(multiscale)) {
    if (isTRUE(multiscale)) {
      multiscale_mode <- if (is.null(ms_mode)) "pacf_weighted" else match.arg(ms_mode, ms_modes)
    }
  } else {
    multiscale_mode <- match.arg(multiscale, ms_modes)
  }
  if (!is.null(ms_mode) && (!is.logical(multiscale) || isTRUE(multiscale))) {
    multiscale_mode <- match.arg(ms_mode, ms_modes)
  }

  n <- nrow(resid)
  Rsets <- if (is.null(runs)) list(seq_len(n)) else split(seq_len(n), as.integer(runs))
  run_mats <- lapply(Rsets, function(idx) resid[idx, , drop = FALSE])

  est_run <- function(mat) {
    if (method == "ar") {
      n_eff <- nrow(mat)
      if (n_eff <= 1L) {
        return(list(phi = numeric(0), theta = numeric(0), order = c(p = 0L, q = 0L)))
      }
      p_cap <- min(as.integer(p_max), n_eff - 1L)
      gamma <- .run_avg_acvf(mat, p_cap)
      best_phi <- numeric(0)
      best_order <- c(p = 0L, q = 0L)
      n_eff_log <- log(n_eff)
      sigma0 <- pmax(gamma[1], 1e-12)
      best_bic <- 2 * n_eff * log(sigma0) + n_eff_log
      if (p_cap >= 1L) {
        for (pp in seq_len(p_cap)) {
          gamma_pp <- gamma[seq_len(pp + 1L)]
          yw <- yw_from_acvf_fast(gamma_pp, pp)
          sigma2 <- pmax(yw$sigma2, 1e-12)
          bic <- 2 * n_eff * log(sigma2) + (pp + 1L) * n_eff_log
          if (bic < best_bic) {
            best_bic <- bic
            best_phi <- yw$phi
            best_order <- c(p = pp, q = 0L)
          }
        }
      }
      list(phi = best_phi, theta = numeric(0), order = best_order)
    } else {
      y_mean <- rowMeans(mat)
      pp <- if (identical(p, "auto")) min(2L, p_max) else as.integer(p)
      qq <- as.integer(q)
      hr_arma(y_mean, p = pp, q = qq, iter = as.integer(hr_iter), step1 = step1)
    }
  }

  if (pooling == "parcel") {
    if (!identical(method, "ar")) stop("Parcel pooling currently supports method = 'ar' only")
    stopifnot(!is.null(parcels))
    parcels <- as.integer(parcels)
    stopifnot(length(parcels) == ncol(resid))

    run_starts0 <- .full_run_starts(runs, censor = NULL, n = n)

    estimator <- function(y) .estimate_ar_series(y, p_max)
    M_fine <- .parcel_means(resid, parcels)

    target <- if (is.null(p_target)) p_max else min(as.integer(p_target), p_max)

    if (is.null(parcel_sets)) {
      est_f <- .ms_estimate_scale(M_fine, estimator, run_starts0)
      if (is.null(multiscale_mode) || target == 0L) {
        phi_parcel <- est_f$phi
      } else if (identical(multiscale_mode, "pacf_weighted")) {
        shrink <- 0.6
        kap_mat <- vapply(est_f$phi, function(phi) .ms_pad(ar_to_pacf(phi), target), numeric(target))
        avg_kap <- if (target > 0L) rowMeans(kap_mat, na.rm = TRUE) else numeric(0)
        avg_kap <- pmin(pmax(avg_kap, -0.99), 0.99)
        phi_parcel <- lapply(est_f$phi, function(phi) {
          kap_f <- .ms_pad(ar_to_pacf(phi), target)
          kap_mix <- (1 - shrink) * kap_f + shrink * avg_kap
          pacf_to_ar(pmin(pmax(kap_mix, -0.99), 0.99))
        })
      } else {
        shrink <- 0.6
        acvf_mat <- vapply(est_f$acvf, function(g) .ms_pad(g, target + 1L), numeric(target + 1L))
        avg_g <- rowMeans(acvf_mat, na.rm = TRUE)
        phi_parcel <- lapply(est_f$acvf, function(g) {
          g_pad <- .ms_pad(g, target + 1L)
          g_mix <- (1 - shrink) * g_pad + shrink * avg_g
          yw <- yw_from_acvf_fast(g_mix, target)
          enforce_stationary_ar(yw$phi)
        })
      }
    } else {
      required_keys <- c("coarse", "medium", "fine")
      stopifnot(all(required_keys %in% names(parcel_sets)))
      parcels_coarse <- as.integer(parcel_sets$coarse)
      parcels_medium <- as.integer(parcel_sets$medium)
      parcels_fine <- as.integer(parcel_sets$fine)
      stopifnot(length(parcels_coarse) == ncol(resid))
      stopifnot(length(parcels_medium) == ncol(resid))
      stopifnot(all(parcels_fine == parcels))

      M_coarse <- .parcel_means(resid, parcels_coarse)
      M_medium <- .parcel_means(resid, parcels_medium)
      est_c <- .ms_estimate_scale(M_coarse, estimator, run_starts0)
      est_m <- .ms_estimate_scale(M_medium, estimator, run_starts0)
      est_f <- .ms_estimate_scale(M_fine, estimator, run_starts0)

      parents <- .ms_parent_maps(parcels_fine, parcels_medium, parcels_coarse)
      sizes <- list(
        n_t = nrow(resid),
        n_runs = if (is.null(runs)) 1L else length(unique(as.integer(runs))),
        beta = beta,
        coarse = as.list(table(parcels_coarse)),
        medium = as.list(table(parcels_medium)),
        fine = as.list(table(parcels_fine))
      )
      disp_list <- list(
        coarse = .ms_dispersion(resid, parcels_coarse),
        medium = .ms_dispersion(resid, parcels_medium),
        fine = .ms_dispersion(resid, parcels_fine)
      )
      if (is.null(multiscale_mode)) {
        phi_parcel <- est_f$phi
      } else {
        phi_parcel <- .ms_combine_to_fine(
          phi_by_coarse = est_c$phi,
          phi_by_medium = est_m$phi,
          phi_by_fine   = est_f$phi,
          acvf_by_coarse = if (identical(multiscale_mode, "acvf_pooled")) est_c$acvf else NULL,
          acvf_by_medium = if (identical(multiscale_mode, "acvf_pooled")) est_m$acvf else NULL,
          acvf_by_fine   = if (identical(multiscale_mode, "acvf_pooled")) est_f$acvf else NULL,
          parents = parents,
          sizes = sizes,
          disp_list = disp_list,
          p_target = target,
          mode = multiscale_mode
        )
      }
    }

    if (is.null(multiscale_mode) && !is.null(p_target) && target > 0L) {
      phi_parcel <- mapply(function(phi, g) {
        g_pad <- .ms_pad(g, target + 1L)
        yw <- yw_from_acvf_fast(g_pad, target)
        enforce_stationary_ar(yw$phi)
      }, phi_parcel, est_f$acvf, SIMPLIFY = FALSE)
    }

    theta_parcel <- setNames(vector("list", length(phi_parcel)), names(phi_parcel))
    order_vec <- c(p = max(vapply(phi_parcel, length, 0L)), q = 0L)

    parcel_ids <- names(phi_parcel)
    if (is.null(parcel_ids)) parcel_ids <- as.character(sort(unique(parcels)))
    return(new_whiten_plan(
      phi = NULL,
      theta = NULL,
      order = order_vec,
      runs = runs,
      exact_first = (exact_first == "ar1"),
      method = method,
      pooling = "parcel",
      parcels = parcels,
      parcel_ids = parcel_ids,
      phi_by_parcel = phi_parcel,
      theta_by_parcel = theta_parcel
    ))
  }

  estimates <- lapply(run_mats, est_run)

  if (pooling == "global") {
    lens <- vapply(Rsets, length, 0L)
    pmax_len <- max(vapply(estimates, function(e) length(e$phi), 0L))
    qmax_len <- max(vapply(estimates, function(e) length(e$theta), 0L))
    Phi <- matrix(0, length(estimates), pmax_len)
    Th <- matrix(0, length(estimates), qmax_len)
    for (i in seq_along(estimates)) {
      if (length(estimates[[i]]$phi))   Phi[i, seq_along(estimates[[i]]$phi)]   <- estimates[[i]]$phi
      if (length(estimates[[i]]$theta)) Th[i, seq_along(estimates[[i]]$theta)] <- estimates[[i]]$theta
    }
    w <- lens / sum(lens)
    phi_list <- list(as.numeric(drop(crossprod(w, Phi))))
    theta_list <- list(as.numeric(drop(crossprod(w, Th))))
  } else {
    phi_list <- lapply(estimates, `[[`, "phi")
    theta_list <- lapply(estimates, `[[`, "theta")
  }

  order_vec <- if (pooling == "global") {
    c(p = length(phi_list[[1]]), q = length(theta_list[[1]]))
  } else {
    c(p = max(vapply(phi_list, length, 0L)), q = max(vapply(theta_list, length, 0L)))
  }

  new_whiten_plan(
    phi = phi_list,
    theta = theta_list,
    order = order_vec,
    runs = runs,
    exact_first = (exact_first == "ar1"),
    method = method,
    pooling = pooling
  )
}

#' Apply a whitening plan to design and data matrices
#'
#' @param plan Whitening plan from [fit_noise()].
#' @param X Numeric matrix of predictors (time x regressors).
#' @param Y Numeric matrix of data (time x voxels).
#' @param runs Optional run labels.
#' @param run_starts Optional 0-based run start indices (alternative to `runs`).
#' @param censor Optional indices of censored TRs (1-based); filter resets after gaps.
#' @param parcels Optional parcel labels (length = ncol(Y)) when using parcel plans.
#' @param inplace Modify inputs in place (logical).
#' @param parallel Use OpenMP parallelism if available.
#' @return List with whitened data. Parcel plans return `X_by` per parcel; others return a single `X` matrix.
#' @examples
#' # Create example design matrix and data
#' n_time <- 200
#' n_pred <- 3
#' n_voxels <- 50
#' X <- matrix(rnorm(n_time * n_pred), n_time, n_pred)
#' Y <- X %*% matrix(rnorm(n_pred * n_voxels), n_pred, n_voxels) +
#'      matrix(rnorm(n_time * n_voxels), n_time, n_voxels)
#'
#' # Fit noise model from residuals
#' residuals <- Y - X %*% solve(crossprod(X), crossprod(X, Y))
#' plan <- fit_noise(residuals, method = "ar", p = 2)
#'
#' # Apply whitening
#' whitened <- whiten_apply(plan, X, Y)
#' Xw <- whitened$X
#' Yw <- whitened$Y
#' @export
whiten_apply <- function(plan, X, Y, runs = NULL, run_starts = NULL, censor = NULL, parcels = NULL,
                         inplace = FALSE, parallel = TRUE) {
  stopifnot(inherits(plan, "fmriAR_plan"))
  if (!is.matrix(X)) X <- as.matrix(X)
  if (!is.matrix(Y)) Y <- as.matrix(Y)
  if (anyNA(X) || anyNA(Y)) {
    stop("whiten_apply: NA values detected in X or Y")
  }
  n <- nrow(X)
  stopifnot(nrow(Y) == n)

  if (!is.null(run_starts)) run_starts <- as.integer(run_starts)
  if (is.null(runs) && !is.null(run_starts)) {
    runs <- .runs_from_starts0(run_starts, n)
  }
  if (is.null(runs) && !is.null(plan$runs) && length(plan$runs) == n) runs <- plan$runs
  if (is.null(runs)) runs <- rep_len(1L, n)

  if (identical(plan$pooling, "parcel")) {
    parcels_vec <- plan$parcels
    if (!is.null(parcels)) {
      stopifnot(length(parcels) == ncol(Y))
      parcels_vec <- as.integer(parcels)
    }
    stopifnot(!is.null(parcels_vec))
    stopifnot(length(parcels_vec) == ncol(Y))

    run_starts_vec <- .full_run_starts(runs, censor, n)
    parcel_ids <- if (!is.null(plan$parcel_ids)) plan$parcel_ids else sort(unique(parcels_vec))
    phi_by <- plan$phi_by_parcel
    theta_by <- plan$theta_by_parcel

    Yw <- matrix(NA_real_, n, ncol(Y))
    X_by <- setNames(vector("list", length(parcel_ids)), as.character(parcel_ids))
    X_base <- X

    for (pid in parcel_ids) {
      cols <- which(parcels_vec == pid)
      if (!length(cols)) next
      key <- as.character(pid)
      phi <- phi_by[[key]]
      if (is.null(phi)) phi <- numeric(0)
      theta <- theta_by[[key]]
      if (is.null(theta)) theta <- numeric(0)
      Y_sub <- Y[, cols, drop = FALSE]
      X_sub <- X_base
      out <- arma_whiten_inplace(
        Y = Y_sub,
        X = X_sub,
        phi = phi,
        theta = theta,
        run_starts = run_starts_vec,
        exact_first_ar1 = isTRUE(plan$exact_first),
        parallel = parallel
      )
      Yw[, cols] <- out$Y
      X_by[[key]] <- out$X
    }

    if (inplace) {
      Y[,] <- Yw
      return(invisible(list(X = NULL, X_by = X_by, Y = Y)))
    }
    return(list(X = NULL, X_by = X_by, Y = Yw))
  }

  rsplits <- split(seq_len(n), as.integer(runs))

  censor_by_run <- lapply(rsplits, function(idx) integer(0L))
  if (!is.null(censor)) {
    censor <- as.integer(censor)
    for (ri in seq_along(rsplits)) {
      idx <- rsplits[[ri]]
      c_in <- intersect(censor, idx)
      if (length(c_in)) censor_by_run[[ri]] <- as.integer(c_in - min(idx) + 1L)
    }
  }

  phi_list <- if (length(plan$phi) == 1L) rep(plan$phi, length(rsplits)) else plan$phi
  theta_list <- if (length(plan$theta) == 1L) rep(plan$theta, length(rsplits)) else plan$theta

  Xw_list <- vector("list", length(rsplits))
  Yw_list <- vector("list", length(rsplits))

  for (ri in seq_along(rsplits)) {
    idx <- rsplits[[ri]]
    Xr <- X[idx, , drop = FALSE]
    Yr <- Y[idx, , drop = FALSE]
    rs <- .sub_run_starts(n_run = nrow(Xr), censor_idx_rel = censor_by_run[[ri]])
    out <- arma_whiten_inplace(
      Yr,
      Xr,
      phi = phi_list[[ri]],
      theta = theta_list[[ri]],
      run_starts = rs,
      exact_first_ar1 = isTRUE(plan$exact_first),
      parallel = parallel
    )
    Xw_list[[ri]] <- out$X
    Yw_list[[ri]] <- out$Y
  }

  Xw <- do.call(rbind, Xw_list)
  Yw <- do.call(rbind, Yw_list)

  if (inplace) {
    X[,] <- Xw
    Y[,] <- Yw
    invisible(list(X = X, Y = Y))
  } else {
    list(X = Xw, Y = Yw)
  }
}

#' Fit and apply whitening in one call
#'
#' @param X Design matrix (time x regressors).
#' @param Y Data matrix (time x voxels).
#' @param runs Optional run labels.
#' @param censor Optional censor indices.
#' @param ... Additional parameters passed to [fit_noise()].
#' @return List with whitened `X` and `Y` matrices.
#' @examples
#' # Create example data
#' n_time <- 200
#' n_pred <- 3
#' n_voxels <- 50
#' X <- matrix(rnorm(n_time * n_pred), n_time, n_pred)
#' Y <- X %*% matrix(rnorm(n_pred * n_voxels), n_pred, n_voxels) +
#'      matrix(rnorm(n_time * n_voxels, sd = 2), n_time, n_voxels)
#'
#' # One-step whitening
#' whitened <- whiten(X, Y, method = "ar", p = 2)
#' @export
whiten <- function(X, Y, runs = NULL, censor = NULL, ...) {
  res <- Y - X %*% qr.solve(X, Y)
  plan <- fit_noise(resid = res, runs = runs, ...)
  whiten_apply(plan, X, Y, runs = runs, censor = censor)
}
