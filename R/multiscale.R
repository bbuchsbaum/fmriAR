# Multi-scale spatial pooling helpers -----------------------------------------

.ms_parent_maps <- function(parcels_fine, parcels_medium, parcels_coarse) {
  pf <- as.integer(parcels_fine)
  pm <- as.integer(parcels_medium)
  pc <- as.integer(parcels_coarse)
  stopifnot(length(pf) == length(pm), length(pf) == length(pc))

  fine_ids <- sort(unique(pf))
  parent_medium <- integer(max(fine_ids))
  parent_coarse <- integer(max(fine_ids))

  for (fid in fine_ids) {
    idx <- which(pf == fid)
    mids <- pm[idx]
    cids <- pc[idx]
    mid <- as.integer(names(which.max(table(mids))))
    cid <- as.integer(names(which.max(table(cids))))
    parent_medium[fid] <- mid
    parent_coarse[fid] <- cid
  }

  list(parent_medium = parent_medium, parent_coarse = parent_coarse)
}

.ms_dispersion <- function(resid, parcels) {
  parcels <- as.integer(parcels)
  vvar <- apply(resid, 2L, stats::var)
  # Relative dispersion of voxel variances within each parcel. The raw MAD is in
  # squared data units, so the weights 1 / (1 + disp) -- and with them phi --
  # changed when the data were merely rescaled. Normalising by the median
  # variance of all voxels keeps the weights unit-free while leaving them
  # unchanged for unit-variance data.
  scale <- stats::median(vvar)
  if (!is.finite(scale) || scale <= 0) scale <- 1
  tapply(vvar / scale, parcels, function(z) stats::mad(z, constant = 1))
}

# Order for the Yule-Walker solve in acvf_pooled mode. An explicit p is fitted
# at the pooling target. Under p = "auto", fitting every parcel at p_max
# discards the order selection entirely; use the largest order selected at any
# contributing scale, which is the effective order pacf_weighted mode produces.
.ms_yw_order <- function(p, target, phi_lists) {
  target <- as.integer(target)
  if (!identical(p, "auto")) return(target)
  sel <- unlist(lapply(phi_lists, function(l) vapply(l, length, 0L)))
  min(target, if (length(sel)) max(sel) else 0L)
}

.ms_weights <- function(n_t, n_runs, sizes, disp, beta = 0.5, eps = 1e-8) {
  s <- (n_t * n_runs) * (sizes ^ beta)
  h <- 1 / (1 + pmax(disp, 0))
  w <- s * h
  w[w < eps] <- eps
  w
}

.ms_pad <- function(x, len) {
  if (length(x) >= len) return(x[seq_len(len)])
  c(x, rep(0, len - length(x)))
}

.segment_acvf <- function(y, run_starts0, lag_max, unbiased = FALSE, center = TRUE) {
  segmented_acvf_fast(y, run_starts0, max_lag = lag_max, unbiased = unbiased, center = center)
}

.ms_combine_to_fine <- function(phi_by_coarse, phi_by_medium, phi_by_fine,
                                acvf_by_coarse = NULL, acvf_by_medium = NULL, acvf_by_fine = NULL,
                                parents, sizes, disp_list, p_target,
                                mode = c("pacf_weighted", "acvf_pooled"),
                                kappa_clip = 0.99, yw_order = p_target) {
  mode <- match.arg(mode)
  pids_fine <- sort(as.integer(names(phi_by_fine)))
  out_phi <- setNames(vector("list", length(pids_fine)), as.character(pids_fine))

  for (fid in pids_fine) {
    mid <- parents$parent_medium[fid]
    cid <- parents$parent_coarse[fid]
    key_f <- as.character(fid)
    key_m <- as.character(mid)
    key_c <- as.character(cid)

    size_vec <- c(
      coarse = as.numeric(sizes$coarse[[key_c]]),
      medium = as.numeric(sizes$medium[[key_m]]),
      fine   = as.numeric(sizes$fine[[key_f]])
    )
    disp_vec <- c(
      coarse = as.numeric(disp_list$coarse[[key_c]]),
      medium = as.numeric(disp_list$medium[[key_m]]),
      fine   = as.numeric(disp_list$fine[[key_f]])
    )
    w <- .ms_weights(sizes$n_t, sizes$n_runs, size_vec, disp_vec, beta = sizes$beta)
    w <- w / sum(w)

    if (mode == "pacf_weighted") {
      kap_c <- ar_to_pacf(phi_by_coarse[[key_c]])
      kap_m <- ar_to_pacf(phi_by_medium[[key_m]])
      kap_f <- ar_to_pacf(phi_by_fine[[key_f]])

      kap_c <- .ms_pad(kap_c, p_target)
      kap_m <- .ms_pad(kap_m, p_target)
      kap_f <- .ms_pad(kap_f, p_target)

      kap <- w["coarse"] * kap_c + w["medium"] * kap_m + w["fine"] * kap_f
      kap <- pmin(pmax(kap, -kappa_clip), kappa_clip)
      out_phi[[key_f]] <- pacf_to_ar(kap)
    } else {
      g_c <- acvf_by_coarse[[key_c]]
      g_m <- acvf_by_medium[[key_m]]
      g_f <- acvf_by_fine[[key_f]]

      g_c <- .ms_pad(g_c, p_target + 1L)
      g_m <- .ms_pad(g_m, p_target + 1L)
      g_f <- .ms_pad(g_f, p_target + 1L)

      g <- w["coarse"] * g_c + w["medium"] * g_m + w["fine"] * g_f
      if (yw_order == 0L) {
        out_phi[[key_f]] <- numeric(0)
      } else {
        yw <- yw_from_acvf_fast(g[seq_len(yw_order + 1L)], yw_order)
        out_phi[[key_f]] <- enforce_stationary_ar(yw$phi)
      }
    }
  }

  out_phi
}

.parcel_means <- function(resid, parcels, na.rm = FALSE) {
  parcel_means_fast(resid, parcels, na.rm = na.rm)
}

# Per-parcel estimates at one scale. `R` holds the valid frames (rows) of every
# voxel and `labels` the parcel of each voxel at this scale; each parcel is
# estimated from its voxels' pooled autocovariance.
#
# `lag_max` must cover the pooling target, not just the order selected for this
# scale. Sizing the acvf to the selected order alone leaves .ms_pad() zero-
# filling the remaining lags, and Yule-Walker on a zero-filled autocovariance
# returns explosive coefficients that only the stationarity clamp hides.
.ms_estimate_scale <- function(R, labels, estimator, run_starts0 = NULL, lag_max = 0L,
                               center_id = NULL, correction = NULL) {
  labels <- as.integer(labels)
  ids <- as.character(sort(unique(labels)))
  n <- nrow(R)
  starts0 <- if (is.null(run_starts0)) 0L else as.integer(run_starts0)
  if (!length(starts0)) starts0 <- 0L
  seg_id <- cumsum(seq_len(n) %in% (starts0 + 1L))

  phi_by <- setNames(vector("list", length(ids)), ids)
  acvf_by <- setNames(vector("list", length(ids)), ids)
  sigma2_by <- setNames(vector("list", length(ids)), ids)
  n_eff_by <- setNames(vector("list", length(ids)), ids)
  for (id in ids) {
    Y <- R[, labels == as.integer(id), drop = FALSE]
    fit <- estimator(Y)
    phi_by[[id]] <- fit$phi
    sigma2_by[[id]] <- if (is.null(fit$sigma2)) NA_real_ else fit$sigma2
    lag_id <- max(as.integer(lag_max), fit$order[["p"]], 1L)
    # The mean is estimated per run, and lag products stay inside segments;
    # centring each censoring fragment on its own mean destroys the
    # autocorrelation being measured.
    lag_acc <- if (is.null(correction)) lag_id else
      min(max(lag_id, nrow(correction) - 1L), n - 1L)
    pooled <- .pooled_acvf_segments(Y, seg_id, lag_acc, center_id = center_id)
    acvf_by[[id]] <- .acvf_from_pooled(pooled, order = lag_id, correction = correction)
    # Frames times effective independent voxels: the evidence behind this
    # parcel's estimate, for empirical-Bayes pooling.
    n_eff_by[[id]] <- .bic_n(Y, nrow(Y), groups = center_id)
  }
  list(phi = phi_by, acvf = acvf_by, sigma2 = sigma2_by, n_eff = n_eff_by)
}

# Empirical-Bayes pooling across scales ----------------------------------------
#
# The size-based weights (.ms_weights) depend only on how many voxels each
# scale holds, never on whether the scales actually agree. Three regions of
# 1,500 voxels with AR(1) coefficients 0.20 / 0.45 / 0.70 were each pulled
# 60% of the way to their common mean (0.35 / 0.45 / 0.55), and their
# false-positive rates drifted off 5% in both directions.
#
# Here each child (finer-scale) estimate is shrunk toward its parent's by
#   lambda = tau2 / (tau2 + s2),
# where s2 is the child's sampling variance and tau2 the real spread of
# children around their parents, estimated by moments from the data and
# pooled across parents. Children that genuinely differ (tau2 >> s2) keep
# their own estimates; small, noisy children of a homogeneous parent are
# pulled most of the way in. Applied lag by lag to partial autocorrelations
# (pacf_weighted) or autocorrelations (acvf_pooled).
.ms_eb_shrink <- function(K_child, S2_child, K_parent_of_child) {
  out <- K_child
  for (j in seq_len(nrow(K_child))) {
    dev <- K_child[j, ] - K_parent_of_child[j, ]
    s2 <- pmax(S2_child[j, ], 0)
    tau2 <- max(0, mean(dev^2) - mean(s2))
    lam <- ifelse(tau2 + s2 > 0, tau2 / (tau2 + s2), 1)
    out[j, ] <- lam * K_child[j, ] + (1 - lam) * K_parent_of_child[j, ]
  }
  out
}

# Per-lag feature matrix (lags x parcels) for one scale: partial
# autocorrelations, or autocorrelations from the pooled autocovariance.
.ms_features <- function(est, ids, target, mode) {
  vapply(ids, function(id) {
    if (identical(mode, "pacf_weighted")) {
      .ms_pad(ar_to_pacf(est$phi[[id]]), target)
    } else {
      g <- .ms_pad(est$acvf[[id]], target + 1L)
      if (!is.finite(g[1]) || g[1] <= 0) numeric(target) else g[-1L] / g[1]
    }
  }, numeric(target))
}

# Sampling variance of each feature: (1 - k^2) / N, N = frames x effective
# number of independent voxels in the parcel.
.ms_sampling_var <- function(F, n_eff) {
  N <- matrix(rep(pmax(n_eff, 1), each = nrow(F)), nrow(F))
  (1 - pmin(F^2, 0.99)) / N
}

# Turn pooled features back into AR coefficients for each fine parcel.
.ms_features_to_phi <- function(F, est_fine, ids, mode, yw_order) {
  setNames(lapply(seq_along(ids), function(i) {
    if (identical(mode, "pacf_weighted")) {
      pacf_to_ar(pmin(pmax(F[, i], -0.99), 0.99))
    } else {
      if (yw_order == 0L) return(numeric(0))
      g0 <- est_fine$acvf[[ids[i]]][1]
      g <- .shrink_to_pd(g0 * c(1, F[seq_len(yw_order), i]))
      yw <- yw_from_acvf_fast(g[seq_len(yw_order + 1L)], yw_order)
      enforce_stationary_ar(yw$phi)
    }
  }), ids)
}

# Three-scale empirical-Bayes pooling: medium toward coarse, then fine toward
# the pooled medium. `parents` maps each fine parcel to its medium and coarse
# parents (.ms_parent_maps); est_* come from .ms_estimate_scale.
.ms_combine_eb <- function(est_c, est_m, est_f, parents, target, mode, yw_order) {
  ids_f <- names(est_f$phi); ids_m <- names(est_m$phi); ids_c <- names(est_c$phi)
  Fc <- .ms_features(est_c, ids_c, target, mode)
  Fm <- .ms_features(est_m, ids_m, target, mode)
  Ff <- .ms_features(est_f, ids_f, target, mode)
  Fc <- matrix(Fc, nrow = target); Fm <- matrix(Fm, nrow = target)
  Ff <- matrix(Ff, nrow = target)
  # Each medium parcel's coarse parent, read off its fine children.
  fid <- as.integer(ids_f)
  med_of_f <- parents$parent_medium[fid]
  coarse_of_f <- parents$parent_coarse[fid]
  coarse_of_m <- vapply(as.integer(ids_m), function(m) {
    cc <- coarse_of_f[med_of_f == m]
    if (length(cc)) as.integer(names(which.max(table(cc)))) else NA_integer_
  }, 0L)
  Fm_parent <- Fc[, match(coarse_of_m, as.integer(ids_c)), drop = FALSE]
  Fm_star <- .ms_eb_shrink(Fm, .ms_sampling_var(Fm, unlist(est_m$n_eff[ids_m])),
                           Fm_parent)
  Ff_parent <- Fm_star[, match(med_of_f, as.integer(ids_m)), drop = FALSE]
  Ff_star <- .ms_eb_shrink(Ff, .ms_sampling_var(Ff, unlist(est_f$n_eff[ids_f])),
                           Ff_parent)
  .ms_features_to_phi(Ff_star, est_f, ids_f, mode, yw_order)
}

# Two-scale version, without parcel_sets: each parcel toward the estimate
# pooled over all voxels.
.ms_combine_eb_global <- function(est_g, est_f, target, mode, yw_order) {
  ids_f <- names(est_f$phi)
  Ff <- matrix(.ms_features(est_f, ids_f, target, mode), nrow = target)
  Fg <- matrix(.ms_features(est_g, names(est_g$phi), target, mode), nrow = target)
  Ff_parent <- Fg[, rep(1L, length(ids_f)), drop = FALSE]
  Ff_star <- .ms_eb_shrink(Ff, .ms_sampling_var(Ff, unlist(est_f$n_eff[ids_f])),
                           Ff_parent)
  .ms_features_to_phi(Ff_star, est_f, ids_f, mode, yw_order)
}
