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
  tapply(vvar, parcels, function(z) stats::mad(z, constant = 1))
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
                                kappa_clip = 0.99) {
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
      yw <- yw_from_acvf_fast(g, p_target)
      out_phi[[key_f]] <- enforce_stationary_ar(yw$phi)
    }
  }

  out_phi
}

.parcel_means <- function(resid, parcels, na.rm = FALSE) {
  parcel_means_fast(resid, parcels, na.rm = na.rm)
}

.ms_estimate_scale <- function(M, estimator, run_starts0 = NULL) {
  ids <- colnames(M)
  phi_by <- setNames(vector("list", length(ids)), ids)
  acvf_by <- setNames(vector("list", length(ids)), ids)
  for (id in ids) {
    fit <- estimator(M[, id])
    phi_by[[id]] <- fit$phi
    lag_max <- max(0L, fit$order[["p"]] + 1L)
    acvf_by[[id]] <- .segment_acvf(M[, id], run_starts0, lag_max)
  }
  list(phi = phi_by, acvf = acvf_by)
}
