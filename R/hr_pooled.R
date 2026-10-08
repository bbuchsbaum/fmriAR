# Pooled, segment-aware Hannan-Rissanen ARMA estimation -----------------------
#
# fit_noise() used to estimate ARMA by running Hannan-Rissanen on the
# voxel-MEAN time series. The mean's autocovariance is the average of all voxel
# cross-covariances, so with many voxels it is dominated by whatever is shared
# across them (global signal, physiology) rather than by the voxel-level noise
# being whitened: near-white voxels plus a weak shared slow component returned
# phi = 0.86 against a voxel lag-1 autocorrelation of 0.25. The AR path pools
# per-voxel autocovariances instead, and this does the same for ARMA: every
# voxel is a replicate series sharing (phi, theta), and both regressions of the
# Hannan-Rissanen procedure are solved on cross-products summed over voxels.
#
# Two further defects of the single-series fit are handled here:
#   * Lag products never span a censoring gap or run boundary. Each unit's
#     valid frames are split into contiguous segments and all lags are taken
#     within a segment, so censored data no longer biases theta.
#   * The long-AR residuals are discarded for the first p_big frames of each
#     segment, where the zero pre-sample makes them poor innovation proxies.

# units: list of list(mat = valid-frame residual matrix (rows x voxels),
#                     starts0 = 0-based starts of the contiguous segments)
# Columns are centred per unit here (one mean per run, as on the AR path).
.hr_arma_pooled <- function(units, p, q, iter = 0L, p_big = NULL,
                            ar_bound = 0.99) {
  p <- as.integer(p); q <- as.integer(q); iter <- as.integer(iter)
  stopifnot(p >= 0L, q >= 0L)
  units <- lapply(units, function(u) {
    mat <- u$mat
    if (!is.matrix(mat)) mat <- as.matrix(mat)
    storage.mode(mat) <- "double"
    if (nrow(mat)) mat <- mat - rep(colMeans(mat), each = nrow(mat))
    starts0 <- as.integer(u$starts0)
    if (!length(starts0)) starts0 <- 0L
    seg_len <- diff(c(starts0, nrow(mat)))
    rel <- sequence(seg_len) - 1L
    list(mat = mat, starts0 = starts0, seg_len = seg_len, rel = rel,
         seg_id = rep(seq_along(seg_len), seg_len))
  })
  units <- units[vapply(units, function(u) nrow(u$mat) >= 2L && ncol(u$mat) >= 1L, TRUE)]
  failure <- function(reason) {
    warning("fit_noise: ARMA estimation failed (", reason, "); ",
            "returning a white-noise model", call. = FALSE)
    list(phi = numeric(p), theta = numeric(q), sigma2 = NA_real_,
         order = c(p = p, q = q), p_big = NA_integer_, iterations = iter,
         method = "hr_pooled", ok = FALSE)
  }
  if (!length(units)) return(failure("no valid data"))

  seg_lens <- unlist(lapply(units, `[[`, "seg_len"))
  k <- p + q
  if (k == 0L) {
    s2 <- sum(vapply(units, function(u) sum(u$mat^2), 0)) /
      sum(vapply(units, function(u) length(u$mat), 0))
    return(list(phi = numeric(0), theta = numeric(0), sigma2 = s2,
                order = c(p = 0L, q = 0L), p_big = 0L, iterations = iter,
                method = "hr_pooled", ok = TRUE))
  }

  usable_rows <- function(pb) sum(pmax(0L, seg_lens - max(p, pb + q)))
  if (is.null(p_big)) {
    p_big <- min(max(8L, k + 5L, ceiling(10 * log10(max(seg_lens)))), 40L)
  }
  p_big <- as.integer(min(p_big, max(seg_lens) - 1L))
  # Shrink the long-AR order until enough rows survive the burn-in.
  while (p_big > max(1L, k) && usable_rows(p_big) < 10L * (k + 1L)) p_big <- p_big - 1L
  if (usable_rows(p_big) <= k) return(failure("too few usable timepoints"))

  # Step 1: long AR from the pooled, segment-respecting autocovariance.
  num <- numeric(p_big + 1L); pairs <- numeric(p_big + 1L)
  for (u in units) {
    # Columns are already centred per unit, so go straight to the lag sums.
    pl <- pooled_acvf_seg_cpp(u$mat, u$seg_id, p_big)
    num <- num + pl$num * ncol(u$mat)
    pairs <- pairs + pl$pairs * ncol(u$mat)
  }
  g <- .acvf_from_pooled(list(num = num, pairs = pairs), order = p_big)
  if (length(g) < 2L) return(failure("degenerate autocovariance"))
  pb <- min(p_big, length(g) - 1L)
  phi_big <- yw_from_acvf_fast(g[seq_len(pb + 1L)], pb)$phi

  innovations <- function(u, ph, th) {
    arma_whiten_inplace(u$mat + 0, matrix(0, nrow(u$mat), 0L), ph, th,
                        run_starts = u$starts0, exact_first = FALSE,
                        parallel = TRUE, n_threads = .n_threads())$Y
  }

  regress <- function(E_list) {
    G <- matrix(0, k, k); b <- numeric(k)
    for (i in seq_along(units)) {
      ne <- hr_normal_eq_cpp(units[[i]]$mat, E_list[[i]], units[[i]]$rel,
                             p, q, max(p, pb + q))
      G <- G + ne$G
      b <- b + ne$b
    }
    tryCatch(solve(G, b), error = function(e) NULL)
  }

  E_list <- lapply(units, innovations, ph = phi_big, th = numeric(0))
  phi <- numeric(p); theta <- numeric(q)
  for (it in 0:iter) {
    coef <- regress(E_list)
    if (is.null(coef) || !all(is.finite(coef))) return(failure("singular regression"))
    phi <- coef[seq_len(p)]
    theta <- coef[p + seq_len(q)]
    if (p) {
      phi <- enforce_stationary_ar(phi, ar_bound)
      if (length(phi) != p) phi <- numeric(p)
    }
    if (q) theta <- enforce_invertible_ma(theta)
    if (it < iter) E_list <- lapply(units, innovations, ph = phi, th = theta)
  }

  # Innovation variance at voxel scale, from the exact (stationary-start)
  # whitening of every valid frame.
  ss <- 0; cnt <- 0
  for (u in units) {
    W <- arma_whiten_inplace(u$mat + 0, matrix(0, nrow(u$mat), 0L), phi, theta,
                             run_starts = u$starts0, exact_first = TRUE,
                             parallel = TRUE, n_threads = .n_threads())$Y
    ss <- ss + sum(W^2); cnt <- cnt + length(W)
  }

  list(phi = as.numeric(phi), theta = as.numeric(theta), sigma2 = ss / cnt,
       order = c(p = p, q = q), p_big = as.integer(pb), iterations = iter,
       method = "hr_pooled", ok = TRUE)
}
