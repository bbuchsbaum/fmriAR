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
#
# `p` and `q` may each be a vector of candidate orders. The order pair is then
# chosen by BIC on the Hannan-Rissanen regression (Hannan & Rissanen, 1982):
# every candidate is a column subset of one regression on the largest lag set,
# fitted on a common set of rows, so selection costs one extra pass over the
# data however many candidates there are. The BIC sample size is frames times
# the effective number of independent voxels (see .effective_voxels()).
.hr_arma_pooled <- function(units, p, q, iter = 0L, p_big = NULL,
                            ar_bound = 0.99) {
  p_grid <- sort(unique(as.integer(p))); q_grid <- sort(unique(as.integer(q)))
  stopifnot(length(p_grid) >= 1L, length(q_grid) >= 1L,
            all(p_grid >= 0L), all(q_grid >= 0L))
  p <- max(p_grid); q <- max(q_grid); iter <- as.integer(iter)
  selecting <- length(p_grid) > 1L || length(q_grid) > 1L
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
  selection <- NULL
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

  if (selecting) {
    # One regression on the full lag set over common rows, then BIC for every
    # (p, q) subset from its normal equations: RSS = y'y - b_s' coef_s.
    start_all <- max(p, pb + q)
    G <- matrix(0, k, k); b <- numeric(k); yy <- 0; rows <- 0
    for (i in seq_along(units)) {
      ne <- hr_normal_eq_cpp(units[[i]]$mat, E_list[[i]], units[[i]]$rel,
                             p, q, start_all)
      G <- G + ne$G; b <- b + ne$b; yy <- yy + ne$yy; rows <- rows + ne$rows
    }
    V <- ncol(units[[1]]$mat)
    if (rows < 2) return(failure("too few usable timepoints"))
    # Voxels are replicate series, but not independent ones: count them at
    # their effective number. Counting frames only (one voxel's worth of
    # evidence) under-fitted shared slow noise -- leaving residual
    # autocorrelation -- while counting every voxel as independent
    # over-fitted when they share fluctuations.
    n_bic <- rows * .effective_voxels(units)
    cand <- expand.grid(p = p_grid, q = q_grid)
    cand$bic <- vapply(seq_len(nrow(cand)), function(i) {
      cols <- c(seq_len(cand$p[i]), p + seq_len(cand$q[i]))
      rss <- yy
      if (length(cols)) {
        cf <- tryCatch(solve(G[cols, cols, drop = FALSE], b[cols]),
                       error = function(e) NULL)
        if (is.null(cf) || !all(is.finite(cf))) return(Inf)
        rss <- yy - sum(b[cols] * cf)
      }
      s2 <- rss / (rows * V)
      if (!is.finite(s2) || s2 <= 0) return(Inf)
      n_bic * log(s2) + (cand$p[i] + cand$q[i]) * log(n_bic)
    }, 0)
    best <- which.min(cand$bic)
    selection <- cand
    p <- cand$p[best]; q <- cand$q[best]; k <- p + q
    if (k == 0L) {
      s2 <- sum(vapply(units, function(u) sum(u$mat^2), 0)) /
        sum(vapply(units, function(u) length(u$mat), 0))
      return(list(phi = numeric(0), theta = numeric(0), sigma2 = s2,
                  order = c(p = 0L, q = 0L), p_big = as.integer(pb),
                  iterations = iter, method = "hr_pooled", ok = TRUE,
                  selection = selection))
    }
  }

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
       method = "hr_pooled", ok = TRUE, selection = selection)
}

# Effective number of independent voxels (Kish design effect). Voxels that
# share fluctuations carry less independent evidence about the noise model
# than their count suggests: with mean inter-voxel correlation rbar,
# V_eff = V / (1 + (V - 1) rbar). rbar follows from the variance of the voxel
# mean, Var(mean) = mean(Var) * (1 + (V - 1) rbar) / V, so no V x V correlation
# matrix is formed. Bounded to [1, V].
.effective_voxels <- function(units) {
  V <- ncol(units[[1]]$mat)
  if (V <= 1L) return(1)
  vbar <- 0; vmean <- 0; nn <- 0
  for (u in units) {
    m <- u$mat
    vbar <- vbar + sum(colSums(m * m)) / V
    rm <- rowMeans(m)
    vmean <- vmean + sum(rm * rm)
    nn <- nn + nrow(m)
  }
  if (!(vbar > 0)) return(1)
  ratio <- V * vmean / vbar           # = 1 + (V - 1) * rbar
  rbar <- min(max((ratio - 1) / (V - 1), 0), 1)
  max(1, min(V, V / (1 + (V - 1) * rbar)))
}
