#' Run- and censor-aware noise autocovariance
#'
#' Estimates the autocovariance of fMRI noise directly, without fitting an
#' AR model first. Lag products never cross a run boundary or a censoring gap,
#' and the mean is removed per run rather than per fragment, so scrubbed data
#' does not have its autocorrelation destroyed by the centering.
#'
#' This is the estimator [fit_noise()] uses internally, exposed so consumers can
#' get the covariance itself. For a run-stationary process with autocovariance
#' `gamma`, the covariance of any linear contrast `L` factors as
#' `sum_h gamma_h * L T_h L'` with `T_h` the lag-`h` indicator, so `gamma` plus
#' a design is enough to compute exact contrast variances without refitting.
#'
#' Unlike [acorr_diagnostics()], which returns normalized autocorrelation for
#' inspection, this returns covariances on the scale of the data.
#'
#' @param resid Numeric matrix of residuals (timepoints x voxels).
#' @param runs Optional run labels, length `nrow(resid)`. Each label must occupy
#'   one contiguous block and may not be missing.
#' @param censor Optional indices of censored timepoints (1-based).
#' @param max_lag Highest lag to return.
#' @param pooling `"global"` for one autocovariance over all runs, `"run"` for
#'   one per run, `"parcel"` for one per parcel.
#' @param parcels Parcel labels (length `ncol(resid)`), required when
#'   `pooling = "parcel"`.
#' @param design Optional design matrix whose projection produced `resid`. When
#'   supplied the residual bias is corrected; see [acvf_bias_matrix()].
#' @param correction_max_lag Lag budget used to estimate and undo residual
#'   projection bias when `design` is supplied. This is independent of
#'   `max_lag`: correction needs enough tail information even when only a few
#'   output lags are requested.
#' @return An object of class `fmriAR_acvf`: a list with
#'   \itemize{
#'     \item `acvf`: named list of autocovariance vectors, lags 0 upward.
#'     \item `pairs`: matching list of the pair count behind each lag, so a
#'       consumer can tell an estimate backed by thousands of pairs from one
#'       backed by three.
#'     \item `n_segments`, `segment_lengths`: the segmentation actually used.
#'     \item `max_lag`, `pooling`, `corrected`.
#'   }
#'   Lags with no surviving pairs are dropped rather than reported as zero.
#' @seealso [fit_noise()], [acvf_bias_matrix()]
#' @examples
#' set.seed(1)
#' r <- matrix(rnorm(400 * 5), 400, 5)
#' a <- noise_acvf(r, runs = rep(1:2, each = 200), max_lag = 5, pooling = "run")
#' a$acvf[["1"]]
#' @export
noise_acvf <- function(resid, runs = NULL, censor = NULL, max_lag = 20L,
                       pooling = c("global", "run", "parcel"), parcels = NULL,
                       design = NULL, correction_max_lag = 25L) {
  pooling <- match.arg(pooling)
  if (!is.matrix(resid)) resid <- as.matrix(resid)
  storage.mode(resid) <- "double"
  n <- nrow(resid)
  if (n < 2L) stop("'resid' needs at least two timepoints")
  if (length(max_lag) != 1L || !is.finite(max_lag) || max_lag < 0) {
    stop("'max_lag' must be one non-negative finite number")
  }
  if (any(!is.finite(resid))) stop("'resid' contains NA, NaN, or Inf")
  max_lag <- as.integer(max_lag)

  if (!is.null(censor)) {
    if (is.logical(censor)) {
      stopifnot(length(censor) == n)
      censor <- which(censor)
    }
    censor <- sort(unique(as.integer(censor)))
    censor <- censor[censor >= 1L & censor <= n]
    if (!length(censor)) censor <- NULL
  }

  corr_by_run <- if (is.null(design)) {
    NULL
  } else {
    if (length(correction_max_lag) != 1L || !is.finite(correction_max_lag) ||
        correction_max_lag < 1) {
      stop("'correction_max_lag' must be one positive finite number")
    }
    .drop_unusable_corrections(
      .acvf_bias_by_run(design, n, runs = runs, censor = censor,
                        max_lag = as.integer(correction_max_lag)))
  }

  # Units to estimate over. "run" and "global" split by run; "global" then
  # pools the per-run results the way fit_noise does.
  Rsets <- .run_sets(runs, n)
  if (is.null(runs)) names(Rsets) <- "1"

  if (pooling == "parcel") {
    if (is.null(parcels)) stop("'parcels' is required when pooling = 'parcel'")
    parcels <- .parcel_codes(parcels)
    stopifnot(length(parcels) == ncol(resid))
    if (!is.null(design)) {
      stop("residual-bias correction is not yet supported for pooling = 'parcel'")
    }
  }

  one_unit <- function(idx, cols, corr) {
    keep <- if (is.null(censor)) idx else setdiff(idx, censor)
    if (length(keep) < 2L) return(NULL)
    seg_id <- cumsum(c(1L, as.integer(diff(keep) > 1L)))
    lag_budget <- if (is.null(corr)) max_lag else
      max(max_lag, nrow(corr) - 1L)
    pooled <- .pooled_acvf_segments(resid[keep, cols, drop = FALSE], seg_id,
                                    lag_budget)
    result <- .acvf_from_pooled(pooled, order = max_lag, correction = corr,
                                return_status = TRUE)
    g <- result$acvf
    if (!length(g)) return(NULL)
    list(acvf = g, pairs = pooled$pairs[seq_along(g)],
         n_seg = max(seg_id), seg_len = as.integer(table(seg_id)),
         corrected = result$corrected)
  }

  all_cols <- seq_len(ncol(resid))
  if (pooling == "parcel") {
    ids <- sort(unique(parcels))
    # One segmentation for every parcel: segments break at run boundaries and
    # censoring gaps alike, and the mean comes off per run rather than per
    # fragment. Centering a two-frame fragment on its own mean forces its
    # lag-1 product negative by construction.
    seg <- .valid_segments(n, runs = runs, censor = censor)
    if (length(seg$idx) < 2L) stop("no valid segments remain after censoring")
    seg_id <- cumsum(seq_along(seg$idx) %in% (seg$starts0 + 1L))
    units <- setNames(lapply(ids, function(pid) {
      cols <- which(parcels == pid)
      if (!length(cols)) return(NULL)
      pooled <- .pooled_acvf_segments(resid[seg$idx, cols, drop = FALSE], seg_id,
                                      max_lag, center_id = seg$run_id)
      g <- .acvf_from_pooled(pooled, order = max_lag)
      if (!length(g)) return(NULL)
      list(acvf = g, pairs = pooled$pairs[seq_along(g)],
           n_seg = max(seg_id), seg_len = as.integer(table(seg_id)),
           corrected = FALSE)
    }), as.character(ids))
  } else {
    units <- setNames(lapply(seq_along(Rsets), function(i) {
      one_unit(Rsets[[i]], all_cols, if (is.null(corr_by_run)) NULL else corr_by_run[[i]])
    }), names(Rsets))
  }

  units <- units[!vapply(units, is.null, logical(1))]
  if (!length(units)) stop("no valid segments remain after censoring")

  if (pooling == "global" && length(units) > 1L) {
    # Truncate to the shortest before averaging: a zero-padded autocovariance is
    # not a valid covariance, and averaging equal-length PD Toeplitz matrices is
    # PD by convexity while averaging padded ones is not.
    L <- min(vapply(units, function(u) length(u$acvf), 0L))
    w <- vapply(units, function(u) sum(u$seg_len), 0)
    w <- w / sum(w)
    G <- do.call(rbind, lapply(units, function(u) u$acvf[seq_len(L)]))
    P <- do.call(rbind, lapply(units, function(u) u$pairs[seq_len(L)]))
    all_corrected <- all(vapply(units, function(u) isTRUE(u$corrected), logical(1)))
    units <- list(`1` = list(acvf = as.numeric(drop(crossprod(w, G))),
                             pairs = colSums(P),
                             n_seg = sum(vapply(units, function(u) u$n_seg, 0L)),
                             seg_len = unlist(lapply(units, function(u) u$seg_len),
                                              use.names = FALSE),
                             corrected = all_corrected))
  }

  structure(list(
    acvf = lapply(units, `[[`, "acvf"),
    pairs = lapply(units, `[[`, "pairs"),
    n_segments = vapply(units, function(u) as.integer(u$n_seg), 0L),
    segment_lengths = lapply(units, `[[`, "seg_len"),
    max_lag = max_lag,
    pooling = pooling,
    corrected = all(vapply(units, function(u) isTRUE(u$corrected), logical(1)))
  ), class = "fmriAR_acvf")
}

#' @export
print.fmriAR_acvf <- function(x, ...) {
  cat("fmriAR noise autocovariance\n")
  cat("  pooling:   ", x$pooling, if (x$corrected) " (bias-corrected)" else "", "\n", sep = "")
  cat("  units:     ", length(x$acvf), " (", paste(utils::head(names(x$acvf), 6L),
      collapse = ", "), if (length(x$acvf) > 6L) ", ..." else "", ")\n", sep = "")
  cat("  max lag:   ", x$max_lag, "\n", sep = "")
  cat("  segments:  ", paste(x$n_segments, collapse = ", "), "\n", sep = "")
  g <- x$acvf[[1L]]
  cat("  unit 1:    gamma_0 = ", format(g[1], digits = 4),
      ", lags returned = ", length(g) - 1L, "\n", sep = "")
  if (length(g) > 1L) {
    cat("             rho = ", paste(format(g[-1] / g[1], digits = 3)[
      seq_len(min(6L, length(g) - 1L))], collapse = " "),
      if (length(g) > 7L) " ..." else "", "\n", sep = "")
  }
  invisible(x)
}
