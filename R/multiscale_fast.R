#' @keywords internal
parcel_means_fast <- function(resid, parcels, na.rm = FALSE) {
  stopifnot(is.matrix(resid), length(parcels) == ncol(resid))
  parcels <- as.integer(parcels)
  ids <- sort(unique(parcels))
  idx <- match(parcels, ids)
  K <- length(ids)

  use_cpp <- isTRUE(getOption("fmriAR.use_cpp_means", TRUE)) && exists("parcel_means_cpp")
  if (use_cpp) {
    out <- parcel_means_cpp(resid, idx, K = K, na_rm = na.rm)
  } else {
    n <- nrow(resid)
    sums <- matrix(0, n, K)
    if (!na.rm) {
      counts <- tabulate(idx, nbins = K)
      for (j in seq_len(ncol(resid))) {
        k <- idx[j]
        sums[, k] <- sums[, k] + resid[, j]
      }
      out <- sweep(sums, 2L, pmax(counts, 1L), `/`, check.margin = FALSE)
    } else {
      counts <- matrix(0L, nrow = n, ncol = K)
      for (j in seq_len(ncol(resid))) {
        k <- idx[j]
        x <- resid[, j]
        ok <- !is.na(x)
        sums[ok, k] <- sums[ok, k] + x[ok]
        counts[ok, k] <- counts[ok, k] + 1L
      }
      counts[counts == 0L] <- 1L
      out <- sums / counts
    }
  }

  colnames(out) <- as.character(ids)
  out
}

#' @keywords internal
segmented_acvf_fast <- function(y, run_starts0, max_lag, unbiased = FALSE, center = TRUE) {
  stopifnot(is.numeric(y))
  run_starts0 <- as.integer(run_starts0)
  if (length(run_starts0) == 0L || run_starts0[1] != 0L)
    stop("run_starts must be 0-based and start at 0")
  segmented_acvf_cpp(y, run_starts0, as.integer(max_lag), unbiased = unbiased, center = center)
}

#' @keywords internal
yw_from_acvf_fast <- function(gamma, p) {
  stopifnot(is.numeric(gamma), length(gamma) >= p + 1L)
  res <- yw_from_acvf_cpp(gamma, as.integer(p))
  list(phi = as.numeric(res$phi), sigma2 = as.numeric(res$sigma2))
}
