#' Pretty-print an fmriAR whitening plan
#'
#' @param x An object returned by [fit_noise()].
#' @param ... Unused; included for S3 compatibility.
#' @export
print.fmriAR_plan <- function(x, ...) {
  cat("fmriAR whitening plan\n")

  method <- if (!is.null(x$method)) x$method else "ar"
  order <- x$order
  pooling <- if (!is.null(x$pooling)) x$pooling else "global"

  cat("  Method: ", toupper(method), "\n", sep = "")
  if (!is.null(order) && length(order)) {
    p <- order[["p"]]
    q <- order[["q"]]
    cat("  Orders: p = ", if (length(p)) p else 0L,
        ", q = ", if (length(q)) q else 0L, "\n", sep = "")
  }
  cat("  Pooling: ", pooling, "\n", sep = "")

  runs <- x$runs
  if (!is.null(runs)) {
    run_levels <- unique(runs)
    shown <- head(run_levels, 5L)
    cat("  Runs: ", length(run_levels),
        if (length(run_levels)) paste0(" (", paste(as.character(shown), collapse = ", "),
                                       if (length(run_levels) > length(shown)) ", …" else "",
                                       ")") else "",
        "\n", sep = "")
  }

  exact_first <- if (isTRUE(x$exact_first)) "AR(1)" else "none"
  cat("  Exact first-sample scaling: ", exact_first, "\n", sep = "")

  fmt_vec <- function(vec) {
    if (is.null(vec) || length(vec) == 0L) return("(none)")
    paste(format(signif(vec, 3), trim = TRUE), collapse = ", ")
  }

  align_lists <- function(lhs, rhs) {
    len_lhs <- length(lhs)
    len_rhs <- length(rhs)
    if (len_lhs < len_rhs) {
      lhs <- c(lhs, vector("list", len_rhs - len_lhs))
    } else if (len_rhs < len_lhs) {
      rhs <- c(rhs, vector("list", len_lhs - len_rhs))
    }
    list(lhs, rhs)
  }

  cat("  Coefficients:\n")
  if (identical(pooling, "parcel")) {
    phi_by <- x$phi_by_parcel
    theta_by <- x$theta_by_parcel
    parcel_ids <- x$parcel_ids
    if (is.null(parcel_ids) && !is.null(phi_by)) parcel_ids <- names(phi_by)
    if (is.null(parcel_ids) && !is.null(x$parcels)) parcel_ids <- sort(unique(x$parcels))
    n_parcels <- length(parcel_ids)
    cat("    Parcel-level coefficients for ", n_parcels, " parcels\n", sep = "")
    if (n_parcels > 0L && length(phi_by)) {
      show_idx <- seq_len(min(5L, n_parcels))
      cat("    Showing first ", length(show_idx), ":\n", sep = "")
      phi_names <- names(phi_by)
      theta_names <- names(theta_by)
      parcel_names <- names(parcel_ids)
      for (idx in show_idx) {
        pid <- parcel_ids[[idx]]
        key <- NULL
        if (!is.null(phi_names) && length(phi_names) >= idx && nzchar(phi_names[[idx]])) key <- phi_names[[idx]]
        if (is.null(key) && !is.null(parcel_names) && length(parcel_names) >= idx && nzchar(parcel_names[[idx]])) key <- parcel_names[[idx]]
        if (is.null(key)) key <- as.character(pid)
        phi_vec <- NULL
        if (!is.null(phi_by) && length(phi_by)) {
          if (!is.null(phi_by[[key]])) phi_vec <- phi_by[[key]] else if (idx <= length(phi_by)) phi_vec <- phi_by[[idx]]
        }
        if (is.null(phi_vec)) phi_vec <- numeric(0L)
        theta_vec <- NULL
        if (!is.null(theta_by) && length(theta_by)) {
          if (!is.null(theta_by[[key]])) theta_vec <- theta_by[[key]] else if (idx <= length(theta_by)) theta_vec <- theta_by[[idx]]
        }
        if (is.null(theta_vec)) theta_vec <- numeric(0L)
        line <- paste0("      Parcel ", pid, ": phi = ", fmt_vec(phi_vec))
        if (length(theta_vec)) {
          line <- paste0(line, "; theta = ", fmt_vec(theta_vec))
        }
        cat(line, "\n", sep = "")
      }
      if (n_parcels > length(show_idx)) cat("      …\n")
    }
  } else {
    phi_list <- x$phi
    theta_list <- x$theta
    if (is.null(phi_list)) phi_list <- vector("list", length(theta_list))
    if (is.null(theta_list)) theta_list <- vector("list", length(phi_list))
    tmp <- align_lists(phi_list, theta_list)
    phi_list <- tmp[[1]]
    theta_list <- tmp[[2]]

    labels <- names(phi_list)
    if (is.null(labels) || all(labels == "")) {
      if (identical(pooling, "global")) {
        labels <- rep("global", length(phi_list))
      } else {
        labels <- paste0("run", seq_along(phi_list))
      }
    }

    show_idx <- seq_len(min(length(phi_list), 8L))
    for (i in show_idx) {
      phi_vec <- phi_list[[i]]
      theta_vec <- theta_list[[i]]
      line <- paste0("    ", labels[[i]], ": phi = ", fmt_vec(phi_vec))
      if (!is.null(theta_vec) && length(theta_vec)) {
        line <- paste0(line, "; theta = ", fmt_vec(theta_vec))
      }
      cat(line, "\n", sep = "")
    }
    if (length(phi_list) > length(show_idx)) cat("    …\n")
  }

  invisible(x)
}

