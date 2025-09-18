# compat.R -- small interface layer for fmrireg/fmriAR integration

`%||%` <- function(x, y) if (!is.null(x)) x else y

compat_env <- local({

  plan_from_phi <- function(phi, theta = NULL,
                            runs = NULL, parcels = NULL,
                            pooling = c("global","run","parcel"),
                            exact_first = TRUE,
                            method = c("ar","arma")) {
    pooling <- match.arg(pooling)
    method <- match.arg(method)

    if (pooling %in% c("global","run")) {
      phi_list   <- if (is.list(phi)) phi else list(phi)
      theta_list <- if (is.null(theta)) list() else if (is.list(theta)) theta else list(theta)
      order_vec  <- c(p = max(vapply(phi_list, length, 0L)),
                      q = max(vapply(theta_list, length, 0L)))
      return(new_whiten_plan(phi = phi_list,
                             theta = theta_list,
                             order = order_vec,
                             runs = runs,
                             exact_first = isTRUE(exact_first),
                             method = method,
                             pooling = pooling))
    }

    stopifnot(!is.null(parcels))
    parcels <- as.integer(parcels)

    if (pooling == "parcel") {
      stopifnot(is.list(phi))
      theta_list <- if (is.null(theta))
        setNames(vector("list", length(phi)), names(phi))
      else theta
      order_vec <- c(p = max(vapply(phi, length, 0L)),
                     q = max(vapply(theta_list, length, 0L)))
      return(new_whiten_plan(phi = NULL,
                             theta = NULL,
                             order = order_vec,
                             runs = runs,
                             exact_first = isTRUE(exact_first),
                             method = method,
                             pooling = "parcel",
                             parcels = parcels,
                             parcel_ids = sort(unique(parcels)),
                             phi_by_parcel = phi,
                             theta_by_parcel = theta_list))
    }

    stop("Unsupported pooling value")
  }

  whiten_with_phi <- function(X, Y, phi, theta = NULL,
                              runs = NULL, parcels = NULL,
                              pooling = c("global","run","parcel"),
                              exact_first = FALSE,
                              parallel = TRUE) {
    plan <- plan_from_phi(phi = phi,
                          theta = theta,
                          runs = runs,
                          parcels = parcels,
                          pooling = match.arg(pooling),
                          exact_first = exact_first,
                          method = if (is.null(theta) || (is.list(theta) && all(vapply(theta, length, 0L) == 0L))) "ar" else "arma")
    whiten_apply(plan, X, Y, runs = runs, parcels = parcels, parallel = parallel)
  }

  update_plan <- function(prev_plan, resid, p = NULL, q = NULL, p_max = 6L) {
    stopifnot(inherits(prev_plan, "fmriAR_plan"))
    method <- prev_plan$method
    pooling <- prev_plan$pooling
    runs <- prev_plan$runs
    parcels <- prev_plan$parcels
    p_use <- if (!is.null(p)) p else prev_plan$order[["p"]]
    q_use <- if (!is.null(q)) q else prev_plan$order[["q"]]

    fit_noise(resid = resid,
              runs = runs,
              parcels = parcels,
              method = method,
              p = if (!is.null(p_use) && p_use > 0L) as.integer(p_use) else "auto",
              q = if (!is.null(q_use)) as.integer(q_use) else 0L,
              p_max = as.integer(p_max),
              pooling = pooling)
  }

  plan_info <- function(plan) {
    stopifnot(inherits(plan, "fmriAR_plan"))
    phi_mat <- NULL
    if (!is.null(plan$phi_by_parcel)) {
      p_ord <- plan$order[["p"]]
      phi_list <- plan$phi_by_parcel
      if (p_ord > 0L) {
        phi_mat <- vapply(phi_list, function(phi) .ms_pad(phi, p_ord), numeric(p_ord))
      } else {
        phi_mat <- matrix(0, nrow = 0L, ncol = length(phi_list))
      }
    }
    list(method = plan$method,
         pooling = plan$pooling,
         order = plan$order,
         has_runs = !is.null(plan$runs),
         has_parcels = !is.null(plan$parcels),
         n_parcels = if (!is.null(plan$parcels)) length(unique(plan$parcels)) else 0L,
         phi_by_parcel = phi_mat)
  }

  whiteness_score <- function(resid, K = 3L, aggregate = c("mean","median")) {
    aggregate <- match.arg(aggregate)
    if (!is.matrix(resid)) resid <- as.matrix(resid)
    scores <- apply(resid, 2L, function(y) {
      ac <- stats::acf(y, plot = FALSE, lag.max = K, demean = TRUE)$acf[-1L]
      mean(abs(ac))
    })
    if (aggregate == "mean") mean(scores) else stats::median(scores)
  }

  list(plan_from_phi = plan_from_phi,
       whiten_with_phi = whiten_with_phi,
       update_plan = update_plan,
       plan_info = plan_info,
       whiteness_score = whiteness_score)
})

#' fmrireg compatibility interface
#'
#' Stable entry points to help upstream packages reuse fmriAR whitening without
#' rewriting existing pipelines.
#'
#' @return A list environment containing compatibility functions:
#'   \itemize{
#'     \item \code{plan_from_phi}: Create whitening plan from AR coefficients
#'     \item \code{whiten_with_phi}: Apply whitening given AR coefficients
#'     \item \code{update_plan}: Update existing plan with new residuals
#'     \item \code{plan_info}: Extract information from a plan object
#'     \item \code{whiteness_score}: Compute whiteness metric from residuals
#'   }
#'
#' @examples
#' # Create compatibility interface
#' compat_funcs <- compat
#'
#' # Example: Create whitening plan from AR coefficients
#' phi <- c(0.3, 0.1)  # AR(2) coefficients
#' plan <- compat_funcs$plan_from_phi(phi, exact_first = TRUE)
#'
#' # Example: Compute whiteness score
#' resid <- matrix(rnorm(100 * 10), 100, 10)
#' score <- compat_funcs$whiteness_score(resid)
#'
#' @keywords internal
#' @export
compat <- compat_env
