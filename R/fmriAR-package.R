#' @keywords internal
#' @useDynLib fmriAR, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom stats setNames
"_PACKAGE"

#' @title fmriAR: Fast AR/ARMA prewhitening for fMRI
#' @description
#' Estimate AR/ARMA noise models from residuals and apply matched GLS
#' prewhitening to fMRI design and data matrices. Run-aware and censor-aware.
#'
#' @details
#' The fmriAR package provides efficient implementations for:
#' \itemize{
#'   \item AR and ARMA model estimation from fMRI residuals
#'   \item Run-aware and censor-aware whitening transformations
#'   \item Parcel-based parameter pooling
#'   \item Sandwich standard error computation
#' }
#'
#' @seealso
#' Useful links:
#' \itemize{
#'   \item \code{\link{fit_noise}} for noise model estimation
#'   \item \code{\link{whiten_apply}} for applying whitening transformations
#'   \item \code{\link{whiten}} for one-step whitening
#' }
#'
#' @name fmriAR-package
#' @aliases fmriAR
NULL
