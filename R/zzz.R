.onLoad <- function(libname, pkgname) {
  op <- options()
  op.fmriAR <- list(
    fmriAR.max_threads = NA_integer_,
    fmriAR.debug = FALSE,
    fmriAR.use_cpp_hr = TRUE
  )
  toset <- !(names(op.fmriAR) %in% names(op))
  if (any(toset)) options(op.fmriAR[toset])
  invisible()
}
