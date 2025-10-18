# Helpers specific to AFNI-style validation
poly_from_roots <- function(roots) {
  coefs <- 1 + 0i
  for (r in roots) {
    coefs <- c(coefs, 0) - r * c(0, coefs)
  }
  Re(coefs)
}

sim_ar_process <- function(phi, n, burn = max(200L, length(phi) * 20L), sd = 1) {
  total <- n + burn
  series <- as.numeric(stats::arima.sim(list(ar = phi), n = total,
                                        innov = rnorm(total, sd = sd)))
  tail(series, n)
}

test_that("AFNI AR(3) root mapping matches polynomial expansion", {
  a <- 0.6
  r1 <- 0.7
  t1 <- pi / 7
  phi <- fmriAR:::.afni_phi_ar3(a, r1, t1)
  roots <- c(a, r1 * exp(1i * t1), r1 * exp(-1i * t1))
  coefs <- poly_from_roots(roots)
  expect_equal(unname(coefs[-1]), unname(-phi), tolerance = 1e-10)
})

test_that("AFNI AR(5) root mapping matches polynomial expansion", {
  a <- 0.55
  r1 <- 0.6; t1 <- pi / 6
  r2 <- 0.45; t2 <- pi / 3
  phi <- fmriAR:::.afni_phi_ar5(a, r1, t1, r2, t2)
  roots <- c(a,
             r1 * exp(1i * t1), r1 * exp(-1i * t1),
             r2 * exp(1i * t2), r2 * exp(-1i * t2))
  coefs <- poly_from_roots(roots)
  expect_equal(unname(coefs[-1]), unname(-phi), tolerance = 1e-9)
})

test_that("afni_restricted_plan builds global AR(3) plan without MA", {
  set.seed(101)
  n <- 256
  resid <- matrix(rnorm(n * 3), n, 3)
  spec <- list(a = 0.6, r1 = 0.7, t1 = pi / 6)
  plan <- fmriAR:::afni_restricted_plan(resid = resid, runs = NULL, parcels = NULL,
                                        p = 3L, roots = spec, estimate_ma1 = FALSE)
  expect_identical(plan$method, "afni")
  expect_identical(plan$pooling, "global")
  expect_identical(plan$order[["p"]], 3L)
  expect_identical(plan$order[["q"]], 0L)
  expect_equal(plan$phi[[1]], fmriAR:::.afni_phi_ar3(spec$a, spec$r1, spec$t1))
  expect_identical(plan$theta[[1]], numeric(0))
})

test_that("afni_restricted_plan optionally estimates MA(1)", {
  withr::local_seed(202)
  n <- 1024
  spec <- list(a = 0.5, r1 = 0.6, t1 = pi / 5)
  phi <- fmriAR:::.afni_phi_ar3(spec$a, spec$r1, spec$t1)
  base_series <- sim_ar_process(phi, n)
  meas_noise <- rnorm(n, sd = 0.2)
  resid <- cbind(base_series + meas_noise, base_series)
  plan <- fmriAR:::afni_restricted_plan(resid = resid, runs = NULL, parcels = NULL,
                                        p = 3L, roots = spec, estimate_ma1 = TRUE)
  expect_identical(plan$order[["p"]], 3L)
  expect_identical(plan$order[["q"]], 1L)
  expect_length(plan$theta[[1]], 1L)
  expect_lt(abs(plan$theta[[1]]), 1)

  X <- matrix(1, n, 1)
  out <- whiten_apply(plan, X, resid[, 1, drop = FALSE])
  ehat <- drop(out$Y)
  ac <- acf(ehat, plot = FALSE, lag.max = 6, demean = TRUE)$acf[-1L]
  ci <- 1.96 / sqrt(length(ehat))
  expect_true(max(abs(ac)) < 4 * ci)
})

test_that("afni parcel plan whitens each parcel", {
  withr::local_seed(303)
  n <- 1500
  spec1 <- list(a = 0.55, r1 = 0.65, t1 = pi / 8)
  spec2 <- list(a = 0.35, r1 = 0.5, t1 = pi / 4)
  phi1 <- fmriAR:::.afni_phi_ar3(spec1$a, spec1$r1, spec1$t1)
  phi2 <- fmriAR:::.afni_phi_ar3(spec2$a, spec2$r1, spec2$t1)
  y1 <- sim_ar_process(phi1, n)
  y2 <- sim_ar_process(phi1, n)
  y3 <- sim_ar_process(phi2, n)
  y4 <- sim_ar_process(phi2, n)
  Y <- cbind(y1, y2, y3, y4)
  parcels <- c(1L, 1L, 2L, 2L)
  roots <- list("1" = spec1, "2" = spec2)
  plan <- fmriAR:::afni_restricted_plan(resid = Y, runs = NULL, parcels = parcels,
                                        p = 3L, roots = roots, estimate_ma1 = FALSE)
  expect_identical(plan$pooling, "parcel")
  expect_identical(plan$order[["p"]], 3L)
  expect_identical(plan$order[["q"]], 0L)
  expect_equal(sort(names(plan$phi_by_parcel)), c("1", "2"))
  expect_true(all(vapply(plan$phi_by_parcel, length, 0L) == 3L))

  X <- matrix(rnorm(n * 2), n, 2)
  out <- whiten_apply(plan, X, Y, parcels = parcels)
  Yw <- out$Y
  ci <- 1.96 / sqrt(n)
  for (j in seq_len(ncol(Yw))) {
    ac <- acf(Yw[, j], plot = FALSE, lag.max = 4, demean = TRUE)$acf[-1L]
    expect_true(max(abs(ac)) < 4 * ci)
  }
  expect_equal(length(out$X_by), length(unique(parcels)))
  expect_true(all(vapply(out$X_by, nrow, 0L) == n))
  expect_true(all(vapply(out$X_by, ncol, 0L) == ncol(X)))
})

test_that("afni parcel plan falls back to default roots and stores MA estimates", {
  withr::local_seed(505)
  n <- 900
  runs <- rep(1:3, each = n / 3)
  parcels <- rep(1:3, each = 2)
  spec_primary <- list(a = 0.45, r1 = 0.55, t1 = pi / 7)
  phi_primary <- fmriAR:::.afni_phi_ar3(spec_primary$a, spec_primary$r1, spec_primary$t1)

  # build residuals with the same AR across all parcels, plus small measurement noise
  base_series <- sim_ar_process(phi_primary, n)
  Y <- matrix(base_series, n, length(parcels)) + matrix(rnorm(n * length(parcels), sd = 0.1), n)

  # supply roots for only the first parcel; others should reuse this spec
  roots_partial <- list("1" = spec_primary)
  plan <- fmriAR:::afni_restricted_plan(
    resid = Y,
    runs = runs,
    parcels = parcels,
    p = 3L,
    roots = roots_partial,
    estimate_ma1 = TRUE
  )

  expect_identical(plan$order[["q"]], 1L)
  expect_equal(sort(names(plan$phi_by_parcel)), c("1", "2", "3"))
  expect_true(all(vapply(plan$phi_by_parcel, function(phi) all.equal(phi, phi_primary, tolerance = 1e-8) == TRUE, logical(1))))

  thetas <- plan$theta_by_parcel
  expect_true(all(vapply(thetas, length, 0L) == 1L))
  expect_true(all(abs(unlist(thetas)) < 1))

  X <- cbind(1, rnorm(n))
  out <- whiten_apply(plan, X, Y, runs = runs, parcels = parcels)
  ci <- 1.96 / sqrt(n)
  for (j in seq_len(ncol(out$Y))) {
    ac <- acf(out$Y[, j], plot = FALSE, lag.max = 4, demean = TRUE)$acf[-1L]
    expect_true(max(abs(ac)) < 4 * ci)
  }
})

test_that("afni plan with runs applies MA estimation per segment", {
  withr::local_seed(606)
  run_lengths <- c(300, 350, 250)
  runs <- rep(seq_along(run_lengths), run_lengths)
  n <- length(runs)
  spec <- list(a = 0.5, r1 = 0.6, t1 = pi / 6)
  phi <- fmriAR:::.afni_phi_ar3(spec$a, spec$r1, spec$t1)

  # simulate piecewise series with resets between runs
  Y <- numeric(n)
  idx_start <- c(1, cumsum(run_lengths)[-length(run_lengths)] + 1)
  for (i in seq_along(run_lengths)) {
    seg_idx <- seq(idx_start[i], length.out = run_lengths[i])
    seg_vals <- sim_ar_process(phi, run_lengths[i])
    Y[seg_idx] <- seg_vals + rnorm(run_lengths[i], sd = 0.15)
  }
  resid <- matrix(Y, n, 1)

  plan <- fmriAR:::afni_restricted_plan(
    resid = resid,
    runs = runs,
    parcels = NULL,
    p = 3L,
    roots = spec,
    estimate_ma1 = TRUE
  )

  expect_identical(plan$order[["q"]], 1L)
  expect_length(plan$theta[[1]], 1L)

  X <- cbind(1, rnorm(n))
  out <- whiten_apply(plan, X, resid, runs = runs)

  ci <- 1.96 / sqrt(n)
  ac <- acf(drop(out$Y), plot = FALSE, lag.max = 6, demean = TRUE)$acf[-1L]
  expect_true(max(abs(ac)) < 4 * ci)

  # ensure segment starts have finite values (no carryover across runs)
  run_starts <- c(1, idx_start[-1])
  expect_false(any(is.na(out$Y[run_starts, 1])))
})

test_that("afni AR(5) plan whitens matching process", {
  withr::local_seed(404)
  n <- 2500
  spec <- list(a = 0.5, r1 = 0.55, t1 = pi / 7,
               r2 = 0.4, t2 = pi / 3)
  phi <- fmriAR:::.afni_phi_ar5(spec$a, spec$r1, spec$t1, spec$r2, spec$t2)
  Y <- matrix(sim_ar_process(phi, n), n, 1)
  X <- cbind(1, sin(seq_len(n) / 20))
  plan <- fmriAR:::afni_restricted_plan(resid = Y, runs = NULL, parcels = NULL,
                                        p = 5L, roots = spec, estimate_ma1 = FALSE)
  expect_identical(plan$pooling, "global")
  expect_identical(plan$order[["p"]], 5L)
  expect_identical(plan$order[["q"]], 0L)
  expect_equal(plan$phi[[1]], phi)

  out <- whiten_apply(plan, X, Y)
  ehat <- drop(out$Y)
  ac <- acf(ehat, plot = FALSE, lag.max = 8, demean = TRUE)$acf[-1L]
  ci <- 1.96 / sqrt(length(ehat))
  expect_true(max(abs(ac)) < 4 * ci)
})

test_that("afni_restricted_plan validates inputs", {
  withr::local_seed(707)
  resid <- matrix(rnorm(300), 100, 3)
  spec <- list(a = 0.6, r1 = 0.7, t1 = pi / 6)

  expect_error(
    fmriAR:::afni_restricted_plan(
      resid = resid,
      runs = NULL,
      parcels = NULL,
      p = 4L,
      roots = spec,
      estimate_ma1 = FALSE
    ),
    regexp = "p %in% c"
  )

  bad_spec <- list(r1 = 0.7, t1 = pi / 6)
  expect_error(
    fmriAR:::afni_restricted_plan(
      resid = resid,
      runs = NULL,
      parcels = NULL,
      p = 3L,
      roots = bad_spec,
      estimate_ma1 = FALSE
    ),
    regexp = "roots"
  )

  expect_error(
    fmriAR:::afni_restricted_plan(
      resid = resid,
      runs = NULL,
      parcels = c(1L, 2L),
      p = 3L,
      roots = spec,
      estimate_ma1 = FALSE
    ),
    regexp = "not TRUE"
  )
})
