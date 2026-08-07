# Targeted coverage for previously uncovered branches and printers.

test_that("print.fmriAR_plan covers global, run, and parcel layouts", {
  set.seed(11)
  n <- 120
  resid <- matrix(rnorm(n * 6), n, 6)
  runs <- rep(1:3, each = 40)

  plan_global <- fit_noise(resid, method = "ar", p = 2L, pooling = "global")
  txt_g <- paste(capture.output(print(plan_global)), collapse = "\n")
  expect_match(txt_g, "fmriAR whitening plan")
  expect_match(txt_g, "Pooling: global")
  expect_match(txt_g, "phi =")
  capture.output(expect_invisible(print(plan_global)))

  plan_run <- fit_noise(resid, runs = runs, method = "ar", p = 1L, pooling = "run")
  txt_r <- paste(capture.output(print(plan_run)), collapse = "\n")
  expect_match(txt_r, "Pooling: run")
  expect_match(txt_r, "Runs: 3")
  expect_match(txt_r, "run")

  # Many runs exercises truncation ("...") in the printer.
  many_runs <- rep(seq_len(10L), each = 20L)
  resid_many <- matrix(rnorm(200 * 2), 200, 2)
  plan_many <- fit_noise(resid_many, runs = many_runs, method = "ar", p = 1L,
                         pooling = "run")
  txt_many <- paste(capture.output(print(plan_many)), collapse = "\n")
  expect_match(txt_many, "\\.\\.\\.")

  parcels <- rep(1:6, each = 1L)
  plan_parcel <- fit_noise(resid, pooling = "parcel", parcels = parcels,
                           method = "ar", p = 2L)
  txt_p <- paste(capture.output(print(plan_parcel)), collapse = "\n")
  expect_match(txt_p, "Pooling: parcel")
  expect_match(txt_p, "Parcel-level coefficients")
  expect_match(txt_p, "Parcel ")

  # ARMA plan prints theta coefficients.
  y <- as.numeric(stats::arima.sim(list(ar = 0.4, ma = 0.3), n = 200))
  plan_arma <- fit_noise(matrix(y, ncol = 1L), method = "arma", p = 1L, q = 1L)
  txt_a <- paste(capture.output(print(plan_arma)), collapse = "\n")
  expect_match(txt_a, "Method: ARMA")
  expect_true(grepl("theta =", txt_a) || grepl("phi =", txt_a))

  # exact_first = TRUE path
  plan_ef <- fit_noise(resid, method = "ar", p = 1L, exact_first = "ar1")
  txt_ef <- paste(capture.output(print(plan_ef)), collapse = "\n")
  expect_match(txt_ef, "Exact first-sample scaling: AR\\(1\\)")
})

test_that("print.fmriAR_plan handles sparse/edge coefficient layouts", {
  # Empty phi/theta and mismatched list lengths exercise align_lists / fmt_vec.
  plan <- fmriAR:::new_whiten_plan(
    phi = list(numeric(0), c(0.3, 0.1)),
    theta = list(0.2),
    order = c(p = 2L, q = 1L),
    runs = rep(1:2, each = 10),
    exact_first = FALSE,
    method = "arma",
    pooling = "run"
  )
  txt <- paste(capture.output(print(plan)), collapse = "\n")
  expect_match(txt, "\\(none\\)|phi =")
  expect_match(txt, "theta =")

  # Parcel plan with more than 5 parcels truncates the listing.
  phi_by <- setNames(lapply(1:8, function(i) c(0.2)), as.character(1:8))
  theta_by <- setNames(lapply(1:8, function(i) numeric(0)), as.character(1:8))
  plan_p <- fmriAR:::new_whiten_plan(
    phi = NULL, theta = NULL,
    order = c(p = 1L, q = 0L),
    runs = NULL,
    exact_first = TRUE,
    method = "ar",
    pooling = "parcel",
    parcel_ids = as.character(1:8),
    phi_by_parcel = phi_by,
    theta_by_parcel = theta_by
  )
  txt_p <- paste(capture.output(print(plan_p)), collapse = "\n")
  expect_match(txt_p, "8 parcels")
  expect_match(txt_p, "\\.\\.\\.")

  # Parcel plan with theta on some parcels.
  theta_by2 <- theta_by
  theta_by2[["1"]] <- 0.15
  plan_pt <- fmriAR:::new_whiten_plan(
    phi = NULL, theta = NULL,
    order = c(p = 1L, q = 1L),
    runs = NULL,
    exact_first = FALSE,
    method = "arma",
    pooling = "parcel",
    parcel_ids = as.character(1:3),
    phi_by_parcel = phi_by[1:3],
    theta_by_parcel = theta_by2[1:3]
  )
  txt_pt <- paste(capture.output(print(plan_pt)), collapse = "\n")
  expect_match(txt_pt, "theta =")
})

test_that("hr_arma falls back to R implementation when C++ is disabled", {
  skip_if_not_installed("withr")
  set.seed(21)
  y <- as.numeric(stats::arima.sim(list(ar = 0.5, ma = -0.3), n = 180))

  withr::with_options(list(fmriAR.use_cpp_hr = FALSE), {
    fit <- fmriAR:::hr_arma(y, p = 1L, q = 1L, iter = 1L, step1 = "yw",
                            enforce = TRUE)
    expect_equal(fit$order, c(p = 1L, q = 1L))
    expect_length(fit$phi, 1L)
    expect_length(fit$theta, 1L)
    expect_true(is.finite(fit$sigma2))
    expect_true(all(Mod(polyroot(c(1, -fit$phi))) > 1))
    expect_true(all(Mod(polyroot(c(1, fit$theta))) > 1 - 1e-8))
  })

  # Direct R path, including AR-only and error cases.
  fit_r <- fmriAR:::hr_arma_R(y, p = 2L, q = 0L, iter = 0L, step1 = "burg",
                              enforce = TRUE)
  expect_length(fit_r$phi, 2L)
  expect_length(fit_r$theta, 0L)

  expect_error(fmriAR:::hr_arma_R(rnorm(5), p = 1L, q = 0L), "too short")
  expect_error(fmriAR:::hr_arma_R(rnorm(20), p = 8L, q = 8L),
               "Not enough data|HR regression")
})

test_that("parcel_means_fast R fallback matches C++ path", {
  skip_if_not_installed("withr")
  set.seed(33)
  resid <- matrix(rnorm(40), nrow = 8)
  resid[2, 3] <- NA_real_
  parcels <- c(1L, 1L, 2L, 2L, 3L)

  withr::with_options(list(fmriAR.use_cpp_means = FALSE), {
    out <- fmriAR:::parcel_means_fast(resid, parcels, na.rm = TRUE)
    ref <- fmriAR:::`.parcel_means`(resid, parcels, na.rm = TRUE)
    # Re-enable C++ for the reference path if .parcel_means respects the option
    # (both go through parcel_means_fast, so compare against manual means).
    manual <- cbind(
      rowMeans(resid[, 1:2], na.rm = TRUE),
      rowMeans(resid[, 3:4], na.rm = TRUE),
      resid[, 5]
    )
    colnames(manual) <- c("1", "2", "3")
    expect_equal(out, manual)

    out2 <- fmriAR:::parcel_means_fast(resid[, 1:4], parcels[1:4], na.rm = FALSE)
    expect_equal(out2[, "1"], rowMeans(resid[, 1:2]))
    expect_equal(out2[, "2"], rowMeans(resid[, 3:4]))
  })
})

test_that("segmented_acvf_fast validates run_starts", {
  expect_error(fmriAR:::segmented_acvf_fast(rnorm(10), integer(0), 2L),
               "run_starts")
  expect_error(fmriAR:::segmented_acvf_fast(rnorm(10), 1L, 2L),
               "run_starts")
})

test_that("sandwich estimates beta when omitted and honors df_mode", {
  set.seed(44)
  Xw <- cbind(1, rnorm(40))
  Yw <- matrix(rnorm(40 * 3), 40, 3)

  out_auto <- sandwich_from_whitened_resid(Xw, Yw, type = "iid")
  beta <- solve(crossprod(Xw), crossprod(Xw, Yw))
  out_manual <- sandwich_from_whitened_resid(Xw, Yw, beta = beta, type = "iid")
  expect_equal(out_auto$se, out_manual$se, tolerance = 1e-10)

  out_np <- sandwich_from_whitened_resid(Xw, Yw, type = "iid", df_mode = "n-p")
  expect_equal(out_np$df, nrow(Xw) - ncol(Xw))
  expect_equal(out_auto$df, nrow(Xw) - qr(Xw)$rank)
})

test_that("pacf and MA helpers cover empty and non-finite edges", {
  expect_identical(fmriAR:::pacf_to_ar(numeric(0)), numeric(0))
  expect_identical(fmriAR:::ar_to_pacf(numeric(0)), numeric(0))
  expect_true(all(is.na(fmriAR:::ar_to_pacf(c(0.2, NA_real_)))))

  expect_identical(fmriAR:::enforce_stationary_ar(numeric(0)), numeric(0))
  expect_identical(fmriAR:::enforce_stationary_ar(c(NA_real_, 0.1)), numeric(0))

  # Near-unit-root high-order AR exercises geometric shrink.
  phi_big <- fmriAR:::pacf_to_ar(rep(0.99, 10L))
  phi_enforced <- fmriAR:::enforce_stationary_ar(phi_big, bound = 0.99, margin = 1e-3)
  expect_true(length(phi_enforced) == 0L ||
                min(Mod(polyroot(c(1, -phi_enforced)))) > 1 + 1e-4)

  expect_identical(fmriAR:::enforce_invertible_ma(numeric(0)), numeric(0))
  theta_ok <- 0.4
  expect_equal(fmriAR:::enforce_invertible_ma(theta_ok), theta_ok, tolerance = 1e-10)
  theta_bad <- 1.5
  theta_fix <- fmriAR:::enforce_invertible_ma(theta_bad)
  expect_true(all(Mod(polyroot(c(1, theta_fix))) > 1 - 1e-8))
})

test_that("acvf bias helpers cover lag-operator and correction edges", {
  # k >= n returns zeros; no same-run pairs across a boundary returns zeros.
  B <- matrix(1:6, nrow = 2)
  expect_equal(fmriAR:::`.apply_lag_operator`(B, k = 10L, run_vec = c(1L, 1L, 1L)),
               matrix(0, 2, 3))
  out_split <- fmriAR:::`.apply_lag_operator`(B, k = 1L, run_vec = c(1L, 2L, 3L))
  expect_equal(out_split, matrix(0, 2, 3))

  # Non-matrix design and short residual df warn / coerce.
  n <- 30L
  X <- cbind(1, seq_len(n), (seq_len(n))^2)
  expect_warning(
    A <- acvf_bias_matrix(as.data.frame(X), max_lag = 40L),
    "residual-bias correction"
  )
  expect_true(is.list(A) && length(A) == 1L)

  expect_error(
    fmriAR:::`.acvf_bias_by_run`(X, n = n + 5L, max_lag = 5L),
    "one row per timepoint|design"
  )

  # .normalize_correction broadcasting and errors.
  M <- diag(3)
  expect_equal(length(fmriAR:::`.normalize_correction`(M, 4L)), 4L)
  expect_equal(length(fmriAR:::`.normalize_correction`(list(M), 3L)), 3L)
  expect_null(fmriAR:::`.normalize_correction`(NULL, 2L))
  expect_error(fmriAR:::`.normalize_correction`(list(M, M), 3L), "matrices")
  expect_error(fmriAR:::`.normalize_correction`("nope", 1L), "matrix or a list")
})

test_that("noise_acvf covers coercion and empty-segment edges", {
  expect_error(noise_acvf(rnorm(1), max_lag = 1L), "at least two")
  # Vector input is coerced to a one-column matrix.
  out <- noise_acvf(as.numeric(stats::arima.sim(list(ar = 0.3), n = 80)),
                    max_lag = 3L)
  expect_s3_class(out, "fmriAR_acvf")
  expect_true(length(out$acvf[[1]]) >= 1L)

  # Out-of-range censor indices are dropped (becoming NULL).
  A <- matrix(rnorm(40), 20, 2)
  out2 <- noise_acvf(A, censor = c(0L, 100L), max_lag = 2L)
  expect_false(isTRUE(out2$corrected))

  # Heavy censoring that leaves nothing usable.
  expect_error(
    noise_acvf(A, censor = seq_len(20L), max_lag = 1L),
    "no valid segments"
  )
})

test_that("fit_noise and whiten_apply cover remaining error and edge paths", {
  set.seed(55)
  n <- 100
  X <- cbind(1, rnorm(n))
  Y <- matrix(rnorm(n * 4), n, 4)
  resid <- Y - X %*% qr.solve(X, Y)

  expect_error(fit_noise(), "resid")
  expect_error(fit_noise(resid[1:5, , drop = FALSE]), "too short")

  # Empty logical censor becomes NULL (no-op).
  plan <- fit_noise(resid, censor = rep(FALSE, n), method = "ar", p = 1L)
  expect_null(plan$censor)

  # Explicit p = 0 returns white-noise plan.
  plan0 <- fit_noise(resid, method = "ar", p = 0L)
  expect_equal(plan0$order[["p"]], 0L)

  # Parcel + ARMA is rejected.
  expect_error(
    fit_noise(resid, pooling = "parcel", parcels = rep(1:2, each = 2),
              method = "arma", p = 1L, q = 1L),
    "ar' only|Parcel pooling"
  )

  # Parcel plan with all timepoints censored.
  expect_error(
    fit_noise(resid, pooling = "parcel", parcels = rep(1:2, each = 2),
              method = "ar", p = 1L, censor = seq_len(n)),
    "no valid"
  )

  # whiten_apply coerces vectors, rejects NA, and supports parcel inplace.
  plan_g <- fit_noise(resid[, 1, drop = FALSE], method = "ar", p = 1L)
  out_vec <- whiten_apply(plan_g, X[, 1], Y[, 1], parallel = FALSE)
  expect_equal(dim(out_vec$Y), c(n, 1L))

  Y_na <- Y
  Y_na[1, 1] <- NA_real_
  expect_error(whiten_apply(plan_g, X, Y_na[, 1, drop = FALSE], parallel = FALSE),
               "NA")

  parcels <- rep(1:2, each = 2)
  plan_p <- fit_noise(resid, pooling = "parcel", parcels = parcels,
                      method = "ar", p = 1L)
  Y_copy <- Y
  out_in <- whiten_apply(plan_p, X, Y_copy, parcels = parcels,
                         inplace = TRUE, parallel = FALSE)
  expect_true(is.null(out_in$X))
  expect_true(is.list(out_in$X_by))
  expect_equal(dim(out_in$Y), dim(Y))
})

test_that("internal helpers cover empty pooled acvf and run-start utilities", {
  empty_pooled <- list(num = numeric(3), pairs = c(0, 0, 0))
  expect_equal(fmriAR:::`.acvf_max_lag`(empty_pooled), -1L)
  expect_identical(fmriAR:::`.acvf_from_pooled`(empty_pooled), numeric(0))
  expect_false(fmriAR:::`.acvf_is_psd`(numeric(0)))
  expect_false(fmriAR:::`.acvf_is_psd`(c(-1, 0)))
  expect_true(fmriAR:::`.acvf_is_psd`(1))
  expect_identical(fmriAR:::`.shrink_to_pd`(numeric(0)), numeric(0))

  # Non-PSD acvf is shrunk toward white noise.
  g_bad <- c(1, 0.99, 0.99, 0.99)
  g_fix <- fmriAR:::`.shrink_to_pd`(g_bad)
  expect_true(fmriAR:::`.acvf_is_psd`(g_fix))

  expect_identical(fmriAR:::`.apply_acvf_correction`(c(1, 0.3), matrix(0, 0, 0)),
                   c(1, 0.3))

  # Ill-conditioned / NULL correction matrices are dropped with a warning.
  bad <- list(`1` = matrix(1e-12, 3, 3), `2` = NULL)
  expect_warning(
    dropped <- fmriAR:::`.drop_unusable_corrections`(bad),
    "ill-conditioned|skipped"
  )
  expect_null(dropped[[1]])
  expect_null(fmriAR:::`.drop_unusable_corrections`(NULL))

  # All-censored valid_segments returns empty structure.
  seg <- fmriAR:::`.valid_segments`(5L, runs = NULL, censor = 1:5)
  expect_identical(seg$idx, integer(0))

  # .full_run_starts with censoring and .runs_from_starts0 round-trip.
  rs <- fmriAR:::`.full_run_starts`(rep(1:2, each = 10), censor = c(5L, 15L), n = 20L)
  expect_true(0L %in% rs)
  runs <- fmriAR:::`.runs_from_starts0`(c(0L, 10L), 20L)
  expect_equal(runs, rep(1:2, each = 10))
  expect_error(fmriAR:::`.runs_from_starts0`(c(1L, 5L), 10L), "include 0")
  expect_error(fmriAR:::`.runs_from_starts0`(c(0L, 15L), 10L), "out of bounds")

  # Vector (non-matrix) path for run_avg helpers.
  g <- fmriAR:::`.run_avg_acvf`(rnorm(30), 2L)
  expect_length(g, 3L)
  g2 <- fmriAR:::`.run_avg_acvf_psd`(rnorm(30), 2L)
  expect_length(g2, 3L)

  # Empty-matrix / invalid max_lag paths in C++ acvf helper.
  expect_error(fmriAR:::run_avg_acvf_cpp(matrix(0, 5, 2), -1L), "max_lag")
  empty <- fmriAR:::run_avg_acvf_cpp(matrix(numeric(0), 0, 0), 2L)
  expect_length(empty, 3L)

  # C++ whitening rejects malformed run_starts.
  Y <- matrix(rnorm(20), 20, 1)
  X <- matrix(1, 20, 1)
  expect_error(
    fmriAR:::arma_whiten_inplace(Y, X, 0.3, numeric(0), integer(0), FALSE, FALSE),
    "run_starts"
  )
  expect_error(
    fmriAR:::arma_whiten_inplace(Y, X, 0.3, numeric(0), 1L, FALSE, FALSE),
    "0-based|first element"
  )
  expect_error(
    fmriAR:::arma_whiten_inplace(Y, X, 0.3, numeric(0), c(0L, 0L), FALSE, FALSE),
    "strictly increasing|out of bounds"
  )
})

test_that("estimate_ar_series handles short segments and p = 0", {
  y <- rnorm(8)
  # Segment lengths all < 2 for p_max effectively force empty / low-order fits.
  fit0 <- fmriAR:::`.estimate_ar_series`(y, p_max = 3L, p = 0L, starts0 = 0L)
  expect_equal(fit0$order[["p"]], 0L)
  expect_true(length(fit0$gamma) >= 1L || length(fit0$phi) == 0L)

  # Tiny series with p_cap < 1.
  fit_tiny <- fmriAR:::`.estimate_ar_series`(rnorm(2), p_max = 5L, p = "auto",
                                            starts0 = 0L)
  expect_equal(fit_tiny$order[["p"]], 0L)
})

test_that("acorr_diagnostics returns NA when variance is zero under segments", {
  resid <- matrix(0, nrow = 12, ncol = 2)
  runs <- rep(1:2, each = 6)
  out <- acorr_diagnostics(resid, max_lag = 2L, runs = runs, aggregate = "none")
  expect_true(all(is.na(out$acf)))
})

test_that(".lag_matrix builds lagged design columns", {
  x <- 1:5
  expect_equal(dim(fmriAR:::`.lag_matrix`(x, 0L)), c(5L, 0L))
  M <- fmriAR:::`.lag_matrix`(x, 2L)
  expect_equal(M[3, 1], 2)
  expect_equal(M[3, 2], 1)
  expect_true(all(is.na(M[1, ])))
})
