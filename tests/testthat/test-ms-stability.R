
test_that("Multiscale pooling stabilizes AR coefficients across train folds", {
  skip_if_no_fmriAR()
  skip_if_not_installed("abind")

  h <- make_hierarchy(n_coarse = 4L, medium_per_coarse = 3L, fine_per_medium = 3L, vox_per_fine = 4L)
  sim_all <- simulate_hier_ar2(h, n_train_per_run = 140L, n_test = 140L, runs_train = 3L, seed = 789)
  parcels <- sim_all$parcels_fine

  # Build 3 folds: leave-one-run-out across the 3 training runs
  n_per_run <- 140L
  total_train <- 3L * n_per_run
  run_idx <- rep(1:3, each = n_per_run)

  get_plan <- function(keep_runs, multiscale) {
    idx <- run_idx %in% keep_runs
    Y_tr <- sim_all$Y_train[idx, , drop = FALSE]
    X_tr <- sim_all$X_train[idx, , drop = FALSE]
    fmriAR::fit_noise(Y = Y_tr, X = X_tr, parcels = parcels,
                      pooling = "parcel", multiscale = multiscale, ms_mode = "acvf_pooled", p_target = 2L)
  }

  # Collect phi_by_parcel across folds
  phi_fine_list <- list()
  phi_ms_list   <- list()
  folds <- list(c(2,3), c(1,3), c(1,2))
  for (f in seq_along(folds)) {
    plan_f <- get_plan(folds[[f]], multiscale = FALSE)
    plan_m <- get_plan(folds[[f]], multiscale = TRUE)
    phi_fine_list[[f]] <- plan_phi_by_parcel(plan_f)
    phi_ms_list[[f]]   <- plan_phi_by_parcel(plan_m)
  }

  # Stack and compute SD across folds per parcel and lag
  phi_fine_arr <- abind::abind(phi_fine_list, along = 3L)  # p x K x F
  phi_ms_arr   <- abind::abind(phi_ms_list,   along = 3L)

  sd_fine <- apply(phi_fine_arr, c(1,2), stats::sd, na.rm = TRUE)
  sd_ms   <- apply(phi_ms_arr,   c(1,2), stats::sd, na.rm = TRUE)

  # Compare median SD across parcels and lags
  med_fine <- median(sd_fine)
  med_ms   <- median(sd_ms)

  # Expect at least ~10% reduction in dispersion
  expect_lte(med_ms, 0.9 * med_fine)
})
