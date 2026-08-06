
# Multiscale pooling buys you accuracy when the fine-scale estimates are noisy:
# short runs and few voxels per parcel. Earlier versions of this file asserted a
# fixed 0.5-nat NLL improvement, which passed only because the non-multiscale
# baseline returned all-zero coefficients and therefore did no whitening at all.
# Both arms are now compared against the per-parcel truth that generated the data.

test_that("multiscale pooling estimates parcel AR coefficients more accurately", {
  skip_if_no_fmriAR()

  h <- make_hierarchy(n_coarse = 4L, medium_per_coarse = 3L,
                      fine_per_medium = 3L, vox_per_fine = 2L)
  sim <- simulate_hier_ar2(h, n_train_per_run = 60L, n_test = 200L,
                           runs_train = 2L, seed = 456)
  parcels <- sim$parcels_fine

  plan_fine <- fmriAR::fit_noise(Y = sim$Y_train, X = sim$X_train, parcels = parcels,
                                 pooling = "parcel", multiscale = FALSE, p_target = 2L,
                                 p = 2L)
  plan_ms <- fmriAR::fit_noise(Y = sim$Y_train, X = sim$X_train, parcels = parcels,
                               pooling = "parcel", multiscale = TRUE,
                               ms_mode = "acvf_pooled", p_target = 2L, p = 2L)

  rmse_fine <- parcel_phi_rmse(plan_fine, sim$phi_fine)
  rmse_ms   <- parcel_phi_rmse(plan_ms,   sim$phi_fine)

  # Both arms must be doing real estimation, not returning zeros.
  expect_gt(max(abs(unlist(plan_fine$phi_by_parcel))), 0.05)
  expect_gt(max(abs(unlist(plan_ms$phi_by_parcel))), 0.05)

  # Pooling across scales should beat the fine-scale-only estimate.
  expect_lt(rmse_ms, rmse_fine)
})

test_that("multiscale and fine-scale plans both genuinely whiten held-out data", {
  skip_if_no_fmriAR()

  h <- make_hierarchy(n_coarse = 4L, medium_per_coarse = 3L,
                      fine_per_medium = 3L, vox_per_fine = 2L)
  sim <- simulate_hier_ar2(h, n_train_per_run = 60L, n_test = 200L,
                           runs_train = 2L, seed = 456)
  parcels <- sim$parcels_fine

  plan_ms <- fmriAR::fit_noise(Y = sim$Y_train, X = sim$X_train, parcels = parcels,
                               pooling = "parcel", multiscale = TRUE,
                               ms_mode = "acvf_pooled", p_target = 2L, p = 2L)
  w_ms <- fmriAR::whiten_apply(plan_ms, X = sim$X_test, Y = sim$Y_test,
                               run_starts = sim$run_starts_test0, parcels = parcels)

  # Whitening must change the data and must reduce residual autocorrelation.
  expect_false(isTRUE(all.equal(w_ms$Y, sim$Y_test)))
  ac <- function(M) mean(abs(apply(M, 2L, function(y)
    stats::acf(y, lag.max = 2L, plot = FALSE, demean = TRUE)$acf[-1L])))
  expect_lt(ac(w_ms$Y), ac(sim$Y_test))
  expect_lt(ac(w_ms$Y), 0.12)
})
