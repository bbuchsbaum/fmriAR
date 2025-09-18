
test_that("Multiscale pooling reduces held-out predictive NLL of innovations", {
  skip_if_no_fmriAR()

  h <- make_hierarchy(n_coarse = 4L, medium_per_coarse = 3L, fine_per_medium = 3L, vox_per_fine = 4L)
  sim <- simulate_hier_ar2(h, n_train_per_run = 160L, n_test = 160L, runs_train = 2L, seed = 456)
  parcels <- sim$parcels_fine

  plan_fine <- fmriAR::fit_noise(Y = sim$Y_train, X = sim$X_train, parcels = parcels,
                                 pooling = "parcel", multiscale = FALSE, p_target = 2L)
  plan_ms   <- fmriAR::fit_noise(Y = sim$Y_train, X = sim$X_train, parcels = parcels,
                                 pooling = "parcel", multiscale = TRUE, ms_mode = "acvf_pooled", p_target = 2L)

  w_fine <- fmriAR::whiten_apply(plan_fine, X = sim$X_test, Y = sim$Y_test, run_starts = sim$run_starts_test0)
  w_ms   <- fmriAR::whiten_apply(plan_ms,   X = sim$X_test, Y = sim$Y_test, run_starts = sim$run_starts_test0)

  nll_fine <- mean(series_nll(w_fine$Y))
  nll_ms   <- mean(series_nll(w_ms$Y))

  # Expect lower (better) average NLL under multiscale
  expect_lt(nll_ms, nll_fine - 0.5)
})
