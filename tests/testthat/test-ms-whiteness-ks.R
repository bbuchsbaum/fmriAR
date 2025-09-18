
test_that("Multiscale pooling improves whiteness (KS to Uniform) on held-out run", {
  skip_if_no_fmriAR()

  # Build hierarchy and simulate
  h <- make_hierarchy(n_coarse = 4L, medium_per_coarse = 3L, fine_per_medium = 3L, vox_per_fine = 4L)
  sim <- simulate_hier_ar2(h, n_train_per_run = 160L, n_test = 160L, runs_train = 2L, seed = 123)
  parcels <- sim$parcels_fine

  # Fit noise models on TRAIN
  plan_fine <- fmriAR::fit_noise(Y = sim$Y_train, X = sim$X_train, parcels = parcels,
                                 pooling = "parcel", multiscale = FALSE, p_target = 2L)
  plan_ms   <- fmriAR::fit_noise(Y = sim$Y_train, X = sim$X_train, parcels = parcels,
                                 pooling = "parcel", multiscale = TRUE, ms_mode = "acvf_pooled", p_target = 2L)

  # Apply on TEST only
  w_fine <- fmriAR::whiten_apply(plan_fine, X = sim$X_test, Y = sim$Y_test, run_starts = sim$run_starts_test0)
  w_ms   <- fmriAR::whiten_apply(plan_ms,   X = sim$X_test, Y = sim$Y_test, run_starts = sim$run_starts_test0)

  # Whitened innovations
  E_fine <- w_fine$Y
  E_ms   <- w_ms$Y

  # p-values & KS
  p_fine <- lb_pvals(E_fine, lag = 10L)
  p_ms   <- lb_pvals(E_ms,   lag = 10L)
  KS_fine <- ks_to_uniform(p_fine)
  KS_ms   <- ks_to_uniform(p_ms)

  # Fraction of small p-values
  frac_fine <- mean(p_fine <= 0.05, na.rm = TRUE)
  frac_ms   <- mean(p_ms   <= 0.05, na.rm = TRUE)

  # Expect multiscale to be closer to Uniform (smaller KS) and lower rejection rate
  expect_true(KS_fine - KS_ms >= 0.02 || frac_fine - frac_ms >= 0.10)
})
