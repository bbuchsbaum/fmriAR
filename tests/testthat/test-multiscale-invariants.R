test_that("multiscale combiner is identity when all scales agree", {
  phi_ref <- fmriAR:::pacf_to_ar(c(0.25, -0.1))
  gamma_ref <- c(1.0, 0.45, 0.12)
  parents <- list(parent_medium = c(1L), parent_coarse = c(1L))
  sizes <- list(
    n_t = 120L,
    n_runs = 2L,
    beta = 0.5,
    coarse = list(`1` = 12L),
    medium = list(`1` = 12L),
    fine = list(`1` = 12L)
  )
  disp_list <- list(
    coarse = list(`1` = 0.2),
    medium = list(`1` = 0.2),
    fine = list(`1` = 0.2)
  )

  pacf_out <- fmriAR:::`.ms_combine_to_fine`(
    phi_by_coarse = list(`1` = phi_ref),
    phi_by_medium = list(`1` = phi_ref),
    phi_by_fine = list(`1` = phi_ref),
    parents = parents,
    sizes = sizes,
    disp_list = disp_list,
    p_target = 2L,
    mode = "pacf_weighted"
  )
  acvf_out <- fmriAR:::`.ms_combine_to_fine`(
    phi_by_coarse = list(`1` = phi_ref),
    phi_by_medium = list(`1` = phi_ref),
    phi_by_fine = list(`1` = phi_ref),
    acvf_by_coarse = list(`1` = gamma_ref),
    acvf_by_medium = list(`1` = gamma_ref),
    acvf_by_fine = list(`1` = gamma_ref),
    parents = parents,
    sizes = sizes,
    disp_list = disp_list,
    p_target = 2L,
    mode = "acvf_pooled"
  )

  expect_equal(pacf_out[["1"]], phi_ref, tolerance = 1e-12)
  expect_equal(
    acvf_out[["1"]],
    fmriAR:::enforce_stationary_ar(fmriAR:::yw_from_acvf_fast(gamma_ref, 2L)$phi),
    tolerance = 1e-12
  )
})

test_that("multiscale parcel fitting is invariant to voxel permutation", {
  set.seed(903)
  n <- 220
  v <- 24
  parcels_coarse <- rep(1:2, each = 12)
  parcels_medium <- rep(1:4, each = 6)
  parcels_fine <- rep(1:6, each = 4)

  X <- cbind(1, rnorm(n), seq(-1, 1, length.out = n))
  base_phi <- c(0.25, 0.45, 0.6, 0.35, 0.5, 0.7)
  Y <- matrix(0, nrow = n, ncol = v)
  for (j in seq_len(v)) {
    eps <- as.numeric(stats::arima.sim(list(ar = base_phi[parcels_fine[j]]), n = n))
    Y[, j] <- X %*% c(0.1, -0.2, 0.3) + eps
  }
  resid <- Y - X %*% qr.solve(X, Y)

  plan_a <- fit_noise(
    resid,
    method = "ar",
    p = 2L,
    pooling = "parcel",
    parcels = parcels_fine,
    parcel_sets = list(coarse = parcels_coarse, medium = parcels_medium, fine = parcels_fine),
    multiscale = TRUE,
    ms_mode = "acvf_pooled"
  )
  perm <- sample.int(v)
  plan_b <- fit_noise(
    resid[, perm, drop = FALSE],
    method = "ar",
    p = 2L,
    pooling = "parcel",
    parcels = parcels_fine[perm],
    parcel_sets = list(
      coarse = parcels_coarse[perm],
      medium = parcels_medium[perm],
      fine = parcels_fine[perm]
    ),
    multiscale = TRUE,
    ms_mode = "acvf_pooled"
  )

  flatten_phi <- function(phi_by) unlist(phi_by[sort(names(phi_by))], use.names = FALSE)
  expect_equal(flatten_phi(plan_a$phi_by_parcel), flatten_phi(plan_b$phi_by_parcel), tolerance = 1e-10)

  out_a <- whiten_apply(plan_a, X, Y, parcels = parcels_fine, parallel = FALSE)
  out_b <- whiten_apply(plan_b, X, Y[, perm, drop = FALSE], parcels = parcels_fine[perm], parallel = FALSE)

  expect_equal(out_a$Y, out_b$Y[, order(perm), drop = FALSE], tolerance = 1e-10)
})
