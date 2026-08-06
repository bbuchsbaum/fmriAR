test_that("fit_noise honors fixed AR order requests", {
  set.seed(901)
  n <- 800
  y <- as.numeric(stats::arima.sim(list(ar = 0.7), n = n))
  resid <- matrix(y, ncol = 1L)

  plan_p1 <- fit_noise(resid, method = "ar", p = 1L, p_max = 5L)
  plan_p3 <- fit_noise(resid, method = "ar", p = 3L, p_max = 5L)

  expect_equal(plan_p1$order[["p"]], 1L)
  expect_equal(length(plan_p1$phi[[1]]), 1L)
  expect_equal(plan_p3$order[["p"]], 3L)
  expect_equal(length(plan_p3$phi[[1]]), 3L)
})

test_that("parcel AR fitting honors fixed order in plain and multiscale modes", {
  set.seed(902)
  n <- 240
  v <- 16
  parcels_fine <- rep(1:4, each = 4)
  parcels_medium <- rep(1:2, each = 8)
  parcels_coarse <- rep(1L, v)
  phi_by_voxel <- c(0.25, 0.35, 0.55, 0.65)[parcels_fine]

  resid <- matrix(0, nrow = n, ncol = v)
  for (j in seq_len(v)) {
    resid[, j] <- as.numeric(stats::arima.sim(list(ar = phi_by_voxel[j]), n = n))
  }

  plan_plain <- fit_noise(
    resid,
    method = "ar",
    p = 2L,
    p_max = 5L,
    multiscale = FALSE,
    pooling = "parcel",
    parcels = parcels_fine
  )
  plan_ms <- fit_noise(
    resid,
    method = "ar",
    p = 2L,
    p_max = 5L,
    pooling = "parcel",
    parcels = parcels_fine,
    parcel_sets = list(coarse = parcels_coarse, medium = parcels_medium, fine = parcels_fine),
    multiscale = TRUE,
    ms_mode = "acvf_pooled"
  )

  expect_equal(plan_plain$order[["p"]], 2L)
  expect_true(all(vapply(plan_plain$phi_by_parcel, length, integer(1)) == 2L))
  expect_equal(plan_ms$order[["p"]], 2L)
  expect_true(all(vapply(plan_ms$phi_by_parcel, length, integer(1)) == 2L))
})
