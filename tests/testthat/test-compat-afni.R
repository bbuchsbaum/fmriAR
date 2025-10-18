test_that("compat exposes AFNI restricted plan builder", {
  withr::local_seed(808)
  n <- 200
  resid <- matrix(rnorm(n * 3), n, 3)
  runs <- rep(1:2, each = n / 2)
  roots <- list(a = 0.6, r1 = 0.7, t1 = pi / 6)

  # via compat
  plan1 <- compat$afni_restricted_plan(resid, runs = runs, p = 3L,
                                       roots = roots, estimate_ma1 = FALSE)
  expect_s3_class(plan1, "fmriAR_plan")
  expect_identical(plan1$method, "afni")
  expect_true(plan1$order[["p"]] %in% c(3L, 5L))

  # reference via internal helper
  plan2 <- fmriAR:::afni_restricted_plan(resid, runs = runs, p = 3L,
                                         roots = roots, estimate_ma1 = FALSE)

  # core fields match
  expect_identical(plan1$method, plan2$method)
  expect_identical(plan1$pooling, plan2$pooling)
  expect_identical(plan1$order, plan2$order)
  expect_equal(plan1$phi, plan2$phi)
  expect_equal(plan1$theta, plan2$theta)

  # works with whiten_apply
  X <- cbind(1, rnorm(n))
  out <- whiten_apply(plan1, X, resid, runs = runs)
  expect_equal(dim(out$Y), dim(resid))
  expect_equal(dim(out$X), dim(X))
})

