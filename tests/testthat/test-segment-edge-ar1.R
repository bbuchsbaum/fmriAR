
test_that("Segment starts with exact_first_ar1 yield correctly scaled first innovations", {
  skip_if_not_installed("fmriAR")
  S <- 100L
  L <- 50L
  phi <- 0.8
  sigma <- 1.3

  sim <- simulate_ar1_runs(S = S, L = L, phi = phi, sigma = sigma, seed = 99L)
  y <- sim$y
  rs0 <- sim$run_starts0

  yw <- do_whiten_Y(y, phi = phi, theta = numeric(), run_starts0 = rs0,
                    exact_first_ar1 = TRUE, parallel = FALSE)

  first_idx <- as.integer(rs0 + 1L)
  first_innov <- yw[first_idx]

  # Theoretical variance at segment starts with exact_first_ar1 is sigma^2
  v_emp <- var(first_innov)
  expect_equal(v_emp, sigma^2, tolerance = 0.2)

  # Also check interior innovations have roughly the same variance
  interior <- setdiff(seq_along(yw), first_idx)
  expect_equal(var(yw[interior]), sigma^2, tolerance = 0.2)
})
