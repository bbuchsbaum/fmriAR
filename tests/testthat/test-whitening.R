test_that("AR(1) whitening reduces lag-1 autocorrelation", {
  set.seed(1)
  n <- 600
  phi <- 0.6
  e <- rnorm(n)
  y <- as.numeric(stats::filter(e, filter = phi, method = "recursive"))
  Y <- cbind(y)
  X <- cbind(1)

  out <- arma_whiten_inplace(
    Y = Y,
    X = X,
    phi = phi,
    theta = numeric(0),
    run_starts = 0L,
    exact_first_ar1 = TRUE,
    parallel = FALSE
  )
  yw <- drop(out$Y)
  acf_vals <- stats::acf(yw, plot = FALSE, lag.max = 10, demean = TRUE)$acf[-1L]
  expect_true(all(abs(acf_vals[1:3]) < 0.15))
})

test_that("PACF <-> AR round-trips", {
  kap <- c(0.2, -0.1, 0.3)
  phi <- fmriAR:::pacf_to_ar(kap)
  kap2 <- fmriAR:::ar_to_pacf(phi)
  expect_equal(unname(kap2), unname(kap), tolerance = 1e-8)
})

test_that("run-specific parameters are applied per segment", {
  set.seed(11)
  n1 <- 80
  n2 <- 90
  phi1 <- 0.5
  phi2 <- -0.3
  y1 <- as.numeric(stats::filter(rnorm(n1), filter = phi1, method = "recursive"))
  y2 <- as.numeric(stats::filter(rnorm(n2), filter = phi2, method = "recursive"))
  Y <- rbind(matrix(y1, ncol = 1L), matrix(y2, ncol = 1L))
  X <- rbind(matrix(1, n1, 1L), matrix(1, n2, 1L))
  runs <- c(rep(1L, n1), rep(2L, n2))

  plan <- fmriAR:::new_whiten_plan(
    phi = list(phi1, phi2),
    theta = list(numeric(0), numeric(0)),
    order = c(p = 1L, q = 0L),
    runs = runs,
    exact_first = FALSE,
    method = "ar",
    pooling = "run"
  )

  out <- whiten_apply(plan, X, Y, parallel = FALSE)

  manual1 <- arma_whiten_inplace(matrix(y1, ncol = 1L), matrix(1, n1, 1L),
                                 phi = phi1, theta = numeric(0),
                                 run_starts = 0L, exact_first_ar1 = FALSE,
                                 parallel = FALSE)
  manual2 <- arma_whiten_inplace(matrix(y2, ncol = 1L), matrix(1, n2, 1L),
                                 phi = phi2, theta = numeric(0),
                                 run_starts = 0L, exact_first_ar1 = FALSE,
                                 parallel = FALSE)

  expect_equal(out$Y[seq_len(n1), , drop = FALSE], manual1$Y)
  expect_equal(out$Y[n1 + seq_len(n2), , drop = FALSE], manual2$Y)
})

test_that("censor gaps reset the whitening recursions", {
  set.seed(21)
  n <- 150
  phi <- 0.4
  y <- as.numeric(stats::filter(rnorm(n), filter = phi, method = "recursive"))
  Y <- matrix(y, ncol = 1L)
  X <- matrix(1, n, 1L)
  runs <- rep(1L, n)
  censor <- c(60L, 120L)

  plan <- fmriAR:::new_whiten_plan(
    phi = list(phi),
    theta = list(numeric(0)),
    order = c(p = 1L, q = 0L),
    runs = runs,
    exact_first = FALSE,
    method = "ar",
    pooling = "global"
  )

  out_plan <- whiten_apply(plan, X, Y, runs = runs, censor = censor, parallel = FALSE)

  run_starts <- as.integer(c(0L, censor))
  out_manual <- arma_whiten_inplace(Y, X,
                                    phi = phi, theta = numeric(0),
                                    run_starts = run_starts,
                                    exact_first_ar1 = FALSE,
                                    parallel = FALSE)

  expect_equal(out_plan$Y, out_manual$Y)
  expect_equal(out_plan$X, out_manual$X)
})

test_that("ARMA(1,1) whitening yields near-white innovations", {
  set.seed(31)
  n <- 800
  phi <- 0.3
  theta <- -0.25
  y <- as.numeric(stats::arima.sim(model = list(ar = phi, ma = theta), n = n))
  Y <- matrix(y, ncol = 1L)
  X <- matrix(1, n, 1L)

  out <- arma_whiten_inplace(Y, X,
                             phi = phi, theta = theta,
                             run_starts = 0L,
                             exact_first_ar1 = FALSE,
                             parallel = FALSE)
  innovations <- drop(out$Y)
  acf_vals <- stats::acf(innovations, lag.max = 10, plot = FALSE, demean = TRUE)$acf[-1L]
  expect_true(all(abs(acf_vals[1:5]) < 0.1))
})
