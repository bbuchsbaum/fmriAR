# Correctness of AR estimation against known ground truth.
#
# Every test here asserts recovery of coefficients that generated the data, or a
# structural invariant, rather than comparing two fmriAR outputs to each other.
# Two shipped bugs were invisible to structural tests: parcel pooling returned
# all-zero coefficients (a silent no-op), and censoring drove the global path to
# non-stationary coefficients that amplified noise instead of whitening it.

ar_sim <- function(n, phi, nvox = 1L, sd = 1, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  burn <- 200L
  p <- length(phi)
  out <- matrix(0, n, nvox)
  for (j in seq_len(nvox)) {
    e <- stats::rnorm(n + burn, sd = sd)
    y <- numeric(n + burn)
    for (t in seq_len(n + burn)) {
      lagsum <- 0
      for (k in seq_len(p)) if (t > k) lagsum <- lagsum + phi[k] * y[t - k]
      y[t] <- lagsum + e[t]
    }
    out[, j] <- y[(burn + 1L):(burn + n)]
  }
  out
}

plan_phi <- function(plan, parcel = 1L) {
  if (!is.null(plan$phi_by_parcel)) return(plan$phi_by_parcel[[parcel]])
  if (is.list(plan$phi)) return(plan$phi[[1L]])
  plan$phi
}

min_root <- function(phi) {
  phi <- phi[seq_len(max(which(c(TRUE, phi != 0)) ) - 1L)]
  if (!length(phi)) return(Inf)
  min(Mod(polyroot(c(1, -phi))))
}

lag1 <- function(M) {
  mean(apply(M, 2L, function(y)
    stats::acf(y, lag.max = 1L, plot = FALSE, demean = TRUE)$acf[2L]))
}


# --- Recovery of known coefficients ------------------------------------------

test_that("every pooling mode recovers a known AR(1) coefficient", {
  # Run each mode with a fixed order AND with p = "auto". The auto path is where
  # order selection lives, and it is the path that used to collapse to p = 0.
  resid <- ar_sim(400L, 0.5, nvox = 40L, seed = 101)
  parcels <- rep(1:4, length.out = 40L)

  for (pspec in list(1L, "auto")) {
    plans <- list(
      global = fmriAR::fit_noise(resid, pooling = "global", method = "ar",
                                 p = pspec, p_max = 4L),
      run    = fmriAR::fit_noise(resid, runs = rep(1:2, each = 200L),
                                 pooling = "run", method = "ar",
                                 p = pspec, p_max = 4L),
      parcel = fmriAR::fit_noise(resid, parcels = parcels, pooling = "parcel",
                                 method = "ar", p = pspec, p_max = 4L),
      pacf   = fmriAR::fit_noise(resid, parcels = parcels, pooling = "parcel",
                                 method = "ar", p = pspec, p_max = 4L,
                                 multiscale = TRUE, ms_mode = "pacf_weighted",
                                 p_target = 1L),
      acvf   = fmriAR::fit_noise(resid, parcels = parcels, pooling = "parcel",
                                 method = "ar", p = pspec, p_max = 4L,
                                 multiscale = TRUE, ms_mode = "acvf_pooled",
                                 p_target = 1L)
    )
    label <- if (identical(pspec, "auto")) "auto" else "fixed"
    for (nm in names(plans)) {
      phi <- plan_phi(plans[[nm]])
      expect_equal(phi[1], 0.5, tolerance = 0.08,
                   info = paste("pooling mode:", nm, "| p:", label))
    }
  }
})

test_that("every pooling mode recovers known AR(2) coefficients", {
  # Parcel modes fit a single parcel-mean series, so they carry single-series
  # sampling error of roughly 1/sqrt(n) per coefficient, unlike global pooling
  # which averages the autocovariance over every voxel. Averaging over seeds
  # separates estimator bias (what this test is for) from that sampling noise.
  truth <- c(0.6, 0.25)
  modes <- c("global", "parcel", "pacf", "acvf")
  acc <- setNames(lapply(modes, function(x) matrix(NA_real_, 5L, 2L)), modes)

  for (s in 1:5) {
    resid <- ar_sim(400L, truth, nvox = 40L, seed = 200 + s)
    parcels <- rep(1:4, length.out = 40L)
    plans <- list(
      global = fmriAR::fit_noise(resid, pooling = "global", method = "ar", p = 2L),
      parcel = fmriAR::fit_noise(resid, parcels = parcels, pooling = "parcel",
                                 method = "ar", p = 2L),
      pacf   = fmriAR::fit_noise(resid, parcels = parcels, pooling = "parcel",
                                 method = "ar", p = 2L, multiscale = TRUE,
                                 ms_mode = "pacf_weighted", p_target = 2L),
      acvf   = fmriAR::fit_noise(resid, parcels = parcels, pooling = "parcel",
                                 method = "ar", p = 2L, multiscale = TRUE,
                                 ms_mode = "acvf_pooled", p_target = 2L)
    )
    for (nm in modes) acc[[nm]][s, ] <- plan_phi(plans[[nm]])[seq_len(2)]
  }

  for (nm in modes) {
    expect_equal(colMeans(acc[[nm]]), truth, tolerance = 0.08,
                 info = paste("pooling mode:", nm))
  }
})

test_that("p = 'auto' selects the true order and recovers its coefficients", {
  # Regression guard: order selection in the parcel path scored an innovation
  # sequence built with a feedback filter instead of an FIR one, inflating the
  # variance at every order so that p = 0 always won.
  for (pooling in c("global", "parcel")) {
    resid <- ar_sim(500L, c(0.6, 0.25), nvox = 30L, seed = 303)
    args <- list(resid = resid, method = "ar", p = "auto", p_max = 6L,
                 pooling = pooling)
    if (pooling == "parcel") args$parcels <- rep(1:3, length.out = 30L)
    plan <- do.call(fmriAR::fit_noise, args)

    expect_equal(plan$order[["p"]], 2L, info = pooling)
    phi <- plan_phi(plan)
    expect_equal(phi[seq_len(2)], c(0.6, 0.25), tolerance = 0.10, info = pooling)
  }
})

test_that("auto order selection returns p = 0 on genuine white noise", {
  set.seed(404)
  resid <- matrix(stats::rnorm(400L * 30L), 400L, 30L)
  plan <- fmriAR::fit_noise(resid, pooling = "global", method = "ar",
                            p = "auto", p_max = 6L)
  expect_lte(plan$order[["p"]], 1L)
  expect_lt(max(abs(plan_phi(plan)), 0), 0.12)
})


# --- No silent no-ops ---------------------------------------------------------

test_that("a fitted plan never silently returns the data unchanged", {
  # Regression guard: pooling = "parcel" with the default p = "auto" produced
  # all-zero coefficients while reporting order p = p_max, so whiten_apply
  # returned X and Y bit-identical to their inputs.
  resid <- ar_sim(300L, 0.6, nvox = 24L, seed = 505)
  parcels <- rep(1:4, length.out = 24L)
  X <- cbind(1, stats::rnorm(300L))

  for (ms in list(NULL, "pacf_weighted", "acvf_pooled")) {
    args <- list(resid = resid, parcels = parcels, pooling = "parcel",
                 method = "ar", p = "auto", p_max = 4L)
    if (!is.null(ms)) { args$multiscale <- TRUE; args$ms_mode <- ms }
    plan <- do.call(fmriAR::fit_noise, args)
    label <- if (is.null(ms)) "no multiscale" else ms

    expect_gt(max(abs(unlist(plan$phi_by_parcel))), 0.1, label = label)

    out <- fmriAR::whiten_apply(plan, X, resid, parcels = parcels)
    expect_false(identical(out$Y, resid), label = label)
    expect_gt(max(abs(out$Y - resid)), 1e-8, label = label)
  }
})

test_that("reported order matches the number of non-trivial coefficients", {
  resid <- ar_sim(300L, 0.5, nvox = 20L, seed = 606)
  plan <- fmriAR::fit_noise(resid, pooling = "global", method = "ar",
                            p = "auto", p_max = 5L)
  expect_equal(length(plan_phi(plan)), plan$order[["p"]])
})

test_that("whitening measurably reduces autocorrelation", {
  resid <- ar_sim(400L, 0.6, nvox = 30L, seed = 707)
  X <- cbind(1, stats::rnorm(400L))
  plan <- fmriAR::fit_noise(resid, pooling = "global", method = "ar", p = 1L)
  out <- fmriAR::whiten_apply(plan, X, resid)

  expect_gt(lag1(resid), 0.4)
  expect_lt(abs(lag1(out$Y)), 0.08)
})


# --- Censoring ----------------------------------------------------------------

test_that("censoring does not attenuate the AR estimate", {
  # Regression guard: each scrubbing fragment was centered on its own mean,
  # which destroys the autocorrelation it is meant to measure. phi fell from
  # 0.6 to roughly 0 as censoring reached 40%.
  truth <- 0.6
  for (frac in c(0, 0.2, 0.3, 0.4)) {
    est <- vapply(1:5, function(s) {
      resid <- ar_sim(250L, truth, nvox = 30L, seed = 800 + s)
      set.seed(900 + s)
      cens <- if (frac == 0) NULL else sort(sample.int(250L, round(frac * 250L)))
      plan <- fmriAR::fit_noise(resid, censor = cens, pooling = "global",
                                method = "ar", p = 1L)
      plan_phi(plan)[1]
    }, numeric(1))
    expect_equal(mean(est), truth, tolerance = 0.12,
                 info = sprintf("censor fraction %.2f", frac))
  }
})

test_that("censoring never yields non-stationary coefficients", {
  # Regression guard: at 25-50% scrubbing the global path returned coefficients
  # whose characteristic roots fell inside the unit circle, and whiten_apply
  # then amplified variance instead of reducing it.
  for (frac in c(0.2, 0.3, 0.4, 0.5)) {
    for (s in 1:4) {
      resid <- ar_sim(200L, 0.6, nvox = 25L, seed = 1000 + s)
      set.seed(1100 + s)
      cens <- sort(sample.int(200L, round(frac * 200L)))
      plan <- fmriAR::fit_noise(resid, censor = cens, pooling = "global",
                                method = "ar", p = "auto", p_max = 6L)
      expect_gt(min_root(plan_phi(plan)), 1,
                label = sprintf("min|root| at censor %.2f seed %d", frac, s))
    }
  }
})

test_that("whitening a censored fit reduces rather than amplifies variance", {
  resid <- ar_sim(200L, 0.6, nvox = 25L, seed = 1201)
  set.seed(1202)
  cens <- sort(sample.int(200L, 60L))
  X <- cbind(1, stats::rnorm(200L))
  plan <- fmriAR::fit_noise(resid, censor = cens, pooling = "global",
                            method = "ar", p = "auto", p_max = 6L)
  out <- fmriAR::whiten_apply(plan, X, resid)
  expect_lt(stats::var(as.numeric(out$Y)), stats::var(as.numeric(resid)))
})

test_that("censored frames do not leak into the estimate", {
  # Frames that are censored may hold arbitrary garbage; excluding them must
  # make the fit identical to fitting the clean series with the same gaps.
  resid <- ar_sim(300L, 0.5, nvox = 20L, seed = 1301)
  cens <- seq(10L, 290L, by = 10L)
  spoiled <- resid
  spoiled[cens, ] <- spoiled[cens, ] + 500

  a <- fmriAR::fit_noise(resid,   censor = cens, pooling = "global", method = "ar", p = 2L)
  b <- fmriAR::fit_noise(spoiled, censor = cens, pooling = "global", method = "ar", p = 2L)
  expect_equal(plan_phi(a), plan_phi(b), tolerance = 1e-8)
})

test_that("parcel pooling honours the censor argument", {
  # Regression guard: the parcel path passed censor = NULL internally, so its
  # estimates were bit-identical with and without censoring.
  resid <- ar_sim(300L, 0.5, nvox = 20L, seed = 1401)
  parcels <- rep(1:4, length.out = 20L)
  cens <- seq(10L, 290L, by = 10L)
  spoiled <- resid
  spoiled[cens, ] <- spoiled[cens, ] + 500

  uncensored <- fmriAR::fit_noise(spoiled, parcels = parcels, pooling = "parcel",
                                  method = "ar", p = 1L)
  censored <- fmriAR::fit_noise(spoiled, parcels = parcels, pooling = "parcel",
                                method = "ar", p = 1L, censor = cens)

  expect_false(isTRUE(all.equal(plan_phi(uncensored), plan_phi(censored))))
  expect_equal(plan_phi(censored)[1], 0.5, tolerance = 0.12)
})


# --- Run and segment boundaries ----------------------------------------------

test_that("per-run mean offsets do not contaminate the estimate", {
  # Regression guard: the parcel path estimated on the concatenated series, so
  # a between-run level shift was read as near-perfect autocorrelation
  # (phi ~ 0.99 on data whose true phi was 0.4, and on white noise).
  set.seed(1501)
  n <- 400L; nv <- 40L
  runs <- rep(1:4, each = 100L)
  white <- matrix(stats::rnorm(n * nv), n, nv)
  offsets <- rep(stats::rnorm(4L, sd = 10), each = 100L)
  shifted <- white + offsets

  for (pooling in c("global", "run", "parcel")) {
    args <- list(resid = shifted, runs = runs, pooling = pooling,
                 method = "ar", p = 1L)
    if (pooling == "parcel") args$parcels <- rep(1:4, length.out = nv)
    plan <- do.call(fmriAR::fit_noise, args)
    expect_lt(abs(plan_phi(plan)[1]), 0.15, label = paste("offset leak,", pooling))
  }
})

test_that("lag products never span a run boundary", {
  # Alternating run means make any boundary-spanning product strongly negative,
  # so a fit that ignores boundaries is pulled well away from the truth. The
  # offsets must leave the estimate unchanged, which is the sharpest form of
  # this check: compare each fit against the same data without offsets.
  n <- 480L; nv <- 30L
  runs <- rep(1:6, each = 80L)
  offs <- rep(rep(c(-8, 8), 3L), each = 80L)

  for (pooling in c("global", "run", "parcel")) {
    for (s in 1:3) {
      base <- ar_sim(n, 0.5, nvox = nv, seed = 1600 + s)
      args <- list(runs = runs, pooling = pooling, method = "ar", p = 1L)
      if (pooling == "parcel") args$parcels <- rep(1:3, length.out = nv)
      clean <- do.call(fmriAR::fit_noise, c(list(resid = base), args))
      shifted <- do.call(fmriAR::fit_noise, c(list(resid = base + offs), args))
      expect_equal(plan_phi(shifted), plan_phi(clean), tolerance = 1e-6,
                   info = paste("boundary,", pooling, "seed", s))
    }
  }
})

test_that("parcel and global fits agree on average", {
  # They are different estimators -- global pools the autocovariance over all
  # voxels, parcel fits the parcel-mean series -- so they agree in expectation
  # rather than seed by seed.
  truth <- c(0.55, 0.2)
  g <- matrix(NA_real_, 6L, 2L)
  p1 <- matrix(NA_real_, 6L, 2L)
  for (s in 1:6) {
    resid <- ar_sim(400L, truth, nvox = 30L, seed = 1700 + s)
    g[s, ] <- plan_phi(fmriAR::fit_noise(resid, pooling = "global",
                                         method = "ar", p = 2L))
    p1[s, ] <- plan_phi(fmriAR::fit_noise(resid, parcels = rep(1L, 30L),
                                          pooling = "parcel", method = "ar", p = 2L))
  }
  expect_equal(colMeans(g), truth, tolerance = 0.08)
  expect_equal(colMeans(p1), truth, tolerance = 0.08)
  expect_equal(colMeans(p1), colMeans(g), tolerance = 0.08)
})


# --- Internal estimator invariants -------------------------------------------

test_that("the PSD projection actually repairs non-PSD autocovariance", {
  # Feed .acvf_from_pooled() sequences that are provably NOT positive
  # semi-definite, so the projection has to do work. An earlier version of this
  # test only used well-behaved data, where the projection never ran at all and
  # deleting it left the test passing.
  psd_min_eig <- function(g) {
    min(eigen(stats::toeplitz(g), symmetric = TRUE, only.values = TRUE)$values)
  }
  bad <- list(
    rho_gt_1   = list(num = c(1, 1.4), pairs = c(10, 9)),
    alternating = list(num = c(1, -0.99, 0.98, -0.97), pairs = c(50, 40, 30, 5)),
    long_tail  = list(num = c(1, 0.2, 0.1, 0.9), pairs = c(100, 90, 80, 3)),
    spike      = list(num = c(1, 0.1, 1.3, 0.1), pairs = c(60, 50, 4, 30))
  )
  ran <- 0L
  for (nm in names(bad)) {
    p <- bad[[nm]]
    g_unb <- p$num / p$pairs
    if (psd_min_eig(g_unb) >= 0) next          # only count genuinely bad inputs
    ran <- ran + 1L
    g <- fmriAR:::.acvf_from_pooled(p)
    expect_gte(psd_min_eig(g), -1e-8 * max(1, abs(g[1])), label = nm)
    expect_true(all(is.finite(g)), label = nm)
  }
  expect_gt(ran, 0L)   # the projection must have been exercised
})

test_that("Yule-Walker on the projected autocovariance stays stationary", {
  for (r1 in c(0.9, 0.999, 1.05, 1.5)) {
    g <- fmriAR:::.acvf_from_pooled(list(num = c(1, r1), pairs = c(100, 90)))
    yw <- yw_from_acvf_fast(g[1:2], 1L)
    expect_lt(abs(yw$phi[1]), 1, label = paste("rho1 =", r1))
  }
})

test_that("multiscale autocovariance covers the pooling target", {
  # Regression guard: the acvf was sized to the order selected at this scale,
  # so .ms_pad() zero-filled the remaining lags and Yule-Walker on the padded
  # vector produced coefficients pinned at the stationarity boundary.
  M <- ar_sim(300L, 0.5, nvox = 3L, seed = 2901)
  colnames(M) <- as.character(seq_len(ncol(M)))
  est1 <- fmriAR:::.ms_estimate_scale(
    M, function(y) list(phi = 0.5, order = c(p = 1L, q = 0L)),
    run_starts0 = 0L, lag_max = 5L)
  for (g in est1$acvf) {
    expect_gte(length(g), 6L)            # lags 0..5, never truncated to p + 1
    expect_true(all(is.finite(g)))
    expect_gt(g[1], 0)
    expect_false(all(g[-1] == 0))        # must not be zero-filled
  }
})

test_that("parcel censoring is handled at every multiscale entry point", {
  # Regression guard: .ms_estimate_scale built its acvf with a helper that
  # centres each segment separately, so censoring re-introduced the per-fragment
  # centering defect on the p_target and acvf_pooled paths only, driving phi
  # from +0.7 to -0.07. And parcel_sets combined with censor crashed outright.
  resid <- ar_sim(300L, 0.7, nvox = 36L, seed = 3001)
  cens <- which(rep(c(FALSE, FALSE, TRUE), length.out = 300L))
  pf <- rep(1:6, each = 6); pm <- rep(1:3, each = 12); pc <- rep(1:2, each = 18)

  variants <- list(
    plain    = list(),
    p_target = list(p_target = 1L),
    acvf     = list(multiscale = TRUE, ms_mode = "acvf_pooled", p_target = 1L),
    pacf     = list(multiscale = TRUE, ms_mode = "pacf_weighted", p_target = 1L),
    sets     = list(parcel_sets = list(coarse = pc, medium = pm, fine = pf),
                    multiscale = TRUE, ms_mode = "acvf_pooled", p_target = 1L)
  )
  for (nm in names(variants)) {
    args <- c(list(resid = resid, parcels = pf, pooling = "parcel",
                   method = "ar", p = 1L, censor = cens), variants[[nm]])
    plan <- do.call(fmriAR::fit_noise, args)
    est <- mean(vapply(plan$phi_by_parcel, function(z) z[1], numeric(1)))
    expect_equal(est, 0.7, tolerance = 0.15, info = nm)
  }
})

test_that("short fragments of a run do not force a correlation of -1", {
  # Heavy scrubbing leaves mostly 2- and 3-frame fragments. Centering each on
  # its own mean makes every lag-1 product negative by construction -- a 2-frame
  # fragment gives -(x1 - x2)^2 / 4 whatever the data -- so pooling them drove
  # rho1 to -1. The mean belongs to the run, not the fragment.
  y <- as.numeric(ar_sim(400L, 0.7, nvox = 1L, seed = 2101))
  keep <- which(rep(c(TRUE, TRUE, FALSE), length.out = 400L))  # 2-frame fragments
  seg <- cumsum(c(1L, as.integer(diff(keep) > 1L)))
  pooled <- fmriAR:::.pooled_acvf_segments(matrix(y[keep], ncol = 1L), seg, 1L,
                                           center_id = rep(1L, length(keep)))
  g <- fmriAR:::.acvf_from_pooled(pooled)
  expect_gt(g[2] / g[1], 0.4)   # true rho1 is 0.7; must not collapse to -1
})

test_that("lag products are confined to their segment", {
  # Two segments of a constant each: within-segment products are zero after
  # centering, so a boundary-spanning product is the only way to get signal.
  m <- matrix(c(rep(0, 5L), rep(10, 5L)), ncol = 1L)
  seg <- c(rep(1L, 5L), rep(2L, 5L))
  pooled <- fmriAR:::.pooled_acvf_segments(m, seg, 2L, center_id = seg)
  expect_equal(pooled$num[2], 0, tolerance = 1e-12)
  expect_equal(pooled$pairs[2], 8)   # 4 within each of the two segments
})

test_that("valid segments break at both run changes and censored frames", {
  seg <- fmriAR:::.valid_segments(10L, runs = c(1,1,1,1,1,2,2,2,2,2),
                                  censor = c(3L, 8L))
  expect_equal(seg$idx, c(1L, 2L, 4L, 5L, 6L, 7L, 9L, 10L))
  # segment starts: index 1, after the censored 3, at the run change, after 8
  expect_equal(seg$starts0, c(0L, 2L, 4L, 6L))
  expect_equal(seg$run_id, c(1, 1, 1, 1, 2, 2, 2, 2))
})

# --- whiten_apply argument safety --------------------------------------------

test_that("whiten_apply does not mutate the caller's X or Y", {
  # Regression guard: the parcel branch handed the caller's X straight to an
  # in-place C++ routine, so calling whiten_apply destroyed the design matrix.
  resid <- ar_sim(200L, 0.6, nvox = 20L, seed = 2401)
  parcels <- rep(1:4, length.out = 20L)
  X <- cbind(1, as.numeric(seq_len(200L) > 100L), stats::rnorm(200L))
  X0 <- X + 0
  Y0 <- resid + 0

  for (pooling in c("global", "parcel")) {
    args <- list(resid = resid, method = "ar", p = 1L, pooling = pooling)
    if (pooling == "parcel") args$parcels <- parcels
    plan <- do.call(fmriAR::fit_noise, args)
    fmriAR::whiten_apply(plan, X, resid,
                         parcels = if (pooling == "parcel") parcels else NULL)
    expect_equal(X, X0, info = paste("X preserved,", pooling))
    expect_equal(resid, Y0, info = paste("Y preserved,", pooling))
  }
})

test_that("each parcel's design is whitened with that parcel's coefficients", {
  # Regression guard: all X_by entries aliased a single matrix that had been
  # filtered once per parcel in sequence, so none matched any parcel's filter.
  resid <- ar_sim(200L, 0.6, nvox = 20L, seed = 2501)
  parcels <- rep(1:4, length.out = 20L)
  X <- cbind(1, stats::rnorm(200L))
  X0 <- X + 0
  plan <- fmriAR::fit_noise(resid, parcels = parcels, pooling = "parcel",
                            method = "ar", p = 1L)
  out <- fmriAR::whiten_apply(plan, X, resid, parcels = parcels)

  phis <- vapply(plan$phi_by_parcel, function(z) z[1], numeric(1))
  expect_gt(stats::sd(phis), 1e-4)   # parcels must genuinely differ

  n <- nrow(X0)
  for (k in seq_along(plan$parcel_ids)) {
    phi <- phis[[k]]
    manual <- X0
    manual[2:n, ] <- X0[2:n, , drop = FALSE] - phi * X0[1:(n - 1), , drop = FALSE]
    manual[1, ] <- X0[1, ] * sqrt(1 - phi^2)
    expect_equal(out$X_by[[k]], manual, tolerance = 1e-8,
                 info = paste("parcel", k))
  }
})

test_that("whiten_apply returns rows in input order for any run labelling", {
  # Regression guard: results were reassembled with rbind() over split(), which
  # orders by sorted run label, so a run labelling that is not ascending in time
  # silently returned a row permutation of the correct answer.
  resid <- ar_sim(120L, 0.6, nvox = 6L, seed = 2601)
  X <- cbind(1, stats::rnorm(120L))
  plan <- fmriAR::fit_noise(resid, pooling = "global", method = "ar", p = 1L)

  ref <- fmriAR::whiten_apply(plan, X, resid, run_starts = c(0L, 60L))
  for (lab in list(c(2L, 1L), c(5L, 3L), c(10L, 1L))) {
    out <- fmriAR::whiten_apply(plan, X, resid, runs = rep(lab, each = 60L))
    expect_equal(out$Y, ref$Y, tolerance = 1e-12,
                 info = paste("labels", paste(lab, collapse = "/")))
    expect_equal(out$X, ref$X, tolerance = 1e-12)
  }
})

test_that("enforce_stationary_ar returns roots strictly outside the unit circle", {
  # Clamping reflection coefficients bounds each |kappa| < 1 but at high order
  # leaves the roots on the unit circle: PACF 0.99 repeated ten times gave
  # coefficients above 45 with min|root| equal to 1 to machine precision.
  # `> 1` alone does not bite: clamping the reflection coefficients leaves
  # min|root| at 1.0000000000, which still compares greater than 1 in floating
  # point. Assert the margin the function actually promises.
  margin <- 1e-6
  cases <- list(rep(0.99, 10), rep(0.99, 6), rep(0.9, 3), c(0.5, 0.3),
                c(1.9, -0.99), rep(0.999, 4), rep(0.95, 8))
  checked <- 0L
  for (v in cases) {
    out <- fmriAR:::enforce_stationary_ar(v, 0.99)
    if (!length(out)) next          # refusing to produce a filter is acceptable
    checked <- checked + 1L
    expect_true(all(is.finite(out)))
    expect_gt(min(Mod(polyroot(c(1, -out)))), 1 + margin / 2)
  }
  expect_gt(checked, 0L)

  # The clamp alone is not sufficient, which is what makes the shrink load-bearing.
  clamp_only <- fmriAR:::pacf_to_ar(pmin(pmax(fmriAR:::ar_to_pacf(rep(0.99, 6)), -0.99), 0.99))
  expect_lt(min(Mod(polyroot(c(1, -clamp_only)))), 1 + margin / 2)
})

test_that("degenerate order requests fail cleanly rather than opaquely", {
  # p_max at or near n drove Levinson to a NaN prediction error, surfacing as
  # "missing value where TRUE/FALSE needed" from an internal comparison.
  set.seed(2701)
  for (n in c(20L, 40L, 60L)) {
    for (pm in c(n - 2L, n - 1L, n)) {
      resid <- matrix(stats::rnorm(n * 4L), n, 4L)
      plan <- fmriAR::fit_noise(resid, pooling = "global", method = "ar",
                                p = "auto", p_max = pm)
      expect_s3_class(plan, "fmriAR_plan")
      expect_gt(min_root(plan_phi(plan)), 1)
    }
  }
})

test_that("order selection is bounded by sample size but an explicit p is not", {
  set.seed(2801)
  short <- matrix(stats::rnorm(11L * 5L), 11L, 5L)   # white noise, 11 frames
  auto <- fmriAR::fit_noise(short, pooling = "global", method = "ar",
                            p = "auto", p_max = 8L)
  expect_lte(auto$order[["p"]], 2L)

  X <- cbind(1, stats::rnorm(11L))
  w <- fmriAR::whiten_apply(auto, X, short)
  expect_lt(stats::sd(w$Y) / stats::sd(short), 1.05)   # must not amplify

  fixed <- fmriAR::fit_noise(matrix(stats::rnorm(25L * 5L), 25L, 5L),
                             pooling = "global", method = "ar", p = 4L)
  expect_equal(length(plan_phi(fixed)), 4L)
})

test_that("whiten_apply rejects malformed runs instead of corrupting output", {
  # Regression guard: a short `runs` was silently recycled by split(), and an NA
  # left that row unwritten -- NA in both X and Y, with no error.
  resid <- ar_sim(20L, 0.5, nvox = 3L, seed = 3201)
  X <- cbind(1, stats::rnorm(20L))
  plan <- fmriAR::fit_noise(resid, pooling = "global", method = "ar", p = 1L)

  bad_na <- rep(1:2, each = 10L); bad_na[5] <- NA
  expect_error(fmriAR::whiten_apply(plan, X, resid, runs = bad_na), "NA")
  expect_error(fmriAR::whiten_apply(plan, X, resid, runs = rep(1L, 5L)), "length")

  # Character run labels are accepted and treated as grouping labels.
  out <- fmriAR::whiten_apply(plan, X, resid,
                              runs = rep(c("a", "b"), each = 10L))
  expect_false(anyNA(out$Y))
  expect_equal(dim(out$Y), dim(resid))
})

test_that("compat plans built without an MA part are usable", {
  # Regression guard: theta = NULL, the documented default, produced
  # order q = -Inf (with a max()-over-nothing warning) and an empty theta list,
  # so whiten_apply() failed with "subscript out of bounds". The example in
  # ?compat builds exactly this plan.
  X <- cbind(1, stats::rnorm(60L))
  Y <- ar_sim(60L, 0.5, nvox = 4L, seed = 3101)

  for (pooling in c("global", "run")) {
    plan <- expect_silent(
      fmriAR:::compat$plan_from_phi(c(0.4, 0.2), pooling = pooling,
                                    runs = if (pooling == "run")
                                      rep(1:2, each = 30L) else NULL)
    )
    expect_equal(plan$order[["p"]], 2L)
    expect_equal(plan$order[["q"]], 0L)
    out <- fmriAR::whiten_apply(plan, X, Y,
                                runs = if (pooling == "run")
                                  rep(1:2, each = 30L) else NULL)
    expect_equal(dim(out$Y), dim(Y))
    expect_false(identical(out$Y, Y))
  }

  # An explicit MA part still reports the right order.
  plan_ma <- fmriAR:::compat$plan_from_phi(c(0.4), theta = c(0.3), method = "arma")
  expect_equal(plan_ma$order[["q"]], 1L)
})

test_that("acorr_diagnostics honours its runs argument", {
  # Regression guard: `runs` was documented and accepted but never referenced,
  # so the result was identical with and without it.
  set.seed(2201)
  n <- 400L; nv <- 20L
  runs <- rep(1:4, each = 100L)
  white <- matrix(stats::rnorm(n * nv), n, nv)
  shifted <- white + rep(stats::rnorm(4L, sd = 10), each = 100L)

  ignored <- fmriAR::acorr_diagnostics(shifted, max_lag = 5L)
  aware <- fmriAR::acorr_diagnostics(shifted, runs = runs, max_lag = 5L)

  expect_false(isTRUE(all.equal(ignored$acf, aware$acf)))
  # Between-run offsets fake strong autocorrelation; run-aware sees white noise.
  expect_gt(ignored$acf[1], 0.5)
  expect_lt(abs(aware$acf[1]), 0.15)
})

test_that("acorr_diagnostics recovers a known autocorrelation", {
  resid <- ar_sim(500L, 0.6, nvox = 20L, seed = 2301)
  out <- fmriAR::acorr_diagnostics(resid, max_lag = 3L)
  expect_equal(out$acf[1], 0.6, tolerance = 0.12)
  runs <- rep(1:2, each = 250L)
  out_r <- fmriAR::acorr_diagnostics(resid, runs = runs, max_lag = 3L)
  expect_equal(out_r$acf[1], 0.6, tolerance = 0.12)
})

test_that("stationarity is enforced on every returned plan", {
  # Sweep the awkward regimes together: near-unit-root data, short series with a
  # high requested order, heavy censoring, and multiple runs. Every returned
  # plan must be stationary, on every pooling mode.
  set.seed(1901)
  for (s in 1:6) {
    n <- sample(60:140, 1L)
    cases <- list(
      walk  = matrix(cumsum(stats::rnorm(n * 6L)), n, 6L),
      ar95  = ar_sim(n, 0.95, nvox = 6L, seed = 1900 + s),
      short = matrix(stats::rnorm(n * 6L), n, 6L)
    )
    for (nm in names(cases)) {
      resid <- cases[[nm]]
      runs <- rep(1:2, length.out = n)
      runs <- sort(runs)
      cens <- sort(sample.int(n, floor(0.3 * n)))
      for (pooling in c("global", "run", "parcel")) {
        args <- list(resid = resid, method = "ar", p = "auto", p_max = 6L,
                     pooling = pooling, runs = runs, censor = cens)
        if (pooling == "parcel") args$parcels <- rep(1:2, length.out = 6L)
        plan <- try(do.call(fmriAR::fit_noise, args), silent = TRUE)
        if (inherits(plan, "try-error")) next
        if (!is.null(plan$phi_by_parcel)) {
          for (ph in plan$phi_by_parcel) {
            expect_gt(min_root(ph), 1,
                      label = sprintf("%s/%s seed %d", nm, pooling, s))
          }
        } else {
          expect_gt(min_root(plan_phi(plan)), 1,
                    label = sprintf("%s/%s seed %d", nm, pooling, s))
        }
      }
    }
  }
})

test_that("ARMA warns rather than silently biasing across censoring gaps", {
  # Hannan-Rissanen runs on the surviving frames spliced together, so its
  # regressions span censoring gaps and bias the coefficients. Documented and
  # warned about rather than left silent.
  resid <- sapply(1:10, function(i) as.numeric(stats::arima.sim(list(ar = 0.5, ma = 0.4), 300L)))
  cens <- seq(10L, 290L, by = 10L)

  expect_warning(
    fmriAR::fit_noise(resid, method = "arma", p = 1L, q = 1L,
                      pooling = "global", censor = cens),
    "censoring gaps"
  )
  # No censoring, no warning.
  expect_silent(
    fmriAR::fit_noise(resid, method = "arma", p = 1L, q = 1L, pooling = "global")
  )
})


# --- Noise scale on the plan --------------------------------------------------

test_that("the plan carries the noise scale, not just its shape", {
  # Regression guard: the plan recorded only correlation structure, so two
  # datasets differing 100-fold in variance produced bit-identical plans and
  # the magnitude of the noise was unrecoverable.
  A <- ar_sim(400L, 0.5, nvox = 20L, seed = 3301)
  B <- A * 10

  pa <- fmriAR::fit_noise(A, pooling = "global", method = "ar", p = 1L)
  pb <- fmriAR::fit_noise(B, pooling = "global", method = "ar", p = 1L)

  expect_false(isTRUE(all.equal(unclass(pa), unclass(pb))))
  expect_equal(plan_phi(pa), plan_phi(pb), tolerance = 1e-8)   # same shape
  expect_equal(pb$sigma2[[1]] / pa$sigma2[[1]], 100, tolerance = 1e-6)
  expect_equal(pb$gamma[[1]][1] / pa$gamma[[1]][1], 100, tolerance = 1e-6)
})

test_that("sigma2 and gamma agree with the data they were fitted from", {
  phi <- 0.6
  A <- ar_sim(4000L, phi, nvox = 30L, seed = 3401)
  plan <- fmriAR::fit_noise(A, pooling = "global", method = "ar", p = 1L)

  # gamma[1] is the marginal variance; sigma2 the innovation variance.
  expect_equal(plan$gamma[[1]][1], mean(apply(A, 2L, stats::var)), tolerance = 0.05)
  expect_equal(plan$sigma2[[1]], plan$gamma[[1]][1] * (1 - phi^2), tolerance = 0.05)
  # gamma[h] / gamma[0] should follow phi^h.
  expect_equal(plan$gamma[[1]][2] / plan$gamma[[1]][1], phi, tolerance = 0.05)
  expect_equal(plan$gamma[[1]][3] / plan$gamma[[1]][1], phi^2, tolerance = 0.05)
})

test_that("every pooling mode reports a scale keyed like its coefficients", {
  A <- ar_sim(300L, 0.5, nvox = 20L, seed = 3501)
  parcels <- rep(1:4, length.out = 20L)

  g <- fmriAR::fit_noise(A, pooling = "global", method = "ar", p = 1L)
  expect_length(g$gamma, 1L)
  expect_length(g$sigma2, 1L)
  expect_true(is.finite(g$sigma2[[1]]) && g$sigma2[[1]] > 0)

  r <- fmriAR::fit_noise(A, runs = rep(1:3, each = 100L), pooling = "run",
                         method = "ar", p = 1L)
  expect_length(r$gamma, 3L)
  expect_length(r$sigma2, 3L)
  expect_true(all(vapply(r$sigma2, function(z) is.finite(z) && z > 0, logical(1))))

  p <- fmriAR::fit_noise(A, parcels = parcels, pooling = "parcel",
                         method = "ar", p = 1L)
  expect_equal(names(p$sigma2_by_parcel), names(p$phi_by_parcel))
  expect_equal(names(p$gamma_by_parcel), names(p$phi_by_parcel))
  expect_true(all(vapply(p$sigma2_by_parcel,
                         function(z) is.finite(z) && z > 0, logical(1))))
})

test_that("the reported gamma reconstructs a valid noise covariance", {
  A <- ar_sim(600L, c(0.6, 0.2), nvox = 20L, seed = 3601)
  plan <- fmriAR::fit_noise(A, pooling = "global", method = "ar", p = 3L)
  g <- plan$gamma[[1]]

  Sigma <- stats::toeplitz(g)
  expect_gt(min(eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values), 0)
  expect_equal(nrow(Sigma), length(g))
})
