# Accuracy, calibration and timing benchmark for fmriAR.
#
# Uses only the exported API so the same script runs against any installed
# version, which is how a change is checked for accuracy regressions:
#
#   Rscript tools/validation/accuracy_benchmark.R <lib_path_or_-> <out.rds>
#
# Pass "-" to use the default library. Each scenario has a fixed seed, so two
# versions see identical data.

args <- commandArgs(trailingOnly = TRUE)
lib <- if (length(args) >= 1L && args[1] != "-") args[1] else NULL
out_file <- if (length(args) >= 2L) args[2] else NULL
suppressPackageStartupMessages(library(fmriAR, lib.loc = lib))
`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x

sim_arma <- function(n, v, phi = numeric(0), theta = numeric(0), burn = 200L) {
  vapply(seq_len(v), function(j) {
    as.numeric(stats::arima.sim(list(ar = phi, ma = theta), n, n.start = burn))
  }, numeric(n))
}

# Mean absolute per-voxel autocorrelation over lags 1..K, computed within runs.
whiteness <- function(R, runs, K = 3L) {
  vals <- vapply(seq_len(ncol(R)), function(j) {
    acc <- numeric(K); den <- 0
    for (r in unique(runs)) {
      y <- R[runs == r, j]; y <- y - mean(y)
      den <- den + sum(y * y)
      for (k in seq_len(K)) acc[k] <- acc[k] + sum(y[-seq_len(k)] * y[seq_len(length(y) - k)])
    }
    mean(abs(acc / den))
  }, 0)
  mean(vals)
}

res <- list()
add <- function(scenario, metric, value) {
  res[[length(res) + 1L]] <<- data.frame(scenario = scenario, metric = metric,
                                         value = as.numeric(value))
}

# S1: AR(1) order selection and estimation, global pooling --------------------
set.seed(101)
est <- t(replicate(100, {
  runs <- rep(1:3, each = 150)
  R <- sim_arma(450, 30, 0.5)
  pl <- fit_noise(R, runs = runs, p = "auto", p_max = 6)
  c(phi1 = pl$phi[[1]][1] %||% 0, p = length(pl$phi[[1]]))
}))
add("S1 AR(1) auto, global", "phi1 bias", mean(est[, 1]) - 0.5)
add("S1 AR(1) auto, global", "phi1 RMSE", sqrt(mean((est[, 1] - 0.5)^2)))
add("S1 AR(1) auto, global", "P(order>1)", mean(est[, 2] > 1))

# S1b: single-voxel order selection (where the BIC penalty matters most)
set.seed(102)
pp <- replicate(300, length(fit_noise(sim_arma(200, 1, 0.5), p = "auto", p_max = 6)$phi[[1]]))
add("S1b AR(1) auto, 1 voxel", "P(order>1)", mean(pp > 1))
set.seed(103)
pp2 <- replicate(300, length(fit_noise(sim_arma(200, 1, c(0.5, 0.25)), p = "auto", p_max = 6)$phi[[1]]))
add("S1b AR(2) auto, 1 voxel", "P(order>=2) [power]", mean(pp2 >= 2))

# S2: AR(2) explicit order, unequal runs, global pooling ----------------------
set.seed(201)
truth <- c(0.5, 0.2)
est <- t(replicate(100, {
  lens <- c(120, 200, 80)
  runs <- rep(seq_along(lens), lens)
  R <- do.call(rbind, lapply(lens, function(L) sim_arma(L, 30, truth)))
  fit_noise(R, runs = runs, p = 2)$phi[[1]]
}))
add("S2 AR(2) p=2, unequal runs", "phi bias (L1)", sum(abs(colMeans(est) - truth)))
add("S2 AR(2) p=2, unequal runs", "phi RMSE", sqrt(mean(sweep(est, 2, truth)^2)))

# S3: AR(1) with 15% censoring ------------------------------------------------
set.seed(301)
est <- replicate(100, {
  R <- sim_arma(300, 30, 0.6)
  cen <- sort(sample(300, 45))
  fit_noise(R, censor = cen, p = 1)$phi[[1]]
})
add("S3 AR(1) 15% censor", "phi bias", mean(est) - 0.6)
add("S3 AR(1) 15% censor", "phi RMSE", sqrt(mean((est - 0.6)^2)))

# S4: ARMA(1,1), independent voxels -------------------------------------------
arma_scn <- function(label, seed, censor_frac = 0, shared = 0) {
  set.seed(seed)
  est <- t(replicate(60, {
    R <- sim_arma(240, 30, 0.5, 0.4)
    if (shared > 0) R <- R + shared * as.numeric(stats::filter(rnorm(240), 0.95, method = "recursive"))
    cen <- if (censor_frac > 0) sort(sample(240, round(censor_frac * 240))) else NULL
    pl <- suppressWarnings(fit_noise(R, censor = cen, method = "arma", p = 1, q = 1))
    w <- whiten_apply(pl, matrix(1, 240, 1), R, censor = cen)$Y
    c(pl$phi[[1]], pl$theta[[1]], whiteness(w, rep(1, 240)))
  }))
  add(label, "phi bias", mean(est[, 1]) - 0.5)
  add(label, "theta bias", mean(est[, 2]) - 0.4)
  add(label, "phi RMSE", sqrt(mean((est[, 1] - 0.5)^2)))
  add(label, "theta RMSE", sqrt(mean((est[, 2] - 0.4)^2)))
  add(label, "whitened mean|acf1..3|", mean(est[, 3]))
}
arma_scn("S4 ARMA(1,1)", 401)
arma_scn("S5 ARMA(1,1) 20% censor", 501, censor_frac = 0.2)
arma_scn("S6 ARMA(1,1) + shared signal", 601, shared = 0.4)

# S6b: shared signal, orders chosen automatically. The noise is ARMA(1,1) plus
# a shared AR(0.95), i.e. not ARMA(1,1), and the single shared realisation
# keeps even the empirical voxel autocorrelation ~0.14 from the theoretical
# one, so the attainable target is in-sample: how white the voxels end up and
# how closely the fitted model tracks the voxels' own autocorrelation.
set.seed(602)
s6b <- t(replicate(30, {
  R <- sim_arma(240, 30, 0.5, 0.4) +
    0.4 * as.numeric(stats::filter(rnorm(240), 0.95, method = "recursive"))
  pl <- tryCatch(fit_noise(R, method = "arma", p = "auto", q = "auto"),
                 error = function(e) NULL)
  if (is.null(pl)) return(c(NA, NA, NA))
  ac_emp <- rowMeans(apply(R, 2, function(y) stats::acf(y, 10, plot = FALSE)$acf[-1]))
  r <- stats::ARMAacf(pl$phi[[1]], pl$theta[[1]], 10)[-1]
  w <- whiten_apply(pl, matrix(1, 240, 1), R)$Y
  c(max(abs(r - ac_emp)), whiteness(w, rep(1, 240)),
    length(pl$phi[[1]]) + length(pl$theta[[1]]))
}))
add("S6b ARMA auto orders + shared signal", "max|acf_fit - acf_voxels| lags 1-10", mean(s6b[, 1]))
add("S6b ARMA auto orders + shared signal", "whitened mean|acf1..3|", mean(s6b[, 2]))
add("S6b ARMA auto orders + shared signal", "p+q (mean)", mean(s6b[, 3]))

# S6c: ARMA order recovery (truth (1,1), (0,0), (1,0)).
for (cfg in list(list("S6c ARMA auto, truth (1,1)", 0.5, 0.4, c(1, 1)),
                 list("S6c ARMA auto, truth (0,0)", numeric(0), numeric(0), c(0, 0)),
                 list("S6c ARMA auto, truth (1,0)", 0.5, numeric(0), c(1, 0)))) {
  set.seed(603)
  ok <- replicate(30, {
    R <- sim_arma(240, 30, cfg[[2]], cfg[[3]])
    pl <- tryCatch(fit_noise(R, method = "arma", p = "auto", q = "auto"),
                   error = function(e) NULL)
    if (is.null(pl)) return(NA)
    all(c(length(pl$phi[[1]]), length(pl$theta[[1]])) == cfg[[4]])
  })
  add(cfg[[1]], "P(correct order)", mean(ok))
}

# S7: GLS efficiency and calibration with short runs --------------------------
set.seed(701)
nr <- 4; L <- 60; n <- nr * L
runs <- rep(seq_len(nr), each = L)
task <- rep(rep(c(0, 1), each = 10), length.out = L)
X <- cbind(task = rep(task, nr), model.matrix(~ factor(runs) - 1))
tt <- bb <- NULL
for (rep_i in 1:150) {
  R <- do.call(rbind, lapply(seq_len(nr), function(i) sim_arma(L, 20, c(0.6, 0.2))))
  res0 <- R - X %*% qr.solve(X, R)
  pl <- fit_noise(res0, runs = runs, p = 2)
  w <- whiten_apply(pl, X, R, runs = runs)
  s <- sandwich_from_whitened_resid(w$X, w$Y)
  b <- qr.solve(w$X, w$Y)[1, ]
  tt <- c(tt, b / s$se[1, ]); bb <- c(bb, b)
}
add("S7 GLS AR(2) 4x60 null", "type-I rate (|t|>1.96)", mean(abs(tt) > 1.96))
add("S7 GLS AR(2) 4x60 null", "var(beta_hat) x1000", 1000 * var(bb))

# S8: whiteness at segment starts, known AR(2) -------------------------------
set.seed(801)
Y <- sim_arma(40, 4000, c(1.2, -0.5))
pl <- compat$plan_from_phi(c(1.2, -0.5), exact_first = TRUE)
w <- whiten_apply(pl, matrix(1, 40, 1), Y)$Y
v <- apply(w, 1, var)
add("S8 AR(2) start variance", "var t=1", v[1])
add("S8 AR(2) start variance", "var t=2", v[2])
add("S8 AR(2) start variance", "var t=3..40 (mean)", mean(v[3:40]))

# S9: multiscale parcel pooling ------------------------------------------------
set.seed(901)
ps <- list(coarse = rep(1:2, each = 64), medium = rep(1:4, each = 32), fine = rep(1:16, each = 8))
phi_true <- rep(seq(0.2, 0.6, length.out = 16), each = 8)
err <- replicate(20, {
  R <- vapply(seq_len(128), function(j) as.numeric(stats::arima.sim(list(ar = phi_true[j]), 200)), numeric(200))
  R <- R * rep(runif(128, 0.5, 2), each = 200)
  pl <- fit_noise(R, pooling = "parcel", parcels = ps$fine, parcel_sets = ps, p = 1,
                  multiscale = "pacf_weighted")
  est <- vapply(pl$phi_by_parcel, function(x) x[1] %||% 0, 0)
  sqrt(mean((est - unique(phi_true))^2))
})
add("S9 multiscale parcels", "phi RMSE", mean(err))

# S9b: multiscale acvf_pooled with automatic order (true order 1)
set.seed(902)
ms_b <- t(replicate(20, {
  R <- vapply(seq_len(128), function(j) as.numeric(stats::arima.sim(list(ar = phi_true[j]), 200)), numeric(200))
  pl <- fit_noise(R, pooling = "parcel", parcels = ps$fine, parcel_sets = ps, p = "auto",
                  p_max = 6, multiscale = "acvf_pooled")
  est <- vapply(pl$phi_by_parcel, function(x) x[1] %||% 0, 0)
  w <- whiten_apply(pl, matrix(1, 200, 1), R, parcels = ps$fine)$Y
  c(sqrt(mean((est - unique(phi_true))^2)), pl$order[["p"]], whiteness(w, rep(1, 200)))
}))
add("S9b multiscale acvf_pooled auto", "phi1 RMSE", mean(ms_b[, 1]))
add("S9b multiscale acvf_pooled auto", "fitted order (mean)", mean(ms_b[, 2]))
add("S9b multiscale acvf_pooled auto", "whitened mean|acf1..3|", mean(ms_b[, 3]))

# S13: plain parcel pooling (no multiscale), with and without a slow
# component shared by the voxels of each parcel. The target is the voxels'
# own autocorrelation: phi1 against the mean per-voxel lag-1 autocorrelation,
# and residual autocorrelation after whitening.
parcel_scn <- function(label, seed, shared) {
  set.seed(seed)
  r <- t(replicate(20, {
    pf <- rep(1:8, each = 16)
    phi_p <- seq(0.2, 0.6, length.out = 8)
    R <- vapply(seq_along(pf), function(j)
      as.numeric(stats::arima.sim(list(ar = phi_p[pf[j]]), 200)), numeric(200))
    if (shared > 0) {
      for (k in 1:8) {
        g <- as.numeric(stats::filter(rnorm(200), 0.9, method = "recursive"))
        R[, pf == k] <- R[, pf == k] + shared * g
      }
    }
    pl <- fit_noise(R, pooling = "parcel", parcels = pf, p = 1)
    est <- vapply(pl$phi_by_parcel, function(x) x[1] %||% 0, 0)
    vox_rho <- vapply(1:8, function(k) mean(apply(R[, pf == k], 2, function(y)
      stats::acf(y, 1, plot = FALSE)$acf[2])), 0)
    w <- whiten_apply(pl, matrix(1, 200, 1), R, parcels = pf)$Y
    c(sqrt(mean((est - vox_rho)^2)), sqrt(mean((est - phi_p)^2)), whiteness(w, rep(1, 200)))
  }))
  add(label, "RMSE(phi1 - voxel lag-1 acf)", mean(r[, 1]))
  add(label, "RMSE(phi1 - generating phi)", mean(r[, 2]))
  add(label, "whitened mean|acf1..3|", mean(r[, 3]))
}
parcel_scn("S13a parcel p=1", 1301, shared = 0)
parcel_scn("S13b parcel p=1 + shared slow", 1302, shared = 0.5)

# S14: residual-bias correction with drift regressors, with/without censoring
dct_basis <- function(L, k) sapply(1:k, function(j) cos(pi * j * (seq_len(L) - 0.5) / L))
for (cf in c(0, 0.1)) {
  set.seed(1400 + 10 * cf)
  n <- 300; runs <- rep(1:2, each = 150)
  Xd <- cbind(model.matrix(~ factor(runs) - 1),
              rbind(cbind(dct_basis(150, 6), matrix(0, 150, 6)),
                    cbind(matrix(0, 150, 6), dct_basis(150, 6))),
              rep(rep(c(0, 1), each = 10), length.out = n))
  est <- t(replicate(30, {
    E <- sim_arma(n, 30, 0.4)
    R <- E - Xd %*% qr.solve(Xd, E)
    cen <- if (cf > 0) sort(sample(n, cf * n)) else NULL
    g <- suppressWarnings(fit_noise(R, runs = runs, p = 1, censor = cen, design = Xd))$phi[[1]]
    pc <- tryCatch(suppressWarnings(fit_noise(R, runs = runs, p = 1, censor = cen, design = Xd,
                                              pooling = "parcel", parcels = rep(1:3, each = 10))),
                   error = function(e) NULL)
    c(g, if (is.null(pc)) NA else mean(unlist(pc$phi_by_parcel)))
  }))
  lab <- sprintf("S14 design correction, %d%% censor", round(100 * cf))
  add(lab, "global phi RMSE (truth 0.4)", sqrt(mean((est[, 1] - 0.4)^2)))
  add(lab, "parcel phi RMSE (truth 0.4)", sqrt(mean((est[, 2] - 0.4)^2)))
}
add("timing", "acvf_bias_matrix n=800, 25 lags [s]", {
  Xb <- cbind(1, poly(seq_len(800), 3), matrix(rnorm(800 * 6), 800))
  min(vapply(1:3, function(i) system.time(acvf_bias_matrix(Xb, runs = rep(1:2, each = 400), max_lag = 25))[["elapsed"]], 0))
})

# S11: global ARMA, two runs, different run lengths
set.seed(1101)
est <- t(replicate(40, {
  runs <- rep(1:2, c(150, 250))
  R <- rbind(sim_arma(150, 20, 0.5, 0.4), sim_arma(250, 20, 0.5, 0.4))
  pl <- fit_noise(R, runs = runs, method = "arma", p = 1, q = 1)
  c(pl$phi[[1]], pl$theta[[1]])
}))
add("S11 ARMA(1,1) global, 2 runs", "phi RMSE", sqrt(mean((est[, 1] - 0.5)^2)))
add("S11 ARMA(1,1) global, 2 runs", "theta RMSE", sqrt(mean((est[, 2] - 0.4)^2)))

# S12: GLS with ARMA(1,1) noise and short runs (exact ARMA start-up)
set.seed(1201)
nr <- 6; L <- 40; n <- nr * L
runs <- rep(seq_len(nr), each = L)
task <- rep(rep(c(0, 1), each = 8), length.out = L)
X <- cbind(task = rep(task, nr), model.matrix(~ factor(runs) - 1))
tt <- bb <- NULL
for (rep_i in 1:150) {
  R <- do.call(rbind, lapply(seq_len(nr), function(i) sim_arma(L, 20, 0.7, 0.5)))
  res0 <- R - X %*% qr.solve(X, R)
  pl <- fit_noise(res0, runs = runs, method = "arma", p = 1, q = 1)
  w <- whiten_apply(pl, X, R, runs = runs)
  s <- sandwich_from_whitened_resid(w$X, w$Y)
  b <- qr.solve(w$X, w$Y)[1, ]
  tt <- c(tt, b / s$se[1, ]); bb <- c(bb, b)
}
add("S12 GLS ARMA(1,1) 6x40 null", "type-I rate (|t|>1.96)", mean(abs(tt) > 1.96))
add("S12 GLS ARMA(1,1) 6x40 null", "var(beta_hat) x1000", 1000 * var(bb))

# S10: whitened-residual whiteness, AR auto, run pooling -----------------------
set.seed(1001)
wv <- replicate(30, {
  runs <- rep(1:2, each = 150)
  R <- sim_arma(300, 30, c(0.4, 0.2), 0.3)
  pl <- fit_noise(R, runs = runs, pooling = "run", p = "auto")
  whiteness(whiten_apply(pl, matrix(1, 300, 1), R, runs = runs)$Y, runs)
})
add("S10 AR auto on ARMA noise", "whitened mean|acf1..3|", mean(wv))

# Timing -----------------------------------------------------------------------
set.seed(1)
n <- 400; V <- 20000
Xw <- matrix(rnorm(n * 10), n, 10)
Yw <- matrix(rnorm(n * V), n, V)
pl <- compat$plan_from_phi(c(0.4, 0.2), exact_first = TRUE)
tm <- function(expr) {
  e <- substitute(expr); env <- parent.frame()
  min(vapply(1:3, function(i) system.time(eval(e, env))[["elapsed"]], 0))
}
add("timing", "whiten_apply 400x20000 AR(2) [s]", tm(whiten_apply(pl, Xw, Yw)))
add("timing", "fit_noise 400x20000 auto [s]", tm(fit_noise(Yw, p = "auto")))
add("timing", "sandwich hc0 400x20000 [s]", tm(sandwich_from_whitened_resid(Xw, Yw, type = "hc0")))
add("timing", "fit_noise ARMA(1,1) 400x20000 [s]", tm(fit_noise(Yw, method = "arma", p = 1, q = 1)))

out <- do.call(rbind, res)
print(out, row.names = FALSE, digits = 4)
if (!is.null(out_file)) saveRDS(out, out_file)
