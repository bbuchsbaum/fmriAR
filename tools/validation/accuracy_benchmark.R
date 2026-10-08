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
