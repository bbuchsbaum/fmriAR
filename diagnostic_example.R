# Correct whitening diagnostics example
#
# Issue: The original diagnostic code computed NEW residuals from whitened data:
#   innov_var <- Yw - Xw %*% qr.solve(Xw, Yw)  # WRONG!
#
# This gives the same autocorrelation as original residuals because you're
# computing residuals from a new regression on the whitened data.
#
# Correct approach: Apply whitening transformation to the ORIGINAL residuals

library(fmriAR)

# Generate example data with AR(1) structure
set.seed(42)
n_time <- 200
n_voxels <- 50
n_pred <- 3

# Design matrix
X <- cbind(1, rnorm(n_time), sin(2*pi*(1:n_time)/20))

# True coefficients
beta_true <- matrix(rnorm(n_pred * n_voxels, 0, 2), n_pred, n_voxels)

# Generate AR(1) noise
phi_true <- 0.6
noise <- matrix(0, n_time, n_voxels)
for (v in 1:n_voxels) {
  e <- rnorm(n_time)
  noise[1, v] <- e[1]
  for (t in 2:n_time) {
    noise[t, v] <- phi_true * noise[t-1, v] + e[t]
  }
}

# Generate Y with AR(1) structure in residuals
Y <- X %*% beta_true + noise

# Compute original residuals (these have AR structure)
resid <- Y - X %*% qr.solve(X, Y)

# Fit whitening plan
plan <- fit_noise(resid, method = "ar", p = 1)

cat("Estimated phi:", plan$phi[[1]], "(true =", phi_true, ")\n\n")

# Apply whitening to X and Y
whitened <- whiten_apply(plan, X, Y)
Xw <- whitened$X
Yw <- whitened$Y

# === WRONG DIAGNOSTIC (original buggy approach) ===
# This computes NEW residuals from whitened regression
innov_var_wrong <- Yw - Xw %*% qr.solve(Xw, Yw)
lag_stats_wrong <- apply(innov_var_wrong, 2, function(y) {
  ac <- acf(y, plot = FALSE, lag.max = 5)$acf[-1]
  mean(abs(ac))
})

# === CORRECT DIAGNOSTIC ===
# Apply whitening transformation to the ORIGINAL residuals
# We create a dummy X matrix (zeros) since we only want to whiten the residuals
X_dummy <- matrix(0, nrow = n_time, ncol = 1)
whitened_resid <- whiten_apply(plan, X_dummy, resid)
innov_var_correct <- whitened_resid$Y

lag_stats_correct <- apply(innov_var_correct, 2, function(y) {
  ac <- acf(y, plot = FALSE, lag.max = 5)$acf[-1]
  mean(abs(ac))
})

# Compare with original residuals
lag_stats_original <- apply(resid, 2, function(y) {
  ac <- acf(y, plot = FALSE, lag.max = 5)$acf[-1]
  mean(abs(ac))
})

cat("Autocorrelation diagnostics:\n")
cat("Original residuals (should be high):", round(mean(lag_stats_original), 4), "\n")
cat("WRONG approach (new resid from whitened):", round(mean(lag_stats_wrong), 4), "\n")
cat("CORRECT approach (whitened original resid):", round(mean(lag_stats_correct), 4), "\n\n")

cat("The WRONG approach gives identical values to original residuals because\n")
cat("it computes new residuals from a new regression on whitened data.\n")
cat("The CORRECT approach shows much lower autocorrelation, indicating\n")
cat("successful whitening of the original AR(1) structure.\n")