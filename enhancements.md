Here’s a prioritized roadmap to enhance your fMRI prewhitening utility while keeping the surface area small and the core blazing fast. I group items by impact/effort, include concise implementation sketches, and call out tests you can add so nothing regresses.

⸻

A) High‑impact, low‑risk (keep the API small)

1) Built‑in stability + selection sanity

Why: Fewer edge‑case failures, better generalization.
What:
	•	Enforce AR stationarity and MA invertibility inside HR (root reflection or PACF‑clip), then recompute σ².
	•	Fix BIC parameter count for AR order selection: use k = p + 1 (includes innovation variance).
	•	Offer q="auto" (tiny grid like q ∈ {0,1} with BIC/AICc) when p is small.

Tests: synthetic AR and ARMA grids; ensure selected (p,q) matches truth ≥70% on n∈[300,1000].

⸻

2) Unique‑filter caching

Why: Huge speed win when many parcels share (nearly) identical filters (common after multiscale pooling).
What:
	•	Hash each run’s (phi, theta) to a key (e.g., rounded to 1e‑6 + digest::digest), whiten X once per unique key, map to all parcels using that key.

Sketch (R):

hash_filter <- function(phi, theta, tol = 1e-6) {
  paste0(round(c(phi, NA, theta), 6), collapse = ",")
}
# whiten_apply(): build a map from hash -> column indices, compute Xw per hash, reuse.

Expected: 2–10× less X‑whitening when many parcels collapse to the same plan.

⸻

3) Deterministic threading + user knob

Why: Reproducibility & UX.
What: Export set_threads(n) wrapper over the internal option; pin OMP_NUM_THREADS for whitening calls.

Sketch (R):

set_threads <- function(n) {
  stopifnot(n >= 1L)
  options(fmriAR.max_threads = as.integer(n))
  invisible(n)
}

Tests: the parallel determinism test you already have (exact equality across 1/2/8 threads).

⸻

4) Parcel means: fast path by default

Why: Multiscale estimation time becomes negligible.
What: Use your parcel_means_cpp() everywhere. Keep R fallback behind an option (you already have this).

Tests: verify equality to R implementation on random small inputs; micro‑benchmark.

⸻

5) Whiteness QA helpers (S3)

Why: Users need a one‑liner to validate plans.
What: Add check_whiteness(plan, X, Y, run_starts, lag=10) that returns:
	•	KS distance to Uniform for Ljung–Box p‑values
	•	fraction p≤0.05
	•	lag‑wise ACF energy

Expose plot() method (base or ggplot2) without adding heavy deps by making ggplot2 Suggests.

Tests: your “multiscale improves KS/NLL” tests double as acceptance tests.

⸻

6) Robust SE variants & contrasts

Why: Downstream inference benefits.
What: Extend sandwich_from_whitened_resid() to support type = c("iid","hc0","hc1","hc2","hc3") and a contrast(beta, L) helper.

Sketch (R):

sandwich_from_whitened_resid <- function(Xw, ew, type=c("iid","hc0","hc1","hc2","hc3")) { ... }
contrast <- function(beta, Sigma_beta, L) {
  est <- as.numeric(L %*% beta); var <- L %*% Sigma_beta %*% t(L)
  list(estimate=est, se=sqrt(diag(var)), Sigma=var)
}

Tests: HC3 ≥ HC0 elementwise on heteroskedastic sims; contrasts match closed‑form on simple designs.

⸻

B) Modeling capabilities (medium effort, big value)

7) Kalman innovations (exact ARMA, missing data tolerant)

Why: Exact likelihood and innovations for ARMA; handles censored samples and arbitrary segmenting without ad‑hoc resets.
What: Add a kalman_whiten_inplace() (state‑space for ARMA p,q) used when missing or q>0 with complex segment patterns; fallback to current fast loop for pure AR or dense data.

Implementation notes:
	•	Build companion form for AR(p), augment with MA(q) to innovations form; use steady‑state Kalman when possible; otherwise per‑segment reinit.
	•	Keep identical API; select engine internally.

Tests: missing points per voxel; compare to dense case where NA dropped → identical innovations on overlapping indices.

⸻

8) Local (slowly varying) AR: piecewise‑stationary

Why: Nonstationarities across long runs (drift, vigilance).
What: Optional sliding‑window AR with overlap (e.g., windows of 60–120 TRs, 50% overlap), smooth φ across windows (e.g., quadratic penalty), and whiten using time‑varying filter (change points handled at window edges with small cross‑fade).

API idea: fit_noise(..., local_ar = list(window=90, overlap=0.5, penalty="quad")).

Tests: simulate φ(t) drift; multiscale+local reduces ACF energy vs static AR.

⸻

9) Empirical‑Bayes weight learning for multiscale

Why: Replace hand‑tuned weights with data‑driven shrinkage strength.
What: Work in PACF‑space: use z = atanh(κ) (Gaussianizable); put a normal prior centered at parent z with variance τ²; estimate τ² per lag (and per level) by EB (marginal likelihood or SURE). The posterior mean becomes the shrinkage target.

Benefit: principled control of shrinkage; adapts across datasets.

Tests: on hierarchical sims, EB‑weights should beat fixed weights on held‑out KS/NLL.

⸻

10) Physio‑aware nuisance templates (optional input)

Why: Better residuals when physio traces exist.
What: Provide a tiny helper to ingest fMRIPrep confounds or RETROICOR regressors; nothing fancy—just a convenience build_confounds() that returns an X with recommended sels.

Tests: presence/absence yields improved whiteness on physio‑simulated data.

⸻

C) Performance & scale

11) Streaming whitening (chunked)

Why: Reduce peak memory & enable very long runs.
What: Process Y in chunks by time with overlap L = max(p,q); maintain per‑column state across chunks.

API: whiten_apply(..., stream_by = 2000).

Tests: stream vs full → identical outputs (within machine epsilon).

⸻

12) SIMD‑friendly inner loops

Why: Small but real speedups, especially AR(1–2).
What:
	•	Add #pragma omp simd on the time loop;
	•	Provide specialized kernels for p=1 and p=2 (branch once per column).

Sketch (C++):

if (p==1 && q==0) whiten_ar1(colPtr, n, phi0, segs);
else if (p==2 && q==0) whiten_ar2(colPtr, n, phi0, phi1, segs);
else whiten_arma_generic(...);

Tests: speed micro‑benchmarks; outputs identical.

⸻

13) Memoize filtered design blocks

Why: X often repeats across subjects/runs; avoid recomputation across calls.
What: Option to cache Xw per (design_hash, filter_hash, run) in a user‑provided cache env.

API: whiten_apply(..., cache = new.env()).

Tests: repeated calls reuse cache (instrument counters).

⸻

D) Developer experience & reliability

14) Plan I/O & reproducibility

Why: Shareable, cached plans across sessions.
What: save_plan(plan, file) / load_plan(file) using jsonlite; include version, p/q, φ/θ by parcel, weights, order rule, options.

Tests: round‑trip equality and compatibility across minor versions.

⸻

15) Vignettes & autoplot

Why: Onboarding & trust.
What:
	•	“Prewhitening for fMRI in practice” (why, how, diagnostics).
	•	“Multiscale pooling: when it helps and when it doesn’t.”
	•	“Exact ARMA via Kalman with missing data.”

⸻

16) CI + sanitizers

Why: Catch UB and threading bugs early.
What: GH Actions matrix: Linux/macOS, R release/devel; run ASan/UBSan builds for C++; run testthat with OMP_NUM_THREADS=1,2,8.

⸻

17) Configuration for OpenMP portability

Why: CRAN friendliness.
What: Autodetect OMP; if absent, compile a serial fallback; expose capabilities() helper to show threading status.

⸻

E) Tiny polish (nice to have)
	•	print() / summary() for fmriAR_plan (show (p,q), #unique filters, weights summary, KS from training).
	•	Better errors on invalid run_starts: echo first bad index and expectations (strictly increasing, 0‑based, includes 0).
	•	Progress reporting (optional): use progressr in R; no C++ calls from worker threads.

⸻

Concrete code snippets you can drop in

A robust OLS fallback in HR (prevents singularity issues)

arma::vec coef;
bool ok = arma::solve(coef, Z, ysub, arma::solve_opts::fast);
if (!ok || coef.n_elem != Z.n_cols) {
  ok = arma::solve(coef, Z, ysub); // QR/LAPACK fallback
}
if (!ok) return hr_failure(p, q, p_big, iter);

AR(1) specialized kernel (illustrative)

inline void whiten_ar1(double* y, int n, double phi0,
                       const int* seg_beg, const int* seg_end, int nseg,
                       double scale_first) {
  for (int s=0; s<nseg; ++s) {
    int a = seg_beg[s], b = seg_end[s];
    double ym1 = 0.0, em1 = 0.0; // not used for pure AR
    if (a < b) {
      double e = (y[a] - phi0 * ym1) * scale_first;
      ym1 = y[a];
      y[a] = e;
      for (int t=a+1; t<b; ++t) {
        double st = y[t] - phi0 * ym1;
        ym1 = y[t];
        y[t] = st;
      }
    }
  }
}

set_threads() + capability probe

set_threads <- function(n) {
  stopifnot(n >= 1L)
  options(fmriAR.max_threads = as.integer(n))
  invisible(n)
}
capabilities_fmriAR <- function() {
  list(openmp = fmriAR:::omp_capable(),  # C++ tiny function returning bool
       simd   = TRUE)                    # if you add omp simd, toggle here
}

Plan I/O

save_plan <- function(plan, file) {
  jsonlite::write_json(fmriAR:::compat$plan_info(plan), path = file, auto_unbox = TRUE, pretty = TRUE)
}
load_plan <- function(file) {
  info <- jsonlite::read_json(file, simplifyVector = TRUE)
  fmriAR:::compat$plan_from_phi(info$phi_by_parcel, info$p, info$q, info$options)
}


⸻

Suggested acceptance tests (add to tests/testthat/)
	•	EB multiscale vs fixed: EB weights must not regress KS/NLL relative to fixed weights on 3 synthetic seeds.
	•	Kalman vs generic: For data without missingness, Kalman innovations ≈ generic innovations (1e‑12).
	•	Streaming vs full: Chunked whitening equals full whitening (1e‑12) on random data with multiple segments.
	•	Cache hit rate: With duplicated filters across parcels, number of X‑whitenings is reduced exactly to #unique_filters.

⸻

Summary

You can keep the interface tiny and still make the package feel “complete” by:
	•	hardening estimation (stability + selection),
	•	cutting repeated work (filter caching, streaming),
	•	offering one exact method for tricky cases (Kalman),
	•	giving users first‑class diagnostics and an easy way to reproduce/share plans.

If you want, I can prep a branch layout (feature flags, files to touch) and a minimal bench + QA vignette skeleton to drop into vignettes/.