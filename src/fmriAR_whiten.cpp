#define ARMA_DONT_PRINT_FAST_MATH_WARNING
#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <vector>
#ifdef _OPENMP
  #include <omp.h>
#endif
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

inline void build_segments(const std::vector<int>& run_starts,
                           int n_time,
                           std::vector<int>& seg_beg,
                           std::vector<int>& seg_end) {
  if (run_starts.empty()) {
    Rcpp::stop("run_starts must contain at least one index (0)" );
  }
  if (run_starts.front() != 0) {
    Rcpp::stop("run_starts must be 0-based with first element equal to 0");
  }
  int prev = -1;
  for (int s : run_starts) {
    if (s < 0 || s >= n_time) {
      Rcpp::stop("run_starts value out of bounds");
    }
    if (s <= prev) {
      Rcpp::stop("run_starts must be strictly increasing");
    }
    prev = s;
  }

  const std::size_t S = run_starts.size();
  seg_beg.resize(S);
  seg_end.resize(S);
  for (std::size_t idx = 0; idx < S; ++idx) {
    seg_beg[idx] = run_starts[idx];
    seg_end[idx] = (idx + 1 < S) ? run_starts[idx + 1] : n_time;
  }
}

// Autocovariance of a causal ARMA(p, q) with unit innovation variance,
//   y_t = sum phi_k y_{t-k} + e_t + sum theta_j e_{t-j},
// at lags 0..max_lag, from the exact linear system of Brockwell & Davis (1991,
// sec. 3.3). Returns false when the system is singular or gamma_0 <= 0, which
// signals a non-stationary AR part.
static bool arma_acvf_unit(const std::vector<double>& phi,
                           const std::vector<double>& theta,
                           int max_lag,
                           std::vector<double>& gamma) {
  const int p = static_cast<int>(phi.size());
  const int q = static_cast<int>(theta.size());

  std::vector<double> psi(q + 1, 0.0);
  psi[0] = 1.0;
  for (int j = 1; j <= q; ++j) {
    double s = theta[j - 1];
    for (int k = 1; k <= std::min(j, p); ++k) s += phi[k - 1] * psi[j - k];
    psi[j] = s;
  }
  auto th = [&](int j) { return j == 0 ? 1.0 : theta[j - 1]; };
  auto rhs = [&](int k) {
    double s = 0.0;
    for (int j = k; j <= q; ++j) s += th(j) * psi[j - k];
    return s;
  };

  const int L = std::max(max_lag, p);
  gamma.assign(L + 1, 0.0);
  arma::mat A(p + 1, p + 1, arma::fill::zeros);
  arma::vec b(p + 1);
  for (int k = 0; k <= p; ++k) {
    A(k, k) += 1.0;
    for (int j = 1; j <= p; ++j) A(k, std::abs(k - j)) -= phi[j - 1];
    b[k] = rhs(k);
  }
  arma::vec g;
  bool ok = arma::solve(g, A, b, arma::solve_opts::no_approx);
  if (!ok || !g.is_finite() || g[0] <= 0.0) return false;
  for (int k = 0; k <= p; ++k) gamma[k] = g[k];
  for (int k = p + 1; k <= L; ++k) {
    double s = rhs(k);
    for (int j = 1; j <= p; ++j) s += phi[j - 1] * gamma[k - j];
    gamma[k] = s;
  }
  gamma.resize(max_lag + 1);
  return true;
}

// Exact (stationary) initialisation via Ansley's (1979) transformation.
//
// With m = max(p, q), z_t = y_t for t < m and z_t = phi(B) y_t afterwards. The
// covariance Omega of z is banded with half-bandwidth m, and z = A y with A unit
// lower triangular, so chol(Omega) = A chol(Sigma) and chol(Omega)^-1 z equals
// chol(Sigma)^-1 y: exact GLS whitening of a stationary segment, in O(n m) per
// column. For a pure AR the factor is the identity after row p, so this reduces
// to the usual recursion with only the first p outputs rescaled; for AR(1) that
// is exactly the Prais-Winsten sqrt(1 - phi^2) scaling.
//
// The factor for a shorter segment is the leading block of the factor for a
// longer one, so one factor of the longest segment length serves every segment.
struct ExactFactor {
  int m = 0;
  int n = 0;
  std::vector<double> L;  // row t stores L(t, t - d) at index t * (m + 1) + d
  double at(int t, int d) const { return L[static_cast<std::size_t>(t) * (m + 1) + d]; }
};

static bool build_exact_factor(const std::vector<double>& phi,
                               const std::vector<double>& theta,
                               int n, ExactFactor& F) {
  const int p = static_cast<int>(phi.size());
  const int q = static_cast<int>(theta.size());
  const int m = std::max(p, q);
  F.m = m;
  F.n = n;
  if (m == 0 || n <= 0) return true;

  std::vector<double> gamma;
  if (!arma_acvf_unit(phi, theta, m + p, gamma)) return false;
  auto th = [&](int j) { return j == 0 ? 1.0 : theta[j - 1]; };

  auto omega = [&](int i, int j) {
    // 0-based, i <= j
    const int d = j - i;
    if (j < m) return gamma[d];
    if (i < m) {
      double s = gamma[d];
      for (int k = 1; k <= p; ++k) s -= phi[k - 1] * gamma[std::abs(d - k)];
      return s;
    }
    if (d > q) return 0.0;
    double s = 0.0;
    for (int l = 0; l + d <= q; ++l) s += th(l) * th(l + d);
    return s;
  };

  const int w = m + 1;
  F.L.assign(static_cast<std::size_t>(n) * w, 0.0);
  for (int t = 0; t < n; ++t) {
    const int j0 = std::max(0, t - m);
    for (int j = j0; j <= t; ++j) {
      double s = omega(j, t);
      const int k0 = std::max(j0, j - m);
      for (int k = k0; k < j; ++k) {
        s -= F.L[static_cast<std::size_t>(t) * w + (t - k)] *
             F.L[static_cast<std::size_t>(j) * w + (j - k)];
      }
      if (j == t) {
        if (!(s > 0.0) || !std::isfinite(s)) return false;
        F.L[static_cast<std::size_t>(t) * w] = std::sqrt(s);
      } else {
        F.L[static_cast<std::size_t>(t) * w + (t - j)] =
          s / F.L[static_cast<std::size_t>(j) * w];
      }
    }
  }
  return true;
}

// Conditional (zero pre-sample) ARMA inverse filter for one segment.
static inline void filter_conditional(double* x, int beg, int end,
                                      const std::vector<double>& phi,
                                      const std::vector<double>& theta,
                                      double first_scale,
                                      std::vector<double>& prev_in,
                                      std::vector<double>& prev_out) {
  const int p = static_cast<int>(phi.size());
  const int q = static_cast<int>(theta.size());
  std::fill(prev_in.begin(), prev_in.end(), 0.0);
  std::fill(prev_out.begin(), prev_out.end(), 0.0);
  for (int t = beg; t < end; ++t) {
    const double orig = x[t];
    double st = orig;
    for (int k = 0; k < p; ++k) st -= phi[k] * prev_in[k];
    if (t == beg && first_scale != 1.0) st *= first_scale;
    double zt = st;
    for (int j = 0; j < q; ++j) zt -= theta[j] * prev_out[j];
    if (p > 0) {
      for (int k = p - 1; k > 0; --k) prev_in[k] = prev_in[k - 1];
      prev_in[0] = orig;
    }
    if (q > 0) {
      for (int j = q - 1; j > 0; --j) prev_out[j] = prev_out[j - 1];
      prev_out[0] = zt;
    }
    x[t] = zt;
  }
}

// Exact stationary-start filter for one segment using the banded factor.
static inline void filter_exact(double* x, int beg, int end,
                                const std::vector<double>& phi,
                                const ExactFactor& F,
                                std::vector<double>& y_hist,
                                std::vector<double>& out_hist) {
  const int p = static_cast<int>(phi.size());
  const int m = F.m;
  const int len = end - beg;
  // Ring buffers of the last p inputs and last m outputs, indexed by t mod size.
  const int hp = std::max(1, p);
  const int hm = std::max(1, m);
  for (int t = 0; t < len; ++t) {
    const double yt = x[beg + t];
    double z = yt;
    if (t >= m) {
      for (int k = 1; k <= p; ++k) z -= phi[k - 1] * y_hist[(t - k) % hp];
    }
    const int j0 = std::max(0, t - m);
    for (int j = j0; j < t; ++j) z -= F.at(t, t - j) * out_hist[j % hm];
    const double xt = z / F.at(t, 0);
    if (p > 0) y_hist[t % hp] = yt;
    if (m > 0) out_hist[t % hm] = xt;
    x[beg + t] = xt;
  }
}

template <typename Matrix>
inline void whiten_matrix_arma_impl(Matrix& M,
                                    const arma::vec& phi_in,
                                    const arma::vec& theta_in,
                                    const std::vector<int>& run_starts,
                                    bool exact_first,
                                    bool parallel,
                                    int n_threads) {
  const int n_time = M.nrow();
  const int n_cols = M.ncol();
  const std::vector<double> phi(phi_in.begin(), phi_in.end());
  const std::vector<double> theta(theta_in.begin(), theta_in.end());
  const int p = static_cast<int>(phi.size());
  const int q = static_cast<int>(theta.size());
  if (n_time == 0 || n_cols == 0) return;
  if (p == 0 && q == 0) return;

  std::vector<int> seg_beg;
  std::vector<int> seg_end;
  build_segments(run_starts, n_time, seg_beg, seg_end);
  const int S = static_cast<int>(seg_beg.size());

  ExactFactor F;
  bool use_exact = false;
  if (exact_first) {
    int max_len = 0;
    for (int s = 0; s < S; ++s) max_len = std::max(max_len, seg_end[s] - seg_beg[s]);
    use_exact = build_exact_factor(phi, theta, max_len, F);
  }
  // Fallback when the exact factor is unavailable (non-stationary AR part):
  // the conditional recursion, with the historic AR(1) first-sample scaling.
  double first_scale = 1.0;
  if (exact_first && !use_exact && p == 1 && q == 0) {
    const double s2 = 1.0 - phi[0] * phi[0];
    if (s2 > 0.0) first_scale = std::sqrt(s2);
  }

  const int m = std::max(p, q);
  auto do_column = [&](int col, std::vector<double>& b1, std::vector<double>& b2) {
    double* colPtr = &M(0, col);
    for (int s = 0; s < S; ++s) {
      if (use_exact) {
        filter_exact(colPtr, seg_beg[s], seg_end[s], phi, F, b1, b2);
      } else {
        filter_conditional(colPtr, seg_beg[s], seg_end[s], phi, theta,
                           first_scale, b1, b2);
      }
    }
  };
  const std::size_t n1 = static_cast<std::size_t>(std::max({1, p, m}));
  const std::size_t n2 = static_cast<std::size_t>(std::max({1, q, m}));

#ifdef _OPENMP
  if (parallel && n_cols > 1) {
    // No R API calls (including checkUserInterrupt) inside the parallel
    // region: an R error longjmps / a C++ exception escaping an OpenMP region
    // terminates the process.
    const int nt = n_threads > 0 ? n_threads : omp_get_max_threads();
    #pragma omp parallel num_threads(nt)
    {
      std::vector<double> b1(n1, 0.0), b2(n2, 0.0);
      #pragma omp for schedule(static)
      for (int col = 0; col < n_cols; ++col) do_column(col, b1, b2);
    }
    return;
  }
#else
  (void) parallel;
  (void) n_threads;
#endif
  std::vector<double> b1(n1, 0.0), b2(n2, 0.0);
  for (int col = 0; col < n_cols; ++col) {
    if ((col & 255) == 0) Rcpp::checkUserInterrupt();
    do_column(col, b1, b2);
  }
}

// [[Rcpp::export]]
Rcpp::List arma_whiten_inplace(Rcpp::NumericMatrix Y,
                               Rcpp::NumericMatrix X,
                               const arma::vec& phi,
                               const arma::vec& theta,
                               Rcpp::IntegerVector run_starts,
                               bool exact_first = false,
                               bool parallel = true,
                               int n_threads = 0) {
  std::vector<int> rs(run_starts.begin(), run_starts.end());
  whiten_matrix_arma_impl(Y, phi, theta, rs, exact_first, parallel, n_threads);
  whiten_matrix_arma_impl(X, phi, theta, rs, exact_first, parallel, n_threads);
  return Rcpp::List::create(Rcpp::Named("Y") = Y,
                            Rcpp::Named("X") = X);
}

// [[Rcpp::export]]
void arma_whiten_void(Rcpp::NumericMatrix Y,
                      Rcpp::NumericMatrix X,
                      const arma::vec& phi,
                      const arma::vec& theta,
                      Rcpp::IntegerVector run_starts,
                      bool exact_first = false,
                      bool parallel = true,
                      int n_threads = 0) {
  std::vector<int> rs(run_starts.begin(), run_starts.end());
  whiten_matrix_arma_impl(Y, phi, theta, rs, exact_first, parallel, n_threads);
  whiten_matrix_arma_impl(X, phi, theta, rs, exact_first, parallel, n_threads);
}

// Exact unit-innovation ARMA autocovariance, exposed for tests and for R-side
// covariance construction. Returns NA when the AR part is not stationary.
// [[Rcpp::export]]
Rcpp::NumericVector arma_acvf_cpp(const arma::vec& phi, const arma::vec& theta,
                                  int max_lag) {
  if (max_lag < 0) Rcpp::stop("max_lag must be >= 0");
  std::vector<double> g;
  std::vector<double> ph(phi.begin(), phi.end());
  std::vector<double> th(theta.begin(), theta.end());
  if (!arma_acvf_unit(ph, th, max_lag, g)) {
    Rcpp::NumericVector out(max_lag + 1, NA_REAL);
    return out;
  }
  return Rcpp::NumericVector(g.begin(), g.end());
}
