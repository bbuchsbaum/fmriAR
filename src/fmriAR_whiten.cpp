#define ARMA_DONT_PRINT_FAST_MATH_WARNING
#include <RcppArmadillo.h>
#include <algorithm>
#include <vector>
#ifdef _OPENMP
  #include <omp.h>
#endif
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

auto ar1_first_scale(double phi1) -> double {
  const double s2 = 1.0 - phi1 * phi1;
  return (s2 > 0.0) ? std::sqrt(s2) : 1.0;
}

inline void reset_buffers(std::vector<double>& prev_in,
                          std::vector<double>& prev_out) {
  std::fill(prev_in.begin(), prev_in.end(), 0.0);
  std::fill(prev_out.begin(), prev_out.end(), 0.0);
}

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

template <typename Matrix>
inline void whiten_matrix_arma_impl(Matrix& M,
                                    const arma::vec& phi,
                                    const arma::vec& theta,
                                    const std::vector<int>& run_starts,
                                    bool exact_first_ar1,
                                    bool parallel) {
  const int n_time = M.nrow();
  const int n_cols = M.ncol();
  const int p = static_cast<int>(phi.n_elem);
  const int q = static_cast<int>(theta.n_elem);
  if (p == 0 && q == 0) return;

  std::vector<int> seg_beg;
  std::vector<int> seg_end;
  build_segments(run_starts, n_time, seg_beg, seg_end);
  const int S = static_cast<int>(seg_beg.size());

  const double scale_ar1 = (exact_first_ar1 && p == 1 && q == 0)
                             ? ar1_first_scale(phi[0]) : 1.0;

  #pragma omp parallel for if(parallel) schedule(static)
  for (int col = 0; col < n_cols; ++col) {
#ifdef _OPENMP
    if (omp_get_thread_num() == 0 && (col & 63) == 0) Rcpp::checkUserInterrupt();
#else
    if ((col & 63) == 0) Rcpp::checkUserInterrupt();
#endif
    double* colPtr = &M(0, col);
    std::vector<double> prev_in(std::max(1, p), 0.0);
    std::vector<double> prev_out(std::max(1, q), 0.0);

    for (int s = 0; s < S; ++s) {
      reset_buffers(prev_in, prev_out);

      for (int t = seg_beg[s]; t < seg_end[s]; ++t) {
        const double orig = colPtr[t];

        double st = orig;
        for (int k = 0; k < p; ++k) st -= phi[k] * prev_in[k];
        if (t == seg_beg[s] && scale_ar1 != 1.0) st *= scale_ar1;

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

        colPtr[t] = zt;
      }
    }
  }
}

// [[Rcpp::export]]
Rcpp::List arma_whiten_inplace(Rcpp::NumericMatrix Y,
                               Rcpp::NumericMatrix X,
                               const arma::vec& phi,
                               const arma::vec& theta,
                               Rcpp::IntegerVector run_starts,
                               bool exact_first_ar1 = false,
                               bool parallel = true) {
  std::vector<int> rs(run_starts.begin(), run_starts.end());
  whiten_matrix_arma_impl(Y, phi, theta, rs, exact_first_ar1, parallel);
  whiten_matrix_arma_impl(X, phi, theta, rs, exact_first_ar1, parallel);
  return Rcpp::List::create(Rcpp::Named("Y") = Y,
                            Rcpp::Named("X") = X);
}

// [[Rcpp::export]]
void arma_whiten_void(Rcpp::NumericMatrix Y,
                      Rcpp::NumericMatrix X,
                      const arma::vec& phi,
                      const arma::vec& theta,
                      Rcpp::IntegerVector run_starts,
                      bool exact_first_ar1 = false,
                      bool parallel = true) {
  std::vector<int> rs(run_starts.begin(), run_starts.end());
  whiten_matrix_arma_impl(Y, phi, theta, rs, exact_first_ar1, parallel);
  whiten_matrix_arma_impl(X, phi, theta, rs, exact_first_ar1, parallel);
}
