#include <Rcpp.h>
#include <algorithm>

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix parcel_means_cpp(const NumericMatrix& resid,
                               const IntegerVector& parcels,
                               int K = -1,
                               bool na_rm = false) {
  const int n = resid.nrow();
  const int v = resid.ncol();
  if (parcels.size() != v) stop("parcels length must equal ncol(resid)");

  int Keff = K > 0 ? K : *std::max_element(parcels.begin(), parcels.end());
  if (Keff <= 0) stop("invalid K");

  NumericMatrix sums(n, Keff);

  if (!na_rm) {
    std::vector<int> counts(Keff, 0);
    for (int j = 0; j < v; ++j) {
      int k = parcels[j];
      if (k < 1 || k > Keff) stop("parcel id out of range");
      double* outcol = &sums(0, k - 1);
      const double* col = &resid(0, j);
      for (int i = 0; i < n; ++i) outcol[i] += col[i];
      counts[k - 1] += 1;
    }
    for (int k = 0; k < Keff; ++k) {
      int cnt = counts[k] > 0 ? counts[k] : 1;
      double* outcol = &sums(0, k);
      for (int i = 0; i < n; ++i) outcol[i] /= cnt;
    }
  } else {
    IntegerMatrix counts(n, Keff);
    for (int j = 0; j < v; ++j) {
      int k = parcels[j];
      if (k < 1 || k > Keff) stop("parcel id out of range");
      double* outcol = &sums(0, k - 1);
      int* cntcol = &counts(0, k - 1);
      const double* col = &resid(0, j);
      for (int i = 0; i < n; ++i) {
        double val = col[i];
        if (!NumericVector::is_na(val)) {
          outcol[i] += val;
          cntcol[i] += 1;
        }
      }
    }
    for (int k = 0; k < Keff; ++k) {
      double* outcol = &sums(0, k);
      int* cntcol = &counts(0, k);
      for (int i = 0; i < n; ++i) {
        int cnt = cntcol[i] > 0 ? cntcol[i] : 1;
        outcol[i] /= cnt;
      }
    }
  }

  return sums;
}

// [[Rcpp::export]]
NumericVector run_avg_acvf_cpp(const NumericMatrix& mat, int max_lag) {
  const int n = mat.nrow();
  const int v = mat.ncol();
  if (max_lag < 0) stop("max_lag must be >= 0");
  if (n == 0 || v == 0) {
    return NumericVector(static_cast<R_xlen_t>(std::max(0, max_lag) + 1));
  }

  const int capped_lag = std::max(0, std::min(max_lag, n - 1));
  NumericVector gamma(capped_lag + 1);
  std::vector<double> means(v, 0.0);

  for (int j = 0; j < v; ++j) {
    const double* col = &mat(0, j);
    double sum = 0.0;
    for (int i = 0; i < n; ++i) sum += col[i];
    means[j] = sum / static_cast<double>(n);
  }

  double acc0 = 0.0;
  for (int j = 0; j < v; ++j) {
    const double* col = &mat(0, j);
    const double mu = means[j];
    for (int i = 0; i < n; ++i) {
      const double centered = col[i] - mu;
      acc0 += centered * centered;
    }
  }
  gamma[0] = acc0 / (static_cast<double>(n) * static_cast<double>(v));

  for (int lag = 1; lag <= capped_lag; ++lag) {
    double acc = 0.0;
    const int pairs = n - lag;
    if (pairs <= 0) {
      gamma[lag] = 0.0;
      continue;
    }
    for (int j = 0; j < v; ++j) {
      const double* col = &mat(0, j);
      const double mu = means[j];
      for (int t = lag; t < n; ++t) {
        const double yt = col[t] - mu;
        const double ylag = col[t - lag] - mu;
        acc += yt * ylag;
      }
    }
    gamma[lag] = acc / (static_cast<double>(pairs) * static_cast<double>(v));
  }

  return gamma;
}

// [[Rcpp::export]]
NumericVector segmented_acvf_cpp(const NumericVector& y,
                                 const IntegerVector& run_starts,
                                 int max_lag,
                                 bool unbiased = false,
                                 bool center = true) {
  const int n = y.size();
  if (n == 0) stop("y is empty");
  if (max_lag < 0) stop("max_lag must be >= 0");

  std::vector<int> starts;
  starts.reserve(run_starts.size() + 1);
  for (int idx = 0; idx < run_starts.size(); ++idx) {
    int s = run_starts[idx];
    if (s < 0 || s >= n) stop("run_starts out of bounds");
    if (!starts.empty() && s <= starts.back()) stop("run_starts must be strictly increasing");
    starts.push_back(s);
  }
  if (starts.empty() || starts.front() != 0) stop("run_starts must start at 0");
  starts.push_back(n);

  NumericVector gamma(max_lag + 1);
  std::vector<double> counts(max_lag + 1, 0.0);
  double total_n = 0.0;

  for (size_t seg = 0; seg + 1 < starts.size(); ++seg) {
    const int beg = starts[seg];
    const int end = starts[seg + 1];
    const int len = end - beg;
    if (len <= 0) continue;
    total_n += len;

    double mu = 0.0;
    if (center) {
      for (int i = beg; i < end; ++i) mu += y[i];
      mu /= static_cast<double>(len);
    }

    for (int lag = 0; lag <= max_lag; ++lag) {
      const int cnt = len - lag;
      if (cnt <= 0) break;
      double acc = 0.0;
      for (int t = lag; t < len; ++t) {
        const double yt = y[beg + t] - mu;
        const double ylag = y[beg + t - lag] - mu;
        acc += yt * ylag;
      }
      gamma[lag] += acc;
      counts[lag] += cnt;
    }
  }

  if (unbiased) {
    for (int lag = 0; lag <= max_lag; ++lag) {
      double denom = counts[lag] > 0.0 ? counts[lag] : 1.0;
      gamma[lag] /= denom;
    }
  } else {
    double denom = total_n > 0.0 ? total_n : 1.0;
    for (int lag = 0; lag <= max_lag; ++lag) {
      gamma[lag] /= denom;
    }
  }

  return gamma;
}

// [[Rcpp::export]]
Rcpp::List yw_from_acvf_cpp(const NumericVector& gamma, int p) {
  if (p < 0) stop("p must be >= 0");
  if (gamma.size() < p + 1) stop("gamma must have length >= p+1");

  if (p == 0) {
    return Rcpp::List::create(Rcpp::Named("phi") = NumericVector(0),
                              Rcpp::Named("sigma2") = gamma[0]);
  }

  double E_prev = gamma[0];
  if (!R_finite(E_prev) || E_prev <= 1e-12) {
    NumericVector phi(p);
    phi.fill(0.0);
    return Rcpp::List::create(Rcpp::Named("phi") = phi,
                              Rcpp::Named("sigma2") = 0.0);
  }

  std::vector<double> phi_prev(p, 0.0);
  std::vector<double> phi_cur(p, 0.0);

  for (int m = 1; m <= p; ++m) {
    double num = gamma[m];
    for (int j = 1; j <= m - 1; ++j) num -= phi_prev[j - 1] * gamma[m - j];
    double kappa = (std::abs(E_prev) > 1e-20) ? (num / E_prev) : 0.0;

    for (int j = 1; j <= m - 1; ++j) {
      phi_cur[j - 1] = phi_prev[j - 1] - kappa * phi_prev[(m - j) - 1];
    }
    phi_cur[m - 1] = kappa;

    for (int j = 0; j < m; ++j) phi_prev[j] = phi_cur[j];

    E_prev *= (1.0 - kappa * kappa);
    if (E_prev <= 0.0) E_prev = 1e-12;
  }

  NumericVector phi(p);
  for (int i = 0; i < p; ++i) phi[i] = phi_prev[i];

  return Rcpp::List::create(Rcpp::Named("phi") = phi,
                            Rcpp::Named("sigma2") = E_prev);
}

// Segment-respecting lag-product sums for a (pre-centred) matrix.
// seg_id labels contiguous segments; a pair (t, t - lag) contributes only when
// both lie in the same segment. Returns per-lag sums averaged over columns and
// the pair counts, matching .pooled_acvf_segments() after centring.
// [[Rcpp::export]]
Rcpp::List pooled_acvf_seg_cpp(const NumericMatrix& mat,
                               const IntegerVector& seg_id,
                               int max_lag) {
  const int n = mat.nrow();
  const int v = mat.ncol();
  if (seg_id.size() != n) stop("seg_id must have one entry per row");
  if (max_lag < 0) max_lag = 0;
  std::vector<int> seg_start(n, 0);
  for (int t = 1; t < n; ++t) {
    seg_start[t] = (seg_id[t] == seg_id[t - 1]) ? seg_start[t - 1] : t;
  }
  NumericVector num(max_lag + 1);
  NumericVector pairs(max_lag + 1);
  for (int t = 0; t < n; ++t) {
    const int reach = std::min(max_lag, t - seg_start[t]);
    for (int l = 0; l <= reach; ++l) pairs[l] += 1.0;
  }
  std::vector<double> acc(max_lag + 1, 0.0);
  for (int j = 0; j < v; ++j) {
    const double* x = &mat(0, j);
    for (int t = 0; t < n; ++t) {
      const int reach = std::min(max_lag, t - seg_start[t]);
      const double xt = x[t];
      for (int l = 0; l <= reach; ++l) acc[l] += xt * x[t - l];
    }
  }
  if (v > 0) {
    for (int l = 0; l <= max_lag; ++l) num[l] = acc[l] / static_cast<double>(v);
  }
  return Rcpp::List::create(Rcpp::Named("num") = num, Rcpp::Named("pairs") = pairs);
}

// Hannan-Rissanen normal equations summed over columns. Regressors are lags
// 1..p of Y and lags 1..q of E; rows enter when their within-segment position
// rel[t] >= start (start >= max(p, q) keeps every lag inside the segment).
// [[Rcpp::export]]
Rcpp::List hr_normal_eq_cpp(const NumericMatrix& Y,
                            const NumericMatrix& E,
                            const IntegerVector& rel,
                            int p, int q, int start) {
  const int n = Y.nrow();
  const int v = Y.ncol();
  const int k = p + q;
  if (E.nrow() != n || E.ncol() != v) stop("Y and E must have the same shape");
  if (rel.size() != n) stop("rel must have one entry per row");
  if (start < std::max(p, q)) stop("start must be >= max(p, q)");
  NumericMatrix G(k, k);
  NumericVector b(k);
  std::vector<double> z(k);
  for (int j = 0; j < v; ++j) {
    const double* y = &Y(0, j);
    const double* e = &E(0, j);
    for (int t = 0; t < n; ++t) {
      if (rel[t] < start) continue;
      for (int a = 0; a < p; ++a) z[a] = y[t - a - 1];
      for (int a = 0; a < q; ++a) z[p + a] = e[t - a - 1];
      const double yt = y[t];
      for (int a = 0; a < k; ++a) {
        b[a] += z[a] * yt;
        for (int c = a; c < k; ++c) G(a, c) += z[a] * z[c];
      }
    }
  }
  for (int a = 0; a < k; ++a) for (int c = 0; c < a; ++c) G(a, c) = G(c, a);
  return Rcpp::List::create(Rcpp::Named("G") = G, Rcpp::Named("b") = b);
}
