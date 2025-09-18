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
