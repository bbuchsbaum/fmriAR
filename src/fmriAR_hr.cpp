#include <RcppArmadillo.h>
#include <algorithm>
#include <vector>
#include <complex>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

static inline arma::vec acvf_biased(const arma::vec& x, const int L) {
  const int n = x.n_elem;
  const double mu = arma::mean(x);
  arma::vec c(L + 1, arma::fill::zeros);
  for (int lag = 0; lag <= L; ++lag) {
    double s = 0.0;
    for (int t = lag; t < n; ++t) {
      s += (x[t] - mu) * (x[t - lag] - mu);
    }
    c[lag] = s / static_cast<double>(n);
  }
  return c;
}

static void levinson_durbin(const arma::vec& gamma, const int p,
                            arma::vec& phi_out, arma::vec& kappa_out, double& sigma2_out) {
  const double eps = 1e-12;
  phi_out.set_size(p);
  kappa_out.set_size(p);
  arma::vec phi_prev;
  arma::vec phi_cur;

  double v = std::max(eps, gamma[0]);
  for (int m = 1; m <= p; ++m) {
    double acc = gamma[m];
    for (int k = 1; k <= m - 1; ++k) acc -= phi_prev[k - 1] * gamma[m - k];
    double kappa = acc / std::max(eps, v);
    if (kappa > 0.999999) kappa = 0.999999;
    if (kappa < -0.999999) kappa = -0.999999;
    kappa_out[m - 1] = kappa;

    phi_cur.set_size(m);
    for (int j = 1; j <= m - 1; ++j) {
      phi_cur[j - 1] = phi_prev[j - 1] - kappa * phi_prev[(m - j) - 1];
    }
    phi_cur[m - 1] = kappa;

    v *= (1.0 - kappa * kappa);
    if (v < eps) v = eps;
    phi_prev = phi_cur;
  }
  phi_out = phi_prev;
  sigma2_out = v;
}

static arma::vec arma_innovations_vec(const arma::vec& y,
                                      const arma::vec& phi,
                                      const arma::vec& theta) {
  const int n = y.n_elem;
  const int p = phi.n_elem;
  const int q = theta.n_elem;
  arma::vec e(n, arma::fill::zeros);
  for (int t = 0; t < n; ++t) {
    double st = y[t];
    for (int k = 0; k < p; ++k) {
      const int idx = t - k - 1;
      st -= phi[k] * (idx >= 0 ? y[idx] : 0.0);
    }
    double et = st;
    for (int j = 0; j < q; ++j) {
      const int idx = t - j - 1;
      et -= theta[j] * (idx >= 0 ? e[idx] : 0.0);
    }
    e[t] = et;
  }
  return e;
}

static inline arma::cx_vec poly_from_roots(const arma::cx_vec& roots) {
  const arma::uword n = roots.n_elem;
  arma::cx_vec coeff(1);
  coeff[0] = arma::cx_double(1.0, 0.0);
  for (arma::uword idx = 0; idx < n; ++idx) {
    arma::cx_vec next(coeff.n_elem + 1, arma::fill::zeros);
    for (arma::uword j = 0; j < coeff.n_elem; ++j) {
      next[j] += coeff[j];
      next[j + 1] -= coeff[j] * roots[idx];
    }
    coeff = next;
  }
  return coeff;
}

// e_t = y_t - sum phi_i y_{t-i} - sum theta_j e_{t-j}; enforce AR stability
static inline void enforce_ar_stationarity(arma::vec& phi) {
  const int p = static_cast<int>(phi.n_elem);
  if (p == 0) return;

  arma::cx_vec poly(p + 1, arma::fill::zeros);
  for (int i = 0; i < p; ++i) poly[i] = -phi[p - i - 1];
  poly[p] = 1.0;

  arma::cx_vec roots = arma::roots(poly);
  bool changed = false;
  for (auto& z : roots) {
    const double mag = std::abs(z);
    if (mag <= 1.0) {
      z = (1.0 / std::conj(z)) * 1.0001;
      changed = true;
    }
  }

  if (changed) {
    arma::cx_vec coeff = poly_from_roots(roots);
    const arma::cx_double scale = coeff[coeff.n_elem - 1];
    if (std::abs(scale) < 1e-12) return;
    coeff /= scale;
    for (int i = 0; i < p; ++i) {
      phi[i] = -coeff[p - (i + 1)].real();
    }
  }
}

// ensure MA polynomial 1 + sum theta_j z^j has roots outside unit circle
static inline void enforce_ma_invertibility(arma::vec& theta) {
  const int q = static_cast<int>(theta.n_elem);
  if (q == 0) return;

  arma::cx_vec poly(q + 1, arma::fill::zeros);
  for (int j = 0; j < q; ++j) poly[j] = theta[q - j - 1];
  poly[q] = 1.0;

  arma::cx_vec roots = arma::roots(poly);
  bool changed = false;
  for (auto& z : roots) {
    const double mag = std::abs(z);
    if (mag <= 1.0) {
      z = (1.0 / std::conj(z)) * 1.0001;
      changed = true;
    }
  }

  if (changed) {
    arma::cx_vec coeff = poly_from_roots(roots);
    const arma::cx_double scale = coeff[coeff.n_elem - 1];
    if (std::abs(scale) < 1e-12) return;
    coeff /= scale;
    for (int j = 0; j < q; ++j) {
      theta[j] = coeff[q - (j + 1)].real();
    }
  }
}

static Rcpp::List hr_failure(int p, int q, int p_big, int iter) {
  arma::vec phi_fail = arma::zeros<arma::vec>(p);
  arma::vec theta_fail = arma::zeros<arma::vec>(q);
  return Rcpp::List::create(
    Rcpp::Named("phi")    = phi_fail,
    Rcpp::Named("theta")  = theta_fail,
    Rcpp::Named("sigma2") = NA_REAL,
    Rcpp::Named("p_big")  = p_big,
    Rcpp::Named("iter")   = iter,
    Rcpp::Named("ok")     = false
  );
}

// [[Rcpp::export]]
Rcpp::List hr_arma_fit_cpp(const arma::vec& y_in,
                           int p, int q,
                           int p_big = 0,
                           int iter = 0) {
  if (p < 0 || q < 0) stop("p and q must be >= 0");
  const int n = y_in.n_elem;
  if (n < 10) stop("series too short");

  arma::vec y = y_in - arma::mean(y_in);

  if (p_big <= 0) {
    p_big = std::min(std::max(8, p + q + 5), std::min(std::max(2, n - 2), 40));
  }

  arma::vec gamma = acvf_biased(y, p_big);
  arma::vec phi_big;
  arma::vec kappa_big;
  double v_big = 0.0;
  levinson_durbin(gamma, p_big, phi_big, kappa_big, v_big);

  arma::vec ehat = arma_innovations_vec(y, phi_big, arma::vec());

  const int mlag = std::max(p, q);
  if (n - mlag <= p + q) return hr_failure(p, q, p_big, iter);

  if ((p + q) == 0) {
    const double sigma2 = arma::mean(arma::square(y));
    arma::vec empty_phi = arma::zeros<arma::vec>(0);
    return Rcpp::List::create(
      Rcpp::Named("phi")    = empty_phi,
      Rcpp::Named("theta")  = empty_phi,
      Rcpp::Named("sigma2") = sigma2,
      Rcpp::Named("p_big")  = p_big,
      Rcpp::Named("iter")   = iter,
      Rcpp::Named("ok")     = true
    );
  }

  arma::vec phi(p, arma::fill::zeros);
  arma::vec theta(q, arma::fill::zeros);

  for (int it = 0; it <= iter; ++it) {
    const int rows = n - mlag;
    const int cols = p + q;
    arma::mat Z(rows, cols, arma::fill::zeros);
    arma::vec ysub = y.subvec(mlag, n - 1);

    for (int i = 1; i <= p; ++i) {
      Z.col(i - 1) = y.subvec(mlag - i, n - 1 - i);
    }
    for (int j = 1; j <= q; ++j) {
      Z.col(p + j - 1) = ehat.subvec(mlag - j, n - 1 - j);
    }

    arma::vec coef;
    bool ok = arma::solve(coef, Z, ysub);
    if (!ok || coef.n_elem != cols) {
      return hr_failure(p, q, p_big, iter);
    }

    if (p > 0)  phi   = coef.subvec(0, p - 1);
    if (q > 0)  theta = coef.subvec(p, p + q - 1);
    enforce_ar_stationarity(phi);
    enforce_ma_invertibility(theta);

    ehat = arma_innovations_vec(y, phi, theta);
  }

  double sigma2 = arma::mean(arma::square(ehat));

  return Rcpp::List::create(
    Rcpp::Named("phi")    = phi,
    Rcpp::Named("theta")  = theta,
    Rcpp::Named("sigma2") = sigma2,
    Rcpp::Named("p_big")  = p_big,
    Rcpp::Named("iter")   = iter,
    Rcpp::Named("ok")     = true
  );
}
