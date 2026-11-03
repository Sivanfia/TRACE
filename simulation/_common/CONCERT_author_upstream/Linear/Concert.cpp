#include <RcppArmadillo.h>
#include <iostream>
#include <stdio.h>
#include <Rcpp/Benchmark/Timer.h>
#include <chrono>
#include <vector>

using namespace Rcpp;
using namespace std;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
vec sample_ind_C(vec n_vec, int num_K);
vec lasso_est(mat X, vec y, double lambda);
double Digamma(double a);

// Optimized, algebraically equivalent implementation of the upstream linear
// CONCERT CAVI updates.  The statistical updates and Gauss-Seidel coordinate
// order are unchanged.  The speedup comes from:
//   1) slicing target/source matrices only once;
//   2) caching column squared norms and X^T y;
//   3) maintaining y - X E(beta) residuals during coordinate sweeps, so
//      X_j^T X_{-j} beta_{-j} is never formed explicitly;
//   4) evaluating variance updates with matrix-vector products instead of
//      repeatedly forming X^T X.
// [[Rcpp::export]]
List Concert(mat X, vec y, vec n_vec,
            vec tau, double q_k, double eta, double q_0,
            double threshold=1e-6, int max_iter=1000){
  const int p = X.n_cols;
  const int K = n_vec.n_elem - 1;
  const vec ind = sample_ind_C(n_vec, K + 1);
  const double a0 = 2.0, b0 = 1.0;

  Timer timer;
  timer.step("start");
  const auto total_start = std::chrono::steady_clock::now();

  // ------------------------------------------------------------------
  // Partition the stacked data exactly once and cache O(p) cross-products.
  // ------------------------------------------------------------------
  const auto precompute_start = std::chrono::steady_clock::now();

  mat X_0 = X.rows((uword)ind(0), (uword)ind(1));
  vec y_0 = y.subvec((uword)ind(0), (uword)ind(1));

  std::vector<mat> X_source(K);
  std::vector<vec> y_source(K);
  std::vector<vec> col_sq_source(K);
  std::vector<vec> Xty_source(K);
  std::vector<double> yty_source(K);

  const vec col_sq_0 = sum(square(X_0), 0).t();
  const double yty_0 = dot(y_0, y_0);

  for (int k = 0; k < K; ++k) {
    const uword first = (uword)ind(k + 1) + 1;
    const uword last  = (uword)ind(k + 2);
    X_source[k] = X.rows(first, last);
    y_source[k] = y.subvec(first, last);
    col_sq_source[k] = sum(square(X_source[k]), 0).t();
    Xty_source[k] = X_source[k].t() * y_source[k];
    yty_source[k] = dot(y_source[k], y_source[k]);
  }

  const auto precompute_end = std::chrono::steady_clock::now();
  timer.step("precompute");

  // ------------------------------------------------------------------
  // Initialization: identical lasso initializers and variance moments.
  // ------------------------------------------------------------------
  const auto initialization_start = std::chrono::steady_clock::now();

  vec m_beta0, gamma_Z(p, fill::zeros), V_beta0(p, fill::ones);
  mat m_betak(p, K), gamma_I(p, K, fill::zeros), V_betak(p, K, fill::ones);
  vec E_inv_sigmasq(K + 1, fill::zeros), E_log_sigmasq(K + 1, fill::zeros);
  vec tmp_vec;
  uvec index;

  m_beta0 = lasso_est(X_0, y_0, 0.1);
  index = find(m_beta0 != 0);
  if (index.n_elem > 0) gamma_Z(index).ones();
  E_inv_sigmasq(0) = 1.0 / var(y_0);
  E_log_sigmasq(0) = log(var(y_0));

  for (int k = 0; k < K; ++k) {
    m_betak.col(k) = lasso_est(X_source[k], y_source[k], 0.1);
    index = find(m_betak.col(k) != 0);
    tmp_vec.zeros(p);
    if (index.n_elem > 0) tmp_vec(index).ones();
    gamma_I.col(k) = tmp_vec;
    E_inv_sigmasq(k + 1) = 1.0 / var(y_source[k]);
    E_log_sigmasq(k + 1) = log(var(y_source[k]));
  }

  const auto initialization_end = std::chrono::steady_clock::now();
  timer.step("initialization");

  // ------------------------------------------------------------------
  // Residual caches.  These are kept exactly consistent with the current
  // posterior-mean coefficient combinations after every coordinate update.
  // ------------------------------------------------------------------
  vec coef_target = gamma_Z % m_beta0;
  vec residual_target = y_0 - X_0 * coef_target;

  std::vector<vec> coef_source(K);
  std::vector<vec> residual_source(K);
  for (int k = 0; k < K; ++k) {
    const vec gamma_k = gamma_I.col(k);
    coef_source[k] = gamma_k % gamma_Z % m_beta0 +
      (1.0 - gamma_k) % m_betak.col(k);
    residual_source[k] = y_source[k] - X_source[k] * coef_source[k];
  }

  // ------------------------------------------------------------------
  // CAVI: same update order as the upstream implementation.
  // ------------------------------------------------------------------
  const auto cavi_start = std::chrono::steady_clock::now();

  double diff_m = 10.0;
  double sum_var, sum_mean, logit;
  vec old_m, a1_sigmasq(K + 1), b1_sigmasq(K + 1);
  uvec order_beta0, order_betak;
  int iter_num = 0;

  while (diff_m > threshold && iter_num < max_iter) {
    old_m = m_beta0;
    ++iter_num;

    // Refresh once per sweep to prevent accumulated floating-point drift in
    // the rank-one residual updates.
    residual_target = y_0 - X_0 * coef_target;
    for (int k = 0; k < K; ++k) {
      residual_source[k] = y_source[k] - X_source[k] * coef_source[k];
    }

    // beta0: m_beta0, V_beta0; Z: gamma_Z
    order_beta0 = sort_index(abs(m_beta0), "descending");
    for (int j = 0; j < p; ++j) {
      const uword j_ord = order_beta0(j);

      // X_j^T(y - X_{-j} b_{-j}) = X_j^T residual + ||X_j||^2 b_j.
      const double partial_target =
        dot(X_0.col(j_ord), residual_target) +
        col_sq_0(j_ord) * coef_target(j_ord);

      sum_var = E_inv_sigmasq(0) *
        (col_sq_0(j_ord) + pow(eta, -2));
      sum_mean = E_inv_sigmasq(0) * partial_target;

      for (int k = 0; k < K; ++k) {
        const double gamma_kj = gamma_I(j_ord, k);
        const double partial_source =
          dot(X_source[k].col(j_ord), residual_source[k]) +
          col_sq_source[k](j_ord) * coef_source[k](j_ord);

        sum_var += E_inv_sigmasq(k + 1) *
          (col_sq_source[k](j_ord) * gamma_kj +
           pow(tau(k), -2) * (1.0 - gamma_kj));

        sum_mean += E_inv_sigmasq(k + 1) *
          (gamma_kj * partial_source +
           (1.0 - gamma_kj) * pow(tau(k), -2) * m_betak(j_ord, k));
      }

      const double old_coef_target = coef_target(j_ord);

      V_beta0(j_ord) = 1.0 / sum_var;
      m_beta0(j_ord) = sum_mean / sum_var;

      logit = pow(m_beta0(j_ord), 2) / (2.0 * V_beta0(j_ord)) -
        0.5 * E_log_sigmasq(0) +
        log(q_0 * sqrt(V_beta0(j_ord)) / ((1.0 - q_0) * eta));
      gamma_Z(j_ord) = 1.0 / (1.0 + exp(-logit));

      coef_target(j_ord) = gamma_Z(j_ord) * m_beta0(j_ord);
      const double delta_target = coef_target(j_ord) - old_coef_target;
      if (delta_target != 0.0) {
        residual_target -= X_0.col(j_ord) * delta_target;
      }

      for (int k = 0; k < K; ++k) {
        const double gamma_kj = gamma_I(j_ord, k);
        const double delta_source = gamma_kj * delta_target;
        coef_source[k](j_ord) += delta_source;
        if (delta_source != 0.0) {
          residual_source[k] -= X_source[k].col(j_ord) * delta_source;
        }
      }
    }

    // Refresh source residuals before the source-coordinate sweeps.  This is
    // only O(sum_k n_k p) and keeps the residual representation numerically
    // aligned with the explicit upstream inner products.
    for (int k = 0; k < K; ++k) {
      residual_source[k] = y_source[k] - X_source[k] * coef_source[k];
    }

    // betak: m_betak, V_betak; I_betak: gamma_I
    for (int k = 0; k < K; ++k) {
      order_betak = sort_index(abs(m_betak.col(k)), "descending");

      for (int j = 0; j < p; ++j) {
        const uword j_ord = order_betak(j);
        const double old_coef = coef_source[k](j_ord);

        // Same excluded-coordinate inner product used in both the m_betak and
        // gamma_I updates.  The updated m_betak(j) does not affect it because
        // coordinate j is excluded in the upstream expression.
        const double partial_source =
          dot(X_source[k].col(j_ord), residual_source[k]) +
          col_sq_source[k](j_ord) * old_coef;
        const double cross_excluding_j =
          Xty_source[k](j_ord) - partial_source;

        V_betak(j_ord, k) = 1.0 /
          (E_inv_sigmasq(k + 1) *
           (col_sq_source[k](j_ord) + pow(tau(k), -2)));

        m_betak(j_ord, k) = E_inv_sigmasq(k + 1) *
          (partial_source +
           pow(tau(k), -2) * gamma_Z(j_ord) * m_beta0(j_ord)) *
          V_betak(j_ord, k);

        logit = -pow(m_betak(j_ord, k), 2) /
          (2.0 * V_betak(j_ord, k)) +
          0.5 * E_log_sigmasq(k + 1) -
          log((1.0 - q_k) * sqrt(V_betak(j_ord, k)) /
              (q_k * tau(k))) +
          E_inv_sigmasq(k + 1) * gamma_Z(j_ord) *
          (Xty_source[k](j_ord) * m_beta0(j_ord) -
           0.5 * (col_sq_source[k](j_ord) - pow(tau(k), -2)) *
             (pow(m_beta0(j_ord), 2) + V_beta0(j_ord)) -
           m_beta0(j_ord) * cross_excluding_j);

        gamma_I(j_ord, k) = 1.0 / (1.0 + exp(-logit));

        coef_source[k](j_ord) =
          gamma_I(j_ord, k) * gamma_Z(j_ord) * m_beta0(j_ord) +
          (1.0 - gamma_I(j_ord, k)) * m_betak(j_ord, k);
        const double delta_source = coef_source[k](j_ord) - old_coef;
        if (delta_source != 0.0) {
          residual_source[k] -= X_source[k].col(j_ord) * delta_source;
        }
      }
    }

    // sigmasq_0: same expression as upstream, evaluated without X^T X.
    a1_sigmasq(0) = 0.5 * (n_vec(0) + sum(gamma_Z)) + a0;
    const vec fitted_target = X_0 * coef_target;
    b1_sigmasq(0) = 0.5 *
      (yty_0 - 2.0 * dot(y_0, fitted_target) + dot(fitted_target, fitted_target) +
       dot(coef_target,
           col_sq_0 % ((1.0 - gamma_Z) % m_beta0)) +
       dot(col_sq_0, gamma_Z % V_beta0) +
       sum(gamma_Z % (m_beta0 % m_beta0 + V_beta0)) * pow(eta, -2)) + b0;
    E_inv_sigmasq(0) = a1_sigmasq(0) / b1_sigmasq(0);
    E_log_sigmasq(0) = log(b1_sigmasq(0)) - Digamma(a1_sigmasq(0));

    // sigmasq_k: same upstream terms, but quadratic forms are computed as
    // ||Xv||^2 and cross-products of Xv rather than by constructing X^T X.
    for (int k = 0; k < K; ++k) {
      const vec gamma_k = gamma_I.col(k);
      const vec m_k = m_betak.col(k);
      const vec u = (1.0 - gamma_k) % m_k;
      const vec v = gamma_k % gamma_Z % m_beta0;
      const vec Xu = X_source[k] * u;
      const vec Xv = X_source[k] * v;
      const vec fitted_source = Xu + Xv;
      const vec &d = col_sq_source[k];

      const double offdiag_cross =
        dot(Xu, Xv) - dot(d, u % v);

      a1_sigmasq(k + 1) =
        0.5 * (n_vec(k + 1) + sum(gamma_k)) + a0;

      b1_sigmasq(k + 1) =
        0.5 * (yty_source[k] - 2.0 * dot(y_source[k], fitted_source)) +
        0.5 *
          (dot(Xu, Xu) +
           dot(u, d % (gamma_k % m_k)) +
           dot(d, (1.0 - gamma_k) % V_betak.col(k))) +
        offdiag_cross +
        0.5 *
          (dot(Xv, Xv) +
           dot(v, d % (gamma_k % (1.0 - gamma_Z) % m_beta0)) +
           dot(d, gamma_k % gamma_k % gamma_Z % V_beta0)) +
        0.5 * dot(gamma_k % (1.0 - gamma_k),
                  d % (gamma_Z % (m_beta0 % m_beta0 + V_beta0))) +
        0.5 * pow(tau(k), -2) *
          sum((1.0 - gamma_k) %
              (m_k % m_k + V_betak.col(k) +
               (1.0 - gamma_Z) %
                 (m_beta0 % m_beta0 + V_beta0 - 2.0 * m_k % m_beta0))) +
        b0;

      E_inv_sigmasq(k + 1) = a1_sigmasq(k + 1) / b1_sigmasq(k + 1);
      E_log_sigmasq(k + 1) =
        log(b1_sigmasq(k + 1)) - Digamma(a1_sigmasq(k + 1));
    }

    diff_m = accu(square(m_beta0 - old_m));
  }

  const auto cavi_end = std::chrono::steady_clock::now();
  timer.step("CAVI");
  const auto total_end = cavi_end;

  const bool converged = (diff_m <= threshold);
  const bool hit_max_iter = (!converged && iter_num >= max_iter);

  const double precompute_time_sec =
    std::chrono::duration<double>(precompute_end - precompute_start).count();
  const double initialization_time_sec =
    std::chrono::duration<double>(initialization_end - initialization_start).count();
  const double cavi_time_sec =
    std::chrono::duration<double>(cavi_end - cavi_start).count();
  const double cpp_total_time_sec =
    std::chrono::duration<double>(total_end - total_start).count();

  return List::create(
    _["m_beta0"] = m_beta0,
    _["gamma_Z"] = gamma_Z,
    _["m_betak"] = m_betak,
    _["gamma_I"] = gamma_I,
    _["V_beta0"] = V_beta0,
    _["V_betak"] = V_betak,
    _["timer"] = timer,
    _["n_iter"] = iter_num,
    _["final_diff_m"] = diff_m,
    _["converged"] = converged,
    _["hit_max_iter"] = hit_max_iter,
    _["precompute_time_sec"] = precompute_time_sec,
    _["initialization_time_sec"] = initialization_time_sec,
    _["cavi_time_sec"] = cavi_time_sec,
    _["cpp_total_time_sec"] = cpp_total_time_sec,
    _["optimized_equivalent_backend"] = true
  );
};


// sample_ind
vec sample_ind_C(vec n_vec, int num_K) {
  // The indices of different source data sets
  vec indices(num_K+1);

  indices(0)=0;
  for (int k = 0; k < num_K; ++k) {
    if (k == 0) {
      indices(k+1) = n_vec(k)-1;
    } else if (k == 1) {
      indices(k+1) = sum(n_vec.subvec(0,1))-1;
    } else {
      indices(k+1) = sum(n_vec.subvec(0,k))-1;
    }
  }
  return indices;
};

// lasso_est
vec lasso_est(mat X, vec y, double lambda){
  Function lasso_est_R("lasso_est_R", Environment::global_env());
  NumericVector tmp = lasso_est_R(X, y, lambda);
  vec res = tmp;
  return res;
};

// Digamma
double Digamma(double a){
  return R::digamma(a);
};


// [[Rcpp::export]]
List solo(mat X, vec y, double eta, double q_0,
                  double threshold=1e-6, double max_iter=1000){
  int p=X.n_cols, n0=X.n_rows;
  double a0=2, b0=1;

  Timer timer;
  timer.step("start");

  //*************************************************************
  vec m_beta0, gamma_Z(p,fill::zeros), V_beta0(p,fill::ones);
  double E_inv_sigmasq, E_log_sigmasq;
  vec tmp_vec;
  uvec index;
  //*************************************************************
  // Initialization
  m_beta0=lasso_est(X, y, 0.1);  
  index=find(m_beta0!=0);
  gamma_Z(index)=linspace(1,1,index.n_elem);
  E_inv_sigmasq=1.0/var(y);
  E_log_sigmasq=log(var(y));

  timer.step("initialization");

  //*************************************************************
  double diff_m, sum_var, sum_mean, logit;
  vec old_m;
  double a1_sigmasq, b1_sigmasq;
  uvec order_beta0, index_beta0;
  //*************************************************************
  // CAVI
  diff_m=10;
  int iter_num = 0;
  
  while (diff_m > threshold){
    old_m = m_beta0;
    iter_num = iter_num+1;
    if (iter_num > max_iter){
      break;
    };
    
    // beta0: m_beta0, V_beta0
    // Z: gamma_Z
    order_beta0 = sort_index(abs(m_beta0), "descending");
    for (int j = 0; j < p; ++j){
      uword j_ord = order_beta0(j);
      index_beta0 = find(linspace(0,p-1,p) != j_ord);
      sum_var = E_inv_sigmasq*(sum(X.col(j_ord).t()*X.col(j_ord))+pow(eta,-2));
      sum_mean = E_inv_sigmasq*(sum(y.t()*X.col(j_ord))-
        sum(X.col(j_ord).t()*X.cols(index_beta0)*(gamma_Z(index_beta0)%m_beta0(index_beta0))));
      V_beta0(j_ord) = 1/sum_var;
      m_beta0(j_ord) = sum_mean/sum_var;

      logit = pow(m_beta0(j_ord),2)/(2*V_beta0(j_ord))-
        0.5*E_log_sigmasq+log(q_0*sqrt(V_beta0(j_ord))/((1-q_0)*eta));
      gamma_Z(j_ord) = 1/(1+exp(-logit));
    };
    
    // sigmasq_0
    a1_sigmasq = 0.5*(n0+sum(gamma_Z))+a0;
    b1_sigmasq = 0.5*(sum(y.t()*y)-2*sum(y.t()*X*(gamma_Z%m_beta0))+
      sum((gamma_Z%m_beta0).t()*X.t()*X*(gamma_Z%m_beta0))+
      sum((gamma_Z%m_beta0).t()*diagmat(X.t()*X)*((1-gamma_Z)%m_beta0))+
      sum(diagmat(X.t()*X*diagmat(gamma_Z%V_beta0))*linspace(1,1,p))+
      sum(gamma_Z.t()*(m_beta0%m_beta0+V_beta0))*pow(eta,-2))+b0;
    E_inv_sigmasq=a1_sigmasq/b1_sigmasq;
    E_log_sigmasq=log(b1_sigmasq)-Digamma(a1_sigmasq);

    diff_m = 0;
    for (int j = 0; j < p; ++j){
      diff_m = diff_m+pow(m_beta0(j)-old_m(j),2);
    };
  };

  timer.step("CAVI");

  return List::create(
    _["m_beta0"] = m_beta0,
    _["gamma_Z"] = gamma_Z,
    _["V_beta0"] = V_beta0,
    _["timer"] = timer
  );
};




