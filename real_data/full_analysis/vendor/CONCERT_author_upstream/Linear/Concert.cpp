#include <RcppArmadillo.h>
#include <iostream>
#include <stdio.h>
#include <Rcpp/Benchmark/Timer.h>

using namespace Rcpp;
using namespace std;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
vec sample_ind_C(vec n_vec, int num_K);
vec lasso_est(mat X, vec y, double lambda);
double Digamma(double a);

// [[Rcpp::export]]
List Concert(mat X, vec y, vec n_vec,
            vec tau, double q_k, double eta, double q_0,
            double threshold=1e-6, int max_iter=1000){
	int p=X.n_cols, K=n_vec.n_elem-1;
  vec ind=sample_ind_C(n_vec, K+1);
  double a0=2, b0=1;
	
	Timer timer;
	timer.step("start");
	
  //*************************************************************
	vec m_beta0, gamma_Z(p,fill::zeros), V_beta0(p,fill::ones);
  mat m_betak(p,K), gamma_I(p,K,fill::zeros), V_betak(p,K,fill::ones);
  vec E_inv_sigmasq(K+1,fill::zeros), E_log_sigmasq(K+1,fill::zeros);
  vec y_0, y_k;
  mat X_0, X_k;
	vec tmp_vec;
  uvec index;
	//*************************************************************
	// Initialization
  X_0 = X.rows(ind(0),ind(1)); y_0 = y.subvec(ind(0),ind(1));
  m_beta0=lasso_est(X_0, y_0, 0.1);  
  index=find(m_beta0!=0);
  gamma_Z(index)=linspace(1,1,index.n_elem);
  E_inv_sigmasq(0)=1.0/var(y_0);
  E_log_sigmasq(0)=log(var(y_0));

  for (int k=1; k<K+1; ++k){
    X_k = X.rows(ind(k)+1,ind(k+1)); y_k = y.subvec(ind(k)+1,ind(k+1));
    m_betak.col(k-1)=lasso_est(X_k, y_k, 0.1);
    index=find(m_betak.col(k-1)!=0);
    tmp_vec=linspace(0,0,p); tmp_vec(index)=linspace(1,1,index.n_elem);
    gamma_I.col(k-1)=tmp_vec;
    E_inv_sigmasq(k)=1.0/var(y_k);
    E_log_sigmasq(k)=log(var(y_k));
  };
  

  timer.step("initialization");

  //*************************************************************
  double diff_m, sum_var, sum_mean, logit;
  vec old_m, gamma_I_k, m_betak_k, a1_sigmasq(K+1), b1_sigmasq(K+1);
  mat X2, diag_X, offdiag_IX;
  uvec order_beta0, order_betak, index_beta0, index_betak;
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
      sum_var = E_inv_sigmasq(0)*(sum(X_0.col(j_ord).t()*X_0.col(j_ord))+pow(eta,-2));
      sum_mean = E_inv_sigmasq(0)*(sum(y_0.t()*X_0.col(j_ord))-
        sum(X_0.col(j_ord).t()*X_0.cols(index_beta0)*(gamma_Z(index_beta0)%m_beta0(index_beta0))));
      for (int k=1; k<K+1; ++k){
        X_k = X.rows(ind(k)+1,ind(k+1)); y_k = y.subvec(ind(k)+1,ind(k+1));
        gamma_I_k = gamma_I.col(k-1); m_betak_k = m_betak.col(k-1);
        sum_var = sum_var+
          E_inv_sigmasq(k)*(sum(X_k.col(j_ord).t()*X_k.col(j_ord))*gamma_I_k(j_ord)+
          pow(tau(k-1),-2)*(1-gamma_I_k(j_ord)));
        sum_mean = sum_mean+
          E_inv_sigmasq(k)*(gamma_I_k(j_ord)*(sum(y_k.t()*X_k.col(j_ord))-
                              sum(X_k.col(j_ord).t()*X_k.cols(index_beta0)*
                                (gamma_I_k(index_beta0)%gamma_Z(index_beta0)%m_beta0(index_beta0)+
                                (1-gamma_I_k(index_beta0))%m_betak_k(index_beta0))))+
          (1-gamma_I_k(j_ord))*pow(tau(k-1),-2)*m_betak_k(j_ord));
      };
      V_beta0(j_ord) = 1/sum_var;
      m_beta0(j_ord) = sum_mean/sum_var;

      logit = pow(m_beta0(j_ord),2)/(2*V_beta0(j_ord))-0.5*E_log_sigmasq(0)+log(q_0*sqrt(V_beta0(j_ord))/((1-q_0)*eta));
      gamma_Z(j_ord) = 1/(1+exp(-logit));
    };
    
    
    // betak: m_betak, V_betak
    // I_betak: gamma_I_k
    for (int k=1; k<K+1; ++k){
      X_k = X.rows(ind(k)+1,ind(k+1)); y_k = y.subvec(ind(k)+1,ind(k+1));
      m_betak_k = m_betak.col(k-1);
      order_betak = sort_index(abs(m_betak_k), "descending");
      for (int j = 0; j < p; ++j){
        uword j_ord = order_betak(j);
        index_betak = find(linspace(0,p-1,p) != j_ord);
        gamma_I_k = gamma_I.col(k-1);
        V_betak(j_ord,k-1) = 1/(E_inv_sigmasq(k)*(sum(X_k.col(j_ord).t()*X_k.col(j_ord))+pow(tau(k-1),-2)));
        m_betak(j_ord,k-1) = E_inv_sigmasq(k)*(sum(y_k.t()*X_k.col(j_ord))-
                              sum(X_k.col(j_ord).t()*X_k.cols(index_betak)*
                                (gamma_I_k(index_betak)%gamma_Z(index_betak)%m_beta0(index_betak)+
                                  (1-gamma_I_k(index_betak))%m_betak_k(index_betak)))+
                              pow(tau(k-1),-2)*gamma_Z(j_ord)*m_beta0(j_ord))*V_betak(j_ord,k-1);
        
        m_betak_k = m_betak.col(k-1);
        logit = -pow(m_betak_k(j_ord),2)/(2*V_betak(j_ord,k-1))+0.5*E_log_sigmasq(k)-log((1-q_k)*sqrt(V_betak(j_ord,k-1))/(q_k*tau(k-1)))+
          E_inv_sigmasq(k)*gamma_Z(j_ord)*(sum(y_k.t()*X_k.col(j_ord))*m_beta0(j_ord)-
                          0.5*(sum(X_k.col(j_ord).t()*X_k.col(j_ord))-pow(tau(k-1),-2))*
                                (pow(m_beta0(j_ord),2)+V_beta0(j_ord))-
                          m_beta0(j_ord)*sum(X_k.col(j_ord).t()*X_k.cols(index_betak)*
                                            (gamma_I_k(index_betak)%gamma_Z(index_betak)%m_beta0(index_betak)+
                                              (1-gamma_I_k(index_betak))%m_betak_k(index_betak))));
        gamma_I(j_ord,k-1) = 1/(1+exp(-logit));
      };
    };

    // sigmasq
    // sigmasq_0
    a1_sigmasq(0) = 0.5*(n_vec(0)+sum(gamma_Z))+a0;
    b1_sigmasq(0) = 0.5*(sum(y_0.t()*y_0)-2*sum(y_0.t()*X_0*(gamma_Z%m_beta0))+
      sum((gamma_Z%m_beta0).t()*X_0.t()*X_0*(gamma_Z%m_beta0))+
      sum((gamma_Z%m_beta0).t()*diagmat(X_0.t()*X_0)*((1-gamma_Z)%m_beta0))+
      sum(diagmat(X_0.t()*X_0*diagmat(gamma_Z%V_beta0))*linspace(1,1,p))+
      sum(gamma_Z.t()*(m_beta0%m_beta0+V_beta0))*pow(eta,-2))+b0;
    E_inv_sigmasq(0)=a1_sigmasq(0)/b1_sigmasq(0);
    E_log_sigmasq(0)=log(b1_sigmasq(0))-Digamma(a1_sigmasq(0));

    // sigmasq_k
    for (int k=1; k<K+1; ++k){
      X_k = X.rows(ind(k)+1,ind(k+1)); y_k = y.subvec(ind(k)+1,ind(k+1));
      gamma_I_k = gamma_I.col(k-1); m_betak_k = m_betak.col(k-1);
      X2 = X_k.t()*X_k; diag_X = diagmat(X2);
      offdiag_IX = diagmat(1-gamma_I_k)*X2*diagmat(gamma_I_k); offdiag_IX -= diagmat(offdiag_IX);
      a1_sigmasq(k) = 0.5*(n_vec(k)+sum(gamma_I_k))+a0;
      b1_sigmasq(k) = 0.5*(sum(y_k.t()*y_k)-2*sum(y_k.t()*X_k*((1-gamma_I_k)%m_betak_k+gamma_I_k%gamma_Z%m_beta0)))+
        0.5*(sum(((1-gamma_I_k)%m_betak_k).t()*X2*((1-gamma_I_k)%m_betak_k))+
             sum(((1-gamma_I_k)%m_betak_k).t()*diag_X*(gamma_I_k%m_betak_k))+
             sum(diagmat(X_k.t()*X_k*diagmat((1-gamma_I_k)%V_betak.col(k-1)))*linspace(1,1,p)))+
        sum(m_betak_k.t()*offdiag_IX*(gamma_Z%m_beta0))+
        0.5*(sum((gamma_I_k%gamma_Z%m_beta0).t()*X2*(gamma_I_k%gamma_Z%m_beta0))+
             sum((gamma_I_k%gamma_Z%m_beta0).t()*diag_X*(gamma_I_k%(1-gamma_Z)%m_beta0))+
             sum(diagmat(X_k.t()*X_k*diagmat(gamma_I_k%gamma_I_k%gamma_Z%V_beta0))*linspace(1,1,p)))+
        0.5*(sum((gamma_I_k%(1-gamma_I_k)).t()*diag_X*(gamma_Z%(m_beta0%m_beta0+V_beta0))))+
        0.5*(pow(tau(k-1),-2)*sum((1-gamma_I_k).t()*(m_betak_k%m_betak_k+V_betak.col(k-1)+
                                   (1-gamma_Z)%(m_beta0%m_beta0+V_beta0-2*m_betak_k%m_beta0))))+b0;
      E_inv_sigmasq(k)=a1_sigmasq(k)/b1_sigmasq(k);
      E_log_sigmasq(k)=log(b1_sigmasq(k))-Digamma(a1_sigmasq(k));
    }

    diff_m = 0;
    for (int j = 0; j < p; ++j){
      diff_m = diff_m+pow(m_beta0(j)-old_m(j),2);
    };

  };
	
  timer.step("CAVI");
	
	return List::create(
		_["m_beta0"] = m_beta0,
		_["gamma_Z"] = gamma_Z,
    _["m_betak"] = m_betak,
    _["gamma_I"] = gamma_I,
    _["V_beta0"] = V_beta0,
    _["V_betak"] = V_betak,
    _["timer"] = timer
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




