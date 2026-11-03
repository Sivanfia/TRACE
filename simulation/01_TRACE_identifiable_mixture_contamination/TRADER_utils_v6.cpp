// TRADER residual backend v5
// Exact Gibbs kernel with observation-space residual maintenance.
// Main beta/sigma2 sweep complexity: O(N p) rather than O(p^2).
// [[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
#include <RcppEigen.h>
#include <string>
#include <algorithm>


using namespace std;
using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::export]]
double test(double shape, double scale)
{
  return(R::rgamma(shape, scale));
}

// [[Rcpp::export]]
std::string TRADER_backend_version()
{
  return("residual_np_v5");
}

// [[Rcpp::export]]
List TRADER_mcmc(
    S4 TRADER,
    const Eigen::MatrixXd xtx,
    const Eigen::VectorXd xty,
    const double yty,
    List sample,
    const Eigen::SparseMatrix<double> LD_pairs,
    const string tau_method,
    const double tau2_lb_inv,
    const double tau2_ub_inv,
    const unsigned int n_chain,
    const NumericVector temp,
    const unsigned int n_mcmc, 
    const unsigned int thin, 
    const unsigned int seed)
{
  Function set_seed("set.seed");
  set_seed(seed);
  
  // Data & hyper-parameters
  const int p = TRADER.slot("p");
  const int N = TRADER.slot("n");
  const int K = TRADER.slot("K");
  
  const double np_half = (N+p) * 0.5;
  const double p1_half = (p+1) * 0.5;
  
  const Eigen::VectorXd alpha = TRADER.slot("alpha");
  const Eigen::MatrixXd w = TRADER.slot("w");
  const Eigen::MatrixXd x = TRADER.slot("x");
  const Eigen::VectorXd y = TRADER.slot("y");
  const double v = TRADER.slot("v");

  // The original interface supplies yty for compatibility.  The residual
  // backend evaluates ||y - X beta||^2 directly, so yty is no longer needed
  // in the main Gibbs sweep.
  (void)yty;
  
  // Initial samples
  int i = 0, ii;
  int j, j1, j2, k, l;
  int n;
  
  MatrixXd beta = sample[0];
  MatrixXd lambda2 = sample[1];
  MatrixXd eta = sample[2];
  VectorXd res = sample[3];
  VectorXd sigma2 = sample[4];
  VectorXd tau2 = sample[5];
  
  MatrixXd keep_beta = MatrixXd::Zero(p, n_mcmc);
  MatrixXd keep_lambda2 = MatrixXd::Zero(p, n_mcmc);
  MatrixXd keep_eta = MatrixXd::Zero(K, n_mcmc);
  VectorXd keep_sigma2 = VectorXd::Zero(n_mcmc);
  VectorXd keep_tau2 = VectorXd::Zero(n_mcmc);
  
  MatrixXd mu;
  if(K > 0)
  {
    mu = w * eta;
  } else {
    mu = MatrixXd::Zero(p, n_chain);
  }

  // Observation-space residuals.  For chain n,
  //   residual.col(n) = y - x * beta.col(n).
  // Maintaining this matrix changes the beta sweep from O(p^2) to O(N p):
  //   x_j' y - sum_{l != j} x_j' x_l beta_l
  // = x_j' residual + (x_j' x_j) beta_j.
  MatrixXd residual = MatrixXd::Zero(N, n_chain);
  for(n = 0; n < static_cast<int>(n_chain); ++n)
  {
    residual.col(n).noalias() = y - x * beta.col(n);
  }

  // Periodically rebuild the residual from beta to eliminate accumulated
  // floating-point drift.  This is deterministic bookkeeping and does not
  // change the MCMC transition kernel.
  const unsigned int residual_refresh_every = 100;
  unsigned int residual_refresh_count = 0;
  double max_residual_refresh_error = 0;

  // mcmc
  double u_j, tsigma2_j, prob;
  double a_k, b_k, rss = 0, rss_1 = 0, rss_2 = 0, ch;
  double beta_old, beta_new, beta1_old, beta2_old, beta1_new, beta2_new;
  double r, u, lb, ub, Fub, Flb, p_12 = 0, p_21;
  VectorXd ch_v(p), ch_e(K), ch_obs(N);
  Function order("order");
  IntegerVector order_beta;
  
  while(i < n_mcmc)
  {
    for(ii = 0; ii < thin; ++ii)
    {
      for(n = 0; n < n_chain; ++n)
      {
        for(j = 0; j < p; ++j)
        {
          //1. lambda2
          u = R::runif(0, lambda2(j,n) / (1+lambda2(j,n)));
          ub = (1-u) / u;
          r = pow(beta(j,n) - mu(j,n), 2) / 2 / sigma2(n) / tau2(n) / temp(n);
          if(r == 0)
          {
            lambda2(j,n) = 1e-10;
          } else {
            u = R::runif(0, 1);
            Fub = min(expm1(-r*ub), -1e-4);
            lambda2(j,n) = -r / log1p(u * Fub);
          }
          
          //2. beta: exact residual-based coordinate update
          tsigma2_j = 1 / (xtx(j,j) + 1 / lambda2(j,n) / tau2(n));

          beta_old = beta(j,n);
          u_j = x.col(j).dot(residual.col(n)) + xtx(j,j) * beta_old;
          if(mu(j,n) != 0)
          {
            u_j += mu(j,n) / lambda2(j,n) / tau2(n);
          }

          beta_new = R::rnorm(
            tsigma2_j * u_j,
            sqrt(tsigma2_j * sigma2(n) * temp(n))
          );
          beta(j,n) = beta_new;
          residual.col(n) -= x.col(j) * (beta_new - beta_old);
        }
        
        //3. eta
        if(K > 0)
        {
          for(k = 0; k < K; ++k)
          {
            res(n) += eta(k,n);
            
            a_k = b_k = 0;
            for(j = 0; j < p; ++j)
            {
              a_k += pow(w(j,k),2) / lambda2(j,n);
              rss = -beta(j,n);
              for(l = 0; l < K; ++l)
              {
                if(l != k)
                {
                  rss += eta(l,n) * w(j,l);
                }
              }
              rss *= w(j,k) / lambda2(j,n);
              b_k += rss;
            }
            b_k /= a_k;
            // cout << "a_k = " << a_k << ", b_k = " << b_k << endl;
            
            u = R::runif(0, exp(-a_k * pow(eta(k,n) + b_k,2) / 2 / sigma2(n) / 
              tau2(n) / temp(n)));
            ub = sqrt(-2 * sigma2(n) * tau2(n) * temp(n) * log(u) / a_k);
            lb = -ub;
            ub -= b_k;
            lb -= b_k;
            ub = min(ub / res(n), 1.0);
            lb = max(lb / res(n), 0.0);
            // cout << "lb = " << lb << ", ub = " << ub << endl;
            
            Flb = R::pbeta(lb, alpha(k), alpha(K), true, false);
            Fub = R::pbeta(ub, alpha(k), alpha(K), true, false);
            // cout << "Flb = " << Flb << ", Fub = " << Fub << endl;
            
            u = R::runif(Flb, Fub);
            eta(k,n) = R::qbeta(u, alpha(k), alpha(K), true, false) * res(n);
            
            res(n) -= eta(k,n);
          }
          mu.col(n) = w * eta.col(n);
        }
        
        //4. sigma2.  This equals
        // y'y + beta'X'X beta - 2 beta'X'y, but costs O(N).
        rss_1 = residual.col(n).squaredNorm();
        rss_2 = 0;
        for(j = 0; j < p; ++j)
        {
          rss_2 += pow(beta(j,n) - mu(j,n),2) / lambda2(j,n);
        }
        rss_2 /= tau2(n);
        // rss_2 = max(rss_2, 1e-4);
        
        rss = (rss_1 + rss_2) / temp(n);
        
        sigma2(n) = 1 / R::rgamma(np_half + v, 1 / (rss / 2 + v));
        
        //5. tau2
        if(tau_method == "Cauchy")
        {
          rss = 0;
          for(j = 0; j < p; ++j)
          {
            rss += pow(beta(j,n) - mu(j,n),2) / lambda2(j,n);
          }
          rss /= (2 * sigma2(n) * temp(n));
          // rss = max(rss, 5e-5);
          
          u = R::runif(0, tau2(n) / (1+tau2(n)));
          ub = min((1-u) / u, tau2_lb_inv);
          Fub = R::pgamma(ub, p1_half, 1 / rss, true, false);
          Flb = R::pgamma(tau2_ub_inv, p1_half, 1 / rss, true, false);
          
          u = R::runif(Flb, Fub);
          tau2(n) = 1 / R::qgamma(u, p1_half, 1 / rss, true, false);
        }
      }
    }
    
    keep_beta.col(i) = beta.col(0);
    keep_lambda2.col(i) = lambda2.col(0);
    keep_eta.col(i) = eta.col(0);
    keep_sigma2(i) = sigma2(0);
    keep_tau2(i) = tau2(0);
    
    ++i;
    
    // 6. Switch 2 variables after every 10 steps.
    if(((i + 1) % 20 == 10) && (i < n_mcmc))
    {
      for(n = 0; n < n_chain; ++n)
      {
        order_beta = order(beta.col(n).array().abs(), _["decreasing"] = true);
        order_beta = order_beta - 1;
        
        k = 50 + R::rgeom(0.3);
        if(k >= p) k = p-1;
        
        for(j = 0; j < k; ++j)
        {
          j1 = order_beta(j); //index of variable 1
          if(LD_pairs.col(j1).nonZeros() >= 1)
          {
            r = R::runif(0,1);
            u = 0;
            for(SparseMatrix<double>::InnerIterator it(LD_pairs,j1); it; ++it)
            {
              u += it.value();
              if(u > r)
              {
                j2 = it.row();
                p_12 = it.value();
                break;
              }
            }
            p_21 = LD_pairs.coeff(j1,j2);
            
            rss = 0;
            if(xtx(j1,j2) > 0)
            {
              for(l = 0; l < p; ++l)
              {
                if((l != j1) && (l != j2))
                {
                  rss += (xtx(l,j1) - xtx(l,j2)) * beta(l,n);
                }
              }
              rss *= (beta(j1,n) - beta(j2,n));
              
              rss += (xtx(j1,j1) - xtx(j2,j2)) * 
                (pow(beta(j1,n),2) - pow(beta(j2,n),2)) / 2 +
                (xty(j1) - xty(j2)) * (beta(j2,n) - beta(j1,n));
              rss += ((2*beta(j1,n) - mu(j1,n) - mu(j2,n)) * 
                (mu(j2,n) - mu(j1,n)) / lambda2(j1,n) +
                (2*beta(j2,n) - mu(j1,n) - mu(j2,n)) * (mu(j1,n) - mu(j2,n)) / 
                lambda2(j2,n)) / 2 / tau2(n);
            } else {
              for(l = 0; l < p; ++l)
              {
                if((l != j1) && (l != j2))
                {
                  rss += (xtx(l,j1) + xtx(l,j2)) * beta(l,n);
                }
              }
              rss *= (beta(j1,n) + beta(j2,n));
              
              rss += (xtx(j1,j1) - xtx(j2,j2)) * 
                (pow(beta(j1,n),2) - pow(beta(j2,n),2)) / 2 -
                (xty(j1) + xty(j2)) * (beta(j2,n) + beta(j1,n));
              rss -= ((2*beta(j1,n) - mu(j1,n) + mu(j2,n)) * 
                (mu(j2,n) + mu(j1,n)) / lambda2(j1,n) +
                (2*beta(j2,n) + mu(j1,n) - mu(j2,n)) * (mu(j1,n) + mu(j2,n)) / 
                lambda2(j2,n)) / 2 / tau2(n);
            }
            rss /= (sigma2(n) * temp(n));
            
            prob = exp(rss) * p_21 / p_12;
            if((prob > 1) || (R::runif(0,1) < prob))
            {
              beta1_old = beta(j1,n);
              beta2_old = beta(j2,n);
              if(xtx(j1,j2) > 0)
              {
                beta1_new = beta2_old;
                beta2_new = beta1_old;
              } else {
                beta1_new = -beta2_old;
                beta2_new = -beta1_old;
              }

              beta(j1,n) = beta1_new;
              beta(j2,n) = beta2_new;
              residual.col(n) -= x.col(j1) * (beta1_new - beta1_old);
              residual.col(n) -= x.col(j2) * (beta2_new - beta2_old);

              ch = lambda2(j1,n);
              lambda2(j1,n) = lambda2(j2,n);
              lambda2(j2,n) = ch;
            }
          }
        }
        
        //6.3. eta
        if(K > 0)
        {
          for(k = 0; k < K; ++k)
          {
            res(n) += eta(k,n);
            
            a_k = b_k = 0;
            for(j = 0; j < p; ++j)
            {
              a_k += pow(w(j,k),2) / lambda2(j,n);
              rss = -beta(j,n);
              for(l = 0; l < K; ++l)
              {
                if(l != k)
                {
                  rss += eta(l,n) * w(j,l);
                }
              }
              rss *= w(j,k) / lambda2(j,n);
              b_k += rss;
            }
            b_k /= a_k;
            // cout << "a_k = " << a_k << ", b_k = " << b_k << endl;
            
            u = R::runif(0, exp(-a_k * pow(eta(k,n) + b_k,2) / 2 / sigma2(n) / 
              tau2(n) / temp(n)));
            ub = sqrt(-2 * sigma2(n) * tau2(n) * temp(n) * log(u) / a_k);
            lb = -ub;
            ub -= b_k;
            lb -= b_k;
            ub = min(ub / res(n), 1.0);
            lb = max(lb / res(n), 0.0);
            // cout << "lb = " << lb << ", ub = " << ub << endl;
            
            Flb = R::pbeta(lb, alpha(k), alpha(K), true, false);
            Fub = R::pbeta(ub, alpha(k), alpha(K), true, false);
            // cout << "Flb = " << Flb << ", Fub = " << Fub << endl;
            
            u = R::runif(Flb, Fub);
            eta(k,n) = R::qbeta(u, alpha(k), alpha(K), true, false) * res(n);
            
            res(n) -= eta(k,n);
          }
          mu.col(n) = w * eta.col(n);
        }
        
        //6.4. sigma2
        rss_1 = residual.col(n).squaredNorm();
        rss_2 = 0;
        for(j = 0; j < p; ++j)
        {
          rss_2 += pow(beta(j,n) - mu(j,n),2) / lambda2(j,n);
        }
        rss_2 /= tau2(n);
        // rss_2 = max(rss_2, 1e-4);
        
        rss = (rss_1 + rss_2) / temp(n);
        
        sigma2(n) = 1 / R::rgamma(np_half + v, 1 / (rss / 2 + v));
        
        //6.5. tau2
        if(tau_method == "Cauchy")
        {
          rss = 0;
          for(j = 0; j < p; ++j)
          {
            rss += pow(beta(j,n) - mu(j,n),2) / lambda2(j,n);
          }
          rss /= (2 * sigma2(n) * temp(n));
          // rss = max(rss, 5e-5);
          
          u = R::runif(0, tau2(n) / (1+tau2(n)));
          ub = min((1-u) / u, tau2_lb_inv);
          Fub = R::pgamma(ub, p1_half, 1 / rss, true, false);
          Flb = R::pgamma(tau2_ub_inv, p1_half, 1 / rss, true, false);
          
          u = R::runif(Flb, Fub);
          tau2(n) = 1 / R::qgamma(u, p1_half, 1 / rss, true, false);
        }
      }
      
      keep_beta.col(i) = beta.col(0);
      keep_lambda2.col(i) = lambda2.col(0);
      keep_eta.col(i) = eta.col(0);
      keep_sigma2(i) = sigma2(0);
      keep_tau2(i) = tau2(0);
      
      ++i;
    } else if(((i + 1) % 20 == 0) && (i < n_mcmc)) {
      for(n = n_chain-2; n >= 0; --n)
      {
        rss_1 = rss_2 = 0;
        for(j = 0; j < p; ++j)
        {
          rss_1 += pow(beta(j,n) - mu(j,n), 2) / lambda2(j,n);
          rss_2 += pow(beta(j,n+1) - mu(j,n+1), 2) / lambda2(j,n+1);
        }
        rss_1 /= tau2(n);
        rss_2 /= tau2(n+1);
        
        rss_1 += residual.col(n).squaredNorm();
        rss_2 += residual.col(n+1).squaredNorm();
        rss_1 /= sigma2(n);
        rss_2 /= sigma2(n+1);
        
        rss = (rss_1 - rss_2) * (1 / temp(n) - 1 / temp(n+1)) / 2;
        prob = exp(rss);
        if((prob > 1) || (R::runif(0,1) < prob))
        {
          // cout << "Step " << i << ", chain " << n << " and chain " << n+1 <<
          //   " switch, prob = " << prob << "." << endl;
          
          ch_v = beta.col(n);
          beta.col(n) = beta.col(n+1);
          beta.col(n+1) = ch_v;

          ch_obs = residual.col(n);
          residual.col(n) = residual.col(n+1);
          residual.col(n+1) = ch_obs;
          
          ch_v = lambda2.col(n);
          lambda2.col(n) = lambda2.col(n+1);
          lambda2.col(n+1) = ch_v;
          
          if(K > 0)
          {
            ch_e = eta.col(n);
            eta.col(n) = eta.col(n+1);
            eta.col(n+1) = ch_e;
            
            ch = res(n);
            res(n) = res(n+1);
            res(n+1) = ch;
            
            ch_v = mu.col(n);
            mu.col(n) = mu.col(n+1);
            mu.col(n+1) = ch_v;
          }
          
          ch = sigma2(n);
          sigma2(n) = sigma2(n+1);
          sigma2(n+1) = ch;
          
          ch = tau2(n);
          tau2(n) = tau2(n+1);
          tau2(n+1) = ch;
        }
      }
      
      keep_beta.col(i) = beta.col(0);
      keep_lambda2.col(i) = lambda2.col(0);
      keep_eta.col(i) = eta.col(0);
      keep_sigma2(i) = sigma2(0);
      keep_tau2(i) = tau2(0);
      
      ++i;
    }

    if((i % residual_refresh_every) == 0)
    {
      for(n = 0; n < static_cast<int>(n_chain); ++n)
      {
        VectorXd exact_residual = y - x * beta.col(n);
        const double refresh_error =
          (residual.col(n) - exact_residual).cwiseAbs().maxCoeff();
        max_residual_refresh_error = std::max(
          max_residual_refresh_error, refresh_error
        );
        residual.col(n) = exact_residual;
      }
      ++residual_refresh_count;
    }
  }

  // Final consistency refresh and diagnostic.
  for(n = 0; n < static_cast<int>(n_chain); ++n)
  {
    VectorXd exact_residual = y - x * beta.col(n);
    const double refresh_error =
      (residual.col(n) - exact_residual).cwiseAbs().maxCoeff();
    max_residual_refresh_error = std::max(
      max_residual_refresh_error, refresh_error
    );
    residual.col(n) = exact_residual;
  }
  ++residual_refresh_count;
  
  if(tau_method == "constant")
  {
    return(List::create(Named("beta") = keep_beta,
                        Named("lambda2") = keep_lambda2,
                        Named("eta") = keep_eta,
                        Named("sigma2") = keep_sigma2,
                        Named("last_sample") = List::create(
                          Named("beta") = beta,
                          Named("lambda2") = lambda2,
                          Named("eta") = eta,
                          Named("res") = res,
                          Named("sigma2") = sigma2),
                        Named("backend") = "residual_np_v5",
                        Named("residual_refresh_count") = residual_refresh_count,
                        Named("max_residual_refresh_error") = max_residual_refresh_error));
  } else {
    return(List::create(Named("beta") = keep_beta,
                        Named("lambda2") = keep_lambda2,
                        Named("eta") = keep_eta,
                        Named("sigma2") = keep_sigma2,
                        Named("tau2") = keep_tau2,
                        Named("last_sample") = List::create(
                          Named("beta") = beta,
                          Named("lambda2") = lambda2,
                          Named("eta") = eta,
                          Named("res") = res,
                          Named("sigma2") = sigma2,
                          Named("tau2") = tau2),
                        Named("backend") = "residual_np_v5",
                        Named("residual_refresh_count") = residual_refresh_count,
                        Named("max_residual_refresh_error") = max_residual_refresh_error));
  }
}
