library(Rcpp)
library(here)

sourceCpp(here("Linear","Concert.cpp"))
#' Load the Concert() function from Rcpp
#'
#' @description
#' This function wraps the compiled C++ implementation of the CONCERT algorithm
#' through Rcpp and provides an R interface.
#'
#' @param X Predictor matrix combining the target and all sources by rows.
#' @param Y Response vector combining the target and all sources.
#' @param n_vec A list specifying sample sizes for the target and each source.
#' @param tau Slab standard deviation in the conditional spike-and-slab prior for \code{beta_k}.
#' @param q_k Prior transferable probability for each source.
#' @param eta Slab standard deviation in the spike-and-slab prior for \code{beta_0}.
#' @param q_0 Prior inclusion probability for the target.
#' @param threshold Convergence threshold for CAVI.
#' @param max_iter Maximum number of iterations for CAVI.
#'
#' @return A list containing the following elements:
#' \itemize{
#'   \item \code{m_beta0}: Estimated slab mean for \code{beta_0}.
#'   \item \code{gamma_Z}: Estimated variational posterior inclusion probability for the target.
#'   \item \code{m_betak}: Estimated slab mean for \code{beta_k}.
#'   \item \code{gamma_I}: Estimated variational posterior transferable probability for each source.
#'   \item \code{V_beta0}: Estimated slab variance for \code{beta_0}.
#'   \item \code{V_betak}: Estimated slab variance for \code{beta_k}.
#' }


sample_ind <- function(n_vec, num_K){
  # The indices of different source data sets
  ind <- list()
  for (k in 1:num_K){
    if (k == 1){
      ind[[1]] <- 1:n_vec[1]
    } else{
      ind[[k]] <- (sum(n_vec[1:(k-1)])+1): sum(n_vec[1:k])
    }
  }
  return(ind)
}

eval_pre <- function(est,X_test=NULL,y_test=NULL){
  pred_err<-NA
  if(!is.null(X_test)& !is.null(y_test)){
    pred_err<- sqrt(mean((y_test-X_test%*%est)^2))
  }
  return(list(pred_err=pred_err))
}

perf <- function(theta, theta0){
  Itheta0 <- as.numeric(theta0 != 0)
  Itheta <- as.numeric(theta != 0)
  
  TP <- sum(Itheta[which(Itheta0 == 1)])/sum(Itheta0)
  FP <- sum(Itheta[which(Itheta0 == 0)])/max(1,sum(Itheta))
  
  return(c(TP, FP))
}

lasso_est_R <- function(X, y, lambda){
  beta <- as.numeric(glmnet::glmnet(X, y, lambda=lambda)$beta)
  return(beta)
}

