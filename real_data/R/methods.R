fit_trace_real <- function(dat,cfg) {
  stopifnot(is.null(dat$beta_true),is.null(dat$X_test),is.null(dat$y_test))
  #rho<-min(.1,max(1/dat$n,1/sqrt(dat$p)))
  rho <- cfg$rho_z
  prior<-list(A_z=50*rho,B_z=50*(1-rho),rho_z_mean=rho,
    A_c=50*cfg$rho_c,B_c=50*(1-cfg$rho_c),rho_c_mean=cfg$rho_c)
  fit_borrow_model_fixedhyper(dat,prior,alpha_method='old',init_method='full',
    init_alpha='uniform',init_gamma='theory_zero',sigma_beta_sq=cfg$sigma_beta_sq,
    max_iter=cfg$max_iter,tol=cfg$tol,min_iter=8L,theory_initialization=TRUE,
    target_init_until_converged=FALSE,target_init_steps=cfg$target_init_steps,
    preliminary_alpha_update='huber',huber_tau_rule='fixed',
    borrow_update_schedule='paper_block',alpha_update_strategy='em_qp',
    alpha_update_every=1L,monitor_elbo=TRUE,warn_elbo_drop=FALSE,verbose=FALSE)
}
fit_acvb_real <- function(dat,cfg) {
  # The simultaneous ACVB map can alternate between almost-full/empty models.
  # Relax only after detecting a two-cycle; convergence requires the undamped
  # map residual as well, so a small damped step cannot imply false convergence.
  #rho<-min(.1,max(1/dat$n,1/sqrt(dat$p)))
  rho <- cfg$rho_z
  w<-rep(1,dat$p);previous<-NULL;damped<-FALSE;converged<-FALSE
  for(iter in seq_len(cfg$max_iter)) {
    proposal<-.acvb_update_woodbury(w,dat,cfg$sigma_beta_sq,1,1,qlogis(rho))$omega
    proposal<-clamp(proposal,.Machine$double.eps,1-.Machine$double.eps)
    raw_residual<-max(abs(proposal-w))
    if(isTRUE(cfg$acvb_adaptive_damping)&&!is.null(previous)&&
       raw_residual>1e-4&&max(abs(proposal-previous))<max(1e-6,1e-3*raw_residual))damped<-TRUE
    new<-if(damped).5*w+.5*proposal else proposal
    change<-max(abs(new-w));previous<-w;w<-new
    if(iter>2L&&change<cfg$tol&&raw_residual<cfg$tol){converged<-TRUE;break}
  }
  braw<-w*as.numeric(.acvb_solve_woodbury(w,w*dat$Xty,dat,cfg$sigma_beta_sq))
  selected<-which(w>.5);b<-rep(0,dat$p)
  if(length(selected))b[selected]<-as.numeric(solve_reg(crossprod(dat$X[,selected,drop=FALSE]),
    crossprod(dat$X[,selected,drop=FALSE],dat$y)))
  list(method='ACVB',omega=w,beta_raw=braw,beta_hat=b,pi_prior=rho,
    sigma_beta_sq=cfg$sigma_beta_sq,n_iter=iter,converged=converged,
    convergence_reason=if(converged)'raw_map_and_step_tolerance'else'max_iter_reached',
    adaptive_damping_used=damped,raw_map_residual=raw_residual,step_residual=change)
}
fit_real_method <- function(method,dat,cfg,seed) {
  set.seed(seed)
  if(method=='TRACE')return(fit_trace_real(dat,cfg))
  if(method=='ACVB')return(fit_acvb_real(dat,cfg))
  if(method=='TRADER')return(fit_trader_baseline(dat,source_input='B',
    estimator='spikeslab',profile='current',allow_estimator_fallback=FALSE,
    corr_threshold=cfg$trader_corr_threshold,n_mcmc=cfg$trader_n_mcmc,
    n_burnin=cfg$trader_burnin,thin=cfg$trader_thin,
    sign_threshold=cfg$trader_sign_threshold,seed=seed,verbose=FALSE))
  if(method=='Trans-Lasso') {
    if(cfg$trans_backend=='simulation_ranked')return(fit_trans_lasso_ranked_fallback(dat))
    # The aggregation folds are entirely within outer training donors.
    n0<-dat$n;u<-unique(dat$target_donor);half<-sample(u,floor(length(u)/2))
    I1<-which(dat$target_donor%in%half);I2<-setdiff(seq_len(n0),I1)
    f1<-Trans.lasso(dat$X_all,dat$y_all,dat$n_vec_trans,I.til=I1,l1=TRUE)
    f2<-Trans.lasso(dat$X_all,dat$y_all,dat$n_vec_trans,I.til=I2,l1=TRUE)
    return(list(method='Trans_Lasso_author_twofold',beta_hat=(f1$beta.hat+f2$beta.hat)/2,
      detail=list(f1,f2),converged=NA,n_iter=NA_integer_))
  }
  if(method=='CONCERT')return(fit_concert_real(dat,cfg,seed=seed))
  stop('Unknown method ',method)
}
get_selected_support_real <- function(fit, p) {
  
  ## TRACE / ACVB
  if (!is.null(fit$omega)) {
    return(which(fit$omega > 0.5))
  }
  
  ## TRADER
  if (!is.null(fit$selected)) {
    return(which(as.logical(fit$selected)))
  }
  
  ## Trans-Lasso
  if (!is.null(fit$beta_hat)) {
    return(which(abs(fit$beta_hat) > 1e-8))
  }
  
  integer(0)
}


ridge_refit_beta <- function(X, y, lambda) {
  
  k <- ncol(X)
  
  if (k == 0L) {
    return(numeric(0))
  }
  
  A <- crossprod(X) + diag(lambda, k)
  rhs <- crossprod(X, y)
  
  as.numeric(solve_spd(A, rhs))
}


post_selection_ridge_cv <- function(dat, fit, cfg, seed = 1L) {
  
  selected <- get_selected_support_real(fit, dat$p)
  
  beta_hat <- rep(0, dat$p)
  selected_indicator <- rep(FALSE, dat$p)
  
  if (length(selected) == 0L) {
    return(list(
      beta = beta_hat,
      selected = selected_indicator,
      lambda = NA_real_,
      cv_mse = NA_real_
    ))
  }
  
  selected_indicator[selected] <- TRUE
  
  X <- dat$X[, selected, drop = FALSE]
  y <- dat$y
  
  ## donor-level CV
  donor <- dat$target_donor
  
  if (is.null(donor)) {
    donor <- seq_len(dat$n)
  }
  
  unique_donor <- unique(donor)
  
  Kfold <- min(
    cfg$refit_folds,
    length(unique_donor)
  )
  
  set.seed(seed)
  
  unique_donor <- sample(unique_donor)
  
  donor_fold <- setNames(
    rep(seq_len(Kfold), length.out = length(unique_donor)),
    as.character(unique_donor)
  )
  
  fold_id <- as.integer(
    donor_fold[as.character(donor)]
  )
  
  ## X is standardized, so use lambda = ratio * n
  lambda_grid <-
    cfg$refit_lambda_ratio_grid * nrow(X)
  
  cv_error <- rep(NA_real_, length(lambda_grid))
  
  for (ll in seq_along(lambda_grid)) {
    
    lambda <- lambda_grid[ll]
    
    pred_cv <- rep(NA_real_, length(y))
    
    for (ff in seq_len(Kfold)) {
      
      idx_test <- which(fold_id == ff)
      idx_train <- which(fold_id != ff)
      
      beta_fold <- ridge_refit_beta(
        X[idx_train, , drop = FALSE],
        y[idx_train],
        lambda
      )
      
      pred_cv[idx_test] <-
        as.numeric(
          X[idx_test, , drop = FALSE] %*% beta_fold
        )
    }
    
    cv_error[ll] <-
      mean((y - pred_cv)^2)
  }
  
  best <- which.min(cv_error)
  lambda_best <- lambda_grid[best]
  
  beta_selected <- ridge_refit_beta(
    X,
    y,
    lambda_best
  )
  
  beta_hat[selected] <- beta_selected
  
  list(
    beta = beta_hat,
    selected = selected_indicator,
    lambda = lambda_best,
    cv_mse = cv_error[best]
  )
}
evaluate_real_fit <- function(fit, prepared, cfg, seed = 1L) {
  
  ## =========================================================
  ## Native beta returned by each original method
  ## =========================================================
  beta_native <- fit$beta_hat
  
  stopifnot(
    length(beta_native) == prepared$training$p,
    all(is.finite(beta_native))
  )
  
  ## =========================================================
  ## Common post-selection beta estimator
  ## =========================================================
  refit <- post_selection_ridge_cv(
    dat = prepared$training,
    fit = fit,
    cfg = cfg,
    seed = seed
  )
  
  beta_refit <- refit$beta
  selected <- refit$selected
  
  ## =========================================================
  ## Fair main prediction:
  ## same target-only ridge refit for every method
  ## =========================================================
  yp <- as.numeric(
    prepared$test$X %*% beta_refit
  )
  
  ## Native prediction, retained only as diagnostic
  yp_native <- as.numeric(
    prepared$test$X %*% beta_native
  )
  
  y_std <- prepared$test$y
  s <- prepared$target_scaling
  
  pred <-
    yp * s$y_sd + s$y_mean
  
  pred_native <-
    yp_native * s$y_sd + s$y_mean
  
  actual <-
    prepared$test$y_log2_tpm
  
  null_mse <-
    mean((actual - s$y_mean)^2)
  
  ## =========================================================
  ## Main common-refit metrics
  ## =========================================================
  mse <- mean((pred - actual)^2)
  
  mae <- mean(abs(pred - actual))
  
  r2 <- if (sum((actual - mean(actual))^2) > 0) {
    1 -
      sum((pred - actual)^2) /
      sum((actual - mean(actual))^2)
  } else {
    NA_real_
  }
  
  mse_standardized <-
    mean((yp - y_std)^2)
  
  relative_mse <-
    mse / null_mse
  
  ## =========================================================
  ## Original method-specific beta metrics
  ## =========================================================
  mse_native <-
    mean((pred_native - actual)^2)
  
  mae_native <-
    mean(abs(pred_native - actual))
  
  r2_native <- if (sum((actual - mean(actual))^2) > 0) {
    1 -
      sum((pred_native - actual)^2) /
      sum((actual - mean(actual))^2)
  } else {
    NA_real_
  }
  
  mse_standardized_native <-
    mean((yp_native - y_std)^2)
  
  relative_mse_native <-
    mse_native / null_mse
  
  list(
    
    metrics = list(
      
      ## main fair comparison
      mse = mse,
      mae = mae,
      r2 = r2,
      mse_standardized = mse_standardized,
      relative_mse = relative_mse,
      
      ## native estimator diagnostic
      mse_native = mse_native,
      mae_native = mae_native,
      r2_native = r2_native,
      mse_standardized_native = mse_standardized_native,
      relative_mse_native = relative_mse_native,
      
      null_mse = null_mse,
      model_size = sum(selected),
      
      refit_lambda = refit$lambda
    ),
    
    predictions = data.frame(
      sample = prepared$test$sample,
      donor = prepared$test$donor,
      observed = actual,
      
      ## main prediction
      predicted = pred,
      
      ## original estimator
      predicted_native = pred_native
    ),
    
    selected = selected,
    
    beta_refit = beta_refit,
    beta_native = beta_native,
    
    refit_lambda = refit$lambda
  )
}
