# Pinned author Trans-Lasso, with mathematically equivalent stable aggregation.
# No environment/package discovery, no automatic alternative estimator.
.translasso_dir <- local({
  paths <- Filter(function(x)is.character(x)&&length(x)==1L&&nzchar(x),
    lapply(sys.frames(),function(x)x$ofile))
  if(!length(paths))stop('Source translasso_official.R from its file path.')
  dirname(normalizePath(tail(paths,1L)[[1L]],mustWork=TRUE))
})
.translasso_reference_file <- file.path(.translasso_dir,'TransLasso-author-285009e.R')
.translasso_adapter_file <- file.path(.translasso_dir,'translasso_official.R')
.translasso_expected_md5 <- 'b37cc462f1af08c3b3ee61476e21eef8'
.translasso_author_commit <- '285009edf3e638567474108f823adcc1dc6c7c4a'
.translasso_reference_env <- NULL

.translasso_stable_aggregate <- function(B,X.test,y.test,total.step=10,selection=FALSE) {
  B<-as.matrix(B);p<-nrow(B);K<-ncol(B)
  if(K<1L||any(!is.finite(B)))stop('Trans-Lasso has invalid candidate coefficients.')
  if(all(B==0))return(list(theta=rep(1/K,K),beta=rep(0,p),beta.ew=rep(0,p)))
  colnames(B)<-NULL
  softmax<-function(logw) {
    if(any(!is.finite(logw)))stop('Trans-Lasso aggregation has nonfinite log weights.')
    w<-exp(logw-max(logw));w/sum(w)
  }
  if(selection) {
    khat<-which.min(colSums((y.test-X.test%*%B)^2))
    theta.hat<-rep(0,K);theta.hat[khat]<-1
    beta<-B[,khat];beta.ew<-NULL
  }else{
    theta.hat<-softmax(-colSums((y.test-X.test%*%B)^2)/2)
    theta.old<-theta.hat;beta<-as.numeric(B%*%theta.hat);beta.ew<-beta
    for(ss in seq_len(total.step)) {
      theta.hat<-softmax(-colSums((y.test-X.test%*%B)^2)/2+
        colSums((as.vector(X.test%*%beta)-X.test%*%B)^2)/8)
      beta<-as.numeric(B%*%theta.hat*1/4+3/4*beta)
      if(sum(abs(theta.hat-theta.old))<1e-3)break
      theta.old<-theta.hat
    }
  }
  list(theta=theta.hat,beta=beta,beta.ew=beta.ew)
}

translasso_backend_info <- function() {
  if(!file.exists(.translasso_reference_file))stop('Pinned author source is missing: ',.translasso_reference_file)
  md5<-unname(tools::md5sum(.translasso_reference_file))
  if(!identical(md5,.translasso_expected_md5))stop('Pinned Trans-Lasso source checksum mismatch. Restore the bundled author file.')
  list(backend='author_Q_aggregation_twofold',commit=.translasso_author_commit,
    source_md5=md5,adapter_md5=unname(tools::md5sum(.translasso_adapter_file)),
    numerics='stable_softmax_and_zero_candidate_list',rng_policy='preserve_caller_state')
}

load_translasso_function <- function() {
  translasso_backend_info()
  if(!requireNamespace('glmnet',quietly=TRUE))stop('Trans-Lasso requires glmnet.')
  if(is.null(.translasso_reference_env)) {
    e<-new.env(parent=asNamespace('glmnet'))
    sys.source(.translasso_reference_file,envir=e)
    if(!is.function(e$Trans.lasso)||!is.function(e$agg.fun))stop('Invalid author Trans-Lasso source.')
    # Keep the author's ranking, CV, thresholding, glmnet defaults and Q update.
    # Override only exponent normalization and the all-zero return structure.
    e$agg.fun<-.translasso_stable_aggregate
    .translasso_reference_env <<- e
  }
  .translasso_reference_env$Trans.lasso
}

translasso_preflight <- function(verbose=TRUE) {
  load_translasso_function();info<-translasso_backend_info()
  if(verbose)message('[Trans-Lasso] official author ',info$commit,'; Q-aggregation; strict backend; stable numerics')
  invisible(info)
}

get_beta_from_translasso_output <- function(out) {
  if(!is.list(out)||is.null(out$beta.hat))stop('Author Trans-Lasso did not return beta.hat.')
  as.numeric(out$beta.hat)
}

fit_translasso_official_if_available <- function(dat,split=c('twofold','first_half'),l1=TRUE) {
  split<-match.arg(split)
  if(!isTRUE(l1))stop('This pinned simulation baseline supports the author l1 algorithm only.')
  # Isolate reference CV randomness from every subsequent simulation method.
  had_rng<-exists('.Random.seed',envir=.GlobalEnv,inherits=FALSE)
  rng<-if(had_rng)get('.Random.seed',envir=.GlobalEnv)else NULL
  on.exit({if(had_rng)assign('.Random.seed',rng,envir=.GlobalEnv)
    else if(exists('.Random.seed',envir=.GlobalEnv,inherits=FALSE))rm('.Random.seed',envir=.GlobalEnv)},add=TRUE)
  info<-translasso_backend_info();fun<-load_translasso_function()
  X<-as.matrix(dat$X_all);y<-as.numeric(dat$y_all);nv<-dat$n_vec_trans
  if(length(nv)<2L||any(!is.finite(nv))||any(nv<1)||any(nv!=floor(nv))||
     nrow(X)!=sum(nv)||length(y)!=nrow(X)||ncol(X)!=dat$p||nv[1]!=dat$n||
     any(!is.finite(X))||any(!is.finite(y)))stop('Invalid stacked training data for author Trans-Lasso.')
  n0<-as.integer(nv[1]);half<-floor(n0/2)
  if(half<8L)stop('Each target half needs >=8 observations for the author 8-fold CV.')
  # The author's marginal screening indexes round(n_train/3) coordinates.
  # Reject unsupported tiny-p input rather than silently changing that rule.
  if(dat$p<max(1L,round(ceiling(n0/2)/3)))stop('p is too small for the author marginal-screening rule.')
  I1<-seq_len(half);I2<-seq.int(half+1L,n0)
  f1<-fun(X=X,y=y,n.vec=nv,I.til=I1,l1=TRUE)
  b1<-get_beta_from_translasso_output(f1)
  if(split=='twofold') {
    f2<-fun(X=X,y=y,n.vec=nv,I.til=I2,l1=TRUE)
    b2<-get_beta_from_translasso_output(f2)
    if(length(b2)!=dat$p||any(!is.finite(b2)))stop('Nonfinite or invalid second-fold author coefficients.')
    beta<-(b1+b2)/2;detail<-list(split1=f1,split2=f2)
  }else{beta<-b1;detail<-f1}
  if(length(b1)!=dat$p||length(beta)!=dat$p||any(!is.finite(beta)))stop('Nonfinite or invalid author Trans-Lasso coefficients.')
  list(method='Trans_Lasso_official',beta_hat=beta,detail=detail,
    translasso_backend=if(split=='twofold')info$backend else 'author_Q_aggregation_first_half',
    translasso_commit=info$commit,translasso_source_md5=info$source_md5,
    translasso_adapter_md5=info$adapter_md5,translasso_split=split,
    translasso_numerics=info$numerics,translasso_rng_policy=info$rng_policy)
}

fit_trans_lasso_main <- function(dat) {
  # Deliberately no tryCatch(..., error=NULL) and no fallback estimator.
  fit_translasso_official_if_available(dat,split='twofold',l1=TRUE)
}

translasso_validate_performance <- function(perf) {
  if(!is.data.frame(perf)||!'method'%in%names(perf))stop('Invalid performance table.')
  idx<-perf$method%in%c('Trans-Lasso','trans_lasso')
  if(!any(idx))return(invisible(TRUE))
  required<-c('method_name','translasso_backend','translasso_commit','translasso_source_md5','translasso_adapter_md5')
  if(!all(required%in%names(perf)))stop('Legacy Trans-Lasso results lack official provenance. Use a NEW result directory.')
  info<-translasso_backend_info();z<-perf[idx,,drop=FALSE]
  valid<-z$method_name=='Trans_Lasso_official' & z$translasso_backend==info$backend &
    z$translasso_commit==info$commit & z$translasso_source_md5==info$source_md5 &
    z$translasso_adapter_md5==info$adapter_md5
  if(anyNA(valid)||!all(valid))stop('Refusing to reuse nonofficial or differently versioned Trans-Lasso results. Use a NEW result directory.')
  invisible(TRUE)
}

translasso_preflight()
