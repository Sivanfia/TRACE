# Reconstruct target arrays from the saved donor split and training scaling.
# Source summaries and calibrated precision are reused, never re-estimated.
check_retained_job <- function(b) {
  if (!identical(b$status,'completed')) stop('Retained job is not completed')
  expected<-c('ACVB','TRADER','Trans-Lasso','CONCERT')
  if (!identical(sort(b$metrics$method),sort(expected)) || any(b$metrics$status!='ok'))
    stop('Retained metrics must contain each of the four baseline methods exactly once')
  for(nm in c('metrics','predictions','features')) {
    z<-b[[nm]]
    if(is.null(z)||!nrow(z)||any(!z$method%in%expected))stop('Old TRACE or missing baseline in ',nm)
    if(any(z$target!=b$target)||any(z$rep!=b$rep))stop('Retained target/rep mismatch in ',nm)
  }
  if(any(b$metrics$seed!=b$seed))stop('Retained seed mismatch')
  if(any(b$source_audit$target_test_overlap!=0))stop('Saved split has target-test donor overlap')
  invisible(TRUE)
}
reconstruct_real_inputs <- function(b,data_dir) {
  check_retained_job(b)
  d<-read_gtex(data_dir,b$target);genes<-as.character(b$feature_map$gene_id)
  itr<-match(b$split$train_samples,d$sample)
  ref<-b$predictions[b$predictions$method=='ACVB',,drop=FALSE]
  ite<-match(ref$sample,d$sample)
  if(anyNA(c(itr,ite))||anyDuplicated(c(itr,ite)))stop('Stored samples missing or duplicated')
  if(!identical(as.character(d$donor[itr]),as.character(b$split$train_donors)) ||
     !setequal(d$donor[ite],b$split$test_donors))stop('Stored donor split mismatch')
  if(!setequal(c(itr,ite),seq_along(d$sample)))stop('Stored train/test split is incomplete')
  if(any(!genes%in%colnames(d$X))||!identical(rownames(b$B),genes))stop('Stored feature order mismatch')
  s<-b$target_scaling;Xraw<-d$X[itr,genes,drop=FALSE];yraw<-d$y[itr]
  if(max(abs(colMeans(Xraw)-s$x_mean))>1e-12 || max(abs(apply(Xraw,2,sd)-s$x_sd))>1e-12 ||
     abs(mean(yraw)-s$y_mean)>1e-12 || abs(sd(yraw)-s$y_sd)>1e-12)stop('Stored training scales do not match raw data')
  if(max(abs(d$y[ite]-ref$observed))>1e-12)stop('Stored held-out observations do not match raw data')
  X<-sweep(sweep(Xraw,2,s$x_mean,'-'),2,s$x_sd,'/');y<-(yraw-s$y_mean)/s$y_sd
  ns<-b$precision$ns;n<-nrow(X);p<-ncol(X);K<-ncol(b$B)
  if(any(b$metrics$n_train!=n)||any(b$metrics$p!=p)||any(b$metrics$K!=K)||
     any(abs(b$metrics$ns_calibrated-ns)>1e-12))stop('Stored training dimensions/precision mismatch')
  dat<-list(X=X,y=y,n=n,p=p,K=K,B=b$B,ns=ns,target_donor=d$donor[itr],
    XtX=crossprod(X),Xty=as.numeric(crossprod(X,y)),yTy=sum(y^2),x_norm2=colSums(X^2))
  dat<-configure_trace_source_precision(dat,mode='high_dimensional_logp')
  list(training=dat,feature_map=b$feature_map,source_audit=b$source_audit,source_overlap=b$source_overlap,
    precision=b$precision,target_scaling=s,split=b$split,
    test=list(X=sweep(sweep(d$X[ite,genes,drop=FALSE],2,s$x_mean,'-'),2,s$x_sd,'/'),
      y=(d$y[ite]-s$y_mean)/s$y_sd,y_log2_tpm=d$y[ite],sample=d$sample[ite],donor=d$donor[ite]))
}
check_retained_data_files <- function(data_dir,manifest,targets) {
  for(t in targets){
    name<-paste0(t,'.cns_tpm.tsv.gz');path<-file.path(data_dir,name)
    key<-paste0('data/processed/gtex/',name);expected<-manifest$files[[key]]
    if(!file.exists(path)||is.null(expected)||!identical(unname(tools::md5sum(path)),unname(expected)))
      stop('Raw expression data differs from the retained analysis: ',name)
  }
  invisible(TRUE)
}
