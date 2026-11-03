read_gtex <- function(data_dir,tissue) {
  f<-file.path(data_dir,paste0(tissue,'.cns_tpm.tsv.gz'))
  con<-gzfile(f,open='rt'); on.exit(close(con))
  z<-read.delim(con,check.names=FALSE)
  iy<-which(z$symbol=='SH2D2A')
  stopifnot(length(iy)==1L,!anyDuplicated(z$gene_id))
  m<-as.matrix(z[,-c(1,2)]); storage.mode(m)<-'double'
  if(any(!is.finite(m))||any(m<0)) stop('Invalid TPM values in ',tissue)
  m<-log2(m+1); ids<-colnames(m)
  donor<-vapply(strsplit(ids,'-',fixed=TRUE),function(a)paste(a[1:2],collapse='-'),'')
  if(anyDuplicated(ids)||any(!startsWith(donor,'GTEX-')))stop('Invalid sample IDs')
  x<-t(m[-iy,,drop=FALSE]); colnames(x)<-z$gene_id[-iy]
  list(X=x,y=as.numeric(m[iy,]),sample=ids,donor=donor,
    genes=z[-iy,c('gene_id','symbol')],tissue=tissue)
}
subset_tissue <- function(d,i) {
  d$X<-d$X[i,,drop=FALSE];d$y<-d$y[i];d$sample<-d$sample[i];d$donor<-d$donor[i];d
}
scale_train <- function(X,y) {
  mx<-colMeans(X);sx<-apply(X,2,sd);my<-mean(y);sy<-sd(y)
  if(any(!is.finite(sx)|sx<1e-10)||!is.finite(sy)||sy<1e-10)stop('Constant training variable')
  list(X=sweep(sweep(X,2,mx,'-'),2,sx,'/'),y=(y-my)/sy,
    x_mean=mx,x_sd=sx,y_mean=my,y_sd=sy)
}
donor_folds <- function(donor,k=5L) {
  u<-unique(donor); f<-sample(rep(seq_len(min(k,length(u))),length.out=length(u)))
  f[match(donor,u)]
}
ridge_beta <- function(X,y,lambda,weights=rep(1,nrow(X))) {
  w<-weights/sum(weights); mx<-colSums(X*w);my<-sum(y*w)
  xc<-sweep(X,2,mx,'-');yc<-y-my
  solve_spd(crossprod(xc,xc*w)+diag(lambda,ncol(X)),crossprod(xc,yc*w))[,1]
}
# Joint donor bootstrap uses the same multiplicity in every tissue, preserving
# dependence across sources. Scaling and ridge lambda are held fixed; this is
# conditional variance calibration, not a bias correction or posterior interval.
calibrate_precision <- function(src,B,dj,cfg,seed) {
  nk<-vapply(src,function(s)nrow(s$X),0L)
  set.seed(seed);u<-unique(unlist(lapply(src,`[[`,'donor')))
  boot<-array(NA_real_,c(nrow(B),ncol(B),cfg$source_bootstrap))
  index<-lapply(src,function(s)match(s$donor,u))
  for(b in seq_len(cfg$source_bootstrap)) {
    w<-tabulate(sample.int(length(u),length(u),replace=TRUE),nbins=length(u))
    for(k in seq_along(src))boot[,k,b]<-ridge_beta(src[[k]]$X,src[[k]]$y,
      src[[k]]$lambda,w[index[[k]]])
  }
  v<-apply(boot,c(1,2),var)
  q<-as.numeric(quantile(as.vector(v*dj),cfg$precision_quantile,names=FALSE))
  calibrated<-min(min(nk),1/max(q,1e-12))
  ns<-if(cfg$precision_mode=='min_source_n')min(nk)else calibrated
  list(ns=ns,variance=v,normalized_variance_quantile=q,
    sample_size_cap=min(nk),mode=cfg$precision_mode,n_boot=cfg$source_bootstrap)
}
prepare_real_data <- function(target_data,sources,cfg,seed) {
  set.seed(seed);u<-unique(target_data$donor)
  te<-sample(u,ceiling(cfg$test_fraction*length(u)))
  train<-subset_tissue(target_data,which(!target_data$donor%in%te))
  test<-subset_tissue(target_data,which(target_data$donor%in%te))
  excluded<-if(cfg$donor_policy=='disjoint_target')u else te
  audit<-lapply(sources,function(d){
    keep<-!d$donor%in%excluded
    data.frame(source=d$tissue,n_before=nrow(d$X),n_used=sum(keep),
      target_train_overlap=length(intersect(d$donor[keep],train$donor)),
      target_test_overlap=length(intersect(d$donor[keep],te)),
      eligible=sum(keep)>=cfg$min_source_n)
  });audit<-do.call(rbind,audit)
  sources<-lapply(sources,function(d)subset_tissue(d,which(!d$donor%in%excluded)))
  sources<-sources[audit$eligible]
  if(length(sources)<2L)stop('Fewer than two eligible sources after donor exclusion')
  stopifnot(all(audit$target_test_overlap==0L))
  if(cfg$donor_policy=='disjoint_target')stopifnot(all(audit$target_train_overlap==0L))
  genes<-Reduce(intersect,c(list(colnames(train$X)),lapply(sources,function(d)colnames(d$X))))
  mats<-c(list(train$X[,genes,drop=FALSE]),lapply(sources,function(d)d$X[,genes,drop=FALSE]))
  ok<-Reduce('&',lapply(mats,function(x)apply(x,2,sd)>1e-10)); genes<-genes[ok]
  cr<-as.numeric(cor(train$X[,genes,drop=FALSE],train$y))
  ix<-which(is.finite(cr)&abs(cr)>=cfg$cor_threshold)
  ix<-head(ix[order(-abs(cr[ix]),genes[ix])],cfg$max_features); genes<-genes[ix]
  if(length(genes)<2L)stop('Too few predictors under the prespecified screening rule')
  a<-scale_train(train$X[,genes,drop=FALSE],train$y)
  sf<-lapply(seq_along(sources),function(k){
    d<-sources[[k]]; s<-scale_train(d$X[,genes,drop=FALSE],d$y)
    set.seed(seed+100L+k)
    cv<-glmnet::cv.glmnet(s$X,s$y,alpha=0,intercept=FALSE,standardize=FALSE,
      foldid=donor_folds(d$donor))
    s$lambda<-cv$lambda.min;s$beta<-ridge_beta(s$X,s$y,s$lambda)
    s$donor<-d$donor;s$tissue<-d$tissue;s
  })
  B<-do.call(cbind,lapply(sf,`[[`,'beta'));colnames(B)<-vapply(sf,`[[`,'','tissue');rownames(B)<-genes
  n<-nrow(a$X);p<-ncol(a$X);prec<-calibrate_precision(sf,B,colSums(a$X^2)/n,cfg,seed+500L)
  dat<-list(X=a$X,y=a$y,n=n,p=p,K=ncol(B),B=B,ns=prec$ns,target_donor=train$donor,
    XtX=crossprod(a$X),Xty=as.numeric(crossprod(a$X,a$y)),yTy=sum(a$y^2),x_norm2=colSums(a$X^2),
    X_source=lapply(sf,`[[`,'X'),y_source=lapply(sf,`[[`,'y'),
    n_source_raw=vapply(sf,function(s)nrow(s$X),0L))
  dat$X_all<-do.call(rbind,c(list(dat$X),dat$X_source));dat$y_all<-c(dat$y,unlist(dat$y_source))
  dat$n_vec_trans<-c(n,dat$n_source_raw)
  # Do not set n_vec: legacy code takes its mean and would override calibrated ns.
  dat<-configure_trace_source_precision(dat,mode='high_dimensional_logp')
  overlap<-outer(seq_along(sf),seq_along(sf),Vectorize(function(i,j)length(intersect(sf[[i]]$donor,sf[[j]]$donor))))
  dimnames(overlap)<-list(colnames(B),colnames(B))
  fmap<-train$genes[match(genes,train$genes$gene_id),];fmap$train_correlation<-cr[ix]
  list(training=dat,feature_map=fmap,source_audit=audit,source_overlap=overlap,
    precision=prec,target_scaling=a[c('x_mean','x_sd','y_mean','y_sd')],source_fits=sf,
    split=list(seed=seed,train_donors=unique(train$donor),test_donors=te,train_samples=train$sample),
    test=list(X=sweep(sweep(test$X[,genes,drop=FALSE],2,a$x_mean,'-'),2,a$x_sd,'/'),
      y=(test$y-a$y_mean)/a$y_sd,y_log2_tpm=test$y,sample=test$sample,donor=test$donor))
}
