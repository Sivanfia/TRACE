atomic_rds <- function(x,path) {
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  tmp<-paste0(path,'.',Sys.getpid(),'.tmp');saveRDS(x,tmp)
  if(!file.rename(tmp,path))stop('Atomic write failed: ',path)
}
csv <- function(x,path) write.csv(x,path,row.names=FALSE,na='')
combine_real_job <- function(fresh,baseline,signature) {
  check_retained_job(baseline)
  if(!identical(fresh$status,'completed')||!identical(fresh$signature,signature))stop('New TRACE checkpoint is not completed or has a different signature')
  if(!identical(fresh$target,baseline$target)||fresh$rep!=baseline$rep||fresh$seed!=baseline$seed)
    stop('New and retained target/rep/seed do not match')
  if(nrow(fresh$metrics)!=1L||!identical(as.character(fresh$metrics$method),'TRACE')||fresh$metrics$status!='ok')
    stop('Expected one successful new TRACE result')
  if(fresh$metrics$seed!=baseline$seed||fresh$metrics$target!=baseline$target||fresh$metrics$rep!=baseline$rep)
    stop('New TRACE metric identifiers do not match retained split')
  for(nm in c('predictions','features')) {
    z<-fresh[[nm]]
    if(!nrow(z)||any(z$method!='TRACE')||any(z$target!=baseline$target)||any(z$rep!=baseline$rep))stop('Invalid new TRACE ',nm)
  }
  ref<-baseline$predictions[baseline$predictions$method=='ACVB',,drop=FALSE]
  if(!identical(as.character(fresh$predictions$sample),as.character(ref$sample))||
     !identical(as.character(fresh$predictions$donor),as.character(ref$donor))||
     max(abs(fresh$predictions$observed-ref$observed))>1e-12)stop('New held-out samples/observations differ from baseline')
  if(!identical(as.character(fresh$features$gene_id),as.character(baseline$feature_map$gene_id)))stop('New gene order differs from saved screened genes')
  for(nm in c('n_train','n_test','p','K','ns_calibrated','null_mse'))
    if(max(abs(fresh$metrics[[nm]]-baseline$metrics[[nm]]))>1e-10)stop('New comparison scale differs: ',nm)
  out<-fresh
  for(nm in c('metrics','predictions','features')) {
    if(!identical(names(fresh[[nm]]),names(baseline[[nm]])))stop('Mismatched columns in ',nm)
    out[[nm]]<-rbind(fresh[[nm]],baseline[[nm]])
  }
  for(nm in c('concert_transfer','concert_prior','source_audit','source_overlap','feature_map','split','precision','B','source_scaling','target_scaling'))out[[nm]]<-baseline[[nm]]
  out$retained_signature<-baseline$retained_signature
  out
}
run_trace_job <- function(baseline_path,cfg,data_dir,out_dir,signature) {
  b<-readRDS(baseline_path);check_retained_job(b)
  name<-basename(baseline_path);path<-file.path(out_dir,'trace_refit/jobs',name)
  cp<-file.path(out_dir,'comparison/jobs',name)
  if(file.exists(path)) {
    old<-readRDS(path)
    if(!identical(old$signature,signature))stop('Checkpoint signature mismatch. Use a new output directory: ',path)
    if(identical(old$status,'completed')){
      atomic_rds(combine_real_job(old,b,signature),cp)
      message('Resume: ',name);return(invisible(old))
    }
    archive<-paste0(path,'.failed_',format(Sys.time(),'%Y%m%d_%H%M%S'))
    if(!file.rename(path,archive))stop('Could not archive failed checkpoint')
  }
  message('Start TRACE-EM: ',b$target,' rep=',b$rep)
  start<-proc.time()[3]
  fresh<-tryCatch({
    p<-reconstruct_real_inputs(b,data_dir);prep_sec<-unname(proc.time()[3]-start)
    warn<-character();t0<-proc.time()[3]
    # TRACE was method 1 in the original analysis: retain seed+1 for fitting AND refit CV.
    f<-withCallingHandlers(fit_real_method('TRACE',p$training,cfg,b$seed+1L),
      warning=function(w){warn<<-c(warn,conditionMessage(w));invokeRestart('muffleWarning')})
    elapsed<-unname(proc.time()[3]-t0)
    init_error<-c(omega=max(abs(f$target_initialization$omega-b$reference_initialization$omega)),
      beta=max(abs(f$preliminary_alpha_beta_target-b$reference_initialization$beta)))
    if(any(!is.finite(init_error))||any(init_error>1e-7))stop('Target-only initialization mismatch; retained comparison was not merged')
    e<-evaluate_real_fit(f,p,cfg,seed=b$seed+1L)
    m<-data.frame(target=b$target,rep=b$rep,seed=b$seed,method='TRACE',backend='TRACE',
      n_train=p$training$n,n_test=length(p$test$y),p=p$training$p,K=p$training$K,ns_calibrated=p$precision$ns,
      prep_sec=prep_sec,fit_sec=elapsed,warnings=paste(unique(warn),collapse=' | '),status='ok',error='',
      converged=f$converged%||%NA,n_iter=f$n_iter%||%NA_integer_,max_psrf=NA_real_,stringsAsFactors=FALSE)
    for(nm in names(e$metrics))m[[nm]]<-e$metrics[[nm]]
    m<-m[,names(b$metrics),drop=FALSE]
    pred<-cbind(target=b$target,rep=b$rep,method='TRACE',e$predictions)
    ft<-cbind(target=b$target,rep=b$rep,method='TRACE',p$feature_map,
      beta=e$beta_refit,beta_refit=e$beta_refit,beta_native=f$beta_hat,selected=e$selected,
      omega=f$omega,gamma=f$gamma,joint_inclusion_borrowing=f$omega*f$gamma)
    x<-list(target=b$target,rep=b$rep,seed=b$seed,signature=signature,status='completed',
      metrics=m,predictions=pred,features=ft,
      alpha=data.frame(target=b$target,rep=b$rep,source=colnames(p$training$B),alpha=f$alpha),
      fits=list(TRACE=f),input_reused=TRUE,initialization_max_error=init_error,
      elapsed_sec=unname(proc.time()[3]-start))
    combined<-combine_real_job(x,b,signature)
    atomic_rds(x,path);atomic_rds(combined,cp)
    message(sprintf('  TRACE-EM: MSE=%.6g, converged=%s, %.2fs',m$mse,m$converged,x$elapsed_sec))
    x
  },error=function(e){
    x<-list(target=b$target,rep=b$rep,seed=b$seed,signature=signature,status='error',error=conditionMessage(e))
    atomic_rds(x,path)
    if(file.exists(cp))unlink(cp)
    message('Failed ',name,': ',x$error);x
  })
  invisible(fresh)
}
