atomic_rds <- function(x,path) {
  tmp<-paste0(path,'.',Sys.getpid(),'.tmp');saveRDS(x,tmp)
  if(!file.rename(tmp,path))stop('Atomic write failed: ',path)
}
csv <- function(x,path) write.csv(x,path,row.names=FALSE,na='')
run_real_job <- function(target,rep_id,cfg,data_dir,out_dir,signature,retry_failed=FALSE) {
  path<-file.path(out_dir,'jobs',sprintf('%s_rep%03d.rds',target,rep_id))
  if(file.exists(path)) {
    old<-readRDS(path)
    if(!identical(old$signature,signature))stop('Checkpoint signature mismatch: ',path)
    failed<-old$status=='error'||any(old$metrics$status!='ok')
    if(!failed||!retry_failed){message('Resume: ',basename(path));return(invisible(old))}
    archive<-paste0(path,'.failed_',format(Sys.time(),'%Y%m%d_%H%M%S'))
    if(!file.rename(path,archive))stop('Could not archive failed checkpoint')
  }
  seed<-as.integer(cfg$seed+10000L*match(target,cfg$targets)+100L*rep_id)
  message('Start ',target,' rep=',rep_id)
  start<-proc.time()[3]
  result<-tryCatch({
    all_tissues<-sub('.cns_tpm.tsv.gz','',list.files(data_dir,pattern='\\.cns_tpm\\.tsv\\.gz$'),fixed=TRUE)
    sources<-setdiff(if(cfg$source_mode=='brain5')cfg$targets else all_tissues,target)
    p<-prepare_real_data(read_gtex(data_dir,target),lapply(sources,function(t)read_gtex(data_dir,t)),cfg,seed)
    prep_sec<-unname(proc.time()[3]-start)
    fits<-list();rows<-list();pred<-list();features<-list();alphas<-list()
    concert_transfer<-list();concert_prior<-list()
    for(j in seq_along(METHODS)) {
      method<-METHODS[j];warn<-character();t0<-proc.time()[3]
      f<-tryCatch(withCallingHandlers(fit_real_method(method,p$training,cfg,seed+j),
        warning=function(w){warn<<-c(warn,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e)e)
      elapsed<-unname(proc.time()[3]-t0)
      backend<-if(method=='Trans-Lasso')cfg$trans_backend else if(method=='CONCERT')
        paste0('official_author_',.concert_upstream_commit) else method
      base<-data.frame(target=target,rep=rep_id,seed=seed,method=method,
        backend=backend,
        n_train=p$training$n,n_test=length(p$test$y),p=p$training$p,K=p$training$K,
        ns_calibrated=p$precision$ns,prep_sec=prep_sec,fit_sec=elapsed,
        warnings=paste(unique(warn),collapse=' | '),stringsAsFactors=FALSE)
      if(!inherits(f,'error')) {
        e <- tryCatch(
          evaluate_real_fit(
            f,
            p,
            cfg,
            seed = seed + j
          ),
          error = function(e) e
        )
        if(inherits(e,'error'))f<-e
      }
      if(inherits(f,'error')) {
        base$status<-'error';base$error<-conditionMessage(f);base$converged<-NA;base$n_iter<-NA;base$max_psrf<-NA_real_
        for(nm in c(
          'mse',
          'mae',
          'r2',
          'mse_standardized',
          'relative_mse',
          
          'mse_native',
          'mae_native',
          'r2_native',
          'mse_standardized_native',
          'relative_mse_native',
          
          'null_mse',
          'model_size',
          'refit_lambda'
        )) {
          base[[nm]] <- NA_real_
        }
      }else{
        base$status<-'ok';base$error<-'';base$converged<-f$converged%||%NA
        base$n_iter<-f$n_iter%||%NA_integer_;base$max_psrf<-f$max_psrf%||%NA_real_
        for(nm in names(e$metrics))base[[nm]]<-e$metrics[[nm]]
        pred[[method]]<-cbind(target=target,rep=rep_id,method=method,e$predictions)
        features[[method]] <- cbind(
          target = target,
          rep = rep_id,
          method = method,
          p$feature_map,
          
          ## common refit beta
          beta = e$beta_refit,
          beta_refit = e$beta_refit,
          
          ## original estimator from each method
          beta_native = f$beta_hat,
          
          selected = e$selected,
          
          omega =
            f$omega %||%
            rep(NA_real_, p$training$p),
          
          gamma =
            f$gamma %||%
            rep(NA_real_, p$training$p),
          
          joint_inclusion_borrowing =
            if(method == 'TRACE') {
              f$omega * f$gamma
            } else {
              rep(NA_real_, p$training$p)
            }
        )
        if(method=='TRACE')alphas[[method]]<-data.frame(target=target,rep=rep_id,source=colnames(p$training$B),alpha=f$alpha)
        if(method=='CONCERT') {
          gi<-as.matrix(f$gamma_I)
          ix<-expand.grid(feature_index=seq_len(p$training$p),source_index=seq_len(p$training$K))
          concert_transfer[[method]]<-data.frame(
            target=target,rep=rep_id,
            gene_id=p$feature_map$gene_id[ix$feature_index],
            symbol=p$feature_map$symbol[ix$feature_index],
            source=colnames(p$training$B)[ix$source_index],
            transfer_probability=gi[cbind(ix$feature_index,ix$source_index)],
            stringsAsFactors=FALSE)
          concert_prior[[method]]<-data.frame(
            target=target,rep=rep_id,source=colnames(p$training$B),
            q_k_agreement=as.numeric(f$q_k_agreement_by_source),
            q_k_raw=f$q_k_raw,q_k_used=f$q_k_used,q_0_used=f$q_0_used,
            agreement_tolerance=f$q_k_tolerance,
            upstream_commit=f$upstream_commit,stringsAsFactors=FALSE)
        }
        fits[[method]]<-f
      }
      rows[[method]]<-base
      message(sprintf('  %s: %s, %.2fs, MSE=%s, converged=%s',method,base$status,elapsed,format(base$mse,digits=4),base$converged))
    }
    # Keep preprocessing scales and source coefficients to reproduce predictions,
    # but avoid repeating large raw expression matrices in each checkpoint.
    list(target=target,rep=rep_id,seed=seed,signature=signature,status='completed',
      metrics=do.call(rbind,rows),predictions=do.call(rbind,pred),features=do.call(rbind,features),
      alpha=do.call(rbind,alphas),concert_transfer=do.call(rbind,concert_transfer),
      concert_prior=do.call(rbind,concert_prior),fits=fits,feature_map=p$feature_map,
      split=p$split,source_audit=p$source_audit,source_overlap=p$source_overlap,
      precision=p$precision,target_scaling=p$target_scaling,B=p$training$B,
      source_scaling=lapply(p$source_fits,function(s)s[c('x_mean','x_sd','y_mean','y_sd','lambda','donor','tissue')]),
      elapsed_sec=unname(proc.time()[3]-start))
  },error=function(e)list(target=target,rep=rep_id,signature=signature,status='error',error=conditionMessage(e)))
  atomic_rds(result,path)
  if(result$status=='error')message('Job failed: ',result$error)
  invisible(result)
}
