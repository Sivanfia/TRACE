save_panel <- function(path,draw,width=7,height=5.5) {
  grDevices::pdf(paste0(path,'.pdf'),width=width,height=height,useDingbats=FALSE)
  tryCatch(draw(),finally=grDevices::dev.off())
  grDevices::png(paste0(path,'.png'),width=width,height=height,units='in',res=300)
  tryCatch(draw(),finally=grDevices::dev.off())
}
.bind_job_field <- function(jobs,name) {
  x<-Filter(function(z)!is.null(z)&&NROW(z)>0L,lapply(jobs,function(j)j[[name]]))
  if(length(x))do.call(rbind,x)else data.frame()
}
summarize_real_results <- function(out_dir,signature,methods) {
  jobs<-lapply(list.files(file.path(out_dir,'jobs'),pattern='\\.rds$',full.names=TRUE),readRDS)
  if(!length(jobs))stop('No checkpoints to summarize')
  if(any(!vapply(jobs,function(x)identical(x$signature,signature),logical(1))))
    stop('Mixed configurations/data/code in output directory; use a new --out')
  tab<-file.path(out_dir,'tables');fig<-file.path(out_dir,'figures')
  dir.create(tab,showWarnings=FALSE);dir.create(fig,showWarnings=FALSE)
  csv(do.call(rbind,lapply(jobs,function(x)data.frame(target=x$target,rep=x$rep,status=x$status,
    error=x$error%||%''))),file.path(tab,'job_status.csv'))
  ok<-jobs[vapply(jobs,function(x)x$status=='completed',logical(1))]
  if(!length(ok))return(invisible(NULL))
  m<-.bind_job_field(ok,'metrics')
  pred<-.bind_job_field(ok,'predictions');ft<-.bind_job_field(ok,'features')
  al<-.bind_job_field(ok,'alpha');ct<-.bind_job_field(ok,'concert_transfer')
  cp<-.bind_job_field(ok,'concert_prior')
  csv(m,file.path(tab,'metrics_per_rep.csv'));csv(pred,file.path(tab,'heldout_predictions.csv'))
  csv(ft,file.path(tab,'gene_results_per_rep.csv'));csv(al,file.path(tab,'trace_alpha_per_rep.csv'))
  csv(ct,file.path(tab,'concert_transfer_probability_per_rep.csv'))
  csv(cp,file.path(tab,'concert_prior_per_rep.csv'))
  sa<-do.call(rbind,lapply(ok,function(x)cbind(target=x$target,rep=x$rep,x$source_audit)))
  csv(sa,file.path(tab,'source_donor_audit.csv'))
  summ<-list();ii<-0L
  metrics <- c(
    ## common post-selection ridge metrics
    'mse',
    'mae',
    'r2',
    'relative_mse',
    'mse_standardized',
    
    ## original method-specific beta metrics
    'mse_native',
    'mae_native',
    'r2_native',
    'relative_mse_native',
    'mse_standardized_native',
    
    'model_size',
    'fit_sec'
  )
  for(t in unique(m$target))for(method in methods) {
    z<-m[m$target==t&m$method==method,,drop=FALSE]
    for(metric in metrics) {
      x<-z[[metric]][z$status=='ok'&is.finite(z[[metric]])];ii<-ii+1L
      summ[[ii]]<-data.frame(target=t,method=method,metric=metric,attempted=nrow(z),
        successful=sum(z$status=='ok'),failed=sum(z$status!='ok'),
        convergence_checked=sum(z$status=='ok'&!is.na(z$converged)),
        converged=sum(z$status=='ok'&z$converged,na.rm=TRUE),n=length(x),
        mean=if(length(x))mean(x)else NA_real_,sd=if(length(x)>1L)sd(x)else NA_real_,
        median=if(length(x))median(x)else NA_real_)
    }
  }
  csv(do.call(rbind,summ),file.path(tab,'metrics_summary.csv'))
  # Descriptive paired differences only: repeated random splits are dependent.
  paired<-list();ii<-0L
  for(t in unique(m$target))for(b in setdiff(methods,'TRACE')) {
    a<-m[m$target==t&m$method=='TRACE'&m$status=='ok',c('rep','mse')]
    z<-m[m$target==t&m$method==b&m$status=='ok',c('rep','mse')]
    q<-merge(a,z,by='rep',suffixes=c('_trace','_baseline'));d<-q$mse_trace-q$mse_baseline
    ii<-ii+1L;paired[[ii]]<-data.frame(target=t,baseline=b,n_paired=length(d),
      mean_mse_difference=if(length(d))mean(d)else NA_real_,
      sd_mse_difference=if(length(d)>1)sd(d)else NA_real_,trace_win_fraction=if(length(d))mean(d<0)else NA_real_)
  }
  csv(do.call(rbind,paired),file.path(tab,'paired_mse_differences.csv'))
  stable<-list();ii<-0L
  if(!is.null(ft)&&nrow(ft))for(t in unique(ft$target))for(method in methods) {
    z<-ft[ft$target==t&ft$method==method,,drop=FALSE]
    den<-sum(m$target==t&m$method==method&m$status=='ok')
    for(g in unique(z$gene_id)) {
      x<-z[z$gene_id==g,,drop=FALSE];ii<-ii+1L
      stable[[ii]]<-data.frame(target=t,method=method,gene_id=g,symbol=x$symbol[1],
        successful_reps=den,screening_frequency=nrow(x)/den,selection_frequency=sum(x$selected)/den,
        mean_omega_if_screened=if(all(is.na(x$omega)))NA_real_ else mean(x$omega,na.rm=TRUE),
        mean_gamma_if_screened=if(all(is.na(x$gamma)))NA_real_ else mean(x$gamma,na.rm=TRUE),
        mean_joint_inclusion_borrowing=if(all(is.na(x$joint_inclusion_borrowing)))NA_real_ else sum(x$joint_inclusion_borrowing,na.rm=TRUE)/den)
    }
  }
  st<-if(length(stable))do.call(rbind,stable)else data.frame()
  csv(st,file.path(tab,'gene_stability.csv'))
  colors<-c(TRACE='#D55E00',ACVB='#0072B2',TRADER='#009E73','Trans-Lasso'='#CC79A7',CONCERT='#56B4E9')
  for(t in unique(m$target)) {
    d<-file.path(fig,t);dir.create(d,showWarnings=FALSE)
    for(metric in c('mse','mae','r2','relative_mse','model_size')) {
      z<-m[m$target==t&m$status=='ok'&is.finite(m[[metric]]),,drop=FALSE]
      if(!nrow(z))next
      save_panel(file.path(d,paste0(metric,'_comparison')),function(){
        par(mar=c(5,4,4,1)); groups<-split(z[[metric]],factor(z$method,levels=methods))
        boxplot(groups,col=colors[methods],ylab=metric,main=gsub('_',' ',t),outline=FALSE)
        for(k in seq_along(methods))points(rep(k,length(groups[[k]])),groups[[k]],pch=16,col=adjustcolor('black',.25),cex=.7)
        mtext('Held-out donor splits; all successful fits, including nonconverged fits',side=3,line=.3,cex=.7)
        if(metric=='relative_mse')abline(h=1,lty=2,col='gray40')
      })
    }
    if(!is.null(al)&&nrow(al)) {
      z<-al[al$target==t,,drop=FALSE]
      labels<-c(brain_amygdala='Amygdala',brain_anterior_cingulate_cortex_ba24='ACC (BA24)',
        brain_cortex='Cortex',brain_hippocampus='Hippocampus',
        brain_nucleus_accumbens_basal_ganglia='Nucleus accumbens',brain_substantia_nigra='Substantia nigra')
      zn<-unname(labels[z$source]);zn[is.na(zn)]<-gsub('_',' ',z$source[is.na(zn)])
      if(nrow(z))save_panel(file.path(d,'trace_source_weights'),function(){
        par(mar=c(4,if(length(unique(zn))>6)19 else 10,3,1));boxplot(split(z$alpha,zn),las=1,
          horizontal=TRUE,col='#E69F00',ylim=c(0,1),xlab='TRACE alpha',main=gsub('_',' ',t),cex.axis=.8)
      },width=if(length(unique(zn))>6)10 else 7,height=max(5.5,.26*length(unique(zn))))
    }
    if(nrow(ct)) {
      z<-ct[ct$target==t,,drop=FALSE]
      labels<-c(brain_amygdala='Amygdala',brain_anterior_cingulate_cortex_ba24='ACC (BA24)',
        brain_cortex='Cortex',brain_hippocampus='Hippocampus',
        brain_nucleus_accumbens_basal_ganglia='Nucleus accumbens',brain_substantia_nigra='Substantia nigra')
      zn<-unname(labels[z$source]);zn[is.na(zn)]<-gsub('_',' ',z$source[is.na(zn)])
      if(nrow(z))save_panel(file.path(d,'concert_source_transfer_probability'),function(){
        par(mar=c(4,if(length(unique(zn))>6)19 else 10,3,1));boxplot(split(z$transfer_probability,zn),las=1,
          horizontal=TRUE,col=colors['CONCERT'],ylim=c(0,1),xlab='CONCERT posterior transferable probability',
          main=gsub('_',' ',t),cex.axis=.8)
      },width=if(length(unique(zn))>6)10 else 7,height=max(5.5,.26*length(unique(zn))))
    }
    if(!is.null(st)&&nrow(st)) {
      z<-st[st$target==t&st$method=='TRACE',,drop=FALSE]
      z<-head(z[order(-z$selection_frequency,-z$mean_joint_inclusion_borrowing),],20)
      if(nrow(z))save_panel(file.path(d,'trace_gene_selection_stability'),function(){
        par(mar=c(5,12,3,1));barplot(rev(z$selection_frequency),names.arg=rev(paste(z$symbol,sub('\\..*','',z$gene_id))),
          horiz=TRUE,las=1,col='#D55E00',xlim=c(0,1),xlab='Selection frequency over successful splits',main='TRACE: selected genes',cex.names=.7)
      },height=7)
    }
    if(!is.null(pred)&&nrow(pred))for(method in methods) {
      z<-pred[pred$target==t&pred$method==method,,drop=FALSE]
      if(!nrow(z))next
      z<-aggregate(cbind(observed,predicted)~donor,data=z,FUN=mean)
      save_panel(file.path(d,paste0(gsub('-','_',method),'_observed_predicted')),function(){
        lim<-range(z$observed,z$predicted);par(mar=c(4,4,4,1))
        plot(z$observed,z$predicted,pch=16,col=adjustcolor(colors[method],.7),xlim=lim,ylim=lim,
          xlab='Observed SH2D2A log2(TPM + 1)',ylab='Mean held-out prediction',main=paste(method,gsub('_',' ',t)))
        abline(0,1,lty=2);mtext('Repeated predictions averaged within donor for display only',side=3,line=.3,cex=.7)
      })
    }
  }
  files<-list.files(fig,recursive=TRUE,full.names=TRUE)
  csv(data.frame(file=sub(paste0(out_dir,'/'),'',files,fixed=TRUE)),file.path(tab,'figure_index.csv'))
  message('Tables and individual PDF/PNG figures: ',out_dir)
}
