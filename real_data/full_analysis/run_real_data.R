#!/usr/bin/env Rscript
args<-commandArgs(trailingOnly=TRUE)
arg<-function(n,default=NULL) {i<-match(paste0('--',n),args);if(is.na(i))default else {if(i==length(args))stop('Missing --',n);args[i+1L]}}
this<-sub('^--file=','',commandArgs()[grepl('^--file=',commandArgs())][1]);if(!file.exists(this))this<-gsub('~+~',' ',this,fixed=TRUE);ROOT<-dirname(normalizePath(this))
source(arg('config',file.path(ROOT,'config/default.R')))
cfg<-CONFIG
cfg$reps<-as.integer(arg('reps',cfg$reps));cfg$source_mode<-arg('sources',cfg$source_mode)
cfg$donor_policy<-arg('donor-policy',cfg$donor_policy);cfg$trans_backend<-arg('trans-backend',cfg$trans_backend)
cfg$max_features<-as.integer(arg('max-features',cfg$max_features));cfg$source_bootstrap<-as.integer(arg('bootstrap',cfg$source_bootstrap))
cfg$target_init_steps<-as.integer(arg('target-init',cfg$target_init_steps))
cfg$trader_n_mcmc<-as.integer(arg('mcmc',cfg$trader_n_mcmc));cfg$trader_burnin<-as.integer(arg('burnin',cfg$trader_burnin))
cfg$concert_max_iter<-as.integer(arg('concert-max-iter',cfg$concert_max_iter))
ALL_METHODS<-c('TRACE','ACVB','TRADER','Trans-Lasso','CONCERT')
METHODS<-strsplit(arg('methods',paste(ALL_METHODS,collapse=',')),',',fixed=TRUE)[[1]]
if(!length(METHODS)||any(!METHODS%in%ALL_METHODS)||anyDuplicated(METHODS))
  stop('--methods must be a comma-separated subset of: ',paste(ALL_METHODS,collapse=', '))
cfg$methods<-METHODS
stopifnot(cfg$reps>=1L,cfg$source_bootstrap>=2L,cfg$max_features>=2L,
  cfg$source_mode%in%c('brain5','all48'),cfg$donor_policy%in%c('exclude_test','disjoint_target'),
  cfg$trans_backend%in%c('official','simulation_ranked'),
  cfg$precision_mode%in%c('bootstrap_capped','min_source_n'),
  cfg$trader_burnin>=0,cfg$trader_burnin<cfg$trader_n_mcmc,
  cfg$concert_q0>0,cfg$concert_q0<1,cfg$concert_eta>0,cfg$concert_tau>0,
  cfg$concert_qk_tolerance>=0,cfg$concert_probability_epsilon>0,
  cfg$concert_probability_epsilon<0.5,cfg$concert_threshold>0,cfg$concert_max_iter>=1L)
targets<-strsplit(arg('targets',paste(cfg$targets,collapse=',')),',',fixed=TRUE)[[1]]
stopifnot(all(targets%in%cfg$targets),!anyDuplicated(targets))
data_dir<-normalizePath(arg('data-dir',file.path(ROOT,'data/processed/gtex')),mustWork=TRUE)
out_dir<-arg('out',file.path(ROOT,'results_50rep'))
dir.create(file.path(out_dir,'jobs'),recursive=TRUE,showWarnings=FALSE)
out_dir<-normalizePath(out_dir)
required<-unique(c('glmnet','quadprog',
  if('TRADER'%in%METHODS)c('Rcpp','RcppEigen','spikeslab'),
  if('CONCERT'%in%METHODS)c('Rcpp','RcppArmadillo')))
missing<-required[!vapply(required,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing))stop('Missing packages: ',paste(missing,collapse=', '),'. Run scripts/install_dependencies.R')
source(file.path(ROOT,'vendor/trace_acvb_kernels.R'))
if('TRADER'%in%METHODS)source(file.path(ROOT,'vendor/trader_bridge.R'))
if('Trans-Lasso'%in%METHODS)source(file.path(ROOT,'vendor/TransLasso-functions.R'))
if('CONCERT'%in%METHODS)source(file.path(ROOT,'vendor/concert_bridge.R'))
for(f in c('preprocess.R','methods.R','runner.R','summarize.R'))source(file.path(ROOT,'R',f))
rfiles<-list.files(file.path(ROOT,'R'),full.names=TRUE)
files<-c(rfiles[basename(rfiles)!='summarize.R'],
  list.files(file.path(ROOT,'vendor'),pattern='\\.(R|cpp)$',full.names=TRUE,recursive=TRUE),
  file.path(ROOT,'run_real_data.R'),list.files(data_dir,pattern='\\.cns_tpm\\.tsv\\.gz$',full.names=TRUE))
hashes<-tools::md5sum(files);names(hashes)<-sub(paste0('^',ROOT,'/'),'',files)
# Job selectors and execution location are excluded; fitting options and data are included.
sigcfg<-cfg;sigcfg$reps<-NULL
versions<-c(R=paste(R.version$major,R.version$minor,sep='.'),vapply(required,function(p)as.character(packageVersion(p)),''))
tmp<-tempfile();saveRDS(list(config=sigcfg,hashes=hashes,versions=versions),tmp);signature<-unname(tools::md5sum(tmp));unlink(tmp)
meta<-file.path(out_dir,paste0('run_manifest_',signature,'.rds'))
if(!file.exists(meta))atomic_rds(list(signature=signature,config=cfg,files=hashes,
  command=commandArgs(),session=capture.output(sessionInfo())),meta)
jobs<-expand.grid(rep=seq_len(cfg$reps),target=targets,stringsAsFactors=FALSE)
index<-as.integer(arg('job-index',NA_integer_))
if(!is.na(index)){if(index<1||index>nrow(jobs))stop('job-index out of range');jobs<-jobs[index,,drop=FALSE]}
workers<-as.integer(arg('workers',1L));stopifnot(workers>=1L)
if(arg('summarize-only','false')!='true') {
  fun<-function(i)run_real_job(jobs$target[i],jobs$rep[i],cfg,data_dir,out_dir,signature,
    retry_failed=arg('retry-failed','false')=='true')
  result<-if(workers>1L&&.Platform$OS.type!='windows')parallel::mclapply(seq_len(nrow(jobs)),fun,
    mc.cores=workers,mc.preschedule=FALSE,mc.set.seed=FALSE)else lapply(seq_len(nrow(jobs)),fun)
  if(any(vapply(result,function(x)inherits(x,'try-error'),logical(1))))stop('A worker failed; inspect logs')
}
# Slurm array tasks only write their own checkpoint; aggregate once after the array.
if(is.na(index)) {
  summarize_real_results(out_dir,signature,METHODS)
  atomic_rds(list(fit_signature=signature,plot_code_md5=unname(tools::md5sum(file.path(ROOT,'R/summarize.R'))),
    time=as.character(Sys.time())),file.path(out_dir,'plot_manifest.rds'))
}
if(exists('result')&&any(vapply(result,function(x)x$status=='error'||any(x$metrics$status!='ok'),logical(1))))
  quit(status=2L)
