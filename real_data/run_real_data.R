#!/usr/bin/env Rscript
args<-commandArgs(TRUE)
if(length(args)%%2L)stop('Arguments must be --name value pairs')
keys<-if(length(args))args[seq.int(1,length(args),2L)]else character()
allowed<-c('--reps','--workers','--targets','--out','--summarize-only','--skip-summary')
if(any(!keys%in%allowed)||anyDuplicated(keys))stop('Unknown or repeated option: ',paste(setdiff(keys,allowed),collapse=', '))
arg<-function(n,default=NULL){i<-match(paste0('--',n),args);if(is.na(i))default else args[i+1L]}
entry<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
if(!file.exists(entry))entry<-gsub('~+~',' ',entry,fixed=TRUE)
ROOT<-dirname(normalizePath(entry))
source(file.path(ROOT,'config/default.R'));cfg<-CONFIG
reps<-as.integer(arg('reps',50L));workers<-as.integer(arg('workers',48L))
stopifnot(length(reps)==1L,!is.na(reps),reps>=1L,reps<=50L,!is.na(workers),workers>=1L)
targets<-strsplit(arg('targets',paste(cfg$targets,collapse=',')),',',fixed=TRUE)[[1]]
stopifnot(length(targets)>0L,all(targets%in%cfg$targets),!anyDuplicated(targets))
required<-c('glmnet','quadprog')
missing<-required[!vapply(required,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing))stop('Missing packages: ',paste(missing,collapse=', '),'. Run Rscript scripts/install_dependencies.R')
source(file.path(ROOT,'vendor/trace_acvb_kernels.R'))
for(f in c('preprocess.R','methods.R','retained_inputs.R','refit_runner.R','summarize.R','paper_outputs.R'))source(file.path(ROOT,'R',f))
data_dir<-file.path(ROOT,'data/processed/gtex')
retained_dir<-file.path(ROOT,'retained_results/jobs')
original_manifest<-readRDS(file.path(ROOT,'retained_results/original_manifest.rds'))
check_retained_data_files(data_dir,original_manifest,cfg$targets)
files<-c(file.path(ROOT,'run_real_data.R'),list.files(file.path(ROOT,'R'),pattern='\\.R$',full.names=TRUE),
 file.path(ROOT,'vendor/trace_acvb_kernels.R'),file.path(ROOT,'config/default.R'),
 list.files(file.path(ROOT,'retained_results'),pattern='\\.(rds|csv)$',recursive=TRUE,full.names=TRUE),
 list.files(data_dir,pattern='\\.gz$',full.names=TRUE))
hashes<-tools::md5sum(files);names(hashes)<-substring(files,nchar(ROOT)+2L)
versions<-c(R=paste(R.version$major,R.version$minor,sep='.'),vapply(required,function(p)as.character(packageVersion(p)),''))
tmp<-tempfile();saveRDS(list(config=cfg,files=hashes,versions=versions,mode='TRACE-only-corrected-T0'),tmp)
signature<-unname(tools::md5sum(tmp));unlink(tmp)
out_dir<-arg('out',file.path(ROOT,'results_T0_corrected'))
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE);out_dir<-normalizePath(out_dir)
manifest_path<-file.path(out_dir,'run_manifest.rds')
if(file.exists(manifest_path)){
 if(!identical(readRDS(manifest_path)$signature,signature))stop('Code/data/config/environment changed. Use a new --out directory.')
}else atomic_rds(list(signature=signature,config=cfg,files=hashes,versions=versions,session=capture.output(sessionInfo())),manifest_path)
for(d in c('trace_refit/jobs','comparison/jobs'))dir.create(file.path(out_dir,d),recursive=TRUE,showWarnings=FALSE)
jobs<-expand.grid(rep=seq_len(reps),target=targets,stringsAsFactors=FALSE)
paths<-file.path(retained_dir,sprintf('%s_rep%03d.rds',jobs$target,jobs$rep))
if(any(!file.exists(paths)))stop('Missing retained split: ',paths[which(!file.exists(paths))[1]])
failed<-FALSE
if(arg('summarize-only','false')!='true'){
 fun<-function(i)run_trace_job(paths[i],cfg,data_dir,out_dir,signature)
 result<-if(workers>1L&&.Platform$OS.type!='windows')parallel::mclapply(seq_along(paths),fun,mc.cores=workers,mc.preschedule=FALSE,mc.set.seed=FALSE)else lapply(seq_along(paths),fun)
 failed<-any(vapply(result,function(x)inherits(x,'try-error')||!identical(x$status,'completed'),logical(1)))
}
fresh_paths<-list.files(file.path(out_dir,'trace_refit/jobs'),pattern='\\.rds$',full.names=TRUE)
status<-lapply(fresh_paths,function(f){
 x<-readRDS(f)
 if(!identical(x$signature,signature))stop('Mixed fit signatures: ',f)
 cp<-file.path(out_dir,'comparison/jobs',basename(f))
 if(identical(x$status,'completed'))atomic_rds(combine_real_job(x,readRDS(file.path(retained_dir,basename(f))),signature),cp)
 else if(file.exists(cp))unlink(cp)
 data.frame(target=x$target,rep=x$rep,status=x$status,error=x$error%||%'')
})
if(length(status))csv(do.call(rbind,status),file.path(out_dir,'job_status.csv'))
cp_paths<-list.files(file.path(out_dir,'comparison/jobs'),pattern='\\.rds$',full.names=TRUE)
if(any(!basename(cp_paths)%in%basename(fresh_paths)))stop('Comparison job has no corresponding new TRACE checkpoint')
if(arg('skip-summary','false')!='true'&&length(cp_paths)){
 summarize_real_results(file.path(out_dir,'comparison'),signature,c('TRACE','ACVB','TRADER','Trans-Lasso','CONCERT'))
 write_paper_outputs(file.path(out_dir,'comparison'),cfg)
}
message('Output: ',out_dir)
if(failed)quit(status=2L)
if(!length(cp_paths))stop('No successful new TRACE fits available for comparison')
