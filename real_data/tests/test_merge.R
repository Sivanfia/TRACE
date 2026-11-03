a<-commandArgs(TRUE);root<-normalizePath(a[1]);source(file.path(root,'R/retained_inputs.R'));source(file.path(root,'R/refit_runner.R'))
b<-readRDS(list.files(file.path(root,'retained_results/jobs'),pattern='rds$',full.names=TRUE)[1])
x<-list(target=b$target,rep=b$rep,seed=b$seed,status='completed',signature='test',
  metrics=b$metrics[1,,drop=FALSE],predictions=b$predictions[b$predictions$method=='ACVB',,drop=FALSE],
  features=b$features[b$features$method=='ACVB',,drop=FALSE],alpha=data.frame(),fits=list())
for(nm in c('metrics','predictions','features'))x[[nm]]$method<-'TRACE'
x$metrics$mse<-0.123456
z<-combine_real_job(x,b,'test')
stopifnot(nrow(z$metrics)==5,identical(z$metrics$mse[z$metrics$method=='TRACE'],0.123456))
for(nm in c('metrics','predictions','features')){
 got<-z[[nm]][z[[nm]]$method!='TRACE',,drop=FALSE]
 rownames(got)<-rownames(b[[nm]])
 stopifnot(identical(got,b[[nm]]))
}
stopifnot(identical(z$concert_transfer,b$concert_transfer),identical(z$concert_prior,b$concert_prior))
fails<-function(new=x,old=b){stopifnot(inherits(try(combine_real_job(new,old,'test'),silent=TRUE),'try-error'))}
bad<-x;bad$seed<-bad$seed+1L;fails(bad)
bad<-x;bad$predictions$observed[1]<-bad$predictions$observed[1]+1;fails(bad)
bad<-x;bad$metrics<-rbind(bad$metrics,bad$metrics);fails(bad)
bad<-x;bad$status<-'error';fails(bad)
bad<-b;bad$metrics$method[1]<-'TRACE';fails(old=bad)
cat('PASS: unchanged baseline values, TRACE replacement, pair/observation/status/duplicate rejection\n')
tmp<-tempfile('resume_guard_');dir.create(file.path(tmp,'trace_refit/jobs'),recursive=TRUE)
f<-list.files(file.path(root,'retained_results/jobs'),pattern='rds$',full.names=TRUE)[1]
checkpoint<-file.path(tmp,'trace_refit/jobs',basename(f));bad<-x;bad$signature<-'different_code';saveRDS(bad,checkpoint)
h<-tools::md5sum(checkpoint)
stopifnot(inherits(try(run_trace_job(f,list(),'',tmp,'test'),silent=TRUE),'try-error'),identical(h,tools::md5sum(checkpoint)))
unlink(tmp,recursive=TRUE)
cat('PASS: changed checkpoint fingerprint is rejected without overwriting retained output\n')
