a<-commandArgs(TRUE);root<-normalizePath(a[1]);out<-normalizePath(a[2],mustWork=FALSE)
entry<-file.path(root,'run_real_data.R')
status<-system2(file.path(R.home('bin'),'Rscript'),c(shQuote(entry),'--reps','1','--workers','2','--out',shQuote(out)))
stopifnot(status==0L)
newfiles<-list.files(file.path(out,'trace_refit/jobs'),pattern='rds$',full.names=TRUE)
stopifnot(length(newfiles)==6L)
for(f in newfiles){
 x<-readRDS(f);b<-readRDS(file.path(root,'retained_results/jobs',basename(f)))
 z<-readRDS(file.path(out,'comparison/jobs',basename(f)))
 stopifnot(x$status=='completed',nrow(x$metrics)==1L,all(x$metrics$method=='TRACE'),x$input_reused,
   all(x$initialization_max_error<1e-7),nrow(z$metrics)==5L)
 for(nm in c('metrics','predictions','features')){
  q<-z[[nm]][z[[nm]]$method!='TRACE',,drop=FALSE];rownames(q)<-rownames(b[[nm]])
  stopifnot(identical(q,b[[nm]]))
 }
 stopifnot(identical(z$concert_transfer,b$concert_transfer),identical(z$concert_prior,b$concert_prior))
}
before<-tools::md5sum(newfiles)
status<-system2(file.path(R.home('bin'),'Rscript'),c(shQuote(entry),'--reps','1','--workers','2','--out',shQuote(out),'--skip-summary','true'))
stopifnot(status==0L,identical(before,tools::md5sum(newfiles)))
stopifnot(file.exists(file.path(out,'comparison/tables/gtex_mse_table.tex')),
 file.exists(file.path(out,'comparison/paper_figures/real_data_amygdala_mse.pdf')),
 file.exists(file.path(out,'comparison/paper_figures/real_data_all_tissues_mse.png')))
cat('PASS: six real-tissue fits, untouched competitors, cached inputs, old CV seed, resume, CSV/TeX/PDF/PNG output\n')
