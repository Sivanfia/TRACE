#!/usr/bin/env Rscript
entry<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
if(!file.exists(entry))entry<-gsub('~+~',' ',entry,fixed=TRUE)
root<-dirname(dirname(normalizePath(entry)))
scripts<-unlist(lapply(c('scripts','simulation','real_data'),function(d)
 list.files(file.path(root,d),pattern='[.]R$',full.names=TRUE,recursive=TRUE)))
scripts<-scripts[!grepl('/(published_results|retained_results)/',scripts)]
for(f in scripts)parse(f)
cat('PASS: parsed',length(scripts),'R files\n')
folders<-c('01_TRACE_identifiable_mixture_contamination','02_TRACE_support_swap_contamination',
 '03_TRACE_mosaic_contamination','04_TRACE_complementary_sources_fixed_noise')
cells<-c(48L,20L,48L,21L);total<-0L
for(i in seq_along(folders)) {
 p<-file.path(root,'simulation/published_results',folders[i],'comparison/performance_all_scenarios.csv')
 d<-read.csv(p,stringsAsFactors=FALSE)
 stopifnot(!anyDuplicated(d[c('scenario_index','rep','method')]),
  nrow(d)==cells[i]*50L*6L,all(table(d$scenario_index,d$method)==50L),
  all(is.finite(d$f1)),all(d$f1>=0&d$f1<=1),all(is.finite(d$pred_mse)))
 baseline<-read.csv(file.path(root,'simulation/retained_results',substr(folders[i],1,2),'performance_baseline.csv'))
 stopifnot(!any(baseline$method%in%c('TRACE-EM','Oracle TRACE')),all(table(baseline$scenario_index,baseline$method)==50L))
 total<-total+nrow(d)
}
stopifnot(total==41100L)
cat('PASS: 137 simulation cells, 50 replications, six methods, 41100 records\n')
p<-file.path(root,'real_data/published_results/comparison/tables/metrics_per_rep.csv')
m<-read.csv(p,stringsAsFactors=FALSE)
stopifnot(nrow(m)==1500L,!anyDuplicated(m[c('target','rep','method')]),
 all(table(m$target,m$method)==50L),all(m$status=='ok'),all(is.finite(m$mse)))
files<-list.files(file.path(root,'real_data/published_results/comparison/jobs'),pattern='[.]rds$',full.names=TRUE)
stopifnot(length(files)==300L)
jobs<-do.call(rbind,lapply(files,function(f){x<-readRDS(f);stopifnot(x$status=='completed');x$metrics}))
key<-function(d)paste(d$target,d$rep,d$method,sep=':')
ix<-match(key(m),key(jobs));stopifnot(!anyNA(ix))
for(nm in c('mse','model_size','null_mse'))stopifnot(isTRUE(all.equal(m[[nm]],jobs[[nm]][ix],tolerance=1e-12)))
cat('PASS: all 1500 real-data records match the 300 saved five-method checkpoints\n')
expected<-read.csv(file.path(root,'provenance/manuscript_sha256.csv'),stringsAsFactors=FALSE)
stopifnot(all(file.exists(file.path(root,expected$file))))
cat('PASS: current manuscript snapshots and figure sources are present\n')
