#!/usr/bin/env Rscript
# Aggregate saved replicate-level results, then reproduce the current paper.
entry<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
if(!file.exists(entry))entry<-gsub('~+~',' ',entry,fixed=TRUE)
root<-dirname(dirname(normalizePath(entry)))
args<-commandArgs(TRUE)
if(length(args)%%2L)stop('Options must be --name value pairs')
allowed<-c('--sim-root','--real-root','--out')
if(length(args)&&any(!args[seq(1,length(args),2)]%in%allowed))stop('Unknown option')
arg<-function(n,d){i<-match(paste0('--',n),args);if(is.na(i))d else args[i+1L]}
simroot<-normalizePath(arg('sim-root',file.path(root,'simulation/published_results')))
realroot<-normalizePath(arg('real-root',file.path(root,'real_data/published_results')))
if(dir.exists(file.path(realroot,'comparison/tables')))realroot<-file.path(realroot,'comparison')
out<-arg('out',file.path(root,'outputs/reproduced'))
dir.create(out,recursive=TRUE,showWarnings=FALSE);out<-normalizePath(out)
if(out%in%c(root,file.path(root,'manuscript'),simroot,realroot))stop('Use a separate output directory')
strict<-identical(simroot,normalizePath(file.path(root,'simulation/published_results')))&&
 identical(realroot,normalizePath(file.path(root,'real_data/published_results/comparison')))
for(d in c('data','figures','tables','checks'))dir.create(file.path(out,d),showWarnings=FALSE)
plotfiles<-list.files(file.path(root,'plotting'),pattern='[.]R$',full.names=TRUE)
stopifnot(all(file.copy(plotfiles,file.path(out,'figures'),overwrite=TRUE)))
folders<-c('01_TRACE_identifiable_mixture_contamination','02_TRACE_support_swap_contamination',
 '03_TRACE_mosaic_contamination','04_TRACE_complementary_sources_fixed_noise')
methods<-c('TRACE-EM','Oracle TRACE','ACVB','Trans-Lasso','TRADER','CONCERT')
checks<-list();summaries<-list()
for(folder in folders) {
 id<-substr(folder,1,2)
 e<-new.env(parent=.GlobalEnv)
 sys.source(file.path(root,'simulation/_common/simulation_50rep_runner.R'),e)
 if(id=='02')sys.source(file.path(root,'simulation',folder,'simulation_50rep_runner.R'),e)
 perf<-read.csv(file.path(simroot,folder,'comparison/performance_all_scenarios.csv'),stringsAsFactors=FALSE)
 stopifnot(!anyDuplicated(perf[c('scenario_index','rep','seed','method')]),setequal(unique(perf$method),methods))
 counts<-table(perf$scenario_index,perf$method)
 if(any(counts!=counts[,1]))stop('Unpaired method counts in Setting ',id)
 requested<-max(counts)
 summary<-e$sim_summarize(perf,requested)
 write.csv(summary,file.path(out,'data',paste0('setting',id,'_summary.csv')),row.names=FALSE,na='NA')
 summaries[[id]]<-summary
 if(strict) {
  expected<-read.csv(file.path(root,'manuscript/data',paste0('setting',id,'_summary.csv')),stringsAsFactors=FALSE)
  key<-function(x)paste(x$scenario_index,x$method,x$metric,sep=':')
  ix<-match(key(expected),key(summary));stopifnot(!anyNA(ix))
  for(nm in c('mean','sd','se','n_eff','n_completed'))
   stopifnot(isTRUE(all.equal(expected[[nm]],summary[[nm]][ix],tolerance=1e-9,check.attributes=FALSE)))
 }
 checks[[id]]<-data.frame(component=paste0('Setting ',id),cells=nrow(counts),replications=requested,
  methods=ncol(counts),rows=nrow(perf),matches_paper=if(strict)TRUE else NA)
 cat('Aggregated Setting',id,':',nrow(perf),'replicate-method records\n')
}
interval<-function(x,digits=3L){
 x<-x[is.finite(x)];stopifnot(length(x)>0L)
 a<-min(x);b<-max(x);fmt<-paste0('%.',digits,'f')
 if(abs(a-b)<1e-10)sprintf(fmt,a)else paste(sprintf(fmt,a),sprintf(fmt,b),sep='--')
}
get_metric<-function(s,m,metric,filter=rep(TRUE,nrow(s)))
 s$mean[s$method==m&s$metric==metric&filter]
table2<-list();tex2<-character()
for(design in c('I','III'))for(m in methods) {
 d<-summaries[['02']];flt<-d$paper_setting_grid==design
 vals<-vapply(c('precision','recall','pred_mse'),function(a)interval(get_metric(d,m,a,flt)),'')
 table2[[length(table2)+1L]]<-data.frame(Design=design,Method=m,Precision=vals[1],Recall=vals[2],Prediction_MSE=vals[3])
 tex2<-c(tex2,paste0(paste(c(design,m,vals),collapse=' & '),' \\\\'))
}
write.csv(do.call(rbind,table2),file.path(out,'tables/setting02_additional_metrics.csv'),row.names=FALSE)
writeLines(tex2,file.path(out,'tables/setting02_additional_metrics_rows.tex'))
d<-summaries[['04']];table4<-list();tex4<-character()
for(s in c(10,15,30)) {
 a<-d[d$nzc_grid==s&d$method=='TRACE-EM'&d$metric=='f1',]
 b<-d[d$nzc_grid==s&d$method=='Oracle TRACE'&d$metric=='f1',]
 ix<-match(a$common_bad_frac_grid,b$common_bad_frac_grid);stopifnot(!anyNA(ix))
 vals<-c(interval(get_metric(d,'TRACE-EM','alpha_l1_error',d$nzc_grid==s)),
  interval(get_metric(d,'TRACE-EM','borrowing_f1',d$nzc_grid==s)),sprintf('%.4f',max(abs(a$mean-b$mean[ix]))))
 table4[[length(table4)+1L]]<-data.frame(s=s,Alpha_L1=vals[1],Transfer_F1=vals[2],F1_difference=vals[3])
 tex4<-c(tex4,paste0(paste(c(s,vals),collapse=' & '),' \\\\'))
}
write.csv(do.call(rbind,table4),file.path(out,'tables/setting04_diagnostics.csv'),row.names=FALSE)
writeLines(tex4,file.path(out,'tables/setting04_diagnostics_rows.tex'))

m<-read.csv(file.path(realroot,'tables/metrics_per_rep.csv'),stringsAsFactors=FALSE)
stopifnot(!anyDuplicated(m[c('target','rep','method')]),all(m$status=='ok'))
counts<-table(m$target,m$method);stopifnot(ncol(counts)==5L,all(counts==counts[,1]))
if(strict) {
 expected<-read.csv(file.path(root,'manuscript/data/real_data_metrics_per_rep.csv'),stringsAsFactors=FALSE)
 key<-function(x)paste(x$target,x$rep,x$method,sep=':')
 ix<-match(key(expected),key(m));stopifnot(!anyNA(ix))
 for(nm in c('mse','model_size','null_mse'))stopifnot(isTRUE(all.equal(expected[[nm]],m[[nm]][ix],tolerance=1e-12)))
}
write.csv(m,file.path(out,'data/real_data_metrics_per_rep.csv'),row.names=FALSE,na='')
csv<-function(x,p)write.csv(x,p,row.names=FALSE,na='')
source(file.path(root,'real_data/R/summarize.R'))
source(file.path(root,'plotting/plot_real_data_figures.R'))
realout<-file.path(out,'real_data');dir.create(file.path(realout,'tables'),recursive=TRUE,showWarnings=FALSE)
csv(m,file.path(realout,'tables/metrics_per_rep.csv'))
write_paper_outputs(realout,list(targets=names(.tissue_labels)))
for(n in c('gtex_mse_table.csv','gtex_mse_table.tex','model_size_table.csv','model_size_table.tex'))
 stopifnot(file.copy(file.path(realout,'tables',n),file.path(out,'tables',n),overwrite=TRUE))
if(strict)for(n in c('gtex_mse_table','model_size_table')) {
 expected<-read.csv(file.path(root,'manuscript/data',paste0('real_data_',n,'.csv')),check.names=FALSE)
 got<-read.csv(file.path(out,'tables',paste0(n,'.csv')),check.names=FALSE)
 stopifnot(isTRUE(all.equal(expected,got,check.attributes=FALSE)))
}
checks[['real']]<-data.frame(component='Real data',cells=nrow(counts),replications=max(counts),methods=5L,
 rows=nrow(m),matches_paper=if(strict)TRUE else NA)
write.csv(do.call(rbind,checks),file.path(out,'checks/result_counts_and_agreement.csv'),row.names=FALSE)

# Plotting occurs in a separate R process so each plotting script resolves
# its own --file path correctly. Only derived outputs are overwritten.
status<-system2(file.path(R.home('bin'),'Rscript'),shQuote(file.path(out,'figures/reproduce_all_figures.R')))
if(status!=0L)stop('Figure reproduction failed with exit status ',status)
writeLines(capture.output(sessionInfo()),file.path(out,'checks/sessionInfo.txt'))
cat('Reproduced figures, tables and checks:',out,'\n')
