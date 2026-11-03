entry<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
if(!file.exists(entry))entry<-gsub('~+~',' ',entry,fixed=TRUE)
root<-dirname(normalizePath(entry))
args<-commandArgs(TRUE)
base<-if(length(args))normalizePath(args[1],mustWork=TRUE) else file.path(root,'published_results')
source(file.path(root,'_common/simulation_50rep_runner.R'))
source(file.path(root,'_common/paper_figures.R'))
folders<-c('01_TRACE_identifiable_mixture_contamination','02_TRACE_support_swap_contamination',
 '03_TRACE_mosaic_contamination','04_TRACE_complementary_sources_fixed_noise')
dir.create(file.path(base,'paper_figures'),recursive=TRUE,showWarnings=FALSE)
for(folder in folders){
 id<-substr(folder,1,2);comparison<-file.path(base,folder,'comparison')
 path<-file.path(comparison,'summary_long_all_scenarios.csv')
 if(!file.exists(path)){message('No completed combined results: ',folder);next}
 source(file.path(root,'_common/plot_results.R'))
 if(id=='02')source(file.path(root,folder,'plot_results.R'))
 summary<-read.csv(path,stringsAsFactors=FALSE)
 sim_make_plots(summary,comparison,id)
 paper_make_figures(summary,comparison,id)
 files<-list.files(file.path(comparison,'paper_figures'),pattern='[.](pdf|png)$',full.names=TRUE)
 file.copy(files,file.path(base,'paper_figures'),overwrite=TRUE)
}
cat('Paper figures: ',file.path(base,'paper_figures'),'\n',sep='')
