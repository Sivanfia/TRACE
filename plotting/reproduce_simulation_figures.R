# Rebuild main simulation figures from the corrected summaries kept in this project.
entry<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
if(!file.exists(entry))entry<-gsub('~+~',' ',entry,fixed=TRUE)
figure_dir<-dirname(normalizePath(entry));project<-dirname(figure_dir)
if (basename(figure_dir) == "plotting") project <- file.path(project, "manuscript")
source(file.path(figure_dir,'plot_paper_figures.R'))
for(id in c('01','02','03','04')) {
 d<-read.csv(file.path(project,'data',paste0('setting',id,'_summary.csv')),stringsAsFactors=FALSE)
 stopifnot(all(d$n_completed==50),all(d$n_requested==50))
 tmp<-tempfile();dir.create(tmp)
 paper_make_figures(d,tmp,id)
 for(ext in c('pdf','png'))stopifnot(file.copy(file.path(tmp,'paper_figures',paste0('setting',id,'_main_f1.',ext)),file.path(figure_dir,paste0('setting',id,'_main_f1.',ext)),overwrite=TRUE))
 unlink(tmp,recursive=TRUE)
}
