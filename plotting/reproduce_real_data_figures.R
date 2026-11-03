# Rebuild real-data MSE displays from the corrected per-split metric table.
entry<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
if(!file.exists(entry))entry<-gsub('~+~',' ',entry,fixed=TRUE)
figure_dir<-dirname(normalizePath(entry));project<-dirname(figure_dir)
if (basename(figure_dir) == "plotting") project <- file.path(project, "manuscript")
source(file.path(figure_dir,'plot_real_data_figures.R'))
csv<-function(x,p)write.csv(x,p,row.names=FALSE,na='')
save_panel<-function(path,draw,width=7,height=5.5){
 pdf(paste0(path,'.pdf'),width=width,height=height,useDingbats=FALSE)
 tryCatch(draw(),finally=dev.off())
 png(paste0(path,'.png'),width=width,height=height,units='in',res=300,
     type=if(Sys.info()[['sysname']]=='Darwin')'quartz'else'cairo')
 tryCatch(draw(),finally=dev.off())
}
tmp<-tempfile();dir.create(file.path(tmp,'tables'),recursive=TRUE)
file.copy(file.path(project,'data/real_data_metrics_per_rep.csv'),file.path(tmp,'tables/metrics_per_rep.csv'))
write_paper_outputs(tmp,list(targets=names(.tissue_labels)))
for(f in list.files(file.path(tmp,'paper_figures'),full.names=TRUE))file.copy(f,file.path(figure_dir,basename(f)),overwrite=TRUE)
unlink(tmp,recursive=TRUE)
