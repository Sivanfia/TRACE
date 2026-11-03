# Compact manuscript outputs; all inference used the unchanged donor-level protocol.
.tissue_labels <- c(brain_amygdala='Amygdala',brain_anterior_cingulate_cortex_ba24='Anterior cingulate cortex',
 brain_cortex='Cortex',brain_hippocampus='Hippocampus',brain_nucleus_accumbens_basal_ganglia='Nucleus accumbens',brain_substantia_nigra='Substantia nigra')
write_paper_outputs <- function(out_dir,cfg) {
 tab<-file.path(out_dir,'tables');fig<-file.path(out_dir,'paper_figures');dir.create(fig,showWarnings=FALSE)
 m<-read.csv(file.path(tab,'metrics_per_rep.csv'),stringsAsFactors=FALSE,check.names=FALSE)
 methods<-c('TRACE','ACVB','TRADER','Trans-Lasso','CONCERT')
 labels<-c('TRACE-EM','ACVB','TRADER','Trans-Lasso','CONCERT')
 targets<-cfg$targets[cfg$targets%in%m$target]
 # Baselines are kept only on the same successfully recomputed donor splits.
 if(any(duplicated(m[,c('target','rep','method')])))stop('Duplicate target/rep/method in paper tables')
 counts<-table(m$target,m$method)
 if(any(counts!=counts[,1]))stop('Unequal paired method counts in paper tables')
 for(metric in c('mse','model_size')){
  result<-data.frame(Tissue=unname(.tissue_labels[targets]),N=as.integer(counts[targets,1]),check.names=FALSE)
  tex<-c('\\begin{table}[htbp]','\\centering',
   if(metric=='mse')'\\caption{Held-out MSE (mean and standard deviation) across repeated donor splits.}'else'\\caption{Number of selected variables (mean and standard deviation).}',
   if(metric=='mse')'\\label{tab:gtex-mse}'else'\\label{tab:gtex-model-size}',
   '\\begin{tabular}{lrrrrr}','\\toprule','Target tissue & TRACE-EM & ACVB & TRADER & Trans-Lasso & CONCERT \\\\','\\midrule')
  for(t in targets){
   vals<-lapply(methods,function(a)m[[metric]][m$target==t&m$method==a&m$status=='ok'])
   means<-vapply(vals,mean,0.0);sds<-vapply(vals,function(x)if(length(x)>1)sd(x)else NA_real_,0.0)
   text<-if(metric=='mse')sprintf('%.4f (%.4f)',means,sds)else sprintf('%.1f (%.1f)',means,sds)
   for(i in seq_along(methods))result[match(t,targets),labels[i]]<-text[i]
   if(metric=='mse')text[which.min(means)]<-paste0('\\textbf{',text[which.min(means)],'}')
   tex<-c(tex,paste0(unname(.tissue_labels[t]),' & ',paste(text,collapse=' & '),' \\\\'))
  }
  tex<-c(tex,'\\bottomrule','\\end{tabular}','\\end{table}')
  csv(result,file.path(tab,if(metric=='mse')'gtex_mse_table.csv'else'model_size_table.csv'))
  writeLines(tex,file.path(tab,if(metric=='mse')'gtex_mse_table.tex'else'model_size_table.tex'))
 }
 colors<-c('#D55E00','#0072B2','#009E73','#CC79A7','#56B4E9')
 draw_one<-function(t,small=FALSE){
  z<-m[m$target==t&m$status=='ok',,drop=FALSE]
  groups<-lapply(methods,function(a)z$mse[z$method==a]);names(groups)<-labels
  par(mar=if(small)c(3.6,3.2,2.0,.6)else c(4,4,2,.7),mgp=c(1.8,.55,0),cex.axis=if(small).65 else .85)
  boxplot(groups,col=colors,ylab='Test MSE',main=unname(.tissue_labels[t]),outline=FALSE,ylim=range(z$mse,finite=TRUE),
   cex.main=if(small).9 else 1,cex.lab=if(small).85 else 1)
  for(i in seq_along(groups))points(rep(i,length(groups[[i]])),groups[[i]],pch=16,cex=.45,col=adjustcolor('black',.2))
 }
 for(t in targets){
  fname<-switch(t,brain_amygdala='real_data_amygdala_mse',brain_hippocampus='real_data_hippocampus_mse',paste0(t,'_mse'))
  save_panel(file.path(fig,fname),function()draw_one(t),width=6.2,height=3.6)
 }
 save_panel(file.path(fig,'real_data_all_tissues_mse'),function(){par(mfrow=c(2,3));for(t in targets)draw_one(t,TRUE)},width=9,height=5.8)
 origin<-data.frame(target=m$target,rep=m$rep,method=m$method,
  result_origin=ifelse(m$method=='TRACE','recomputed_T0','retained_original'))
 csv(origin,file.path(tab,'result_provenance.csv'))
 message('Manuscript tables: ',tab,'; figures: ',fig)
}
