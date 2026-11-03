# Reconstruct raw inputs and fit one target using the full-analysis entry.
entry<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
if(!file.exists(entry))entry<-gsub('~+~',' ',entry,fixed=TRUE)
root<-dirname(dirname(normalizePath(entry)))
out<-tempfile('TRACE_full_real_entry_')
status<-system2(file.path(R.home('bin'),'Rscript'),c(
 shQuote(file.path(root,'full_analysis/run_real_data.R')),
 '--targets','brain_amygdala','--reps','1','--workers','1','--methods','TRACE',
 '--data-dir',shQuote(file.path(root,'data/processed/gtex')),'--out',shQuote(out)))
stopifnot(status==0L)
x<-readRDS(file.path(out,'jobs/brain_amygdala_rep001.rds'))
stopifnot(x$status=='completed',nrow(x$metrics)==1L,x$metrics$status=='ok',
 x$metrics$method=='TRACE',is.finite(x$metrics$mse),file.exists(file.path(out,'tables/metrics_per_rep.csv')))
cat('PASS: full real-data entry reconstructs raw inputs and fits TRACE.\n')
