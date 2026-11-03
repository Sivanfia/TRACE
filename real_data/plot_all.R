#!/usr/bin/env Rscript
entry<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
if(!file.exists(entry))entry<-gsub('~+~',' ',entry,fixed=TRUE)
root<-dirname(dirname(normalizePath(entry)))
a<-commandArgs(TRUE)
opts<-if(length(a))c('--real-root',shQuote(a[1]))else character()
status<-system2(file.path(R.home('bin'),'Rscript'),c(shQuote(file.path(root,'scripts/reproduce_results.R')),opts))
if(status!=0L)quit(status=status)
