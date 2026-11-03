#!/usr/bin/env Rscript
.entry <- local({
  frames <- Filter(function(x)is.character(x)&&length(x)==1L&&nzchar(x),lapply(sys.frames(),function(x)x$ofile))
  if(length(frames))tail(frames,1)[[1]] else sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
})
if(!file.exists(.entry)).entry<-gsub('~+~',' ',.entry,fixed=TRUE)
.project <- dirname(normalizePath(.entry,mustWork=TRUE))
source(file.path(.project,'..','_common','start_setting.R'),local=.GlobalEnv)
bundle_start_setting(.project)
