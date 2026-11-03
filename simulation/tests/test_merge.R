entry<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
root<-dirname(dirname(normalizePath(entry)))
source(file.path(root,'_common/simulation_50rep_runner.R'))
source(file.path(root,'_common/retained_results.R'))
expect_error<-function(expr){ok<-FALSE;tryCatch(force(expr),error=function(e)ok<<-TRUE);stopifnot(ok)}
old<-data.frame(scenario_index=c(1,1,1,1),rep=c(1,1,2,2),seed=c(101,101,102,102),
                method=c('ACVB','TRADER','ACVB','TRADER'),f1=c(.6,.8,.4,.7))
fresh<-data.frame(scenario_index=c(1,1),rep=c(1,1),seed=c(101,101),
                  method=c('TRACE-EM','Oracle TRACE'),f1=c(.9,1))
z<-retained_combine(fresh,old)
stopifnot(nrow(z)==4,all(z$rep==1),identical(z$f1[z$method=='ACVB'],.6),
          identical(z$f1[z$method=='TRADER'],.8),all(z$result_origin[z$method%in%c('ACVB','TRADER')]=='retained_baseline'))
expect_error(retained_combine(fresh[1,,drop=FALSE],old))
expect_error(retained_combine(rbind(fresh,fresh[1,,drop=FALSE]),old))
bad<-old;bad$seed<-999;expect_error(retained_combine(fresh,bad))
bad<-old;bad$method[1]<-'TRACE-EM';expect_error(retained_combine(fresh,bad))
allfresh<-rbind(fresh,data.frame(scenario_index=1,rep=1,seed=101,method='ACVB',f1=.65))
z<-retained_combine(allfresh,old)
stopifnot(sum(z$method=='ACVB')==1,z$f1[z$method=='ACVB']==.65)
cat('PASS: comparator retention, seed pairing, no old TRACE, no duplicates, and full-mode replacement.\n')
