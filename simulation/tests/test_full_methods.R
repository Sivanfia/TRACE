# Small six-method fitting check. These temporary parameters are not used
# for the publication results or full rerun configuration.
entry<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
if(!file.exists(entry))entry<-gsub('~+~',' ',entry,fixed=TRUE)
root<-dirname(dirname(normalizePath(entry)))
options(device=function(...)grDevices::pdf(file=NULL))
Sys.setenv(SIM_CONFIG_ONLY='TRUE',SIM_MODE='all',SIM_SOURCE_SS_BURNIN='10',SIM_SOURCE_SS_SAMPLES='10')
project<-file.path(root,'01_TRACE_identifiable_mixture_contamination')
source(file.path(project,'run_analysis_50rep.R'),local=.GlobalEnv)
Sys.setenv(BORROW_CORE_PATH=file.path(root,'_common/borrow_functions_core.R'))
source(file.path(project,'functions.R'),local=.GlobalEnv)
source(file.path(root,'_common/concert_simulation.R'),local=.GlobalEnv)
DGP_BASE_ARGS$n<-30L;DGP_BASE_ARGS$p<-40L;DGP_BASE_ARGS$n_test<-40L
ALGO_ARGS$max_iter_borrow<-8L;ALGO_ARGS$verbose<-FALSE
ALGO_ARGS$trader_n_mcmc<-40L;ALGO_ARGS$trader_n_burnin<-20L
ALGO_ARGS$trader_verbose<-FALSE;ALGO_ARGS$concert_max_iter<-8L
out<-tempfile('TRACE_six_method_check_');dir.create(out)
ans<-do.call(run_one_elbo,c(DGP_BASE_ARGS,
 list(seed=5481L,nzc=4L,rho=.3,snr=4,common_bad_frac=.25,common_bad_count=1L),
 ALGO_ARGS,list(save_dir=out)))
stopifnot(nrow(ans$performance)==6L,
 setequal(ans$performance$method,c('TRACE-EM','Oracle TRACE','ACVB','Trans-Lasso','TRADER','CONCERT')),
 all(is.finite(ans$performance$f1)),all(is.finite(ans$performance$pred_mse)))
cat('PASS: all six simulation methods fit the same temporary simulation data.\n')
