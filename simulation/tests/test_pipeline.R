entry<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
root<-dirname(dirname(normalizePath(entry)))
options(device=function(...)grDevices::pdf(file=NULL))
arg<-commandArgs(TRUE);id<-if(length(arg))arg[1] else '01'
design<-if(length(arg)>1)arg[2] else 'I'
folders<-c('01'='01_TRACE_identifiable_mixture_contamination','02'='02_TRACE_support_swap_contamination',
           '03'='03_TRACE_mosaic_contamination','04'='04_TRACE_complementary_sources_fixed_noise')
Sys.setenv(SIM_CONFIG_ONLY='TRUE',SIM_MODE='trace',SIM_NREP='1',SIM_NCORES='1',
           SIM_SOURCE_SS_BURNIN='10',SIM_SOURCE_SS_SAMPLES='10')
PROJECT_DIR<-file.path(root,folders[[id]])
source(file.path(PROJECT_DIR,'run_analysis_50rep.R'),local=.GlobalEnv)
Sys.setenv(BORROW_CORE_PATH=file.path(root,'_common/borrow_functions_core.R'))
source(file.path(PROJECT_DIR,'functions.R'),local=.GlobalEnv)
if(id!='02')source(file.path(root,'_common/concert_simulation.R'),local=.GlobalEnv)
source(file.path(root,'_common/retained_results.R'))
source(file.path(root,'_common/paper_figures.R'))
output<-tempfile(paste0('TRACE_pipeline_',id,'_'));dir.create(output)
fixture<-file.path(output,'baseline');dir.create(fixture)
g<-if(id=='02')scenario_grid[scenario_grid$paper_setting==design & scenario_grid$Ka==1,,drop=FALSE] else scenario_grid[1,,drop=FALSE]
if(id!='02') {
 DGP_BASE_ARGS$n<-30L;DGP_BASE_ARGS$p<-40L;DGP_BASE_ARGS$n_test<-40L
 g$nzc<-4L;g$common_bad_frac<-.25;g$common_bad_count<-1L
} else DGP_BASE_ARGS$source_summary_method<-'oracle'
ALGO_ARGS$max_iter_borrow<-1L;ALGO_ARGS$verbose<-FALSE;ALGO_ARGS$warn_elbo_drop<-FALSE
seed<-BASE_SEED+1000L*as.integer(g$paired_seed_index)
call_args<-list(seed=seed,nzc=as.integer(g$nzc),rho=g$rho,snr=g$snr,common_bad_frac=g$common_bad_frac)
if(id=='02')call_args<-c(call_args,list(n=g$n,p=g$p,K=g$K,n_source_raw=rep(g$n_source,g$K),
 paper_setting=as.character(g$paper_setting),Ka=g$Ka,h=g$h,sigma_fixed=g$sigma_fixed))
if(id=='04')call_args<-c(call_args,list(sigma_fixed=g$sigma_fixed,shared_spur_count=g$shared_spur_count))
save<-file.path(output,'reference');dir.create(save)
ref<-do.call(run_one_elbo,c(DGP_BASE_ARGS,call_args,ALGO_ARGS,list(save_dir=save)))
fit<-ref$fits$borrow_cvb_em_alpha
tr<-fit$theory_init_trace;z<-tr[tr$step=='alpha_preliminary',,drop=FALSE]
cols<-grep('^alpha_[0-9]+$',names(z),value=TRUE);cols<-cols[order(as.integer(sub('alpha_','',cols)))]
targets<-list();targets[[as.character(seed)]]<-list(omega=as.numeric(fit$target_initialization$omega),
 beta=as.numeric(fit$preliminary_alpha_beta_target),sigma_beta_sq=fit$sigma_beta_sq)
alphas<-list();alphas[[sprintf('%d:%d',g$scenario_index,1)]]<-as.numeric(z[1,cols])
saveRDS(list(targets=targets,initial_alpha=alphas),file.path(fixture,'initialization_reference.rds'))
baseline<-data.frame(setting=id,scenario_index=g$scenario_index,rep=1,seed=seed,method='ACVB',
 f1=.75,pred_mse=1,precision=.75,recall=.75,model_size=4,tp=3,fp=1,fn=1)
baseline<-sim_add_grid(baseline,g,1L,seed)
write.csv(baseline,file.path(fixture,'performance_baseline.csv'),row.names=FALSE)
if(id!='02') {
 dir.create(file.path(fixture,'source_cache'))
 cache<-list(B=ref$dat$B,standard_errors=NULL,n_source=ref$dat$n_vec,
  method='spikeslab_gnet_original_scale',diagnostics=ref$dat$source_summary_diagnostics)
 saveRDS(cache,file.path(fixture,'source_cache',sprintf('scen%03d_rep001.rds',g$scenario_index)))
}
original_formals<-formals(make_dat)
retained_install(PROJECT_DIR,root,fixture)
stopifnot(identical(formals(make_dat),original_formals))
ALGO_ARGS$max_iter_borrow<-8L
cell<-file.path(output,'cell');dir.create(cell)
ans<-sim_run_rep(1L,g,cell,resume=FALSE)
if(!isTRUE(ans$ok))stop('Pipeline failure: ',ans$error,'; logs: ',cell)
stopifnot(setequal(ans$performance$method,c('TRACE-EM','Oracle TRACE')),
 ans$scenario$retained_target_reference_error<1e-7,ans$scenario$retained_initial_alpha_error<1e-6)
cached<-sim_run_rep(1L,g,cell,resume=TRUE);stopifnot(isTRUE(cached$cached))
raw<-file.path(output,'trace_refit');dir.create(raw)
write.csv(ans$performance,file.path(raw,'performance_all_scenarios.csv'),row.names=FALSE)
compare<-file.path(output,'comparison')
combined<-retained_finalize(raw,compare,id)
stopifnot(nrow(combined)==3L,combined$f1[combined$method=='ACVB']==.75,
 file.exists(file.path(compare,'paper_figures',paste0('setting',id,'_main_f1.pdf'))),
 file.exists(file.path(compare,'paper_figures',paste0('setting',id,'_main_f1.png'))))
bad<-.retained_state$reference$targets[[as.character(seed)]]$beta
.retained_state$reference$targets[[as.character(seed)]]$beta<-bad+1
failed<-sim_run_rep(1L,g,file.path(output,'bad_cell'),resume=FALSE)
stopifnot(!isTRUE(failed$ok),grepl('Input alignment failed',failed$error))
cat('PASS: real TRACE pipeline, retained comparator merge, input guard, resume and figures; setting ',id,' ',design,'\n',sep='')
cat('Validation-only temporary output: ',output,'\n',sep='')
