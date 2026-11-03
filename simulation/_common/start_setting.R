bundle_start_setting <- function(project_dir) {
  assign('PROJECT_DIR',normalizePath(project_dir,mustWork=TRUE),envir=.GlobalEnv)
  assign('BUNDLE_DIR',normalizePath(file.path(PROJECT_DIR,'..'),mustWork=TRUE),envir=.GlobalEnv)
  setwd(PROJECT_DIR)
  source(file.path(BUNDLE_DIR,'_common/simulation_50rep_runner.R'),local=.GlobalEnv)
  source(file.path(BUNDLE_DIR,'_common/plot_results.R'),local=.GlobalEnv)
  source(file.path(PROJECT_DIR,'config.R'),local=.GlobalEnv)
  if(SETTING_ID=='02') {
    source(file.path(PROJECT_DIR,'simulation_50rep_runner.R'),local=.GlobalEnv)
    source(file.path(PROJECT_DIR,'plot_results.R'),local=.GlobalEnv)
  }
  source(file.path(BUNDLE_DIR,'_common/run_mode.R'),local=.GlobalEnv)
  mode<-configure_run_mode()
  if(SETTING_ID=='02' && mode=='all') {
    ALGO_ARGS$run_concert <<- TRUE
    SIM_METHODS <<- c('TRACE-EM','Oracle TRACE','ACVB','Trans-Lasso','TRADER','CONCERT')
  }
  if(sim_env_bool('SIM_CONFIG_ONLY',FALSE)) {
    cat('Setting ',SETTING_ID,': ',nrow(scenario_grid),' cells; mode=',mode,
      '; methods=',paste(SIM_METHODS,collapse=', '),'\n',sep='')
    return(invisible(NULL))
  }
  required<-c('MASS','Matrix','glmnet','Rcpp','RcppArmadillo','RcppEigen','quadprog','spikeslab')
  missing<-required[!vapply(required,requireNamespace,logical(1),quietly=TRUE)]
  if(length(missing))stop('Missing R packages: ',paste(missing,collapse=', '),'. Run install_dependencies.R.')
  if(!capabilities('png'))stop('R requires PNG graphics support.')
  if(Sys.info()[['sysname']]!='Darwin'&&!capabilities('cairo'))stop('Headless plotting requires Cairo support.')
  Sys.setenv(BORROW_CORE_PATH=file.path(BUNDLE_DIR,'_common/borrow_functions_core.R'))
  source(file.path(PROJECT_DIR,'functions.R'),local=.GlobalEnv)
  source(file.path(BUNDLE_DIR,'_common/concert_simulation.R'),local=.GlobalEnv)
  if(isTRUE(ALGO_ARGS$run_trans_lasso))translasso_preflight()
  raw<-Sys.getenv('SIM_RESULT_ROOT',file.path(BUNDLE_DIR,paste0('results_T0_corrected_',mode),basename(PROJECT_DIR),'trace_refit'))
  comparison<-Sys.getenv('SIM_COMPARISON_ROOT',file.path(dirname(raw),'comparison'))
  Sys.setenv(SIM_RESULT_ROOT=raw)
  source(file.path(BUNDLE_DIR,'_common/retained_results.R'),local=.GlobalEnv)
  source(file.path(BUNDLE_DIR,'_common/paper_figures.R'),local=.GlobalEnv)
  fresh_inputs <- sim_env_bool('SIM_FRESH_INPUTS',FALSE)
  if(fresh_inputs && mode!='all')stop('SIM_FRESH_INPUTS requires all-method mode.')
  if(!fresh_inputs) {
    retained_install(PROJECT_DIR,BUNDLE_DIR)
    retained_record_fingerprint(.retained_state$directory,comparison)
  }
  sim_run(PROJECT_DIR,BUNDLE_DIR)
  if(!fresh_inputs) {
    retained_finalize(raw,comparison,SETTING_ID)
  } else {
    # A full rerun regenerates source summaries and fits all six methods.
    # No preserved method rows or input-reference checks enter this route.
    perf<-read.csv(file.path(raw,'performance_all_scenarios.csv'),stringsAsFactors=FALSE)
    if(!nrow(perf))stop('No successful full-method replications.')
    perf$result_origin<-'full_rerun'
    dir.create(comparison,recursive=TRUE,showWarnings=FALSE)
    sim_csv(perf,file.path(comparison,'performance_all_scenarios.csv'))
    summary<-sim_summarize(perf,sim_env_int('SIM_NREP',50L))
    sim_summary_tables(summary,comparison,suffix='_all_scenarios')
    sim_make_plots(summary,comparison,SETTING_ID)
    paper_make_figures(summary,comparison,SETTING_ID)
    sim_csv(scenario_grid,file.path(comparison,'scenario_grid.csv'))
  }
  invisible(NULL)
}
