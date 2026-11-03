entry<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
root<-dirname(dirname(normalizePath(entry)))
Sys.setenv(SIM_CONFIG_ONLY='TRUE',SIM_MODE='trace')
folders<-c('01_TRACE_identifiable_mixture_contamination','02_TRACE_support_swap_contamination',
           '03_TRACE_mosaic_contamination','04_TRACE_complementary_sources_fixed_noise')
sizes<-integer(0)
for(folder in folders){
 PROJECT_DIR<-file.path(root,folder)
 source(file.path(PROJECT_DIR,'run_analysis_50rep.R'),local=.GlobalEnv)
 stopifnot(identical(SIM_METHODS,c('TRACE-EM','Oracle TRACE')),
           isTRUE(ALGO_ARGS$run_em_alpha),isTRUE(ALGO_ARGS$run_oracle_alpha),
           !isTRUE(ALGO_ARGS$run_acvb),!isTRUE(ALGO_ARGS$run_trader),
           !isTRUE(ALGO_ARGS$run_trans_lasso),!isTRUE(ALGO_ARGS$run_concert))
 sizes<-c(sizes,nrow(scenario_grid))
 saved<-readRDS(file.path(root,'retained_results',substr(folder,1,2),'original_manifest.rds'))
 switches<-c('run_acvb','run_trans_lasso','run_trader','run_concert','run_old_alpha',
 'run_current_alpha_every5','run_fixed_alpha','run_oracle_trans_simple','run_naive_trans_simple','run_alpha_grid_profile')
 for(key in setdiff(names(saved$algorithm),switches))stopifnot(isTRUE(all.equal(ALGO_ARGS[[key]],saved$algorithm[[key]])))
 stopifnot(isTRUE(all.equal(DGP_BASE_ARGS,saved$dgp)),isTRUE(all.equal(scenario_grid,saved$grid)),BASE_SEED==saved$base_seed)
 Sys.setenv(SIM_MODE='all')
 source(file.path(PROJECT_DIR,'run_analysis_50rep.R'),local=.GlobalEnv)
 expected<-c('TRACE-EM','Oracle TRACE','ACVB','Trans-Lasso','TRADER','CONCERT')
 stopifnot(identical(SIM_METHODS,expected),isTRUE(ALGO_ARGS$run_acvb),
           isTRUE(ALGO_ARGS$run_trader),isTRUE(ALGO_ARGS$run_trans_lasso))
 Sys.setenv(SIM_MODE='trace')
}
stopifnot(identical(sizes,c(48L,20L,48L,21L)))
cat('PASS: default TRACE-only routing and all four original grids.\n')
