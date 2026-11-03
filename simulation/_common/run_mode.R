configure_run_mode <- function() {
  mode <- tolower(Sys.getenv('SIM_MODE','trace'))
  if (!mode %in% c('trace','all')) stop('SIM_MODE must be trace or all.')
  if (mode=='trace') {
    for (key in c('run_acvb','run_trans_lasso','run_trader','run_concert',
                  'run_old_alpha','run_current_alpha_every5','run_fixed_alpha',
                  'run_oracle_trans_simple','run_naive_trans_simple','run_alpha_grid_profile')) {
      ALGO_ARGS[[key]] <<- FALSE
    }
    ALGO_ARGS$run_em_alpha <<- TRUE
    ALGO_ARGS$run_oracle_alpha <<- TRUE
    SIM_METHODS <<- c('TRACE-EM','Oracle TRACE')
  }
  mode
}
