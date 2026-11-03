retained_pair_key <- function(x) paste(x$scenario_index,x$rep,x$seed,sep=':')
retained_row_key <- function(x) paste(retained_pair_key(x),x$method,sep=':')

retained_combine <- function(fresh,baseline) {
  trace_methods<-c('TRACE-EM','Oracle TRACE')
  if (!nrow(fresh)) return(fresh)
  if(any(baseline$method %in% trace_methods))stop('Retained baselines must not contain old TRACE rows.')
  if(anyDuplicated(retained_row_key(fresh))||anyDuplicated(retained_row_key(baseline)))stop('Duplicate result keys.')
  pairs<-unique(retained_pair_key(fresh))
  for(key in pairs) {
    if(!setequal(fresh$method[retained_pair_key(fresh)==key & fresh$method%in%trace_methods],trace_methods))
      stop('A completed replicate must contain both new TRACE methods: ',key)
    if(!key%in%retained_pair_key(baseline))stop('No matching retained baseline for cell/rep/seed: ',key)
  }
  old<-baseline[retained_pair_key(baseline)%in%pairs,,drop=FALSE]
  old<-old[!retained_row_key(old)%in%retained_row_key(fresh),,drop=FALSE]
  fresh$result_origin<-'recomputed_T0'
  old$result_origin<-'retained_baseline'
  z<-sim_bind(list(fresh,old))
  methods<-c('TRACE-EM','Oracle TRACE','ACVB','Trans-Lasso','TRADER','CONCERT')
  z<-z[order(z$scenario_index,z$rep,match(z$method,methods)),,drop=FALSE]
  rownames(z)<-NULL
  z
}

retained_check_vector <- function(got,want,label,tolerance=1e-7) {
  if(length(got)!=length(want)||!length(got)||any(!is.finite(got))||
     max(abs(as.numeric(got)-as.numeric(want)))>tolerance) {
    stop('Input alignment failed for ',label,
      '. The regenerated inputs or initialization differ from the retained experiment. ',
      'Use the original server R/package environment. No retained results will be mixed into this replicate.')
  }
  max(abs(as.numeric(got)-as.numeric(want)))
}

retained_install <- function(project_dir,bundle_dir,baseline_dir=NULL) {
  id<-as.character(SETTING_ID)
  if(is.null(baseline_dir))baseline_dir<-file.path(bundle_dir,'retained_results',id)
  state<-new.env(parent=emptyenv())
  state$id<-id;state$directory<-baseline_dir
  state$baseline<-read.csv(file.path(baseline_dir,'performance_baseline.csv'),stringsAsFactors=FALSE)
  state$reference<-readRDS(file.path(baseline_dir,'initialization_reference.rds'))
  state$context<-NULL
  state$original_rep<-sim_run_rep
  state$original_make_dat<-make_dat
  state$original_run_one<-run_one_elbo
  assign('.retained_state',state,envir=.GlobalEnv)

  if(exists('fit_shared_source_summaries',mode='function',inherits=TRUE)) {
    state$original_summary<-get('fit_shared_source_summaries',mode='function')
    assign('fit_shared_source_summaries',function(X_source,y_source,n_source_raw,seed) {
      path<-Sys.getenv('TRACE_SOURCE_CACHE_FILE','')
      if(!nzchar(path))return(.retained_state$original_summary(X_source,y_source,n_source_raw,seed))
      x<-readRDS(path)
      K<-length(X_source);p<-ncol(X_source[[1]])
      if(!identical(dim(x$B),c(as.integer(p),as.integer(K)))||any(!is.finite(x$B))||
         !all(x$n_source==vapply(X_source,nrow,1L))||
         !identical(x$method,'spikeslab_gnet_original_scale')||
         !all(x$diagnostics$preprocessing_seed==as.integer(seed)+seq_len(K)))
        stop('Saved source summary does not match this replication.')
      x$diagnostics$reused_saved_summary<-TRUE
      list(B=x$B,standard_errors=x$standard_errors,counts=x$n_source,
           diagnostics=x$diagnostics,fit_time_sec=0)
    },envir=.GlobalEnv)
  }

  assign('sim_run_rep',function(rr,g,cell_dir,resume) {
    st<-.retained_state
    previous_device<-getOption('device')
    previous_devices<-grDevices::dev.list()
    options(device=function(...)grDevices::pdf(file=NULL))
    on.exit({
      extra<-setdiff(grDevices::dev.list(),previous_devices)
      for(device in extra)try(grDevices::dev.off(device),silent=TRUE)
      options(device=previous_device)
    },add=TRUE)
    seed<-BASE_SEED+1000L*as.integer(g$paired_seed_index)+rr-1L
    key<-paste(as.integer(g$scenario_index),as.integer(rr),as.integer(seed),sep=':')
    if(!key%in%retained_pair_key(st$baseline))stop('Missing retained replicate: ',key)
    st$context<-list(rep=rr,scenario=as.integer(g$scenario_index),seed=as.integer(seed))
    cache<-file.path(st$directory,'source_cache',sprintf('scen%03d_rep%03d.rds',g$scenario_index,rr))
    if(st$id!='02'&&!file.exists(cache))stop('Missing packaged source summary: ',cache)
    previous<-Sys.getenv('TRACE_SOURCE_CACHE_FILE',unset=NA_character_)
    Sys.setenv(TRACE_SOURCE_CACHE_FILE=if(file.exists(cache))cache else '')
    on.exit({st$context<-NULL;if(is.na(previous))Sys.unsetenv('TRACE_SOURCE_CACHE_FILE') else Sys.setenv(TRACE_SOURCE_CACHE_FILE=previous)},add=TRUE)
    st$original_rep(rr,g,cell_dir,resume)
  },envir=.GlobalEnv)

  make_wrapper<-function(...) {
    call_args<-as.list(match.call())[-1L]
    dat<-do.call(.retained_state$original_make_dat,call_args,envir=parent.frame())
    ctx<-.retained_state$context
    if(is.null(ctx))stop('Missing retained-input context.')
    ref<-.retained_state$reference$targets[[as.character(ctx$seed)]]
    if(is.null(ref))stop('Missing target-input reference.')
    preliminary<-.compute_alpha_preliminary_qp_terms(ref$omega,dat,
      sigma_beta_sq=ref$sigma_beta_sq)$beta_target
    dat$retained_target_reference_error<-retained_check_vector(preliminary,ref$beta,'target coefficients')
    dat$retained_source_summary_reused<-nzchar(Sys.getenv('TRACE_SOURCE_CACHE_FILE',''))
    dat
  }
  formals(make_wrapper)<-formals(state$original_make_dat)
  assign('make_dat',make_wrapper,envir=.GlobalEnv)

  assign('run_one_elbo',function(...) {
    ans<-.retained_state$original_run_one(...)
    ctx<-.retained_state$context
    ref<-.retained_state$reference$targets[[as.character(ctx$seed)]]
    fit<-ans$fits$borrow_cvb_em_alpha
    retained_check_vector(fit$target_initialization$omega,ref$omega,'target inclusion initialization')
    retained_check_vector(fit$preliminary_alpha_beta_target,ref$beta,'auxiliary target coefficients')
    tr<-fit$theory_init_trace
    row<-tr[tr$step=='alpha_preliminary',,drop=FALSE]
    cols<-grep('^alpha_[0-9]+$',names(row),value=TRUE)
    cols<-cols[order(as.integer(sub('alpha_','',cols)))]
    alpha_ref<-.retained_state$reference$initial_alpha[[sprintf('%d:%d',ctx$scenario,ctx$rep)]]
    alpha_error<-retained_check_vector(as.numeric(row[1,cols]),alpha_ref,'source-weight initialization',tolerance=1e-6)
    ans$scenario$retained_target_reference_error<-ans$dat$retained_target_reference_error
    ans$scenario$retained_initial_alpha_error<-alpha_error
    ans$scenario$retained_source_summary_reused<-ans$dat$retained_source_summary_reused
    ans
  },envir=.GlobalEnv)
  invisible(state)
}

retained_record_fingerprint <- function(baseline_dir,comparison_dir) {
  paths<-list.files(baseline_dir,recursive=TRUE,full.names=TRUE)
  hashes<-tools::md5sum(paths)
  names(hashes)<-substring(paths,nchar(baseline_dir)+2L)
  dir.create(comparison_dir,recursive=TRUE,showWarnings=FALSE)
  file<-file.path(comparison_dir,'retained_input_fingerprint.rds')
  if(file.exists(file)) {
    if(!identical(readRDS(file),hashes))stop('Retained baseline inputs changed. Use a new result directory.')
  } else saveRDS(hashes,file)
  invisible(hashes)
}

retained_finalize <- function(raw_dir,comparison_dir,id) {
  path<-file.path(raw_dir,'performance_all_scenarios.csv')
  if(!file.exists(path))stop('No new TRACE result table was produced.')
  fresh<-read.csv(path,stringsAsFactors=FALSE)
  if(!nrow(fresh))stop('No completed new TRACE results. Retained baselines alone are not plotted.')
  combined<-retained_combine(fresh,.retained_state$baseline)
  dir.create(comparison_dir,recursive=TRUE,showWarnings=FALSE)
  sim_csv(combined,file.path(comparison_dir,'performance_all_scenarios.csv'))
  nrep<-sim_env_int('SIM_NREP',50L,maximum=50L)
  summary<-sim_summarize(combined,nrep)
  sim_summary_tables(summary,comparison_dir,suffix='_all_scenarios')
  sim_make_plots(summary,comparison_dir,id)
  paper_make_figures(summary,comparison_dir,id)
  sim_csv(scenario_grid,file.path(comparison_dir,'scenario_grid.csv'))
  records<-as.data.frame(table(combined$method,combined$result_origin),stringsAsFactors=FALSE)
  names(records)<-c('method','result_origin','rows')
  sim_csv(records[records$rows>0,,drop=FALSE],file.path(comparison_dir,'result_provenance.csv'))
  cat('[combined results and figures] ',comparison_dir,'\n',sep='')
  invisible(combined)
}
