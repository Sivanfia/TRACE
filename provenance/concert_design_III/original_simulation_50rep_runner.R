## Setting-02 overrides for the shared resumable 50-rep runner.

sim_cell_name <- function(g) {
  sprintf(
    "scen%03d_TRADER_%s_n%d_p%d_s%d_Ka%02d_h%s",
    g$scenario_index, g$paper_setting, g$n, g$p, g$nzc, g$Ka,
    gsub(".", "p", as.character(g$h), fixed = TRUE)
  )
}

SIM_METHODS <- c("TRACE-EM", "Oracle TRACE", "ACVB", "Trans-Lasso", "TRADER", "CONCERT")
SIM_METRICS <- c(
  "pred_mse", "tp", "fp", "fn", "precision", "recall", "f1", "model_size",
  "support_exact", "support_fdr", "support_fpr", "borrowing_tp", "borrowing_fp",
  "borrowing_fn", "borrowing_precision", "borrowing_recall", "borrowing_f1",
  "borrowing_exact", "transferable_size", "alpha_l1_error",
  "alpha_informative_mass", "alpha_uninformative_mass", "time_sec", "n_iter",
  "converged", "target_init_n_iter", "target_init_size", "target_init_converged",
  "target_init_max_change", "max_param_change_final", "mean_omega_signal",
  "mean_omega_null", "gamma_null_mean", "true_C0_size", "trader_max_psrf"
)

sim_summarize <- function(perf, requested) {
  if (!nrow(perf)) return(data.frame())
  keys <- intersect(c(
    "setting", "scenario_index", "paper_setting_grid", "n_grid", "p_grid",
    "nzc_grid", "K_grid", "n_source_grid", "Ka_grid", "h_grid", "rho_grid",
    "snr_grid", "sigma_fixed_grid", "paired_seed_index", "method"
  ), names(perf))
  groups <- split(seq_len(nrow(perf)), interaction(perf$scenario_index, perf$method, drop = TRUE))
  sim_bind(lapply(groups, function(ii) {
    d <- perf[ii, , drop = FALSE]
    sim_bind(lapply(intersect(SIM_METRICS, names(d)), function(m) {
      x <- as.numeric(d[[m]]); x <- x[is.finite(x)]; n_eff <- length(x)
      z <- d[1, keys, drop = FALSE]
      z$metric <- m; z$n_requested <- requested
      z$n_completed <- length(unique(d$rep)); z$n_eff <- n_eff
      z$mean <- if (n_eff) mean(x) else NA_real_
      z$sd <- if (n_eff > 1L) stats::sd(x) else NA_real_
      z$se <- if (n_eff > 1L) z$sd / sqrt(n_eff) else NA_real_
      z
    }))
  }))
}

sim_run_rep <- function(rr, g, cell_dir, resume) {
  seed <- BASE_SEED + 1000L * as.integer(g$paired_seed_index) + rr - 1L
  rep_dir <- file.path(cell_dir, sprintf("rep%03d_seed%d", rr, seed))
  dir.create(rep_dir, recursive = TRUE, showWarnings = FALSE)
  pfile <- file.path(cell_dir, sprintf("performance_rep%03d.csv", rr))
  sfile <- file.path(cell_dir, sprintf("scenario_rep%03d.csv", rr))
  efile <- file.path(cell_dir, sprintf("error_rep%03d.csv", rr))
  marker <- file.path(rep_dir, "completed.rds")
  cached <- if (resume && file.exists(marker)) tryCatch({
    m <- readRDS(marker)
    stopifnot(identical(m$seed, seed), identical(m$rep, rr), length(m$md5) >= 5L)
    paths <- file.path(cell_dir, names(m$md5))
    stopifnot(all(file.exists(paths)), identical(unname(tools::md5sum(paths)), unname(m$md5)))
    ptab <- read.csv(pfile, stringsAsFactors = FALSE)
    stab <- read.csv(sfile, stringsAsFactors = FALSE)
    sim_validate_rep(ptab, stab, rr, seed)
    list(ok = TRUE, cached = TRUE, performance = ptab, scenario = stab)
  }, error = function(e) NULL) else NULL
  if (!is.null(cached)) return(cached)
  unlink(marker)

  output_level <- sink.number(); message_level <- sink.number(type = "message")
  log <- file(file.path(rep_dir, "run.log"), open = "wt")
  sink(log); sink(log, type = "message")
  on.exit({
    while (sink.number() > output_level) sink()
    if (sink.number(type = "message") != message_level) sink(type = "message")
    close(log)
  }, add = TRUE)

  tryCatch({
    cat("setting=", SETTING_ID, " cell=", g$scenario_index,
        " paper_setting=", g$paper_setting, " Ka=", g$Ka,
        " rep=", rr, " seed=", seed, "\n", sep = "")
    dgp_args <- list(
      seed = seed, n = as.integer(g$n), p = as.integer(g$p),
      nzc = as.integer(g$nzc), rho = as.numeric(g$rho), snr = as.numeric(g$snr),
      common_bad_frac = 0, common_bad_count = 0L,
      paper_setting = as.character(g$paper_setting), Ka = as.integer(g$Ka),
      h = as.numeric(g$h), K = as.integer(g$K),
      n_source_raw = rep(as.integer(g$n_source), as.integer(g$K)),
      sigma_fixed = as.numeric(g$sigma_fixed)
    )
    algo_args <- ALGO_ARGS
    ## The paper caps spike-and-slab at s only in high-dimensional Setting III.
    if (identical(as.character(g$paper_setting), "III")) {
      algo_args$trader_max_var <- as.integer(g$nzc)
    } else {
      algo_args$trader_max_var <- NULL
    }
    ans <- do.call(run_one_elbo, c(DGP_BASE_ARGS, dgp_args, algo_args, list(save_dir = rep_dir)))
    ptab <- sim_add_grid(ans$performance, g, rr, seed)
    stab <- sim_add_grid(ans$scenario, g, rr, seed)
    for (nm in intersect(c(
      "paper_setting", "Ka", "h", "source_summary_method", "noise_mode",
      "sigma_true", "realized_snr", "signal_profile", "beta_min",
      "beta_min_information", "beta_min_margin_ratio", "source_precision_ns_raw",
      "source_precision_ns_eff"
    ), names(stab))) ptab[[nm]] <- stab[[nm]][1]
    ptab$support_exact <- as.integer(ptab$fp == 0 & ptab$fn == 0)
    ptab$support_fdr <- ptab$fp / pmax(1, ptab$model_size)
    ptab$support_fpr <- ptab$fp / (as.integer(g$p) - as.integer(g$nzc))
    ptab$borrowing_exact <- as.integer(ptab$borrowing_fp == 0 & ptab$borrowing_fn == 0)
    sim_validate_rep(ptab, stab, rr, seed)
    sim_csv(ptab, pfile); sim_csv(stab, sfile)
    sim_csv(ptab, file.path(rep_dir, "performance.csv"))
    sim_csv(stab, file.path(rep_dir, "scenario.csv"))
    required <- file.path(rep_dir, c("coordinate_diagnostics.csv", "theory_fit_states.rds"))
    stopifnot(all(file.exists(required)))
    files <- c(pfile, sfile, list.files(rep_dir, pattern = "\\.(csv|rds)$", full.names = TRUE))
    md5 <- tools::md5sum(files); names(md5) <- substring(files, nchar(cell_dir) + 2L)
    sim_rds(list(rep = rr, seed = seed, md5 = md5,
                 completed_at = as.character(Sys.time())), marker)
    unlink(efile)
    list(ok = TRUE, cached = FALSE, performance = ptab, scenario = stab)
  }, error = function(e) {
    msg <- conditionMessage(e); cat("ERROR: ", msg, "\n", sep = "")
    sim_csv(data.frame(setting = SETTING_ID, scenario_index = g$scenario_index,
      rep = rr, seed = seed, error = msg), efile)
    list(ok = FALSE, rep = rr, error = msg)
  })
}

sim_run <- function(project_dir, bundle_dir) {
  nrep <- sim_env_int("SIM_NREP", 50L, maximum = 1000L)
  workers <- min(sim_env_int("SIM_NCORES", 8L), nrep)
  if (nzchar(Sys.getenv("SLURM_CPUS_PER_TASK"))) {
    workers <- min(workers, sim_env_int("SLURM_CPUS_PER_TASK", workers))
  }
  if (.Platform$OS.type == "windows") workers <- 1L
  resume <- sim_env_bool("SIM_RESUME", TRUE)
  result_root <- path.expand(Sys.getenv(
    "SIM_RESULT_ROOT", file.path(project_dir, "results_50rep_TRADER_paper_settings_I_III_with_CONCERT")
  ))
  dir.create(result_root, recursive = TRUE, showWarnings = FALSE)
  result_root <- normalizePath(result_root, mustWork = TRUE)
  stopifnot(
    nrow(scenario_grid) == 20L, !anyDuplicated(scenario_grid$scenario_index),
    all(scenario_grid$Ka %in% 1:10), all(scenario_grid$K == 10L),
    all(scenario_grid$nzc %in% ALLOWED_NZC), all(is.finite(scenario_grid$snr)),
    identical(as.character(scenario_grid$paper_setting), rep(c("I", "III"), each = 10L))
  )
  selected <- scenario_grid$scenario_index
  filter <- Sys.getenv("SIM_SCENARIOS", "")
  if (nzchar(filter)) {
    selected <- suppressWarnings(as.numeric(strsplit(filter, ",", fixed = TRUE)[[1]]))
    if (anyNA(selected) || any(!selected %in% scenario_grid$scenario_index) || anyDuplicated(selected)) {
      stop("SIM_SCENARIOS must be comma-separated distinct scenario_index values (1-20).")
    }
  }
  code <- c(
    list.files(project_dir, pattern = "\\.(R|cpp)$", full.names = TRUE),
    list.files(file.path(bundle_dir, "_common"), pattern = "\\.(R|cpp)$",
      full.names = TRUE, recursive = TRUE),
    Sys.getenv("BORROW_CORE_PATH")
  )
  code <- sort(unique(code[nzchar(code)]))
  hashes <- tools::md5sum(code)
  names(hashes) <- ifelse(startsWith(code, paste0(bundle_dir, "/")),
    substring(code, nchar(bundle_dir) + 2L), paste0("external/", basename(code)))
  packages <- c("MASS", "Matrix", "glmnet", "Rcpp", "RcppArmadillo", "RcppEigen", "quadprog", "spikeslab")
  runtime <- list(
    R = as.character(getRversion()), RNG = RNGkind(),
    packages = setNames(vapply(packages, function(x) as.character(utils::packageVersion(x)), ""), packages)
  )
  config <- list(
    schema = 3L, translasso = translasso_backend_info(), setting = SETTING_ID,
    dgp = DGP_BASE_ARGS, algorithm = ALGO_ARGS, grid = scenario_grid,
    base_seed = BASE_SEED, code_md5 = hashes, runtime = runtime
  )
  manifest <- file.path(result_root, "run_manifest.rds")
  if (file.exists(manifest)) {
    if (!identical(readRDS(manifest), config)) {
      stop("Settings, code or R/package versions changed. Use a NEW SIM_RESULT_ROOT.")
    }
  } else {
    if (length(list.files(result_root, all.files = FALSE))) {
      stop("Result root is nonempty without a run manifest. Choose a new SIM_RESULT_ROOT.")
    }
    sim_rds(config, manifest); dput(config, file = file.path(result_root, "run_config.txt"))
  }
  sim_csv(scenario_grid, file.path(result_root, "scenario_grid.csv"))
  sim_csv(data.frame(nrep = nrep, workers = workers,
    selected_cells = paste(selected, collapse = ","),
    started_at = as.character(Sys.time())), file.path(result_root, "last_invocation.csv"))
  writeLines(capture.output(sessionInfo()), file.path(result_root, "sessionInfo.txt"))
  cat(sprintf("[setting %s] %d selected cells x %d reps; %d workers\n",
              SETTING_ID, length(selected), nrep, workers))
  cat("[results] ", result_root, "\n", sep = "")

  all_perf <- all_scen <- statuses <- list()
  for (idx in selected) {
    g <- scenario_grid[scenario_grid$scenario_index == idx, , drop = FALSE]
    cell <- file.path(result_root, sim_cell_name(g)); dir.create(cell, showWarnings = FALSE)
    cat(sprintf("[cell %d/20] paper=%s n=%d p=%d s=%d Ka=%d h=%g SNR=%.5g\n",
      idx, g$paper_setting, g$n, g$p, g$nzc, g$Ka, g$h, g$snr))
    fun <- function(rr) sim_run_rep(as.integer(rr), g, cell, resume)
    results <- if (workers > 1L) parallel::mclapply(
      seq_len(nrep), fun, mc.cores = workers, mc.preschedule = FALSE, mc.set.seed = FALSE
    ) else lapply(seq_len(nrep), fun)
    good <- vapply(results, function(x) is.list(x) && isTRUE(x$ok), logical(1))
    cached_n <- sum(vapply(results[good], function(x) isTRUE(x$cached), logical(1)))
    ptab <- sim_bind(lapply(results[good], `[[`, "performance"))
    stab <- sim_bind(lapply(results[good], `[[`, "scenario"))
    if (nrow(ptab)) ptab <- ptab[order(ptab$rep, match(ptab$method, SIM_METHODS)), , drop = FALSE]
    if (nrow(stab)) stab <- stab[order(stab$rep), , drop = FALSE]
    sim_csv(ptab, file.path(cell, "performance_all_reps.csv"))
    sim_csv(stab, file.path(cell, "scenarios_all_reps.csv"))
    sim_summary_tables(sim_summarize(ptab, nrep), cell)
    all_perf[[as.character(idx)]] <- ptab; all_scen[[as.character(idx)]] <- stab
    statuses[[as.character(idx)]] <- data.frame(
      setting = SETTING_ID, scenario_index = idx, n_requested = nrep,
      n_completed = sum(good), n_resumed = cached_n, n_failed = sum(!good),
      failed_reps = paste(which(!good), collapse = ";"), complete = all(good)
    )
    sim_csv(sim_bind(statuses), file.path(result_root, "run_status.csv"))
    sim_csv(sim_bind(all_perf), file.path(result_root, "performance_all_scenarios.csv"))
    sim_csv(sim_bind(all_scen), file.path(result_root, "scenarios_all_scenarios.csv"))
    cat(sprintf("[cell %d] completed=%d/%d; resumed=%d; failed=%d\n",
                idx, sum(good), nrep, cached_n, sum(!good)))
  }
  summary <- sim_summarize(sim_bind(all_perf), nrep)
  sim_summary_tables(summary, result_root, "_all_scenarios")
  sim_make_plots(summary, result_root, SETTING_ID)
  if (any(!sim_bind(statuses)$complete)) {
    stop("Incomplete repetitions. Inspect run_status.csv and rep logs, then rerun to resume.")
  }
  cat("[complete] Tables and individual PDF/PNG plots: ", result_root, "\n", sep = "")
  invisible(summary)
}
