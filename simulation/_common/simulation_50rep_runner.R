# Replicate-level parallel execution. Scenario cells remain serial.
sim_env_int <- function(name, default, maximum = .Machine$integer.max) {
  raw <- Sys.getenv(name, as.character(default)); x <- suppressWarnings(as.numeric(raw))
  if (length(x) != 1L || !is.finite(x) || x < 1 || x != floor(x) || x > maximum)
    stop(name, " must be a positive integer <= ", maximum)
  as.integer(x)
}
sim_env_bool <- function(name, default) {
  x <- tolower(Sys.getenv(name, if (default) "TRUE" else "FALSE"))
  if (!x %in% c("true", "false", "1", "0", "yes", "no", "t", "f")) stop("Invalid boolean: ", name)
  x %in% c("true", "1", "yes", "t")
}
sim_atomic <- function(path, writer) {
  tmp <- tempfile(".writing_", tmpdir = dirname(path))
  on.exit(unlink(tmp), add = TRUE)
  writer(tmp)
  if (!file.rename(tmp, path)) stop("Cannot commit file: ", path)
}
sim_csv <- function(x, path) sim_atomic(path, function(tmp) write.csv(x, tmp, row.names = FALSE, na = "NA"))
sim_rds <- function(x, path) sim_atomic(path, function(tmp) saveRDS(x, tmp))
sim_read_saved_table <- function(path) {
  header <- readLines(path, n = 2L, warn = FALSE)
  # write.csv(data.frame()) emits only an empty quoted header after all reps fail.
  if (!length(header) || identical(trimws(header), '""')) return(data.frame())
  read.csv(path, stringsAsFactors = FALSE)
}
sim_bind <- function(xs) {
  xs <- Filter(function(x) is.data.frame(x) && nrow(x) > 0, xs)
  if (!length(xs)) return(data.frame())
  cols <- unique(unlist(lapply(xs, names)))
  xs <- lapply(xs, function(x) { for (n in setdiff(cols, names(x))) x[[n]] <- NA; x[, cols, drop = FALSE] })
  out <- do.call(rbind, xs); rownames(out) <- NULL; out
}
sim_add_grid <- function(tab, g, rr, seed) {
  for (nm in names(g)) tab[[if (nm %in% c("scenario_index", "paired_seed_index") || grepl("_grid$", nm)) nm else paste0(nm, "_grid")]] <- g[[nm]][1]
  tab$setting <- SETTING_ID; tab$rep <- rr; tab$seed <- seed; tab
}
sim_cell_name <- function(g) {
  suffix <- if ("shared_spur_count" %in% names(g)) paste0("_spur", g$shared_spur_count) else ""
  paste0(sprintf("scen%03d_nzc%d_rho%s_snr%s_q%02d", g$scenario_index, g$nzc,
    gsub(".", "p", as.character(g$rho), fixed = TRUE),
    gsub(".", "p", as.character(g$snr), fixed = TRUE), round(100 * g$common_bad_frac)), suffix)
}
SIM_METHODS <- c("TRACE-EM", "Oracle TRACE", "ACVB", "Trans-Lasso", "TRADER", "CONCERT")
SIM_METRICS <- c("pred_mse", "tp", "fp", "fn", "precision", "recall", "f1", "model_size",
  "support_exact", "support_fdr", "support_fpr", "borrowing_tp", "borrowing_fp", "borrowing_fn",
  "borrowing_precision", "borrowing_recall", "borrowing_f1", "borrowing_exact", "transferable_size",
  "alpha_l1_error", "time_sec", "n_iter", "converged", "target_init_n_iter", "target_init_size",
  "target_init_converged", "target_init_max_change", "max_param_change_final",
  "mean_omega_signal", "mean_omega_null", "mean_gamma_common_bad", "mean_gamma_transferable_signal",
  "mean_gamma_spurious_null", "gamma_null_mean", "true_C0_size", "trader_max_psrf")
sim_summarize <- function(perf, requested) {
  if (!nrow(perf)) return(data.frame())
  keys <- intersect(c("setting", "scenario_index", "nzc_grid", "rho_grid", "snr_grid",
    "common_bad_frac_grid", "common_bad_count_grid", "shared_spur_count_grid", "sigma_fixed_grid",
    "alpha_identification_boundary_grid", "paired_seed_index", "method"), names(perf))
  groups <- split(seq_len(nrow(perf)), interaction(perf$scenario_index, perf$method, drop = TRUE))
  sim_bind(lapply(groups, function(ii) {
    d <- perf[ii, , drop = FALSE]
    sim_bind(lapply(intersect(SIM_METRICS, names(d)), function(m) {
      x <- as.numeric(d[[m]]); x <- x[is.finite(x)]; n <- length(x)
      z <- d[1, keys, drop = FALSE]
      z$metric <- m; z$n_requested <- requested; z$n_completed <- length(unique(d$rep)); z$n_eff <- n
      z$mean <- if (n) mean(x) else NA_real_
      z$sd <- if (n > 1L) sd(x) else NA_real_
      z$se <- if (n > 1L) z$sd / sqrt(n) else NA_real_
      z
    }))
  }))
}
sim_summary_tables <- function(summary, root, suffix = "") {
  sim_csv(summary, file.path(root, paste0("summary_long", suffix, ".csv")))
  if (!nrow(summary)) return(invisible(NULL))
  ids <- setdiff(names(summary), c("metric", "mean", "sd", "se", "n_eff"))
  for (stat in c("mean", "sd", "se", "n_eff")) {
    z <- reshape(summary[, c(ids, "metric", stat)], idvar = ids, timevar = "metric", direction = "wide")
    names(z) <- sub(paste0("^", stat, "\\."), "", names(z))
    sim_csv(z, file.path(root, paste0("summary_", stat, suffix, ".csv")))
  }
  # One row per cell + metric, methods in parallel columns; counts retained separately.
  ids2 <- setdiff(ids, c("method", "n_completed"))
  for (stat in c("mean", "sd", "se", "n_eff")) {
    z <- reshape(summary[, c(ids2, "metric", "method", stat)], idvar = c(ids2, "metric"), timevar = "method", direction = "wide")
    names(z) <- sub(paste0("^", stat, "\\."), "", names(z))
    sim_csv(z, file.path(root, paste0("wide_", stat, suffix, ".csv")))
  }
}
sim_validate_rep <- function(perf, scen, rr, seed) {
  stopifnot(is.data.frame(perf), nrow(perf) == length(SIM_METHODS),
    setequal(perf$method, SIM_METHODS), !anyDuplicated(perf$method),
    nrow(scen) == 1L, all(perf$rep == rr), all(perf$seed == seed),
    scen$rep == rr, scen$seed == seed,
    all(is.finite(perf$pred_mse)), all(is.finite(perf$f1)))
  translasso_validate_performance(perf)
  if (any(grepl("^(mse_beta|sse_beta|l2_beta|signal_mse|total_mse|null_mse)", names(perf))))
    stop("Coefficient error metrics must not be emitted.")
  invisible(TRUE)
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
    p <- read.csv(pfile, stringsAsFactors = FALSE); s <- read.csv(sfile, stringsAsFactors = FALSE)
    sim_validate_rep(p, s, rr, seed)
    list(ok = TRUE, cached = TRUE, performance = p, scenario = s)
  }, error = function(e) NULL) else NULL
  if (!is.null(cached)) return(cached)
  unlink(marker)
  output_level <- sink.number(); message_level <- sink.number(type = "message")
  log <- file(file.path(rep_dir, "run.log"), open = "wt")
  sink(log); sink(log, type = "message")
  on.exit({ while (sink.number() > output_level) sink();
    if (sink.number(type = "message") != message_level) sink(type = "message"); close(log) }, add = TRUE)
  tryCatch({
    cat("setting=", SETTING_ID, " cell=", g$scenario_index, " rep=", rr, " seed=", seed, "\n")
    set.seed(seed)
    args <- list(seed = seed, nzc = as.integer(g$nzc), rho = as.numeric(g$rho),
      snr = as.numeric(g$snr), common_bad_frac = as.numeric(g$common_bad_frac))
    if ("shared_spur_count" %in% names(g)) args$shared_spur_count <- as.integer(g$shared_spur_count)
    if ("sigma_fixed" %in% names(g)) args$sigma_fixed <- as.numeric(g$sigma_fixed)
    ans <- do.call(run_one_elbo, c(DGP_BASE_ARGS, args, ALGO_ARGS, list(save_dir = rep_dir)))
    # Record the shared input protocol without serializing the large raw datasets.
    actual_counts <- vapply(ans$dat$X_source, nrow, 1L)
    stopifnot(all(actual_counts == DGP_BASE_ARGS$n_source_raw),
      all(ans$dat$n_vec == actual_counts), ans$dat$ns == mean(actual_counts),
      identical(ans$dat$source_summary_method, "spikeslab_gnet_original_scale"),
      identical(ans$dat$source_summary_coefficient_scale, "original_predictors"),
      isTRUE(ans$dat$source_summary_from_raw_samples))
    protocol <- list(source_summary_method = ans$dat$source_summary_method,
      source_summary_coefficient_scale = ans$dat$source_summary_coefficient_scale,
      source_summary_from_raw_samples = TRUE,
      source_sample_size_min = min(actual_counts), source_sample_size_max = max(actual_counts),
      source_summary_fit_time_sec = ans$dat$source_summary_fit_time_sec)
    for (nm in names(protocol)) {
      ans$scenario[[nm]] <- protocol[[nm]]
      ans$performance[[nm]] <- protocol[[nm]]
    }
    sim_csv(ans$dat$source_summary_diagnostics, file.path(rep_dir, "source_summary_diagnostics.csv"))
    sim_rds(list(B = ans$dat$B, standard_errors = ans$dat$source_summary_se,
      n_source = actual_counts, method = ans$dat$source_summary_method,
      from_raw_samples = TRUE), file.path(rep_dir, "source_summaries.rds"))
    p <- sim_add_grid(ans$performance, g, rr, seed)
    s <- sim_add_grid(ans$scenario, g, rr, seed)
    for (nm in intersect(c("noise_mode", "sigma_true", "realized_snr", "signal_profile",
      "theory_signal_magnitude", "beta_min", "beta_max", "beta_min_information",
      "beta_min_log_benchmark", "beta_min_margin_ratio", "common_bad_gap_population",
      "source_weight_identification_eigenvalue", "h6_nonempty_C0", "h7_all_subsets_verified"), names(s))) p[[nm]] <- s[[nm]][1]
    p$support_exact <- as.integer(p$fp == 0 & p$fn == 0)
    p$support_fdr <- p$fp / pmax(1, p$model_size)
    p$support_fpr <- p$fp / (DGP_BASE_ARGS$p - g$nzc)
    p$borrowing_exact <- as.integer(p$borrowing_fp == 0 & p$borrowing_fn == 0)
    sim_validate_rep(p, s, rr, seed)
    sim_csv(p, pfile); sim_csv(s, sfile)
    sim_csv(p, file.path(rep_dir, "performance.csv")); sim_csv(s, file.path(rep_dir, "scenario.csv"))
    required <- file.path(rep_dir, c("coordinate_diagnostics.csv", "theory_fit_states.rds",
      "source_summary_diagnostics.csv", "source_summaries.rds"))
    stopifnot(all(file.exists(required)))
    files <- c(pfile, sfile, list.files(rep_dir, pattern = "\\.(csv|rds)$", full.names = TRUE))
    md5 <- tools::md5sum(files); names(md5) <- substring(files, nchar(cell_dir) + 2L)
    sim_rds(list(rep = rr, seed = seed, md5 = md5, completed_at = as.character(Sys.time())), marker)
    unlink(efile)
    list(ok = TRUE, cached = FALSE, performance = p, scenario = s)
  }, error = function(e) {
    msg <- conditionMessage(e); cat("ERROR: ", msg, "\n")
    sim_csv(data.frame(setting = SETTING_ID, scenario_index = g$scenario_index, rep = rr,
      seed = seed, error = msg), efile)
    list(ok = FALSE, rep = rr, error = msg)
  })
}
sim_run <- function(project_dir, bundle_dir) {
  nrep <- sim_env_int("SIM_NREP", 50L, maximum = 1000L)
  workers <- min(sim_env_int("SIM_NCORES", 4L), nrep)
  if (nzchar(Sys.getenv("SLURM_CPUS_PER_TASK"))) workers <- min(workers, sim_env_int("SLURM_CPUS_PER_TASK", workers))
  if (.Platform$OS.type == "windows") workers <- 1L
  resume <- sim_env_bool("SIM_RESUME", TRUE)
  result_root <- path.expand(Sys.getenv("SIM_RESULT_ROOT", file.path(bundle_dir,
    paste0("results_same_source_", Sys.getenv("SIM_SOURCE_N", "500"), "_spikeslab"), basename(project_dir))))
  dir.create(result_root, recursive = TRUE, showWarnings = FALSE)
  result_root <- normalizePath(result_root, mustWork = TRUE)
  stopifnot(!anyDuplicated(scenario_grid$scenario_index), all(scenario_grid$nzc %in% ALLOWED_NZC),
    DGP_BASE_ARGS$n == 300, DGP_BASE_ARGS$p == 600, all(is.finite(scenario_grid$snr)),
    all(scenario_grid$common_bad_frac >= 0 & scenario_grid$common_bad_frac <= 1))
  selected <- scenario_grid$scenario_index
  filter <- Sys.getenv("SIM_SCENARIOS", "")
  if (nzchar(filter)) {
    selected <- suppressWarnings(as.numeric(strsplit(filter, ",", fixed = TRUE)[[1]]))
    if (anyNA(selected) || any(!selected %in% scenario_grid$scenario_index) || anyDuplicated(selected))
      stop("SIM_SCENARIOS must be comma-separated distinct scenario_index values.")
  }
  code <- c(list.files(project_dir, pattern = "\\.(R|cpp)$", full.names = TRUE),
    list.files(file.path(bundle_dir, "_common"), pattern = "\\.(R|cpp)$",
      full.names = TRUE, recursive = TRUE))
  code <- sort(unique(c(code, Sys.getenv("BORROW_CORE_PATH"))))
  hashes <- tools::md5sum(code)
  names(hashes) <- ifelse(startsWith(code, paste0(bundle_dir, "/")), substring(code, nchar(bundle_dir) + 2L), paste0("external/", basename(code)))
  packages <- c("MASS", "Matrix", "glmnet", "Rcpp", "RcppArmadillo", "RcppEigen", "quadprog", "spikeslab")
  runtime <- list(R = as.character(getRversion()), RNG = RNGkind(),
    packages = setNames(vapply(packages, function(x) as.character(utils::packageVersion(x)), ""), packages))
  config <- list(schema = 4L, source_input_protocol = "shared_raw_data_spikeslab_gnet_original_scale",
    source_preprocessing = list(component = "gnet.scale", coefficient_scale = "original_predictors",
      burnin = sim_env_int("SIM_SOURCE_SS_BURNIN", 500L),
      samples = sim_env_int("SIM_SOURCE_SS_SAMPLES", 500L),
      seed_offset = 610000L,
      max_var_rule = "min(p,max(10,floor(n_k/log(max(p,3)))))"),
    translasso = translasso_backend_info(), setting = SETTING_ID, dgp = DGP_BASE_ARGS, algorithm = ALGO_ARGS,
    grid = scenario_grid, base_seed = BASE_SEED, code_md5 = hashes, runtime = runtime)
  manifest <- file.path(result_root, "run_manifest.rds")
  if (file.exists(manifest)) {
    if (!identical(readRDS(manifest), config)) stop("Settings, code or R/package versions changed. Use a NEW SIM_RESULT_ROOT to avoid mixing results.")
  } else {
    if (length(list.files(result_root, all.files = FALSE))) stop("Result root is nonempty without a run manifest. Choose a new SIM_RESULT_ROOT.")
    sim_rds(config, manifest)
    dput(config, file = file.path(result_root, "run_config.txt"))
  }
  sim_csv(scenario_grid, file.path(result_root, "scenario_grid.csv"))
  sim_csv(data.frame(nrep = nrep, workers = workers, selected_cells = paste(selected, collapse = ","),
    started_at = as.character(Sys.time())), file.path(result_root, "last_invocation.csv"))
  writeLines(capture.output(sessionInfo()), file.path(result_root, "sessionInfo.txt"))
  cat(sprintf("[setting %s] %d cells x %d reps; %d concurrent replicate workers\n", SETTING_ID, length(selected), nrep, workers))
  cat("[results] ", result_root, "\n", sep = "")
  all_perf <- all_scen <- all_summary <- statuses <- list()
  # Retain completed cells when this invocation retries only a selected subset.
  previous_status <- if (file.exists(file.path(result_root, "run_status.csv")))
    read.csv(file.path(result_root, "run_status.csv"), stringsAsFactors = FALSE) else data.frame()
  for (idx in scenario_grid$scenario_index) {
    g <- scenario_grid[scenario_grid$scenario_index == idx, , drop = FALSE]
    cell <- file.path(result_root, sim_cell_name(g))
    pfile <- file.path(cell, "performance_all_reps.csv")
    sfile <- file.path(cell, "scenarios_all_reps.csv")
    statusfile <- file.path(cell, "run_status_cell.csv")
    status <- if (file.exists(statusfile)) read.csv(statusfile, stringsAsFactors = FALSE) else
      if (nrow(previous_status)) previous_status[previous_status$scenario_index == idx, , drop = FALSE] else data.frame()
    if (nrow(status) != 1L || !file.exists(pfile) || !file.exists(sfile)) next
    key <- as.character(idx)
    all_perf[[key]] <- sim_read_saved_table(pfile)
    all_scen[[key]] <- sim_read_saved_table(sfile)
    statuses[[key]] <- status
    all_summary[[key]] <- sim_summarize(all_perf[[key]], status$n_requested)
  }
  for (idx in selected) {
    g <- scenario_grid[scenario_grid$scenario_index == idx, , drop = FALSE]
    cell <- file.path(result_root, sim_cell_name(g)); dir.create(cell, showWarnings = FALSE)
    cat(sprintf("[cell %d/%d] s=%d rho=%g snr=%g q=%g\n", idx, nrow(scenario_grid), g$nzc, g$rho, g$snr, g$common_bad_frac))
    fun <- function(rr) sim_run_rep(as.integer(rr), g, cell, resume)
    results <- if (workers > 1L) parallel::mclapply(seq_len(nrep), fun,
      mc.cores = workers, mc.preschedule = FALSE, mc.set.seed = FALSE) else lapply(seq_len(nrep), fun)
    good <- vapply(results, function(x) is.list(x) && isTRUE(x$ok), logical(1))
    cached <- sum(vapply(results[good], function(x) isTRUE(x$cached), logical(1)))
    p <- sim_bind(lapply(results[good], `[[`, "performance")); s <- sim_bind(lapply(results[good], `[[`, "scenario"))
    if (nrow(p)) p <- p[order(p$rep, match(p$method, SIM_METHODS)), , drop = FALSE]
    if (nrow(s)) s <- s[order(s$rep), , drop = FALSE]
    sim_csv(p, file.path(cell, "performance_all_reps.csv")); sim_csv(s, file.path(cell, "scenarios_all_reps.csv"))
    all_summary[[as.character(idx)]] <- sim_summarize(p, nrep)
    sim_summary_tables(all_summary[[as.character(idx)]], cell)
    all_perf[[as.character(idx)]] <- p; all_scen[[as.character(idx)]] <- s
    statuses[[as.character(idx)]] <- data.frame(setting = SETTING_ID, scenario_index = idx,
      n_requested = nrep, n_completed = sum(good), n_resumed = cached, n_failed = sum(!good),
      failed_reps = paste(which(!good), collapse = ";"), complete = all(good))
    sim_csv(statuses[[as.character(idx)]], file.path(cell, "run_status_cell.csv"))
    sim_csv(sim_bind(statuses), file.path(result_root, "run_status.csv"))
    cat(sprintf("[cell %d] completed=%d/%d; resumed=%d; failed=%d\n", idx, sum(good), nrep, cached, sum(!good)))
    # Refresh combined tables after every cell, including partial runs.
    sim_csv(sim_bind(all_perf), file.path(result_root, "performance_all_scenarios.csv"))
    sim_csv(sim_bind(all_scen), file.path(result_root, "scenarios_all_scenarios.csv"))
  }
  summary <- sim_bind(all_summary)
  sim_summary_tables(summary, result_root, "_all_scenarios")
  sim_make_plots(summary, result_root, SETTING_ID)
  if (any(!sim_bind(statuses)$complete)) stop("Incomplete repetitions. Inspect run_status.csv and rep logs; rerun to retry missing reps.")
  cat("[complete] Tables and PDF/PNG plots written to ", result_root, "\n", sep = "")
  invisible(summary)
}
