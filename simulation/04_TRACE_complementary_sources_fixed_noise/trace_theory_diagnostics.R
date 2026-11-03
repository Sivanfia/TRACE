## Diagnostics for the 2026-09-12 theory-aligned run. Truth is read only here,
## after fitting, and never enters TRACE initialization or its updates.
TRACE_THEORY_VERSION <- "20260913_fixed_target_steps"
.run_one_elbo_theory_base <- run_one_elbo
run_one_elbo <- function(..., save_dir = NULL) {
  ans <- .run_one_elbo_theory_base(..., save_dir = save_dir)
  dat <- ans$dat
  true_active <- dat$beta_true != 0
  population_gap <- as.numeric(dat$W_source_true %*% dat$alpha_true) - dat$beta_true
  true_C0 <- true_active & abs(population_gap) <= 1e-8
  perf <- ans$performance
  perf$theory_version <- TRACE_THEORY_VERSION
  perf$rho_z <- ans$prior_par_borrow$rho_z_mean
  perf$rho_c <- ans$prior_par_borrow$rho_c_mean
  perf$source_ns_raw <- trace_get_raw_ns_eff(dat)
  perf$source_precision_ns <- get_ns_eff(dat)
  perf$true_C0_size <- sum(true_C0)
  perf$target_init_mode <- NA_character_
  perf$target_init_steps_requested <- NA_integer_
  perf$target_phase_complete <- NA
  perf$target_init_converged <- NA
  perf$target_init_n_iter <- NA_integer_
  perf$target_init_max_change <- NA_real_
  perf$target_init_size <- NA_integer_
  perf$target_init_tol <- NA_real_
  perf$target_init_damping_trigger_iter <- NA_integer_
  for (col in c("borrowing_tp", "borrowing_fp", "borrowing_fn", "borrowing_precision",
                "borrowing_recall", "borrowing_f1", "credibility_size", "transferable_size",
                "gamma_null_mean", "null_gamma_distance_from_prior")) perf[[col]] <- NA_real_
  coord <- data.frame(j=seq_len(dat$p), beta_true=dat$beta_true,
    coordinate_type=ifelse(true_active,ifelse(true_C0,"transferable_active","active_nontransferable"),
                          ifelse(seq_len(dat$p) %in% dat$spurious_union,"source_only_null","target_null")),
    source_center_true=as.numeric(dat$W_source_true %*% dat$alpha_true))
  coord$true_C0 <- as.integer(true_C0)
  coord$population_source_gap <- population_gap
  states <- list()
  for (nm in names(ans$fits)) {
    fit <- ans$fits[[nm]]
    states[[nm]] <- fit[intersect(names(fit), c(
      "omega", "gamma", "alpha", "beta_hat", "n_iter", "converged", "convergence_reason",
      "target_initialization", "target_init_converged", "target_init_n_iter", "target_init_max_change", "target_init_size", "max_param_change_final", "sigma_beta_sq", "theory_init_trace", "elbo_trace",
      "fixed_E_logit_rho_z", "fixed_E_logit_rho_c", "preliminary_alpha_beta_target",
      "preliminary_alpha_huber_tau", "preliminary_alpha_huber_kkt_residual"))]
    if (is.null(fit$gamma) || is.null(fit$omega)) next
    row <- if ("method_key" %in% names(perf)) which(perf$method_key == nm) else which(perf$method == nm)
    if (length(row) != 1L) stop("Cannot match diagnostic fit to performance row: ", nm)
    for (field in c("target_init_converged","target_init_n_iter","target_init_max_change","target_init_size"))
      perf[row,field] <- fit[[field]] %||% NA
    if(!is.null(fit$target_initialization)) {
      perf$target_init_mode[row] <- if(isTRUE(fit$target_init_until_converged)) "until_converged" else "fixed_steps"
      perf$target_init_steps_requested[row] <- if(isTRUE(fit$target_init_until_converged)) NA_integer_ else fit$target_init_steps
      perf$target_phase_complete[row] <- TRUE
      perf$target_init_tol[row] <- fit$target_initialization$tol
      perf$target_init_damping_trigger_iter[row] <- fit$target_initialization$damping_trigger_iter
      coord[[paste0(nm,"_omega_target_initialized")]] <- fit$target_initialization$omega
      if(!is.null(save_dir)) write.csv(fit$target_initialization$trace,
        file.path(save_dir,paste0("target_only_trace_",nm,".csv")),row.names=FALSE)
    }
    selected_active <- fit$omega > 0.5
    selected_cred <- fit$gamma > 0.5
    selected_C <- selected_active & selected_cred
    met <- calc_selection_metrics(selected_C, as.numeric(true_C0))
    for (m in c("tp", "fp", "fn", "precision", "recall", "f1")) perf[row,paste0("borrowing_",m)] <- met[[m]]
    perf$credibility_size[row] <- sum(selected_cred)
    perf$transferable_size[row] <- sum(selected_C)
    perf$gamma_null_mean[row] <- mean(fit$gamma[!true_active])
    perf$null_gamma_distance_from_prior[row] <- mean(abs(fit$gamma[!true_active] - perf$rho_c[row]))
    coord[[paste0(nm,"_selected_S")]] <- as.integer(selected_active)
    coord[[paste0(nm,"_selected_Ccred")]] <- as.integer(selected_cred)
    coord[[paste0(nm,"_selected_C")]] <- as.integer(selected_C)
  }
  ans$performance <- perf
  ans$coordinate_diagnostics <- coord
  ans$scenario$theory_version <- TRACE_THEORY_VERSION
  ans$scenario$rho_c <- ans$prior_par_borrow$rho_c_mean
  ans$scenario$true_C0_size <- sum(true_C0)
  ans$scenario$h6_nonempty_C0 <- any(true_C0)
  ans$scenario$h7_all_subsets_verified <- FALSE
  ans$scenario$main_estimator <- "method_returned_beta; TRACE paper model average"
  if (!is.null(save_dir)) {
    write.csv(perf, file.path(save_dir,"performance.csv"), row.names=FALSE)
    write.csv(coord, file.path(save_dir,"coordinate_diagnostics.csv"), row.names=FALSE)
    write.csv(ans$scenario, file.path(save_dir,"scenario.csv"), row.names=FALSE)
    saveRDS(list(version=TRACE_THEORY_VERSION, states=states, prior=ans$prior_par_borrow,
                 source_ns_raw=trace_get_raw_ns_eff(dat), source_precision_ns=get_ns_eff(dat)),
            file.path(save_dir,"theory_fit_states.rds"))
  }
  ans
}
