# Adapter that adds the official linear CONCERT implementation to every
# simulation replicate.  This file is sourced after each setting's functions.R,
# so it wraps the already-installed TRACE/ACVB/Trans-Lasso/TRADER runner.

.concert_adapter_dir <- local({
  ofiles <- Filter(
    function(x) is.character(x) && length(x) == 1L && nzchar(x),
    lapply(sys.frames(), function(x) x$ofile)
  )
  ofile <- if (length(ofiles)) tail(ofiles, 1L)[[1L]] else NULL
  if (is.null(ofile)) getwd() else dirname(normalizePath(ofile))
})

.concert_cpp_file <- file.path(
  .concert_adapter_dir,
  "CONCERT_author_upstream", "Linear", "Concert.cpp"
)
.concert_upstream_commit <- "fea55703c39655a196d3aa5296d1851f0fd53547"
.concert_backend_label <- "official_author_linear_optimized_equivalent"

for (.concert_pkg in c("glmnet", "Rcpp", "RcppArmadillo")) {
  if (!requireNamespace(.concert_pkg, quietly = TRUE)) {
    stop("Package '", .concert_pkg, "' is required for the CONCERT baseline.")
  }
}
if (!file.exists(.concert_cpp_file)) {
  stop("Cannot find the vendored CONCERT C++ source: ", .concert_cpp_file)
}

# Concert.cpp calls this helper from the R global environment for its lasso
# initialization.  The function body matches the upstream Concert.R helper.
lasso_est_R <- function(X, y, lambda) {
  beta <- as.numeric(glmnet::glmnet(X, y, lambda = lambda)$beta)
  beta
}

# Compile once in the parent process.  Forked replicate workers inherit the
# loaded backend, while sourceCpp's cache avoids recompiling on later runs.
Rcpp::sourceCpp(.concert_cpp_file, rebuild = FALSE, verbose = FALSE)
message("[CONCERT] optimized-equivalent author linear backend ", .concert_upstream_commit)

.concert_cv_ridge <- function(X, y, seed) {
  set.seed(as.integer(seed))
  fit <- glmnet::cv.glmnet(
    x = as.matrix(X), y = as.numeric(y), family = "gaussian", alpha = 0
  )
  as.numeric(stats::coef(fit, s = "lambda.min"))[-1L]
}

estimate_concert_qk_author <- function(dat, tolerance = 0.008,
                                       probability_epsilon = 1e-6,
                                       seed = 1L) {
  stopifnot(
    length(dat$X_source) == dat$K,
    length(dat$y_source) == dat$K,
    dat$p == ncol(dat$X)
  )
  tolerance <- as.numeric(tolerance)[1L]
  probability_epsilon <- as.numeric(probability_epsilon)[1L]
  if (!is.finite(tolerance) || tolerance < 0) {
    stop("concert_qk_tolerance must be finite and nonnegative.")
  }
  if (!is.finite(probability_epsilon) || probability_epsilon <= 0 ||
      probability_epsilon >= 0.5) {
    stop("concert_probability_epsilon must lie in (0, 0.5).")
  }

  beta0 <- .concert_cv_ridge(dat$X, dat$y, seed)
  beta_source <- vapply(seq_len(dat$K), function(k) {
    .concert_cv_ridge(dat$X_source[[k]], dat$y_source[[k]], seed + k)
  }, numeric(dat$p))
  if (dat$K == 1L) beta_source <- matrix(beta_source, ncol = 1L)

  agreement <- colMeans(abs(sweep(beta_source, 1L, beta0, "-")) <= tolerance)
  q_raw <- min(agreement)
  q_used <- min(max(q_raw, probability_epsilon), 1 - probability_epsilon)
  list(
    q_raw = q_raw,
    q_used = q_used,
    agreement_by_source = agreement,
    target_ridge_beta = beta0,
    source_ridge_beta = beta_source,
    tolerance = tolerance,
    probability_epsilon = probability_epsilon
  )
}

fit_concert_simulation <- function(dat, seed = 1L,
                                   qk_rule = c("author_ridge", "fixed"),
                                   qk_fixed = 0.5,
                                   qk_tolerance = 0.008,
                                   probability_epsilon = 1e-6,
                                   q0 = 0.1,
                                   tau = 1,
                                   eta = 1,
                                   threshold = 1e-6,
                                   max_iter = 1000L,
                                   selection_threshold = 0.5) {
  qk_rule <- match.arg(qk_rule)
  required <- c("X_all", "y_all", "n_vec_trans", "X_source", "y_source")
  missing <- required[!vapply(required, function(nm) !is.null(dat[[nm]]), logical(1))]
  if (length(missing)) {
    stop("CONCERT requires raw target/source data; missing dat fields: ",
         paste(missing, collapse = ", "))
  }
  stopifnot(
    is.matrix(dat$X_all), ncol(dat$X_all) == dat$p,
    length(dat$y_all) == nrow(dat$X_all),
    length(dat$n_vec_trans) == dat$K + 1L,
    sum(dat$n_vec_trans) == nrow(dat$X_all)
  )

  probability_epsilon <- as.numeric(probability_epsilon)[1L]
  if (!is.finite(probability_epsilon) || probability_epsilon <= 0 ||
      probability_epsilon >= 0.5) {
    stop("concert_probability_epsilon must lie in (0, 0.5).")
  }
  q_info <- if (identical(qk_rule, "author_ridge")) {
    estimate_concert_qk_author(
      dat = dat,
      tolerance = qk_tolerance,
      probability_epsilon = probability_epsilon,
      seed = seed
    )
  } else {
    q_raw <- as.numeric(qk_fixed)[1L]
    if (!is.finite(q_raw) || q_raw < 0 || q_raw > 1) {
      stop("concert_qk_fixed must lie in [0, 1].")
    }
    list(
      q_raw = q_raw,
      q_used = min(max(q_raw, probability_epsilon), 1 - probability_epsilon),
      agreement_by_source = rep(NA_real_, dat$K),
      target_ridge_beta = rep(NA_real_, dat$p),
      source_ridge_beta = matrix(NA_real_, dat$p, dat$K),
      tolerance = qk_tolerance,
      probability_epsilon = probability_epsilon
    )
  }

  q0_used <- min(max(as.numeric(q0)[1L], probability_epsilon),
                 1 - probability_epsilon)
  tau_vec <- rep(as.numeric(tau), length.out = dat$K)
  if (any(!is.finite(tau_vec)) || any(tau_vec <= 0)) {
    stop("concert_tau must be finite and positive.")
  }
  eta <- as.numeric(eta)[1L]
  threshold <- as.numeric(threshold)[1L]
  max_iter <- as.integer(max_iter)[1L]
  selection_threshold <- as.numeric(selection_threshold)[1L]
  if (!is.finite(eta) || eta <= 0 || !is.finite(threshold) || threshold <= 0 ||
      is.na(max_iter) || max_iter < 1L || !is.finite(selection_threshold) ||
      selection_threshold <= 0 || selection_threshold >= 1) {
    stop("Invalid CONCERT eta/threshold/max_iter/selection_threshold.")
  }

  fit <- Concert(
    X = as.matrix(dat$X_all),
    y = as.numeric(dat$y_all),
    n_vec = as.numeric(dat$n_vec_trans),
    tau = tau_vec,
    q_k = q_info$q_used,
    eta = eta,
    q_0 = q0_used,
    threshold = threshold,
    max_iter = max_iter
  )
  required_fit <- c("m_beta0", "gamma_Z", "m_betak", "gamma_I",
                    "V_beta0", "V_betak")
  if (!all(required_fit %in% names(fit))) {
    stop("The official CONCERT backend returned an incomplete result.")
  }
  if (length(fit$m_beta0) != dat$p || length(fit$gamma_Z) != dat$p ||
      !all(dim(fit$gamma_I) == c(dat$p, dat$K))) {
    stop("The official CONCERT backend returned incompatible dimensions.")
  }

  beta_hat <- as.numeric(fit$m_beta0) * as.numeric(fit$gamma_Z)
  if (any(!is.finite(beta_hat)) || any(!is.finite(fit$gamma_Z)) ||
      any(!is.finite(fit$gamma_I))) {
    stop("The official CONCERT backend returned non-finite posterior quantities.")
  }

  fit$method <- "CONCERT_author_linear"
  fit$beta_hat <- beta_hat
  fit$omega <- as.numeric(fit$gamma_Z)
  fit$selected <- fit$omega >= selection_threshold
  fit$q_k_rule <- qk_rule
  fit$q_k_raw <- q_info$q_raw
  fit$q_k_used <- q_info$q_used
  fit$q_k_agreement_by_source <- q_info$agreement_by_source
  fit$q_k_tolerance <- q_info$tolerance
  fit$q_0_used <- q0_used
  fit$tau <- tau_vec
  fit$eta <- eta
  fit$threshold <- threshold
  fit$max_iter <- max_iter
  fit$selection_threshold <- selection_threshold
  fit$upstream_commit <- .concert_upstream_commit
  fit$backend <- .concert_backend_label
  required_diag <- c(
    "n_iter", "final_diff_m", "converged", "hit_max_iter",
    "precompute_time_sec", "initialization_time_sec", "cavi_time_sec",
    "cpp_total_time_sec"
  )
  missing_diag <- required_diag[!required_diag %in% names(fit)]
  if (length(missing_diag)) {
    stop("Optimized CONCERT backend is missing diagnostics: ",
         paste(missing_diag, collapse = ", "))
  }
  fit$n_iter <- as.integer(fit$n_iter)[1L]
  fit$final_diff_m <- as.numeric(fit$final_diff_m)[1L]
  fit$converged <- isTRUE(as.logical(fit$converged)[1L])
  fit$hit_max_iter <- isTRUE(as.logical(fit$hit_max_iter)[1L])
  fit$precompute_time_sec <- as.numeric(fit$precompute_time_sec)[1L]
  fit$initialization_time_sec <- as.numeric(fit$initialization_time_sec)[1L]
  fit$cavi_time_sec <- as.numeric(fit$cavi_time_sec)[1L]
  fit$cpp_total_time_sec <- as.numeric(fit$cpp_total_time_sec)[1L]
  fit
}

make_concert_performance_row <- function(fit, dat, template,
                                         elapsed = NA_real_) {
  if (!exists(".make_na_performance_row", mode = "function", inherits = TRUE) ||
      !exists(".bind_rows_fill", mode = "function", inherits = TRUE)) {
    stop("CONCERT integration must be loaded after trader_bridge.R.")
  }
  row <- .make_na_performance_row(template)
  pred <- if (exists("response_prediction_mse", mode = "function", inherits = TRUE)) {
    response_prediction_mse(fit$beta_hat, dat)
  } else {
    prediction_mse(fit$beta_hat, dat)
  }
  sel <- calc_selection_metrics(
    as.integer(fit$gamma_Z >= fit$selection_threshold), dat$beta_true
  )
  beta_refit <- target_refit_beta(
    dat = dat, fit = fit, thr = fit$selection_threshold, ridge = 1e-6
  )
  pred_refit <- prediction_mse(beta_refit, dat)

  put <- function(nm, value) if (nm %in% names(row)) row[[nm]] <<- value
  put("scenario_id", dat$scenario_id)
  put("method", "CONCERT")
  put("method_name", "CONCERT_author_linear")
  put("estimator_definition", "official_CONCERT_posterior_mean_beta")
  put("pred_mse_raw", pred)
  put("pred_mse_target_refit", pred_refit)
  put("pred_mse_formula", pred)
  put("pred_mse", pred)
  put("pred_mse_definition", PRED_MSE_DEFINITION_V9)
  put("pred_mse_target", PRED_MSE_TARGET_V9)
  put("tp", sel$tp); put("fp", sel$fp); put("fn", sel$fn)
  put("precision", sel$precision); put("recall", sel$recall)
  put("f1", sel$f1); put("model_size", sel$model_size)
  put("mean_omega_signal", mean(fit$gamma_Z[dat$tag_true == 1]))
  put("mean_omega_null", mean(fit$gamma_Z[dat$tag_true == 0]))
  put("criterion_name", "official_CONCERT_CAVI_objective_not_exposed")
  put("criterion_direction", "larger_is_better")
  put("time_sec", elapsed)
  put("n_iter", fit$n_iter)
  put("converged", fit$converged)
  put("convergence_reason", if (fit$converged) {
    "cavi_threshold_reached"
  } else if (fit$hit_max_iter) {
    "max_iter_reached"
  } else {
    "cavi_stopped_without_threshold"
  })

  extras <- data.frame(
    concert_backend = fit$backend,
    concert_upstream_commit = fit$upstream_commit,
    concert_qk_rule = fit$q_k_rule,
    concert_qk_raw = fit$q_k_raw,
    concert_qk_used = fit$q_k_used,
    concert_qk_tolerance = fit$q_k_tolerance,
    concert_q0 = fit$q_0_used,
    concert_tau = paste(signif(fit$tau, 8), collapse = ";"),
    concert_eta = fit$eta,
    concert_threshold = fit$threshold,
    concert_max_iter = fit$max_iter,
    concert_selection_threshold = fit$selection_threshold,
    concert_mean_transfer_probability = mean(fit$gamma_I),
    concert_n_iter = fit$n_iter,
    concert_final_diff_m = fit$final_diff_m,
    concert_converged = fit$converged,
    concert_hit_max_iter = fit$hit_max_iter,
    concert_precompute_time_sec = fit$precompute_time_sec,
    concert_initialization_time_sec = fit$initialization_time_sec,
    concert_cavi_time_sec = fit$cavi_time_sec,
    concert_cpp_total_time_sec = fit$cpp_total_time_sec,
    stringsAsFactors = FALSE
  )
  for (nm in names(extras)) row[[nm]] <- extras[[nm]]
  row
}

.run_one_elbo_without_concert <- run_one_elbo

run_one_elbo <- function(...,
                         run_concert = TRUE,
                         concert_qk_rule = c("author_ridge", "fixed"),
                         concert_qk_fixed = 0.5,
                         concert_qk_tolerance = 0.008,
                         concert_probability_epsilon = 1e-6,
                         concert_q0 = 0.1,
                         concert_tau = 1,
                         concert_eta = 1,
                         concert_threshold = 1e-6,
                         concert_max_iter = 1000L,
                         concert_selection_threshold = 0.5,
                         concert_seed = NULL,
                         save_dir = NULL) {
  concert_qk_rule <- match.arg(concert_qk_rule)
  dots <- list(...)
  ans <- do.call(
    .run_one_elbo_without_concert,
    c(dots, list(save_dir = save_dir))
  )

  ans$scenario$run_concert <- isTRUE(run_concert)
  if (!isTRUE(run_concert)) {
    if (!is.null(save_dir)) {
      write.csv(ans$scenario, file.path(save_dir, "scenario.csv"), row.names = FALSE)
    }
    return(ans)
  }

  seed_base <- dots$seed %||% 123L
  if (is.null(concert_seed)) concert_seed <- as.integer(seed_base) + 900000L

  cat("\n================ fit: CONCERT baseline ================\n")
  timed <- .time_it(fit_concert_simulation(
    dat = ans$dat,
    seed = concert_seed,
    qk_rule = concert_qk_rule,
    qk_fixed = concert_qk_fixed,
    qk_tolerance = concert_qk_tolerance,
    probability_epsilon = concert_probability_epsilon,
    q0 = concert_q0,
    tau = concert_tau,
    eta = concert_eta,
    threshold = concert_threshold,
    max_iter = concert_max_iter,
    selection_threshold = concert_selection_threshold
  ))
  fit <- timed$value
  row <- make_concert_performance_row(
    fit = fit, dat = ans$dat, template = ans$performance,
    elapsed = timed$elapsed
  )
  ans$performance <- .bind_rows_fill(ans$performance, row)
  ans$fits$concert <- fit

  ans$scenario$concert_backend <- fit$backend
  ans$scenario$concert_upstream_commit <- fit$upstream_commit
  ans$scenario$concert_qk_rule <- fit$q_k_rule
  ans$scenario$concert_qk_raw <- fit$q_k_raw
  ans$scenario$concert_qk_used <- fit$q_k_used
  ans$scenario$concert_qk_tolerance <- fit$q_k_tolerance
  ans$scenario$concert_q0 <- fit$q_0_used
  ans$scenario$concert_tau <- paste(signif(fit$tau, 8), collapse = ";")
  ans$scenario$concert_eta <- fit$eta
  ans$scenario$concert_threshold <- fit$threshold
  ans$scenario$concert_max_iter <- fit$max_iter
  ans$scenario$concert_selection_threshold <- fit$selection_threshold
  ans$scenario$concert_seed <- concert_seed
  ans$scenario$concert_n_iter <- fit$n_iter
  ans$scenario$concert_final_diff_m <- fit$final_diff_m
  ans$scenario$concert_converged <- fit$converged
  ans$scenario$concert_hit_max_iter <- fit$hit_max_iter
  ans$scenario$concert_precompute_time_sec <- fit$precompute_time_sec
  ans$scenario$concert_initialization_time_sec <- fit$initialization_time_sec
  ans$scenario$concert_cavi_time_sec <- fit$cavi_time_sec
  ans$scenario$concert_cpp_total_time_sec <- fit$cpp_total_time_sec

  if (!is.null(save_dir)) {
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    write.csv(ans$performance, file.path(save_dir, "performance.csv"), row.names = FALSE)
    write.csv(ans$scenario, file.path(save_dir, "scenario.csv"), row.names = FALSE)
    write.csv(data.frame(
      source = seq_len(ans$dat$K),
      ridge_agreement = fit$q_k_agreement_by_source,
      q_k_raw = fit$q_k_raw,
      q_k_used = fit$q_k_used,
      q_k_rule = fit$q_k_rule,
      q_k_tolerance = fit$q_k_tolerance,
      upstream_commit = fit$upstream_commit,
      backend = fit$backend,
      n_iter = fit$n_iter,
      final_diff_m = fit$final_diff_m,
      converged = fit$converged,
      hit_max_iter = fit$hit_max_iter,
      precompute_time_sec = fit$precompute_time_sec,
      initialization_time_sec = fit$initialization_time_sec,
      cavi_time_sec = fit$cavi_time_sec,
      cpp_total_time_sec = fit$cpp_total_time_sec,
      stringsAsFactors = FALSE
    ), file.path(save_dir, "concert_diagnostics.csv"), row.names = FALSE)
    coordinate <- data.frame(
      j = seq_len(ans$dat$p),
      beta_true = ans$dat$beta_true,
      beta_hat = fit$beta_hat,
      gamma_Z = fit$gamma_Z,
      selected = fit$selected,
      stringsAsFactors = FALSE
    )
    for (k in seq_len(ans$dat$K)) {
      coordinate[[paste0("gamma_I_source", k)]] <- fit$gamma_I[, k]
    }
    write.csv(
      coordinate,
      file.path(save_dir, "concert_coordinate_posterior.csv"),
      row.names = FALSE
    )
  }

  ans
}
