## =========================================================
## Bridge between Borrow-CVB simulation objects and TRADER
## =========================================================

.trader_bridge_dir <- local({
  ofiles <- Filter(function(x) is.character(x) && length(x) == 1L && nzchar(x),
                   lapply(sys.frames(), function(x) x$ofile))
  ofile <- if (length(ofiles)) tail(ofiles, 1L)[[1L]] else NULL
  if (is.null(ofile)) getwd() else dirname(normalizePath(ofile))
})
source(file.path(.trader_bridge_dir, "trader_core.R"))

if (!exists("%||%", mode = "function", inherits = TRUE)) {
  `%||%` <- function(x, y) if (is.null(x)) y else x
}

.trader_elapsed <- function(expr) {
  tm <- system.time(value <- eval.parent(substitute(expr)))
  list(value = value, elapsed = unname(tm[["elapsed"]]))
}

.trader_safe_cor <- function(x, y) {
  x <- as.numeric(x); y <- as.numeric(y)
  if (length(x) != length(y) || length(x) < 2L) return(NA_real_)
  sx <- stats::sd(x); sy <- stats::sd(y)
  if (!is.finite(sx) || !is.finite(sy) || sx <= 0 || sy <= 0) return(NA_real_)
  suppressWarnings(stats::cor(x, y))
}

.trader_spikeslab_fit <- function(x, y,
                                  profile = c("current", "setting_I", "setting_III"),
                                  max_var = NULL) {
  profile <- match.arg(profile)
  if (!requireNamespace("spikeslab", quietly = TRUE)) {
    stop(
      "Package 'spikeslab' is required for the exact preprocessing used in the original TRADER simulations.\n",
      "Install it once with:\n",
      "  install.packages('spikeslab', repos = 'https://cloud.r-project.org', dependencies = TRUE)\n",
      "Then restart R and verify with packageVersion('spikeslab').\n",
      "For a non-exact sensitivity run only, set allow_estimator_fallback = TRUE ",
      "or estimator = 'glmnet'."
    )
  }
  x <- as.matrix(x); y <- as.numeric(y)
  args <- list(y = y, x = x, verbose = FALSE, intercept = FALSE)
  if (profile == "setting_I") {
    args$bigp.smalln <- FALSE
    args$fast <- TRUE
  } else {
    args$bigp.smalln <- TRUE
    if (!is.null(max_var)) args$max.var <- as.integer(max_var)
  }
  fit <- do.call(spikeslab::spikeslab, args)
  beta <- as.numeric(fit$gnet)
  if (length(beta) != ncol(x)) stop("spikeslab returned a coefficient vector of the wrong length.")
  beta[!is.finite(beta)] <- 0
  beta
}

.trader_glmnet_fallback <- function(x, y, seed = 1) {
  if (!requireNamespace("glmnet", quietly = TRUE)) stop("Package 'glmnet' is required for fallback estimation.")
  x <- as.matrix(x); y <- as.numeric(y)
  set.seed(seed)
  nfolds <- min(10L, max(3L, floor(nrow(x) / 10L)))
  fit <- glmnet::cv.glmnet(
    x = x, y = y, alpha = 1, intercept = FALSE,
    standardize = TRUE, nfolds = nfolds
  )
  beta <- as.numeric(stats::coef(fit, s = "lambda.min"))[-1]
  beta[!is.finite(beta)] <- 0
  beta
}

.trader_estimate_beta <- function(x, y,
                                  estimator = c("spikeslab", "glmnet"),
                                  profile = c("current", "setting_I", "setting_III"),
                                  max_var = NULL,
                                  allow_fallback = FALSE,
                                  seed = 1) {
  estimator <- match.arg(estimator)
  profile <- match.arg(profile)
  if (estimator == "glmnet") {
    return(list(beta = .trader_glmnet_fallback(x, y, seed = seed), estimator_used = "glmnet"))
  }
  ans <- tryCatch(
    .trader_spikeslab_fit(x, y, profile = profile, max_var = max_var),
    error = function(e) e
  )
  if (!inherits(ans, "error")) return(list(beta = ans, estimator_used = "spikeslab"))
  if (!isTRUE(allow_fallback)) stop(conditionMessage(ans))
  warning("spikeslab preprocessing failed; using glmnet fallback: ", conditionMessage(ans))
  list(beta = .trader_glmnet_fallback(x, y, seed = seed), estimator_used = "glmnet_fallback")
}

.trader_dimension_only_max_var <- function(dat) {
  min(
    as.integer(dat$p),
    max(10L, as.integer(floor(dat$n / log(max(dat$p, 3L)))))
  )
}

estimate_trader_external_coefficients <- function(
    dat,
    source_input = c("B", "raw"),
    estimator = c("spikeslab", "glmnet"),
    profile = c("current", "setting_I", "setting_III"),
    max_var = NULL,
    allow_fallback = FALSE,
    seed = 1) {
  source_input <- match.arg(source_input)
  estimator <- match.arg(estimator)
  profile <- match.arg(profile)
  set.seed(as.integer(seed))
  ## Final-paper Setting I calls spikeslab without max.var. Setting III passes
  ## max.var=s explicitly from the scenario runner. Other profiles retain the
  ## dimension-only rule, which never reads beta_true.
  if (is.null(max_var) && profile != "setting_I") {
    max_var <- .trader_dimension_only_max_var(dat)
  }
  if (!is.null(max_var)) {
    max_var <- max(1L, min(as.integer(max_var), as.integer(dat$p)))
  }

  target_fit <- .trader_estimate_beta(
    dat$X, dat$y,
    estimator = estimator,
    profile = profile,
    max_var = max_var,
    allow_fallback = allow_fallback,
    seed = seed
  )

  if (source_input == "B") {
    W_hat <- as.matrix(dat$B)
    source_estimators <- rep("provided_B", ncol(W_hat))
  } else {
    if (is.null(dat$X_source) || is.null(dat$y_source)) {
      stop("source_input='raw' requires dat$X_source and dat$y_source.")
    }
    K <- length(dat$X_source)
    W_hat <- matrix(0, dat$p, K)
    source_estimators <- character(K)
    for (k in seq_len(K)) {
      fk <- .trader_estimate_beta(
        dat$X_source[[k]], dat$y_source[[k]],
        estimator = estimator,
        profile = profile,
        max_var = max_var,
        allow_fallback = allow_fallback,
        seed = seed + k
      )
      W_hat[, k] <- fk$beta
      source_estimators[k] <- fk$estimator_used
    }
  }

  list(
    beta_target = target_fit$beta,
    W_hat = W_hat,
    target_estimator = target_fit$estimator_used,
    source_estimators = source_estimators,
    source_input = source_input,
    profile = profile,
    max_var_used = if (is.null(max_var)) NA_integer_ else max_var,
    max_var_rule = if (is.null(max_var)) {
      "official_setting_I_uncapped"
    } else if (profile == "setting_III") {
      "official_setting_III_s"
    } else {
      "dimension_only_floor_n_over_log_p"
    }
  )
}

screen_trader_sources <- function(beta_target, W_hat, corr_threshold = 0.1,
                                  min_column_norm = 1e-10) {
  W_hat <- as.matrix(W_hat)
  K <- ncol(W_hat)
  if (K == 0L) {
    return(list(W_kept = W_hat, correlations = numeric(0), kept = integer(0), dropped = integer(0)))
  }
  correlations <- vapply(seq_len(K), function(k) .trader_safe_cor(W_hat[, k], beta_target), numeric(1))
  norms <- sqrt(colSums(W_hat^2))
  keep <- which(is.finite(correlations) & correlations > corr_threshold &
                  is.finite(norms) & norms > min_column_norm)
  list(
    W_kept = W_hat[, keep, drop = FALSE],
    correlations = correlations,
    kept_correlations = correlations[keep],
    kept = keep,
    dropped = setdiff(seq_len(K), keep),
    norms = norms
  )
}

.trader_selection_metrics <- function(selected, beta_true) {
  selected <- as.logical(selected)
  truth <- beta_true != 0
  tp <- sum(selected & truth)
  fp <- sum(selected & !truth)
  fn <- sum(!selected & truth)
  precision <- if (tp + fp > 0) tp / (tp + fp) else 0
  recall <- if (tp + fn > 0) tp / (tp + fn) else 0
  f1 <- if (precision + recall > 0) 2 * precision * recall / (precision + recall) else 0
  list(tp = tp, fp = fp, fn = fn, precision = precision, recall = recall,
       f1 = f1, model_size = sum(selected), power = recall,
       fdr = if (tp + fp > 0) fp / (tp + fp) else 0)
}

fit_trader_baseline <- function(
    dat,
    source_input = c("B", "raw"),
    estimator = c("spikeslab", "glmnet"),
    profile = c("current", "setting_I", "setting_III"),
    corr_threshold = 0.1,
    max_var = NULL,
    allow_estimator_fallback = FALSE,
    n_mcmc = 10000,
    n_burnin = 5000,
    thin = 2,
    n_chain = 1,
    fold = 1.1,
    v = 0.001,
    tau = "Cauchy",
    tau_lb = NULL,
    tau_ub = 1,
    residual_alpha = 2,
    sign_threshold = 0.95,
    seed = 123,
    external_coefficients = NULL,
    external_preprocess_time_sec = NULL,
    keep_samples = FALSE,
    verbose = TRUE) {
  source_input <- match.arg(source_input)
  estimator <- match.arg(estimator)
  profile <- match.arg(profile)
  if (is.null(tau_lb)) tau_lb <- 1 / dat$p

  if (is.null(external_coefficients)) {
    prep_tm <- .trader_elapsed(estimate_trader_external_coefficients(
      dat = dat,
      source_input = source_input,
      estimator = estimator,
      profile = profile,
      max_var = max_var,
      allow_fallback = allow_estimator_fallback,
      seed = seed
    ))
    ext <- prep_tm$value
  } else {
    ext <- external_coefficients
    if (is.null(ext$beta_target) || is.null(ext$W_hat)) {
      stop("external_coefficients must contain beta_target and W_hat.")
    }
    if (nrow(as.matrix(ext$W_hat)) != dat$p) {
      stop("external_coefficients$W_hat has the wrong number of rows.")
    }
    prep_tm <- list(
      value = ext,
      elapsed = if (is.null(external_preprocess_time_sec)) 0 else as.numeric(external_preprocess_time_sec)
    )
  }
  scr <- screen_trader_sources(ext$beta_target, ext$W_hat, corr_threshold = corr_threshold)

  if (verbose) {
    cat(sprintf("[TRADER] source input=%s | kept %d/%d sources | correlations=(%s)\n",
                source_input, length(scr$kept), ncol(ext$W_hat),
                paste(format(round(scr$correlations, 3), nsmall = 3), collapse = ", ")))
  }

  alpha_dirichlet <- c(pmax(scr$kept_correlations, 1e-4), residual_alpha)
  W_for_trader <- scr$W_kept
  if (ncol(W_for_trader) == 0L) W_for_trader <- matrix(0, dat$p, 0L)

  obj <- CreateTRADERObject(
    y = dat$y,
    x = dat$X,
    w = W_for_trader,
    alpha = alpha_dirichlet,
    v = v,
    tau = tau,
    tau_lb = tau_lb,
    tau_ub = tau_ub
  )

  mcmc_tm <- .trader_elapsed(TRADER_MCMC(
    obj,
    n_mcmc = n_mcmc,
    n_burnin = n_burnin,
    n_chain = n_chain,
    fold = fold,
    thin = thin,
    seed = seed
  ))
  obj <- mcmc_tm$value

  beta_samples <- obj@mcmc_samples$beta
  beta_hat <- rowMeans(beta_samples)
  beta_q <- t(apply(beta_samples, 1, stats::quantile, probs = c(0.025, 0.975)))
  sign_certainty <- apply(beta_samples, 1, function(x) abs(sum(x > 0) - sum(x < 0)) / length(x))
  selected <- sign_certainty >= sign_threshold
  sel <- .trader_selection_metrics(selected, dat$beta_true)

  eta_mean_full <- rep(0, dat$K)
  if (length(scr$kept) > 0L && nrow(obj@mcmc_samples$eta) > 0L) {
    eta_mean_full[scr$kept] <- rowMeans(obj@mcmc_samples$eta)
  }
  psrf <- obj@mcmc_samples$PSRF_beta

  out <- list(
    method = "TRADER",
    beta_hat = beta_hat,
    beta_q025 = beta_q[, 1],
    beta_q975 = beta_q[, 2],
    sign_certainty = sign_certainty,
    selected = selected,
    sel_credible = sel,
    eta_mean = eta_mean_full,
    kept_sources = scr$kept,
    dropped_sources = scr$dropped,
    source_correlations = scr$correlations,
    target_estimate = ext$beta_target,
    W_hat_input = ext$W_hat,
    source_input = source_input,
    estimator_requested = estimator,
    target_estimator = ext$target_estimator,
    source_estimators = ext$source_estimators,
    max_var_used = ext$max_var_used %||% if (is.null(max_var)) {
      .trader_dimension_only_max_var(dat)
    } else {
      as.integer(max_var)
    },
    max_var_rule = ext$max_var_rule %||% if (is.null(max_var)) {
      "dimension_only_floor_n_over_log_p"
    } else {
      "user_supplied_nonoracle"
    },
    profile = profile,
    corr_threshold = corr_threshold,
    sign_threshold = sign_threshold,
    n_mcmc = as.integer(n_mcmc),
    n_burnin = as.integer(n_burnin),
    n_draws = ncol(beta_samples),
    thin = as.integer(thin),
    n_chain = as.integer(n_chain),
    preprocess_time_sec = prep_tm$elapsed,
    mcmc_time_sec = mcmc_tm$elapsed,
    total_time_sec = prep_tm$elapsed + mcmc_tm$elapsed,
    mcmc_backend = obj@mcmc_samples$backend,
    residual_refresh_count = obj@mcmc_samples$residual_refresh_count,
    max_residual_refresh_error = obj@mcmc_samples$max_residual_refresh_error,
    mcmc_setup_time_sec = obj@mcmc_samples$setup_time_sec,
    mcmc_cpp_time_sec = obj@mcmc_samples$cpp_time_sec,
    seconds_per_transition = obj@mcmc_samples$cpp_time_sec /
      max(1, as.numeric(n_mcmc) * as.numeric(thin) * as.numeric(n_chain)),
    max_psrf = if (any(is.finite(psrf))) max(psrf, na.rm = TRUE) else NA_real_,
    mean_psrf = if (any(is.finite(psrf))) mean(psrf, na.rm = TRUE) else NA_real_,
    n_iter = as.integer(n_mcmc),
    converged = if (any(is.finite(psrf))) max(psrf, na.rm = TRUE) < 1.1 else NA,
    convergence_reason = "split-chain PSRF diagnostic",
    final_elbo = NA_real_,
    elbo_trace = data.frame(),
    elbo_drop_count = NA_integer_
  )
  if (isTRUE(keep_samples)) out$beta_samples <- beta_samples
  out
}

.trader_na_value <- function(x) {
  if (is.character(x)) NA_character_
  else if (is.logical(x)) NA
  else if (is.integer(x)) NA_integer_
  else NA_real_
}

.make_na_performance_row <- function(template) {
  if (is.null(template) || nrow(template) == 0L) stop("A nonempty performance template is required.")
  out <- template[1, , drop = FALSE]
  for (nm in names(out)) out[[nm]] <- .trader_na_value(out[[nm]])
  out
}

.bind_rows_fill <- function(x, y) {
  nx <- names(x); ny <- names(y)
  for (nm in setdiff(ny, nx)) x[[nm]] <- .trader_na_value(y[[nm]])
  for (nm in setdiff(nx, ny)) y[[nm]] <- .trader_na_value(x[[nm]])
  y <- y[, names(x), drop = FALSE]
  rbind(x, y)
}

make_trader_performance_row <- function(fit, dat, template,
                                        target_logp_sigma_beta_sq = 100) {
  row <- .make_na_performance_row(template)
  beta_hat <- fit$beta_hat
  ## Use exactly the same observed-response test metric as every other method.
  ## No training-error fallback is allowed because that would change the target.
  pred <- response_prediction_mse(beta_hat, dat)
  sel <- fit$sel_credible

  put <- function(nm, value) if (nm %in% names(row)) row[[nm]] <<- value
  put("scenario_id", dat$scenario_id)
  put("method", "trader")
  put("method_name", "TRADER")

  for (suffix in c("raw", "formula")) {
    put(paste0("pred_mse_", suffix), pred)
  }
  put("pred_mse", pred)
  put("pred_mse_definition", PRED_MSE_DEFINITION_V9)
  put("pred_mse_target", PRED_MSE_TARGET_V9)

  put("tp", sel$tp); put("fp", sel$fp); put("fn", sel$fn)
  put("precision", sel$precision); put("recall", sel$recall)
  put("f1", sel$f1); put("model_size", sel$model_size)
  put("criterion_name", "MCMC_posterior_sign_selection_no_ELBO")
  put("criterion_direction", "larger_is_better")

  if (exists("compute_target_logp_y_z", mode = "function", inherits = TRUE)) {
    lp <- tryCatch(compute_target_logp_y_z(
      z = as.numeric(fit$selected), dat = dat,
      sigma_beta_sq = target_logp_sigma_beta_sq,
      mode = "selected", threshold = 0.5
    )$logp, error = function(e) NA_real_)
    put("criterion_logp", lp)
    put("target_logp_y_z", lp)
  }

  put("time_sec", fit$total_time_sec)
  put("n_iter", fit$n_iter)
  put("final_elbo", NA_real_)
  put("initial_elbo", NA_real_)
  put("total_elbo_gain", NA_real_)
  put("max_elbo_drop", NA_real_)
  put("elbo_drop_count", NA_integer_)
  put("converged", fit$converged)
  put("convergence_reason", fit$convergence_reason)

  extras <- data.frame(
    trader_source_input = fit$source_input,
    trader_profile = fit$profile,
    trader_target_estimator = fit$target_estimator,
    trader_source_estimators = paste(unique(fit$source_estimators), collapse = ";"),
    trader_corr_threshold = fit$corr_threshold,
    trader_sign_threshold = fit$sign_threshold,
    trader_max_var = fit$max_var_used,
    trader_max_var_rule = fit$max_var_rule,
    trader_k_input = dat$K,
    trader_k_kept = length(fit$kept_sources),
    trader_kept_sources = paste(fit$kept_sources, collapse = ";"),
    trader_source_correlations = paste(round(fit$source_correlations, 6), collapse = ";"),
    trader_eta_hat = paste(round(fit$eta_mean, 6), collapse = ";"),
    trader_n_mcmc = fit$n_mcmc,
    trader_n_burnin = fit$n_burnin,
    trader_n_draws = fit$n_draws,
    trader_thin = fit$thin,
    trader_preprocess_time_sec = fit$preprocess_time_sec,
    trader_mcmc_time_sec = fit$mcmc_time_sec,
    trader_mcmc_backend = fit$mcmc_backend,
    trader_mcmc_setup_time_sec = fit$mcmc_setup_time_sec,
    trader_mcmc_cpp_time_sec = fit$mcmc_cpp_time_sec,
    trader_seconds_per_transition = fit$seconds_per_transition,
    trader_residual_refresh_count = fit$residual_refresh_count,
    trader_max_residual_refresh_error = fit$max_residual_refresh_error,
    trader_max_psrf = fit$max_psrf,
    trader_mean_psrf = fit$mean_psrf,
    trader_power = sel$power,
    trader_fdr = sel$fdr,
    stringsAsFactors = FALSE
  )
  for (nm in names(extras)) row[[nm]] <- extras[[nm]]
  row
}
