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
    args$bigp.smalln <- ncol(x) > nrow(x)
    # GTEx can have n > p; spikeslab rejects its high-dimensional path there.
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
  ## Preprocessing capacity must not use beta_true or the true support size.
  if (is.null(max_var)) max_var <- .trader_dimension_only_max_var(dat)
  max_var <- max(1L, min(as.integer(max_var), as.integer(dat$p)))

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
    max_var_used = max_var,
    max_var_rule = "dimension_only_floor_n_over_log_p"
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
