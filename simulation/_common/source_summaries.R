# Fit shared summaries using source observations only, with no population inputs.
validate_shared_source_sizes <- function(n_source_raw, K, p,
                                         ns_summary = NULL,
                                         source_noise_scale = 1) {
  if (length(K) != 1L || !is.finite(K) || K < 1 || K != floor(K) ||
      length(p) != 1L || !is.finite(p) || p < 1 || p != floor(p))
    stop("Source counts require positive integer K and p.")
  if (length(n_source_raw) == 1L) n_source_raw <- rep(n_source_raw, K)
  if (length(n_source_raw) != K || any(!is.finite(n_source_raw)) ||
      any(n_source_raw != floor(n_source_raw)) || any(n_source_raw < 3))
    stop("Source counts must be integers of at least 3, with one count per source.")
  if (!is.null(ns_summary)) {
    if (length(ns_summary) == 1L) ns_summary <- rep(ns_summary, K)
    if (length(ns_summary) != K || any(!is.finite(ns_summary)) ||
        any(ns_summary != n_source_raw))
      stop("ns_summary must equal the actual source counts.")
  }
  if (length(source_noise_scale) != 1L || !is.finite(source_noise_scale) ||
      source_noise_scale != 1)
    stop("source_noise_scale is obsolete: summary errors come from the actual source responses.")
  as.integer(n_source_raw)
}

.source_summary_integer <- function(name, default, minimum = 1L) {
  raw <- Sys.getenv(name, unset = as.character(default))
  value <- suppressWarnings(as.numeric(raw))
  if (length(value) != 1L || !is.finite(value) || value < minimum ||
      value != floor(value) || value > .Machine$integer.max)
    stop(name, " must be an integer >= ", minimum, ".")
  as.integer(value)
}

fit_shared_source_summaries <- function(X_source, y_source, n_source_raw,
                                         seed = 20261012L) {
  started <- proc.time()[["elapsed"]]
  if (!requireNamespace("spikeslab", quietly = TRUE))
    stop("Package 'spikeslab' is required for the shared high-dimensional summaries.")
  K <- length(X_source)
  if (!is.list(X_source) || K < 1L || !is.list(y_source) || length(y_source) != K)
    stop("X_source and y_source must be nonempty lists of equal length.")
  if (!is.matrix(X_source[[1]]) || !is.numeric(X_source[[1]]))
    stop("Source designs must be numeric matrices.")
  p <- ncol(X_source[[1]])
  counts <- validate_shared_source_sizes(n_source_raw, K, p)
  if (length(seed) != 1L || !is.finite(seed) || seed < 1 || seed != floor(seed) ||
      seed > .Machine$integer.max - K)
    stop("Summary seed must be a positive integer with room for per-source seeds.")
  seed <- as.integer(seed)
  burnin <- .source_summary_integer("SIM_SOURCE_SS_BURNIN", 500L)
  samples <- .source_summary_integer("SIM_SOURCE_SS_SAMPLES", 500L)

  # Preserve the DGP's downstream RNG stream and paired designs.
  had_rng <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_rng <- if (had_rng) get(".Random.seed", envir = .GlobalEnv) else NULL
  old_kind <- RNGkind()
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_rng) assign(".Random.seed", old_rng, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)

  B <- matrix(NA_real_, p, K)
  diagnostics <- vector("list", K)
  for (k in seq_len(K)) {
    X <- X_source[[k]]; y <- y_source[[k]]
    if (!is.matrix(X) || !is.numeric(X) || ncol(X) != p ||
        !is.numeric(y) || length(y) != nrow(X) ||
        any(!is.finite(X)) || any(!is.finite(y)))
      stop("Every source must have matching finite numeric X and y with a common predictor dimension.")
    if (nrow(X) != counts[k]) stop("Declared source count differs from the actual observations.")
    if (any(apply(X, 2L, stats::sd) <= 0) || stats::sd(y) <= 0)
      stop("Source predictors and responses must be nonconstant.")

    source_seed <- seed + as.integer(k)
    capacity <- min(p, max(10L, as.integer(floor(nrow(X) / log(max(p, 3L))))))
    high_dimensional <- p >= nrow(X)
    start_k <- proc.time()[["elapsed"]]
    set.seed(source_seed)
    fit <- spikeslab::spikeslab(
      x = X, y = y, intercept = FALSE, verbose = FALSE, fast = TRUE,
      bigp.smalln = high_dimensional, max.var = capacity,
      n.iter1 = burnin, n.iter2 = samples, seed = -source_seed
    )
    # gnet uses standardized-x units; gnet.scale uses original predictor units.
    beta <- as.numeric(fit$gnet.scale)
    if (length(beta) != p || any(!is.finite(beta)))
      stop("spikeslab returned invalid original-scale coefficients; no estimator fallback is used.")
    B[, k] <- beta
    diagnostics[[k]] <- data.frame(
      source = k, n_source = counts[k], p = p,
      method = "spikeslab_gnet_original_scale", coefficient_scale = "original_predictors",
      summary_component = "gnet.scale", bigp_smalln = high_dimensional,
      max_var = capacity, burnin = burnin, samples = samples,
      preprocessing_seed = source_seed, selected_count = sum(beta != 0),
      training_residual_mse = mean((y - as.numeric(X %*% beta))^2),
      coefficient_norm = sqrt(sum(beta^2)),
      package_version = as.character(utils::packageVersion("spikeslab")),
      fit_time_sec = proc.time()[["elapsed"]] - start_k,
      stringsAsFactors = FALSE
    )
  }
  # Regularized summaries do not come with OLS standard errors.
  list(B = B, standard_errors = NULL, counts = counts,
       diagnostics = do.call(rbind, diagnostics),
       fit_time_sec = proc.time()[["elapsed"]] - started)
}
