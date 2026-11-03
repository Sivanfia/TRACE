## Lai et al. (2026), final-paper Simulation Settings I and III.
## Source summaries B are estimated once from the raw source samples and then
## supplied unchanged to both TRACE and TRADER.

.paper_ar1_design <- function(n, p, rho) {
  z <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
  if (p > 1L) {
    innovation_sd <- sqrt(1 - rho^2)
    for (j in 2:p) z[, j] <- rho * z[, j - 1L] + innovation_sd * z[, j]
  }
  z
}

.paper_ar1_signal_variance <- function(s, rho, signal_strength) {
  lag_sum <- sum((s - seq_len(s - 1L)) * rho^seq_len(s - 1L))
  signal_strength^2 * (s + 2 * lag_sum)
}

.paper_spikeslab_summary <- function(x, y, paper_setting, s) {
  if (identical(paper_setting, "III")) {
    fit <- spikeslab::spikeslab(
      y = y, x = x, verbose = FALSE, max.var = s,
      bigp.smalln = TRUE, intercept = FALSE
    )
  } else {
    fit <- spikeslab::spikeslab(
      y = y, x = x, verbose = FALSE, bigp.smalln = FALSE,
      fast = TRUE, intercept = FALSE
    )
  }
  out <- as.numeric(fit$gnet)
  if (length(out) != ncol(x) || any(!is.finite(out))) {
    stop("spikeslab did not return one finite gnet coefficient per covariate.")
  }
  out
}

make_dat <- function(
    seed = 3036,
    n = 200,
    p = 500,
    nzc = 5,
    rho = 0.5,
    snr = NULL,
    common_bad_frac = 0,
    common_bad_count = NULL,
    K = 10,
    n_source_raw = 200,
    n_test = 1000,
    paper_setting = "I",
    Ka = 1,
    h = 5,
    signal_strength = 0.5,
    sigma_fixed = 1,
    noise_mode = "fixed_sigma",
    source_summary_method = c("spikeslab", "oracle")) {
  source_summary_method <- match.arg(source_summary_method)
  seed <- as.integer(seed); n <- as.integer(n); p <- as.integer(p)
  nzc <- as.integer(nzc); K <- as.integer(K); Ka <- as.integer(Ka)
  n_test <- as.integer(n_test)
  if (!paper_setting %in% c("I", "III")) stop("paper_setting must be I or III.")
  expected <- if (paper_setting == "I") c(n = 200, p = 500, s = 5, ns = 200, h = 5) else
    c(n = 100, p = 2000, s = 20, ns = 500, h = 10)
  if (!all(c(n, p, nzc) == unname(expected[c("n", "p", "s")]))) {
    stop("n, p and nzc do not match final-paper Setting ", paper_setting, ".")
  }
  if (K != 10L || Ka < 1L || Ka > K) stop("This DGP requires K=10 and Ka in 1:10.")
  n_source_raw <- rep(as.integer(n_source_raw), length.out = K)
  if (any(n_source_raw != as.integer(expected["ns"]))) {
    stop("Source sample size does not match final-paper Setting ", paper_setting, ".")
  }
  if (!isTRUE(all.equal(as.numeric(h), as.numeric(expected["h"])))) {
    stop("h does not match final-paper Setting ", paper_setting, ".")
  }
  if (!identical(noise_mode, "fixed_sigma") || sigma_fixed <= 0) {
    stop("The paper DGP uses fixed Gaussian noise with positive sigma_fixed.")
  }
  if (abs(rho - 0.5) > 1e-12 || abs(signal_strength - 0.5) > 1e-12) {
    stop("The paper DGP requires target AR(1) rho=0.5 and signal strength 0.5.")
  }

  sig <- seq_len(nzc)
  beta_true <- numeric(p)
  beta_true[sig] <- signal_strength
  signal_var <- .paper_ar1_signal_variance(nzc, rho, signal_strength)
  realized_snr <- signal_var / sigma_fixed^2
  if (!is.null(snr) && is.finite(snr) &&
      abs(snr - realized_snr) > 1e-6 * max(1, realized_snr)) {
    stop("Grid SNR is inconsistent with the paper's fixed beta, AR(1) design and sigma=1.")
  }

  set.seed(seed)
  X <- .paper_ar1_design(n, p, rho)
  y <- as.numeric(X %*% beta_true + stats::rnorm(n, sd = sigma_fixed))
  X_test <- .paper_ar1_design(n_test, p, rho)
  y_test <- as.numeric(X_test %*% beta_true + stats::rnorm(n_test, sd = sigma_fixed))

  W <- matrix(0, nrow = p, ncol = K)
  X_source <- vector("list", K)
  y_source <- vector("list", K)
  source_support <- vector("list", K)
  for (k in seq_len(K)) {
    ## Source-specific streams preserve pairing across Ka: changing whether
    ## source k is informative cannot alter any other source's data.
    set.seed(seed + 100003L + 1000L * k)
    if (k <= Ka) {
      wk <- beta_true + (h / p) * sample(c(-1, 1), p, replace = TRUE)
      source_support[[k]] <- which(wk != 0)
    } else {
      false_pool <- seq.int(2L * nzc + 1L, p)
      strong_false <- sort(unique(c(nzc + seq_len(nzc), sample(false_pool, nzc))))
      wk <- (2 * h / p) * sample(c(-1, 1), p, replace = TRUE)
      wk[strong_false] <- wk[strong_false] + signal_strength
      source_support[[k]] <- strong_false
    }
    set.seed(seed + 200003L + 1000L * k)
    Xk <- matrix(stats::rt(n_source_raw[k] * p, df = 4), nrow = n_source_raw[k], ncol = p)
    yk <- as.numeric(Xk %*% wk + stats::rnorm(n_source_raw[k], sd = sigma_fixed))
    W[, k] <- wk
    X_source[[k]] <- Xk
    y_source[[k]] <- yk
  }

  B <- if (source_summary_method == "oracle") {
    W
  } else {
    vapply(seq_len(K), function(k) {
      set.seed(seed + 300003L + 1000L * k)
      .paper_spikeslab_summary(X_source[[k]], y_source[[k]], paper_setting, nzc)
    }, numeric(p))
  }
  if (is.vector(B)) B <- matrix(B, ncol = K)

  alpha_true <- c(rep(1 / Ka, Ka), rep(0, K - Ka))
  source_only_union <- if (Ka < K) {
    sort(unique(unlist(source_support[(Ka + 1L):K], use.names = FALSE)))
  } else integer(0)
  true_C0 <- beta_true != 0
  gamma_true <- as.numeric(true_C0)
  source_cor <- vapply(seq_len(K), function(k) {
    z <- suppressWarnings(stats::cor(W[, k], beta_true))
    if (is.finite(z)) z else 0
  }, numeric(1))
  source_scale <- sqrt(colSums(W^2) / max(sum(beta_true^2), 1e-12))
  population_center <- as.numeric(W %*% alpha_true)
  center_gap <- population_center - beta_true

  X_all <- do.call(rbind, c(list(X), X_source))
  y_all <- c(y, unlist(y_source, use.names = FALSE))
  scenario_id <- sprintf("setting02_TRADER_%s_Ka%02d", paper_setting, Ka)
  beta_min <- min(abs(beta_true[sig]))
  beta_info <- n * beta_min^2 / sigma_fixed^2
  beta_benchmark <- log(max(n, 2)) + log(max(p, 2))

  list(
    scenario_id = scenario_id,
    scenario_title = paste0("TRADER paper Setting ", paper_setting),
    scenario_description = paste0(
      "Lai et al. (2026) final-paper Setting ", paper_setting,
      ": K=10 sources, Ka=", Ka, " informative sources; target AR(1) Gaussian design, ",
      "source iid t4 designs, beta_j=0.5 on the first s coordinates, and sigma=1."
    ),
    paper_setting = paper_setting,
    Ka = Ka,
    h = h,
    dgp_variant = "lai_2026_trader_final_paper",
    common_bad_mode = "paper_informative_vs_uninformative_sources",
    independent_contamination_factors = FALSE,
    n = n, p = p, K = K,
    ns = mean(n_source_raw), n_vec = n_source_raw,
    X = X, y = y, X_test = X_test, y_test = y_test,
    beta_true = beta_true, tag_true = as.numeric(beta_true != 0),
    gamma_true = gamma_true, true_C0 = true_C0, alpha_true = alpha_true,
    B = B, W_source_true = W,
    X_source = X_source, y_source = y_source,
    Sigma_source = rep(list(NULL), K),
    X_all = X_all, y_all = y_all, n_vec_trans = c(n, n_source_raw),
    n_source_raw = n_source_raw,
    source_type = c(rep("informative", Ka), rep("uninformative", K - Ka)),
    source_summary_method = source_summary_method,
    source_support = source_support,
    signal_layout = "first_s_equal_positive",
    signal_profile = "TRADER_paper_constant_0.5",
    theory_signal_magnitude = signal_strength,
    common_bad_signal_magnitude = NA_real_,
    transferable_signal_magnitude = signal_strength,
    nzc = nzc, rho = rho, noise_mode = noise_mode,
    snr = realized_snr, snr_requested = realized_snr, realized_snr = realized_snr,
    sigma_fixed = sigma_fixed, sigma_true = sigma_fixed,
    beta_min = beta_min, beta_max = max(abs(beta_true[sig])),
    beta_min_information = beta_info,
    beta_min_log_benchmark = beta_benchmark,
    beta_min_margin_ratio = beta_info / beta_benchmark,
    SigX = NULL,
    common_bad = integer(0), common_bad_count = 0L, common_bad_frac = 0,
    common_bad_scale = NA_real_, common_bad_gap_population = NA_real_,
    regular_signal = sig, transferable_count = length(sig),
    source_weight_identification_eigenvalue = NA_real_,
    source_weight_identification_pass = NA,
    alpha_identification_boundary = FALSE,
    theory_advantage_region = Ka > 0L,
    exact_combination_max_error = max(abs(center_gap[sig])),
    source_scale = source_scale,
    source_target_population_correlation = source_cor,
    shared_spurious = source_only_union,
    shared_spurious_count = length(source_only_union),
    shared_spurious_gap_population = if (length(source_only_union)) max(abs(population_center[source_only_union])) else NA_real_,
    spurious_union = source_only_union,
    spur_list = source_support,
    spur2 = if (K >= 2L) source_support[[2L]] else integer(0),
    spur3 = if (K >= 3L) source_support[[3L]] else integer(0),
    common_spurious = integer(0), common_spurious_count = 0L,
    support_swap_spur_ratio = NA_real_, support_swap_size = NA_real_,
    A_trans_oracle = seq_len(Ka), alpha_bias_scale = NA_real_,
    XtX = crossprod(X), Xty = as.numeric(crossprod(X, y)),
    x_norm2 = colSums(X^2), yTy = sum(y^2)
  )
}
