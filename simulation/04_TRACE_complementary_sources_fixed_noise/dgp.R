# Resolve the bundled fitter relative to this file, including source(..., local=...).
.source_summary_dgp_file <- tail(Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile)), 1)[[1]]
source(file.path(dirname(normalizePath(.source_summary_dgp_file, mustWork = TRUE)),
                 "..", "_common", "source_summaries.R"), local = TRUE)
rm(.source_summary_dgp_file)


## =========================================================
## Setting 04: identifiable complementary sources with local contamination
## DGP file
## =========================================================

suppressPackageStartupMessages({
  library(MASS)
})

.safe_prob <- function(x, low = 1e-6, high = 1 - 1e-6) {
  pmin(pmax(x, low), high)
}

.make_ar1_sigma <- function(p, rho) {
  outer(seq_len(p), seq_len(p), function(i, j) rho^abs(i - j))
}

.make_strong_block_beta <- function(p,
                                    nzc = 8,
                                    signal = NULL,
                                    signal_profile = c("theory_flat", "legacy_weak_tail"),
                                    theory_signal_magnitude = 1) {
  signal_profile <- match.arg(signal_profile)
  beta <- rep(0, p)
  if (is.null(signal)) {
    if (identical(signal_profile, "theory_flat")) {
      ## Theory-aligned Setting-01 signal sequence.  Equal magnitudes make
      ## beta_min bounded away from zero uniformly over nzc and prevent q
      ## from changing the signal-strength composition of the common-bad and
      ## transferable sets.  Under the main fixed-SNR design its absolute
      ## magnitude is only a scale convention; the equal-magnitude profile,
      ## rather than a global rescaling, improves the standardized beta-min
      ## margin relative to the legacy 0.50--0.22 weak tail.
      theory_signal_magnitude <- as.numeric(theory_signal_magnitude)[1]
      if (!is.finite(theory_signal_magnitude) || theory_signal_magnitude <= 0) {
        stop("theory_signal_magnitude must be a finite positive number.")
      }
      signal <- theory_signal_magnitude *
        rep(c(1, -1), length.out = nzc)
    } else {
      ## Original Setting-01 stress sequence.  Keep it available as a named
      ## sensitivity design, but do not use it for theory-rate verification:
      ## at nzc=30 its 0.50--0.22 tail confounds contamination with beta-min.
      base <- c(1.50, -1.30, 1.15, -1.05,
                0.95, -0.90, 0.85, -0.80,
                0.70, -0.65, 0.60, -0.55)
      if (nzc <= length(base)) {
        signal <- base[seq_len(nzc)]
      } else {
        tail_mag <- seq(0.50, 0.22, length.out = nzc - length(base))
        tail <- tail_mag * rep(c(1, -1), length.out = length(tail_mag))
        signal <- c(base, tail)
      }
    }
  }
  if (length(signal) < nzc) stop("signal length must be at least nzc.")
  beta[seq_len(nzc)] <- signal[seq_len(nzc)]
  beta
}

.make_z_prior_params_alpha_commonbad <- function(n, p, nzc = 8,
                                                 z_prior = c("dimension_mild", "dimension_mild_dense", "high_dimensional", "truth", "adaptive", "anti_collapse", "n_inverse"),
                                                 min_A_z = 2,
                                                 min_conc_z = 50) {
  z_prior <- match.arg(z_prior)
  truth_sparsity <- nzc / p

  if (z_prior == "high_dimensional") {
    rho_z_mean <- 1 / (p * max(p^(3 / 2), sqrt(n)))
  } else if (z_prior == "dimension_mild") {
    ## A deterministic, non-oracle prior shared by TRACE and ACVB.  Its
    ## expected model size is sqrt(p), capped at 10% of p.  For the Setting-01
    ## design p=600 this gives rho_z=0.0408: materially milder than 1/n while
    ## remaining independent of the unknown number of nonzero coefficients.
    rho_z_mean <- min(0.10, max(1 / max(n, 2), 1 / sqrt(max(p, 2))))
  }else if (z_prior == "dimension_mild_dense") {
    rho_z_mean <- min(
      0.10,
      max(1 / max(n, 2), 1.5 / sqrt(max(p, 2)))
    )
  } else if (z_prior == "truth") {
    rho_z_mean <- truth_sparsity
  } else if (z_prior == "n_inverse") {
    rho_z_mean <- 1 / max(n, 2)
  } else if (z_prior == "anti_collapse") {
    ## Match Setting 03: raise the prior inclusion mass enough to keep weak
    ## signals alive during the first target-only sweep, without allowing an
    ## effectively dense model.
    rho_z_mean <- min(0.30, max(4.0 * truth_sparsity, 1 / n))
  } else {
    if (n > p) {
      rho_z_mean <- min(0.20, max(truth_sparsity, 1 / n))
    } else {
      rho_z_mean <- min(0.12, max(2.0 * truth_sparsity, 1 / n))
    }
  }

  if (z_prior != "high_dimensional") rho_z_mean <- .safe_prob(rho_z_mean)
  ## Match Setting 03: use a less concentrated beta prior so the likelihood
  ## can correct the anti-collapse initialization after the early iterations.
  conc_floor <- if (z_prior == "anti_collapse") 20 else min_conc_z
  conc_z <- max(conc_floor, min_A_z / rho_z_mean)
  A_z <- conc_z * rho_z_mean
  B_z <- conc_z * (1 - rho_z_mean)

  list(
    A_z = A_z,
    B_z = B_z,
    rho_z_mean = rho_z_mean,
    conc_z = conc_z,
    qlogis_rho_z_mean = qlogis(rho_z_mean),
    E_logit_rho_z = qlogis(rho_z_mean),
    z_prior = z_prior,
    z_prior_rule = if (z_prior == "high_dimensional") {
      "1/(p*max(p^(3/2),sqrt(n)))"
    } else if (z_prior == "dimension_mild") {
      "min(0.10, max(1/n, 1/sqrt(p)))"
    } else {
      z_prior
    },
    z_prior_uses_true_nzc = z_prior %in% c("truth", "adaptive", "anti_collapse")
  )
}

make_prior_params_borrow_alpha_commonbad <- function(n, p, nzc = 8,
                                                     z_prior = c("dimension_mild", "dimension_mild_dense", "high_dimensional", "truth", "adaptive", "anti_collapse", "n_inverse"),
                                                     gamma_prior = c("constant_half", "ns_inverse_quarter", "n_inverse", "moderate", "strong", "weak"),
                                                     ns_eff = NULL,
                                                     min_A_z = 2,
                                                     min_conc_z = 50) {
  z_prior <- match.arg(z_prior)
  gamma_prior <- match.arg(gamma_prior)

  zpar <- .make_z_prior_params_alpha_commonbad(
    n = n,
    p = p,
    nzc = nzc,
    z_prior = z_prior,
    min_A_z = min_A_z,
    min_conc_z = min_conc_z
  )

  if (gamma_prior == "constant_half") {
    rho_c_mean <- 0.5
    conc_c <- 50
  } else if (gamma_prior == "ns_inverse_quarter") {
    ## Prior mean for c_j = 1.  Here n^s is the effective source-summary
    ## sample size.  This is rho_c = (n^s)^(-1/4), not an initialization.
    if (is.null(ns_eff)) ns_eff <- n
    rho_c_mean <- as.numeric(ns_eff)^(-1 / 4)
    rho_c_mean <- .safe_prob(rho_c_mean)
    ## Keep the beta prior numerically stable and avoid A_c being too tiny.
    conc_c <- max(50, 2 / rho_c_mean, 2 / (1 - rho_c_mean))
  } else if (gamma_prior == "n_inverse") {
    rho_c_mean <- 1 - 1 / max(n, 2)
    rho_c_mean <- .safe_prob(rho_c_mean)
    conc_c <- max(50, 2 / (1 - rho_c_mean))
  } else if (gamma_prior == "weak") {
    rho_c_mean <- 0.80
    conc_c <- 8
  } else if (gamma_prior == "moderate") {
    rho_c_mean <- 0.88
    conc_c <- 10
  } else {
    rho_c_mean <- 0.95
    conc_c <- 10
  }

  A_c <- conc_c * rho_c_mean
  B_c <- conc_c * (1 - rho_c_mean)

  c(
    zpar,
    list(
      A_c = A_c,
      B_c = B_c,
      rho_c_mean = rho_c_mean,
      conc_c = conc_c,
      E_logit_rho_c = qlogis(rho_c_mean),
      gamma_prior = gamma_prior
    )
  )
}

make_prior_params_vbss_alpha_commonbad <- function(n, p, nzc = 8,
                                                   z_prior = c("dimension_mild", "dimension_mild_dense", "high_dimensional", "truth", "adaptive", "anti_collapse", "n_inverse"),
                                                   min_A_z = 2,
                                                   min_conc_z = 50) {
  z_prior <- match.arg(z_prior)
  .make_z_prior_params_alpha_commonbad(
    n = n,
    p = p,
    nzc = nzc,
    z_prior = z_prior,
    min_A_z = min_A_z,
    min_conc_z = min_conc_z
  )
}

## Backward-compatible names used by the runner.
make_prior_params_borrow_commonbad <- make_prior_params_borrow_alpha_commonbad
make_prior_params_vbss_commonbad <- make_prior_params_vbss_alpha_commonbad

## Active Setting 04 DGP. The configuration determines every simulation cell.
make_dat <- function(
    seed = 3036,
    n = 300,
    p = 600,
    nzc = 30,
    rho = 0.3,
    snr = 2.25,
    noise_mode = c("fixed_sigma", "fixed_snr"),
    sigma_fixed = 3,
    K = 3,
    ns_summary = NULL,
    n_source_raw = rep(500L, 3L),
    n_test = 1000,
    source_mode = "exact_combination_cycle",
    source_multiplier_patterns = rbind(
      c(3, -1, -1),
      c(-0.2, 3.8, -0.2),
      c(0.2, 0.2, 4.2)
    ),
    alpha_true = c(0.50, 0.30, 0.20),
    combination_tolerance = 1e-10,
    common_bad_count = NULL,
    common_bad_frac = NULL,
    common_bad = NULL,
    common_bad_scale = c(-1, -1, -1),
    shared_spur_count = nzc,
    shared_spur_size = 1.0,
    shared_spur_scale = c(1, 1, 1),
    source_noise_scale = 1.0,
    A_trans_oracle = integer(0),
    signal = NULL,
    signal_profile = c("theory_flat", "structured_contrast", "legacy_weak_tail"),
    theory_signal_magnitude = 1,
    common_bad_signal_magnitude = 1.2,
    transferable_signal_magnitude = 0.9
) {
  noise_mode <- match.arg(noise_mode)
  n_source_raw <- validate_shared_source_sizes(n_source_raw, K, p, ns_summary, source_noise_scale)
  signal_profile <- match.arg(signal_profile)
  if (!identical(source_mode, "exact_combination_cycle")) {
    stop("Setting 04 now requires source_mode='exact_combination_cycle'.")
  }
  if (K != 3L) stop("Setting 04 exact-combination patterns require K = 3.")
  if (length(alpha_true) != K) stop("alpha_true must have length K.")
  if (length(common_bad_scale) != K) stop("common_bad_scale must have length K.")
  if (length(shared_spur_scale) != K) stop("shared_spur_scale must have length K.")
  if (length(n_source_raw) == 1L) n_source_raw <- rep(n_source_raw, K)
  if (length(n_source_raw) != K) stop("n_source_raw must have length 1 or K.")
  common_bad_signal_magnitude <- as.numeric(common_bad_signal_magnitude)[1]
  transferable_signal_magnitude <- as.numeric(transferable_signal_magnitude)[1]
  if (!is.finite(common_bad_signal_magnitude) || common_bad_signal_magnitude <= 0 ||
      !is.finite(transferable_signal_magnitude) || transferable_signal_magnitude <= 0) {
    stop("Structured signal magnitudes must be finite and positive.")
  }

  source_multiplier_patterns <- as.matrix(source_multiplier_patterns)
  storage.mode(source_multiplier_patterns) <- "double"
  if (ncol(source_multiplier_patterns) != K || nrow(source_multiplier_patterns) < 2L) {
    stop("source_multiplier_patterns must have K columns and at least two rows.")
  }
  if (any(!is.finite(source_multiplier_patterns))) {
    stop("source_multiplier_patterns must be finite.")
  }
  combination_tolerance <- as.numeric(combination_tolerance)[1]
  if (!is.finite(combination_tolerance) || combination_tolerance <= 0) {
    stop("combination_tolerance must be finite and positive.")
  }
  shared_spur_count <- as.integer(round(as.numeric(shared_spur_count)[1]))
  shared_spur_size <- as.numeric(shared_spur_size)[1]
  if (!is.finite(shared_spur_count) || shared_spur_count < 0L) {
    stop("shared_spur_count must be a non-negative integer.")
  }
  if (!is.finite(shared_spur_size) || shared_spur_size < 0) {
    stop("shared_spur_size must be finite and non-negative.")
  }

  alpha_true <- as.numeric(alpha_true)
  if (any(!is.finite(alpha_true)) || any(alpha_true < 0) || sum(alpha_true) <= 0) {
    stop("alpha_true must be a non-negative simplex vector.")
  }
  alpha_true <- alpha_true / sum(alpha_true)
  pattern_combination_values <- as.numeric(source_multiplier_patterns %*% alpha_true)
  pattern_combination_error <- max(abs(pattern_combination_values - 1))
  if (pattern_combination_error > combination_tolerance) {
    stop(
      "Every source multiplier pattern must satisfy alpha_true' pattern = 1. ",
      "Maximum error is ", signif(pattern_combination_error, 6), "."
    )
  }

  set.seed(seed)
  base_signal_profile <- if (identical(signal_profile, "structured_contrast")) {
    "theory_flat"
  } else {
    signal_profile
  }
  beta_true <- .make_strong_block_beta(
    p = p,
    nzc = nzc,
    signal = signal,
    signal_profile = base_signal_profile,
    theory_signal_magnitude = theory_signal_magnitude
  )
  tag_true <- as.integer(beta_true != 0)
  sig <- which(tag_true == 1L)
  nul <- which(tag_true == 0L)

  ## Common-bad coordinates are nested across fractions under paired seeds.
  if (!is.null(common_bad)) {
    common_bad <- sort(unique(as.integer(common_bad)))
    common_bad <- common_bad[common_bad %in% sig]
  } else {
    if (!is.null(common_bad_count)) {
      n_common_bad <- as.integer(round(as.numeric(common_bad_count)[1]))
    } else {
      common_bad_frac <- as.numeric(common_bad_frac %||% 0)[1]
      if (common_bad_frac > 1) common_bad_frac <- common_bad_frac / 100
      if (!is.finite(common_bad_frac) || common_bad_frac < 0 || common_bad_frac > 1) {
        stop("common_bad_frac must lie in [0,1].")
      }
      n_common_bad <- as.integer(round(length(sig) * common_bad_frac))
    }
    n_common_bad <- max(0L, min(n_common_bad, length(sig)))
    common_bad <- if (n_common_bad > 0L) sig[seq_len(n_common_bad)] else integer(0)
  }
  common_bad <- sort(unique(common_bad[common_bad %in% sig]))
  common_bad_count <- length(common_bad)
  common_bad_frac <- if (length(sig) > 0L) common_bad_count / length(sig) else 0
  regular_signal <- setdiff(sig, common_bad)

  ## Retained sensitivity option only.  The active Setting-04 runner uses
  ## theory_flat so beta_true and target data do not change with q.  This older
  ## structured profile remains available for reproducing development pilots.
  if (identical(signal_profile, "structured_contrast")) {
    if (!is.null(signal)) {
      stop("Do not supply signal when signal_profile='structured_contrast'.")
    }
    active_sign <- rep(c(1, -1), length.out = length(sig))
    names(active_sign) <- sig
    beta_true[common_bad] <-
      common_bad_signal_magnitude * active_sign[as.character(common_bad)]
    beta_true[regular_signal] <-
      transferable_signal_magnitude * active_sign[as.character(regular_signal)]
  }

  SigX <- .make_ar1_sigma(p, rho)
  X <- MASS::mvrnorm(n, rep(0, p), SigX)
  signal_var <- as.numeric(crossprod(beta_true, SigX %*% beta_true))
  snr_requested <- as.numeric(snr)[1]
  if (identical(noise_mode, "fixed_sigma")) {
    sigma_true <- as.numeric(sigma_fixed)[1]
    if (!is.finite(sigma_true) || sigma_true <= 0) {
      stop("sigma_fixed must be finite and positive in fixed_sigma mode.")
    }
  } else {
    if (!is.finite(snr_requested) || snr_requested <= 0) {
      stop("snr must be finite and positive in fixed_snr mode.")
    }
    sigma_true <- sqrt(signal_var / snr_requested)
  }
  realized_snr <- signal_var / sigma_true^2
  y <- as.numeric(X %*% beta_true + rnorm(n, 0, sigma_true))
  X_test <- MASS::mvrnorm(n_test, rep(0, p), SigX)
  y_test <- as.numeric(X_test %*% beta_true + rnorm(n_test, 0, sigma_true))

  ## Pattern membership is tied to the original active-coordinate position,
  ## rather than re-numbered after contamination.  This keeps the population
  ## pattern deterministic and paired across q, even when a finite transferable
  ## subset is not exactly balanced across the three patterns.
  W <- matrix(0, p, K)
  pattern_id <- rep(NA_integer_, p)
  if (length(regular_signal) > 0L) {
    active_position <- match(regular_signal, sig)
    pattern_id[regular_signal] <- ((active_position - 1L) %% nrow(source_multiplier_patterns)) + 1L
    for (ii in seq_along(regular_signal)) {
      j <- regular_signal[ii]
      W[j, ] <- beta_true[j] * source_multiplier_patterns[pattern_id[j], ]
    }
  }
  if (length(common_bad) > 0L) {
    for (k in seq_len(K)) {
      W[common_bad, k] <- common_bad_scale[k] * beta_true[common_bad]
    }
  }

  population_source_center <- as.numeric(W %*% alpha_true)
  exact_combination_max_error <- if (length(regular_signal) > 0L) {
    max(abs(population_source_center[regular_signal] - beta_true[regular_signal]))
  } else {
    0
  }
  if (exact_combination_max_error > combination_tolerance) {
    stop("Transferable rows failed the exact-combination identity.")
  }
  common_bad_gap_population <- if (length(common_bad) > 0L) {
    min(abs(population_source_center[common_bad] - beta_true[common_bad]))
  } else {
    NA_real_
  }

  ## Proposition-1 diagnostic: exact-combination rows vary in two independent
  ## simplex-tangent directions, so alpha is identified even though their
  ## alpha-weighted population center equals beta exactly.
  tangent_basis <- qr.Q(qr(stats::contr.helmert(K)))
  source_weight_identification_eigenvalue <- if (length(regular_signal) == 0L) {
    0
  } else {
    B_transfer <- W[regular_signal, , drop = FALSE]
    Qdiag_transfer <- diag(SigX)[regular_signal]
    identification_matrix <- crossprod(
      B_transfer,
      B_transfer * Qdiag_transfer
    ) / max(length(sig), 1L)
    projected_identification <- crossprod(
      tangent_basis,
      identification_matrix %*% tangent_basis
    )
    max(0, min(eigen(
      (projected_identification + t(projected_identification)) / 2,
      symmetric = TRUE,
      only.values = TRUE
    )$values))
  }
  source_weight_identification_pass <- source_weight_identification_eigenvalue > 1e-8
  alpha_identification_boundary <- length(regular_signal) == 0L
  theory_advantage_region <- !alpha_identification_boundary

  ## Independent omega challenge: a nested target-null block is nonzero in all
  ## sources.  Its size is not a function of common_bad_count.
  set.seed(as.integer(seed) + 300007L)
  null_order <- if (length(nul) > 0L) sample(nul, length(nul), replace = FALSE) else integer(0)
  shared_spur_count <- min(shared_spur_count, length(null_order))
  shared_spurious <- if (shared_spur_count > 0L) {
    sort(null_order[seq_len(shared_spur_count)])
  } else {
    integer(0)
  }
  shared_spurious_sign <- rep(c(1, -1), length.out = shared_spur_count)
  if (shared_spur_count > 0L && shared_spur_size > 0) {
    for (k in seq_len(K)) {
      W[shared_spurious, k] <-
        shared_spur_size * shared_spur_scale[k] * shared_spurious_sign
    }
  }
  population_source_center <- as.numeric(W %*% alpha_true)
  shared_spurious_gap_population <- if (shared_spur_count > 0L) {
    min(abs(population_source_center[shared_spurious]))
  } else {
    NA_real_
  }
  source_target_population_correlation <- vapply(seq_len(K), function(k) {
    source_var_k <- as.numeric(crossprod(W[, k], SigX %*% W[, k]))
    source_target_cov_k <- as.numeric(crossprod(beta_true, SigX %*% W[, k]))
    source_target_cov_k / sqrt(max(signal_var * source_var_k, 1e-14))
  }, numeric(1))

  spur_list <- rep(list(shared_spurious), K)
  spurious_union <- shared_spurious
  gamma_true <- rep(1, p)
  gamma_true[common_bad] <- 0
  gamma_true[spurious_union] <- 0

  X_source <- vector("list", K)
  y_source <- vector("list", K)
  Sigma_source <- vector("list", K)
  set.seed(as.integer(seed) + 200003L)
  for (k in seq_len(K)) {
    Xk <- MASS::mvrnorm(n_source_raw[k], rep(0, p), SigX)
    yk <- as.numeric(Xk %*% W[, k] + rnorm(n_source_raw[k], 0, sigma_true))
    X_source[[k]] <- Xk
    y_source[[k]] <- yk
    Sigma_source[[k]] <- SigX
  }
  source_fit <- fit_shared_source_summaries(X_source, y_source, n_source_raw,
    seed = as.integer(seed) + 610000L)
  B <- source_fit$B
  ns_vec <- source_fit$counts

  X_all <- do.call(rbind, c(list(X), X_source))
  y_all <- c(y, unlist(y_source, use.names = FALSE))
  n_vec_trans <- c(n, n_source_raw)

  scenario_id <- paste0(
    "setting04_exact_combination_nzc", nzc,
    "_rho", gsub("\\.", "p", as.character(rho)),
    if (identical(noise_mode, "fixed_sigma")) {
      paste0("_sigma", gsub("\\.", "p", as.character(sigma_true)))
    } else {
      paste0("_snr", gsub("\\.", "p", as.character(snr_requested)))
    },
    "_cbfrac", round(100 * common_bad_frac),
    "_sharedspur", shared_spur_count
  )

  active_abs <- abs(beta_true[sig])
  beta_min <- min(active_abs)
  beta_max <- max(active_abs)
  beta_min_information <- n * beta_min^2 / sigma_true^2
  beta_min_log_benchmark <- log(max(n, 2)) + log(max(p, 2))
  beta_min_margin_ratio <- beta_min_information / beta_min_log_benchmark
  source_scale <- colMeans(source_multiplier_patterns)

  list(
    scenario_id = scenario_id,
    scenario_title = "Setting 04: globally attractive coordinate-biased sources",
    scenario_description = paste0(
      "Target coefficients are fixed across contamination fractions. Common-bad source effects ",
      "use the configured common_bad_scale; transferable coordinates ",
      "cycle through exact-combination patterns; shared source-only nulls form an independent ",
      "coordinate-selection stress axis. ",
      "Signal profile=", signal_profile, ", noise mode=", noise_mode, "."
    ),
    dgp_variant = "positive_global_coordinate_bias_exact_combination",
    common_bad_mode = "shared_configured_source_scale",
    independent_contamination_factors = TRUE,
    n = n,
    p = p,
    K = K,
    ns = mean(ns_vec),
    n_vec = ns_vec,
    X = X,
    y = y,
    X_test = X_test,
    y_test = y_test,
    beta_true = beta_true,
    tag_true = tag_true,
    gamma_true = gamma_true,
    alpha_true = alpha_true,
    B = B,
    source_summary_method = "spikeslab_gnet_original_scale",
    source_summary_coefficient_scale = "original_predictors",
    source_summary_from_raw_samples = TRUE,
    source_summary_diagnostics = source_fit$diagnostics,
    source_summary_se = source_fit$standard_errors,
    source_summary_fit_time_sec = source_fit$fit_time_sec,
    W_source_true = W,
    X_source = X_source,
    y_source = y_source,
    Sigma_source = Sigma_source,
    X_all = X_all,
    y_all = y_all,
    n_vec_trans = n_vec_trans,
    n_source_raw = n_source_raw,
    source_type = source_mode,
    source_multiplier_patterns = source_multiplier_patterns,
    source_multiplier_pattern_id = pattern_id,
    pattern_combination_values = pattern_combination_values,
    exact_combination_max_error = exact_combination_max_error,
    signal_layout = if (identical(signal_profile, "structured_contrast")) {
      "strong_common_bad_weaker_transferable_alternating_sign"
    } else if (identical(signal_profile, "theory_flat")) {
      "block_equal_magnitude_alternating_sign"
    } else {
      "block_decay_legacy_weak_tail"
    },
    signal_profile = signal_profile,
    theory_signal_magnitude = if (identical(signal_profile, "theory_flat")) {
      theory_signal_magnitude
    } else {
      NA_real_
    },
    common_bad_signal_magnitude = if (identical(signal_profile, "structured_contrast")) {
      common_bad_signal_magnitude
    } else {
      NA_real_
    },
    transferable_signal_magnitude = if (identical(signal_profile, "structured_contrast")) {
      transferable_signal_magnitude
    } else {
      NA_real_
    },
    nzc = nzc,
    rho = rho,
    noise_mode = noise_mode,
    snr = realized_snr,
    snr_requested = snr_requested,
    realized_snr = realized_snr,
    sigma_fixed = if (identical(noise_mode, "fixed_sigma")) sigma_true else NA_real_,
    sigma_true = sigma_true,
    beta_min = beta_min,
    beta_max = beta_max,
    beta_min_information = beta_min_information,
    beta_min_log_benchmark = beta_min_log_benchmark,
    beta_min_margin_ratio = beta_min_margin_ratio,
    SigX = SigX,
    common_bad = common_bad,
    common_bad_count = common_bad_count,
    common_bad_frac = common_bad_frac,
    common_bad_scale = common_bad_scale,
    common_bad_gap_population = common_bad_gap_population,
    regular_signal = regular_signal,
    transferable_count = length(regular_signal),
    source_weight_identification_eigenvalue = source_weight_identification_eigenvalue,
    source_weight_identification_pass = source_weight_identification_pass,
    alpha_identification_boundary = alpha_identification_boundary,
    theory_advantage_region = theory_advantage_region,
    source_scale = source_scale,
    source_target_population_correlation = source_target_population_correlation,
    shared_spurious = shared_spurious,
    shared_spurious_count = shared_spur_count,
    shared_spur_size = shared_spur_size,
    shared_spur_scale = shared_spur_scale,
    shared_spurious_gap_population = shared_spurious_gap_population,
    n_spur_each = rep(shared_spur_count, K),
    spur_size = rep(shared_spur_size, K),
    spur_list = spur_list,
    spurious_union = spurious_union,
    spur2 = if (K >= 2L) spur_list[[2L]] else integer(0),
    spur3 = if (K >= 3L) spur_list[[3L]] else integer(0),
    ## Compatibility fields for existing scenario summaries.
    common_spurious = shared_spurious,
    common_spurious_count = shared_spur_count,
    support_swap_spur_ratio = NA_real_,
    support_swap_size = NA_real_,
    A_trans_oracle = A_trans_oracle,
    alpha_bias_scale = NA_real_,
    XtX = crossprod(X),
    Xty = as.numeric(crossprod(X, y)),
    x_norm2 = colSums(X^2),
    yTy = sum(y^2)
  )
}
