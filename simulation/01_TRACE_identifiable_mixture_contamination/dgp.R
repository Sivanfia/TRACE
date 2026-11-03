# Resolve the bundled fitter relative to this file, including source(..., local=...).
.source_summary_dgp_file <- tail(Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile)), 1)[[1]]
source(file.path(dirname(normalizePath(.source_summary_dgp_file, mustWork = TRUE)),
                 "..", "_common", "source_summaries.R"), local = TRUE)
rm(.source_summary_dgp_file)


## =========================================================
## Setting 01: identifiable source mixture with coordinate contamination
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
                                                 z_prior = c("dimension_mild", "dimension_mild_dense", "high_dimensional", "truth", "adaptive", "anti_collapse", "n_inverse", "dimension_mild_dense"),
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
  } else if (z_prior == "dimension_mild_dense") {
    rho_z_mean <- min(
      0.10,
      max(1 / max(n, 2), 1.5 / sqrt(max(p, 2)))
    )
  }
  else if (z_prior == "truth") {
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
    } else if (z_prior == "dimension_mild_dense") {
      "min(0.10, max(1/n, 1.5/sqrt(p)))"
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

## Main DGP.
make_dat <- function(
    seed = 3036,
    n = 300,
    p = 600,
    nzc = 10,
    rho = 0.3,
    snr = 2.0,
    noise_mode = c("fixed_snr", "fixed_sigma"),
    sigma_fixed = 1,
    K = 3,
    ns_summary = NULL,
    n_source_raw = rep(500L, 3L),
    n_test = 1000,
    source_mode = "alpha_combination",
    source_scale = c(1.00, 0.70, 0.40),
    alpha_true = c(0.50, 0.30, 0.20),
    alpha_bias_scale = 0.25,
    
    ## New: allow either absolute common-bad count or fraction.
    ## If common_bad is supplied explicitly, it has the highest priority.
    ## Else if common_bad_count is supplied, use it.
    ## Else use common_bad_frac.
    common_bad_count = NULL,
    common_bad_frac = NULL,
    common_bad = NULL,
    
    common_bad_scale = c(0, -1, -2),
    n_spur_each = c(0,  0,  0),
    spur_size = c(0, 0, 0),
    source_noise_scale = 1.0,
    A_trans_oracle = integer(0),
    signal = NULL,
    signal_profile = c("theory_flat", "legacy_weak_tail"),
    theory_signal_magnitude = 1
) {
  noise_mode <- match.arg(noise_mode)
  n_source_raw <- validate_shared_source_sizes(n_source_raw, K, p, ns_summary, source_noise_scale)
  signal_profile <- match.arg(signal_profile)
  if (length(source_scale) != K) stop("source_scale must have length K.")
  if (length(alpha_true) != K) stop("alpha_true must have length K.")
  if (length(common_bad_scale) != K) stop("common_bad_scale must have length K.")
  if (length(n_spur_each) != K) stop("n_spur_each must have length K.")
  if (length(spur_size) != K) stop("spur_size must have length K.")
  if (length(n_source_raw) == 1) n_source_raw <- rep(n_source_raw, K)
  if (length(n_source_raw) != K) stop("n_source_raw must have length 1 or K.")
  
  set.seed(seed)
  
  alpha_true <- as.numeric(alpha_true)
  if (sum(alpha_true) <= 0) stop("alpha_true must have positive sum.")
  alpha_true <- alpha_true / sum(alpha_true)
  
  beta_true <- .make_strong_block_beta(
    p = p,
    nzc = nzc,
    signal = signal,
    signal_profile = signal_profile,
    theory_signal_magnitude = theory_signal_magnitude
  )
  tag_true <- as.integer(beta_true != 0)
  
  sig <- which(tag_true == 1)
  nul <- which(tag_true == 0)
  
  ## ------------------------------------------------------------
  ## Generate common-bad coordinates.
  ##
  ## Priority:
  ## 1. If common_bad is explicitly supplied, use it.
  ## 2. Else if common_bad_count is supplied, use the first
  ##    common_bad_count true signals.
  ## 3. Else use common_bad_frac.
  ##
  ## This makes the DGP compatible with both:
  ## - previous fraction-based design
  ## - new fixed-count design
  ## ------------------------------------------------------------
  if (!is.null(common_bad)) {
    common_bad <- sort(unique(as.integer(common_bad)))
    common_bad <- common_bad[common_bad %in% sig]
    n_common_bad <- length(common_bad)
    
  } else {
    if (!is.null(common_bad_count)) {
      n_common_bad <- as.integer(round(common_bad_count))
      n_common_bad <- max(0L, min(n_common_bad, length(sig)))
      
    } else {
      if (is.null(common_bad_frac)) {
        common_bad_frac <- 0
      }
      
      ## Allow both 0.4 and 40 to mean 40%.
      common_bad_frac <- as.numeric(common_bad_frac)
      if (common_bad_frac > 1) {
        common_bad_frac <- common_bad_frac / 100
      }
      common_bad_frac <- min(max(common_bad_frac, 0), 1)
      
      n_common_bad <- as.integer(round(length(sig) * common_bad_frac))
      n_common_bad <- max(0L, min(n_common_bad, length(sig)))
    }
    
    common_bad <- if (n_common_bad > 0) {
      sig[seq_len(n_common_bad)]
    } else {
      integer(0)
    }
  }
  
  common_bad <- sort(unique(common_bad))
  common_bad <- common_bad[common_bad %in% sig]
  n_common_bad <- length(common_bad)
  
  common_bad_count <- n_common_bad
  common_bad_frac <- if (length(sig) > 0) n_common_bad / length(sig) else 0
  
  regular_signal <- setdiff(sig, common_bad)
  
  SigX <- .make_ar1_sigma(p, rho)
  X <- MASS::mvrnorm(n, rep(0, p), SigX)
  
  signal_var <- as.numeric(t(beta_true) %*% SigX %*% beta_true)
  snr_requested <- as.numeric(snr)[1]
  if (identical(noise_mode, "fixed_sigma")) {
    sigma_true <- as.numeric(sigma_fixed)[1]
    if (!is.finite(sigma_true) || sigma_true <= 0) {
      stop("sigma_fixed must be a finite positive number in fixed_sigma mode.")
    }
  } else {
    if (!is.finite(snr_requested) || snr_requested <= 0) {
      stop("snr must be a finite positive number in fixed_snr mode.")
    }
    sigma_true <- sqrt(signal_var / snr_requested)
  }
  realized_snr <- signal_var / sigma_true^2
  y <- as.numeric(X %*% beta_true + rnorm(n, 0, sigma_true))
  
  X_test <- MASS::mvrnorm(n_test, rep(0, p), SigX)
  y_test <- as.numeric(X_test %*% beta_true + rnorm(n_test, 0, sigma_true))
  
  W <- matrix(0, p, K)
  
  if (source_mode == "single_best") {
    for (k in seq_len(K)) {
      W[regular_signal, k] <- source_scale[k] * beta_true[regular_signal]
    }
    
  } else if (source_mode == "alpha_combination") {
    if (K != 3) stop("alpha_combination mode is currently written for K = 3.")
    
    a <- alpha_true
    if (abs(a[3]) < 1e-12) {
      stop("alpha_combination mode requires alpha_true[3] to be nonzero.")
    }
    
    for (ii in seq_along(regular_signal)) {
      j <- regular_signal[ii]
      bj <- beta_true[j]
      
      if (ii %% 2 == 1) {
        d1 <- alpha_bias_scale * bj
        d2 <- -0.15 * bj
      } else {
        d1 <- -0.20 * bj
        d2 <- 0.30 * bj
      }
      
      d3 <- -(a[1] * d1 + a[2] * d2) / a[3]
      W[j, ] <- bj + c(d1, d2, d3)
    }
    
  } else {
    stop("Unknown source_mode: ", source_mode)
  }
  
  ## For common-bad coordinates, all sources are misleading.
  if (length(common_bad) > 0) {
    for (k in seq_len(K)) {
      W[common_bad, k] <- common_bad_scale[k] * beta_true[common_bad]
    }
  }

  population_source_center <- as.numeric(W %*% alpha_true)
  common_bad_gap_population <- if (length(common_bad) > 0) {
    min(abs(population_source_center[common_bad] - beta_true[common_bad]))
  } else {
    NA_real_
  }
  ## Diagnostic for Proposition 1, equations (24)--(25).  This does not alter
  ## the DGP; it records whether transferable source rows identify every
  ## simplex-tangent direction.  In the original single_best Setting-01 the
  ## rows are collinear, so beta can satisfy A10/A11 even when this auxiliary
  ## alpha-identification condition does not.
  source_weight_identification_eigenvalue <- if (K <= 1L) {
    Inf
  } else if (length(regular_signal) == 0L) {
    0
  } else {
    tangent_basis <- qr.Q(qr(stats::contr.helmert(K)))
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
  source_weight_identification_pass <-
    is.infinite(source_weight_identification_eigenvalue) ||
    source_weight_identification_eigenvalue > 1e-8
  
  ## Add source-specific spurious null signals.
  available_null <- nul
  spur_list <- vector("list", K)
  
  for (k in seq_len(K)) {
    nk <- min(n_spur_each[k], length(available_null))
    
    if (nk > 0 && spur_size[k] > 0) {
      spk <- sort(sample(available_null, nk))
      available_null <- setdiff(available_null, spk)
      W[spk, k] <- spur_size[k] * sample(c(-1, 1), nk, replace = TRUE)
      spur_list[[k]] <- spk
    } else {
      spur_list[[k]] <- integer(0)
    }
  }
  
  spurious_union <- sort(unique(unlist(spur_list)))
  
  ## gamma_true = 0 for non-transferable or spurious coordinates.
  gamma_true <- rep(1, p)
  if (length(common_bad) > 0) gamma_true[common_bad] <- 0
  if (length(spurious_union) > 0) gamma_true[spurious_union] <- 0
  
  ## Generate one source dataset per source, shared by all transfer methods.
  X_source <- vector("list", K)
  y_source <- vector("list", K)
  Sigma_source <- vector("list", K)
  
  ## Likewise pair the raw source designs across contamination fractions.
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
    "setting01_identifiable_mixture_nzc", nzc,
    "_rho", gsub("\\.", "p", as.character(rho)),
    if (identical(noise_mode, "fixed_sigma")) {
      paste0("_sigma", gsub("\\.", "p", as.character(sigma_true)))
    } else {
      paste0("_snr", gsub("\\.", "p", as.character(snr_requested)))
    },
    "_cbcount", common_bad_count,
    "_cbfrac", round(100 * common_bad_frac)
  )

  active_abs <- abs(beta_true[sig])
  beta_min <- min(active_abs)
  beta_max <- max(active_abs)
  beta_min_information <- n * beta_min^2 / sigma_true^2
  beta_min_log_benchmark <- log(max(n, 2)) + log(max(p, 2))
  beta_min_margin_ratio <- beta_min_information / beta_min_log_benchmark
  
  list(
    scenario_id = scenario_id,
    scenario_title = "Setting 01: identifiable source mixture with active-coordinate contamination",
    scenario_description = paste0(
      "No source is individually oracle. On transferable coordinates, ",
      "B_j alpha0 equals beta0_j and the alternating source rows identify ",
      "the simplex-tangent directions. On common-bad coordinates the ",
      "weighted source center has a fixed non-vanishing gap. Signal profile=",
      signal_profile, ", noise mode=", noise_mode,
      ". DGP is independent of the Borrow-CVB initialization scheme."
    ),
    dgp_variant = "identifiable_mixture_active_contamination",
    common_bad_mode = "weighted_center_shift",
    
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
    signal_layout = if (identical(signal_profile, "theory_flat")) {
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
    
    source_scale = source_scale,
    
    n_spur_each = n_spur_each,
    spur_size = spur_size,
    spur_list = spur_list,
    spurious_union = spurious_union,
    spur2 = if (K >= 2) spur_list[[2]] else integer(0),
    spur3 = if (K >= 3) spur_list[[3]] else integer(0),
    
    A_trans_oracle = A_trans_oracle,
    alpha_bias_scale = alpha_bias_scale,
    
    XtX = crossprod(X),
    Xty = as.numeric(crossprod(X, y)),
    x_norm2 = colSums(X^2),
    yTy = sum(y^2)
  )
}
