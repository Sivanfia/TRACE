# All methods share actual source observations; default n_s=500 and original-scale spikeslab summaries.
# Current DGP and algorithm settings, preserved from the fixed-five-sweep version.
# Edit this file before starting a NEW result root.
SETTING_ID <- "01"
N_FIXED <- 300L
P_FIXED <- 600L
SOURCE_N <- sim_env_int("SIM_SOURCE_N", 500L)
ALLOWED_NZC <- c(10L, 15L, 30L)
TRADER_N_MCMC <- sim_env_int("SIM_TRADER_N_MCMC", 10000L)
TRADER_N_BURNIN <- sim_env_int("SIM_TRADER_N_BURNIN", as.integer(TRADER_N_MCMC/2))
if(TRADER_N_BURNIN>=TRADER_N_MCMC)stop("TRADER burn-in must be below chain length")
BASE_SEED <- sim_env_int("SIM_BASE_SEED", 3036L)
scenario_grid <- expand.grid(
  nzc = c(10, 15, 30),
  rho = c(0.3, 0.6),
  snr = c(4, 3),
  ## q=0.8 and q=1 are extreme stress tests, not alpha-identification cases.
  common_bad_frac = c(0, 0.2, 0.4, 0.6),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

scenario_grid$scenario_index <- seq_len(nrow(scenario_grid))
scenario_grid$common_bad_count <- as.integer(round(scenario_grid$nzc * scenario_grid$common_bad_frac))
scenario_grid$paired_seed_index <- as.integer(interaction(
  scenario_grid$nzc, scenario_grid$rho, scenario_grid$snr,
  drop = TRUE, lex.order = TRUE
))

## ---------------------------------------------------------
## Fixed DGP options shared by all scenarios.
## You may edit source_scale / contamination strength here.
## ---------------------------------------------------------
DGP_BASE_ARGS <- list(
  n = N_FIXED,
  p = P_FIXED,
  K = 3,
  n_source_raw = rep(SOURCE_N, 3L),
  n_test = 1000,
  source_mode = "alpha_combination",
  source_scale = c(1.00, 0.70, 0.40),
  alpha_true = c(0.50, 0.30, 0.20),
  alpha_bias_scale = 0.25,
  ## alpha0-weighted bad center = -0.70 beta0, hence the population gap is
  ## 1.70 |beta0| and no source preserves the old correct-sign shortcut.
  common_bad_scale = c(0, -1, -2),
  n_spur_each = c(0, 0, 0),
  spur_size = c(0, 0, 0),
  A_trans_oracle = integer(0),
  ## Main fixed-SNR DGP: all active beta_j have equal magnitude and alternating
  ## signs. This removes the 0.50--0.22 weak-tail confounding while preserving
  ## the fixed non-vanishing beta-min condition across the grid.
  signal_profile = "theory_flat",
  theory_signal_magnitude = 1,
  noise_mode = "fixed_snr"
)

## ---------------------------------------------------------
## Algorithm options shared by all scenarios.
## These options are intentionally diagnostic, not a large benchmark.
## ---------------------------------------------------------
ALGO_ARGS <- list(
  ## Finite-sample dimension-only inclusion prior; rho_c remains 0.5.
  z_prior = "dimension_mild_dense",
  gamma_prior = "constant_half",
  ## S1bar_inv = ns * diag(X'X) / (n * log(p)); raw ns is kept for Huber tau.
  source_precision_mode = "high_dimensional_logp",
  source_precision_cap_ratio = 1.0,
  ## Draft Algorithm 1: omega^(0)=1, alpha uniform, gamma^(0)=0.
  init_method_borrow = "full",
  ## Change the initial value of alpha here.
  ## Allowed values follow make_init_alpha(): "uniform", "source1", "random", "oracle".
  init_alpha_borrow = "uniform",

  ## Run a fixed number of target sweeps with gamma=0, then initialize alpha/gamma.
  init_gamma_borrow = "theory_zero",
  warm_start_learned_from_fixed = FALSE,
  theory_initialization = TRUE,
  target_init_until_converged = sim_env_bool("SIM_TARGET_INIT_UNTIL_CONVERGED", FALSE),
  target_init_steps = sim_env_int("SIM_TARGET_INIT_STEPS", 5L),
  target_init_max_iter = 1000L,
  target_init_tol = 1e-6,
  preliminary_alpha_update = "huber",
  huber_tau = NULL,
  huber_tau_rule = "fixed",
  huber_tau_multiplier = 1,

  ## Full-p updates avoid pre-initialization source-dependent screening.
  use_candidate_screening = FALSE,
  candidate_size = 220,
  ## TRUE uses the exact n-dimensional Woodbury backend for paper_block.
  fast_cached_update = TRUE,
  ## Paper block order: update all omega from (gamma^t, alpha^t), then all
  ## gamma from (omega^(t+1), gamma^t, alpha^t), then update alpha.
  borrow_update_schedule = "paper_block",
  warm_start_iter = 20,

  ## Shared target slab scale for TRACE and ACVB.  The fixed value 0.1 is
  ## allowed by Assumption A5 and is calibrated to the O(1) target signals.
  ## The previous value 100 made the slab excessively diffuse and added a
  ## large finite-sample inclusion penalty, especially when fixed total SNR
  ## was spread over nzc=30 coordinates.
  sigma_beta_sq_borrow = 0.1,
  max_iter_borrow = 500,
  tol_borrow = 1e-6,
  ## ELBO is still evaluated and written at every iteration, but it is not a
  ## stopping condition.  TRACE stops from parameter changes only.
  tol_elbo = 1e-7,
  stop_by_elbo = FALSE,
  min_iter = 8,

  ## Two-cycle rescue: eta=0.5 after a detected cycle; raw-map convergence check.
  borrow_adaptive_damping = TRUE,
  borrow_update_damping = 0.5,
  borrow_cycle_detection_tol = 1e-6,
  borrow_cycle_detection_relative_tol = 1e-3,
  borrow_cycle_min_amplitude = NULL,

  ## From the regular cycle onward TRACE uses the ordinary quadratic EM step.
  alpha_update_every_learned = 1,
  alpha_update_strategy_learned = "em_qp",
  alpha_update_strategy_old = "elbo_accept",
  run_old_alpha = FALSE,
  run_em_alpha = TRUE,
  alpha_freeze_iter_learned = 0,
  alpha_grid_step_learned = 0.02,
  alpha_accept_tol = 1e-8,
  run_current_alpha_every5 = FALSE,

  ## Fixed uniform alpha is off by default in this comparison, but can be turned on.
  run_fixed_alpha = FALSE,
  run_oracle_alpha = TRUE,

  ## Baselines. ACVB is target-only VB. Trans-Lasso uses the pinned author
  ## Q-aggregation implementation; errors never trigger another estimator.
  run_acvb = TRUE,
  ## Exact Woodbury/Sherman-Morrison backend: same ACVB coordinate update,
  ## but solves n x n systems and vectorizes the p coordinate comparisons.
  acvb_acceleration = "woodbury",
  run_trans_lasso = TRUE,
  run_oracle_trans_simple = FALSE,
  run_naive_trans_simple = FALSE,

  run_alpha_grid_profile = FALSE,
  alpha_grid_step = 0.02,

  source_prec_ratio = 0.25,
  source_prec_cap = 0.75,
  refit_ridge = 1e-6,
  gamma_ebic = 0.5,

  monitor_elbo = TRUE,
  elbo_print_every = 10,

  ## Include the normalized prior determinant in the collapsed surrogate.
  include_source_normalization = TRUE,

  ## Performance diagnostic for the formula log p(y | z, c, alpha).
  ## "selected" means z_j=1(omega_j>0.5) and c_j=1(gamma_j>0.5),
  ## matching the displayed formula with Z and C as indicator matrices.
  collapsed_logp_mode = "selected",
  collapsed_logp_threshold = 0.5,

  ## TRADER baseline.  Using source_input = "B" gives TRADER exactly the
  ## same external summary coefficient matrix used by TRACE. Change to
  ## "raw" for a sensitivity analysis that re-estimates every source by
  ## spikeslab from X_source/y_source, as in the original TRADER scripts.
  run_trader = TRUE,
  trader_source_input = "B",
  trader_estimator = "spikeslab",
  trader_allow_estimator_fallback = FALSE,
  trader_corr_threshold = 0.10,
  ## NULL invokes the dimension-only floor(n/log(p)) rule.  It never reads
  ## beta_true or nzc.
  trader_max_var = NULL,
  trader_n_mcmc = TRADER_N_MCMC,
  trader_n_burnin = TRADER_N_BURNIN,
  trader_thin = 2,
  trader_n_chain = 1,
  trader_sign_threshold = 0.95,
  trader_keep_samples = FALSE,
  trader_verbose = TRUE,

  ## CONCERT baseline: official author C++ implementation using the full
  ## target/source X and y.  The q_k rule reproduces the author's simulation
  ## scripts: minimum ridge-coefficient agreement over the K sources.
  run_concert = TRUE,
  concert_qk_rule = "author_ridge",
  concert_qk_fixed = 0.5,
  concert_qk_tolerance = 0.008,
  concert_probability_epsilon = 1e-6,
  concert_q0 = 0.1,
  concert_tau = 1,
  concert_eta = 1,
  concert_threshold = 1e-6,
  concert_max_iter = 1000L,
  concert_selection_threshold = 0.5,

  verbose = TRUE
)
