# Adapter for the official linear CONCERT implementation.
#
# Upstream source (vendored without algorithmic edits):
#   https://github.com/RuqianZhang/CONCERT
#   commit fea55703c39655a196d3aa5296d1851f0fd53547

.concert_bridge_dir <- local({
  ofiles <- Filter(function(x) is.character(x) && length(x) == 1L && nzchar(x),
                   lapply(sys.frames(), function(x) x$ofile))
  ofile <- if (length(ofiles)) tail(ofiles, 1L)[[1L]] else NULL
  if (is.null(ofile)) getwd() else dirname(normalizePath(ofile))
})

.concert_cpp_file <- file.path(
  .concert_bridge_dir,
  "CONCERT_author_upstream", "Linear", "Concert.cpp"
)

if (!requireNamespace("Rcpp", quietly = TRUE)) {
  stop("Package 'Rcpp' is required for CONCERT.")
}
if (!requireNamespace("RcppArmadillo", quietly = TRUE)) {
  stop("Package 'RcppArmadillo' is required for CONCERT.")
}
if (!file.exists(.concert_cpp_file)) {
  stop("Cannot find the vendored CONCERT C++ source: ", .concert_cpp_file)
}

# Concert.cpp calls this function from the global R environment for its
# lasso initialization.  This is byte-for-byte the R helper used upstream.
lasso_est_R <- function(X, y, lambda) {
  beta <- as.numeric(glmnet::glmnet(X, y, lambda = lambda)$beta)
  beta
}

# Compile the unmodified upstream C++ source. Rcpp's cache prevents needless
# recompilation in ordinary repeated runs; forked workers inherit the backend.
Rcpp::sourceCpp(.concert_cpp_file, rebuild = FALSE, verbose = FALSE)

.concert_upstream_commit <- "fea55703c39655a196d3aa5296d1851f0fd53547"
message("[CONCERT] official author linear backend ", .concert_upstream_commit)

.concert_cv_ridge <- function(X, y, seed) {
  set.seed(as.integer(seed))
  fit <- glmnet::cv.glmnet(
    x = as.matrix(X), y = as.numeric(y), family = "gaussian", alpha = 0
  )
  as.numeric(fit$glmnet.fit$beta[, fit$lambda == fit$lambda.min, drop = FALSE])
}

estimate_concert_qk_author_gtex <- function(dat, tolerance = 0.008,
                                             probability_epsilon = 1e-6,
                                             seed = 1L) {
  stopifnot(length(dat$X_source) == dat$K, length(dat$y_source) == dat$K)
  beta0 <- .concert_cv_ridge(dat$X, dat$y, seed)
  beta_source <- vapply(seq_len(dat$K), function(k) {
    .concert_cv_ridge(dat$X_source[[k]], dat$y_source[[k]], seed + k)
  }, numeric(dat$p))
  if (dat$K == 1L) beta_source <- matrix(beta_source, ncol = 1L)
  names_source <- colnames(dat$B)
  colnames(beta_source) <- names_source
  agreement <- colMeans(abs(sweep(beta_source, 1L, beta0, "-")) <= tolerance)
  q_raw <- min(agreement)
  q_used <- min(max(q_raw, probability_epsilon), 1 - probability_epsilon)
  list(
    q_raw = q_raw,
    q_used = q_used,
    agreement_by_source = agreement,
    tolerance = tolerance,
    probability_epsilon = probability_epsilon,
    target_ridge_beta = beta0,
    source_ridge_beta = beta_source
  )
}

fit_concert_real <- function(dat, cfg, seed = 1L) {
  stopifnot(
    is.matrix(dat$X_all), length(dat$y_all) == nrow(dat$X_all),
    length(dat$n_vec_trans) == dat$K + 1L,
    sum(dat$n_vec_trans) == nrow(dat$X_all)
  )
  q_info <- estimate_concert_qk_author_gtex(
    dat,
    tolerance = cfg$concert_qk_tolerance,
    probability_epsilon = cfg$concert_probability_epsilon,
    seed = seed
  )
  q0_used <- min(max(cfg$concert_q0, cfg$concert_probability_epsilon),
                 1 - cfg$concert_probability_epsilon)
  fit <- Concert(
    X = dat$X_all,
    y = dat$y_all,
    n_vec = as.numeric(dat$n_vec_trans),
    tau = rep(cfg$concert_tau, dat$K),
    q_k = q_info$q_used,
    eta = cfg$concert_eta,
    q_0 = q0_used,
    threshold = cfg$concert_threshold,
    max_iter = as.integer(cfg$concert_max_iter)
  )
  required <- c("m_beta0", "gamma_Z", "m_betak", "gamma_I", "V_beta0", "V_betak")
  if (!all(required %in% names(fit))) {
    stop("The official CONCERT backend returned an incomplete result.")
  }
  if (length(fit$m_beta0) != dat$p || length(fit$gamma_Z) != dat$p ||
      !identical(dim(fit$gamma_I), c(dat$p, dat$K))) {
    stop("The official CONCERT backend returned incompatible dimensions.")
  }
  beta_hat <- as.numeric(fit$m_beta0) * as.numeric(fit$gamma_Z)
  if (any(!is.finite(beta_hat)) || any(!is.finite(fit$gamma_Z)) ||
      any(!is.finite(fit$gamma_I))) {
    stop("The official CONCERT backend returned non-finite posterior quantities.")
  }
  colnames(fit$gamma_I) <- colnames(dat$B)
  rownames(fit$gamma_I) <- rownames(dat$B)
  fit$method <- "CONCERT_author_linear"
  fit$beta_hat <- beta_hat
  fit$omega <- as.numeric(fit$gamma_Z)
  fit$q_k_raw <- q_info$q_raw
  fit$q_k_used <- q_info$q_used
  fit$q_k_agreement_by_source <- q_info$agreement_by_source
  fit$q_k_tolerance <- q_info$tolerance
  fit$q_0_used <- q0_used
  fit$tau <- rep(cfg$concert_tau, dat$K)
  fit$eta <- cfg$concert_eta
  fit$upstream_commit <- .concert_upstream_commit
  # Upstream does not expose its final CAVI residual or iteration count.
  fit$converged <- NA
  fit$n_iter <- NA_integer_
  fit
}
