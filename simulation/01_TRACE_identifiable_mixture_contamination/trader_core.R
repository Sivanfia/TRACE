## =========================================================
## TRADER core: portable R wrapper around TRADER_utils_v6.cpp
## =========================================================
## The active C++ backend maintains r = y - X beta in observation space.
## This is an algebraically exact implementation of the original coordinate
## Gibbs update, but changes its main sweep from O(p^2) to O(n p).
## TRADER_utils_v6_dense_portable.cpp preserves the previous dense backend.

.trader_core_dir <- local({
  ofiles <- Filter(function(x) is.character(x) && length(x) == 1L && nzchar(x),
                   lapply(sys.frames(), function(x) x$ofile))
  ofile <- if (length(ofiles)) tail(ofiles, 1L)[[1L]] else NULL
  if (is.null(ofile)) getwd() else dirname(normalizePath(ofile))
})

.load_trader_cpp <- function(cpp_file = file.path(.trader_core_dir, "TRADER_utils_v6.cpp"),
                             rebuild = FALSE,
                             verbose = FALSE,
                             required_backend = "residual_np_v5") {
  backend_loaded <- tryCatch({
    if (exists("TRADER_backend_version", mode = "function", inherits = TRUE)) {
      as.character(TRADER_backend_version())
    } else {
      NA_character_
    }
  }, error = function(e) NA_character_)

  has_mcmc_loaded <- exists("TRADER_mcmc", mode = "function", inherits = TRUE)
  already_loaded <- has_mcmc_loaded &&
    isTRUE(identical(backend_loaded, required_backend))
  if (already_loaded && !isTRUE(rebuild)) return(invisible(cpp_file))

  if (!requireNamespace("Rcpp", quietly = TRUE)) stop("Package 'Rcpp' is required for TRADER.")
  if (!requireNamespace("RcppEigen", quietly = TRUE)) stop("Package 'RcppEigen' is required for TRADER.")
  if (!file.exists(cpp_file)) stop("Cannot find TRADER C++ file: ", cpp_file)

  ## Force recompilation when an older dense TRADER_mcmc is already present in
  ## the R session.  Otherwise sourceCpp could leave the old function active.
  rebuild_now <- isTRUE(rebuild) ||
    (has_mcmc_loaded && !isTRUE(identical(backend_loaded, required_backend)))
  Rcpp::sourceCpp(cpp_file, rebuild = rebuild_now, verbose = verbose)

  backend_after <- tryCatch(as.character(TRADER_backend_version()), error = function(e) NA_character_)
  if (!identical(backend_after, required_backend)) {
    stop("Wrong TRADER C++ backend loaded. Expected '", required_backend,
         "' but found '", backend_after, "'. Restart R and source trader_core.R again.")
  }
  invisible(cpp_file)
}

.load_trader_cpp()
.trader_cpp_backend <- as.character(TRADER_backend_version())
message("[TRADER C++ backend] ", .trader_cpp_backend, " (exact residual O(n p) update)")

if (!methods::isClass("TRADER")) {
  methods::setClass(
    "TRADER",
    slots = c(
      n = "numeric",
      p = "numeric",
      K = "numeric",
      y = "numeric",
      x = "matrix",
      w = "matrix",
      alpha = "numeric",
      v = "numeric",
      tau = "ANY",
      tau_lb = "numeric",
      tau_ub = "numeric",
      mcmc_samples = "list",
      mcmc_mean = "list",
      mcmc_quantile = "list"
    )
  )
}

CreateTRADERObject <- function(y, x, w = NULL, alpha = NULL, v = 0,
                               tau = "Cauchy", tau_lb = 0, tau_ub = Inf) {
  x <- as.matrix(x)
  y <- as.numeric(y)
  n <- nrow(x)
  p <- ncol(x)
  if (length(y) != n) stop("length(y) must equal nrow(x).")

  if (is.null(w)) {
    K <- 0L
    w <- matrix(0, nrow = p, ncol = 0L)
    alpha <- 1
  } else {
    w <- as.matrix(w)
    if (nrow(w) != p) stop("nrow(w) must equal ncol(x).")
    K <- ncol(w)
    if (is.null(alpha)) alpha <- rep(1, K + 1L)
  }

  alpha <- as.numeric(alpha)
  if (length(alpha) != K + 1L) {
    stop("TRADER alpha must have length K + 1, including the residual mass parameter.")
  }
  alpha[!is.finite(alpha) | alpha <= 0] <- 1e-4

  methods::new(
    "TRADER",
    n = as.numeric(n),
    p = as.numeric(p),
    K = as.numeric(K),
    y = y,
    x = x,
    w = w,
    alpha = alpha,
    v = as.numeric(v),
    tau = tau,
    tau_lb = as.numeric(tau_lb),
    tau_ub = as.numeric(tau_ub),
    mcmc_samples = list(n_samples = 0, n_burnin = 0),
    mcmc_mean = list(),
    mcmc_quantile = list()
  )
}

TRADER_MCMC <- function(TRADER, n_mcmc = 10000, n_burnin = 5000,
                        n_chain = 1, fold = 1.1, thin = 2, seed = 123,
                        compute_psrf = TRUE) {
  if (!methods::is(TRADER, "TRADER")) stop("TRADER must be a TRADER S4 object.")
  n_mcmc <- as.integer(n_mcmc)
  n_burnin <- as.integer(n_burnin)
  n_chain <- as.integer(n_chain)
  thin <- as.integer(thin)
  if (n_mcmc <= 0) stop("n_mcmc must be positive.")
  if (n_burnin < 0 || n_burnin >= n_mcmc) stop("Require 0 <= n_burnin < n_mcmc.")
  if (n_chain <= 0 || thin <= 0) stop("n_chain and thin must be positive.")

  trader_setup_start <- unname(proc.time()[["elapsed"]])
  p <- as.integer(TRADER@p)
  K <- as.integer(TRADER@K)
  xtx <- crossprod(TRADER@x)
  xty <- as.numeric(crossprod(TRADER@x, TRADER@y))
  yty <- as.numeric(crossprod(TRADER@y))

  d <- sqrt(pmax(diag(xtx), 1e-12))
  LD_pairs <- sweep(sweep(xtx, 1, d, "/"), 2, d, "/")
  diag(LD_pairs) <- 0
  LD_pairs <- LD_pairs^2
  LD_pairs[LD_pairs < 0.16] <- 0
  rs <- rowSums(LD_pairs)
  nz <- rs > 0
  LD_norm <- LD_pairs
  if (any(nz)) LD_norm[nz, ] <- LD_norm[nz, , drop = FALSE] / rs[nz]
  if (any(!nz)) LD_norm[!nz, ] <- 0
  ## Match the uploaded TRADER wrapper exactly: transpose the row-normalized
  ## proposal matrix so C++ column j contains probabilities from variable j.
  LD_pairs <- t(LD_norm)
  if (!requireNamespace("Matrix", quietly = TRUE)) stop("Package 'Matrix' is required for TRADER.")
  LD_pairs <- methods::as(Matrix::Matrix(LD_pairs, sparse = TRUE), "dgCMatrix")

  sample <- list()
  if (identical(TRADER@mcmc_samples$n_samples, 0) || is.null(TRADER@mcmc_samples$last_sample)) {
    sample[[1]] <- matrix(0, p, n_chain)
    sample[[2]] <- matrix(1, p, n_chain)
    sample[[3]] <- matrix(1 / (K + 1), K, n_chain)
    sample[[4]] <- rep(1 / (K + 1), n_chain)
    sample[[5]] <- rep(1, n_chain)
    if (is.numeric(TRADER@tau)) {
      sample[[6]] <- rep(TRADER@tau^2, n_chain)
    } else if (is.infinite(TRADER@tau_ub)) {
      sample[[6]] <- rep(TRADER@tau_lb^2 + 1, n_chain)
    } else {
      sample[[6]] <- rep((TRADER@tau_lb^2 + TRADER@tau_ub^2) / 2, n_chain)
    }
  } else {
    sample <- TRADER@mcmc_samples$last_sample
    if (is.numeric(TRADER@tau)) sample[[6]] <- rep(TRADER@tau^2, n_chain)
  }

  temp <- fold^(seq_len(n_chain) - 1L)
  tau_method <- if (is.numeric(TRADER@tau)) "constant" else as.character(TRADER@tau)
  tau_lb_inv <- if (TRADER@tau_lb <= 0) Inf else 1 / TRADER@tau_lb^2
  tau_ub_inv <- if (is.infinite(TRADER@tau_ub)) 0 else 1 / TRADER@tau_ub^2

  trader_setup_time_sec <- unname(proc.time()[["elapsed"]]) - trader_setup_start
  trader_cpp_tm <- system.time({
    samples <- TRADER_mcmc(
      TRADER, xtx, xty, yty, sample, LD_pairs,
      tau_method, tau_lb_inv, tau_ub_inv,
      n_chain, temp, n_mcmc, thin, as.integer(seed)
    )
  })
  trader_cpp_time_sec <- unname(trader_cpp_tm[["elapsed"]])

  draw_names <- paste0("mcmc_", seq_len(n_mcmc))
  colnames(samples[[1]]) <- draw_names
  colnames(samples[[2]]) <- draw_names
  colnames(samples[[3]]) <- draw_names
  names(samples[[4]]) <- draw_names
  if (!is.numeric(TRADER@tau)) names(samples[[5]]) <- draw_names

  keep <- seq.int(n_burnin + 1L, n_mcmc)
  TRADER@mcmc_samples$beta <- samples[[1]][, keep, drop = FALSE]
  TRADER@mcmc_samples$lambda2 <- samples[[2]][, keep, drop = FALSE]
  TRADER@mcmc_samples$eta <- samples[[3]][, keep, drop = FALSE]
  TRADER@mcmc_samples$sigma2 <- samples[[4]][keep]
  if (!is.numeric(TRADER@tau)) TRADER@mcmc_samples$tau2 <- samples[[5]][keep]
  TRADER@mcmc_samples$last_sample <- samples$last_sample
  TRADER@mcmc_samples$n_burnin <- n_burnin
  TRADER@mcmc_samples$n_samples <- length(keep)
  TRADER@mcmc_samples$backend <- if (!is.null(samples$backend)) {
    as.character(samples$backend)
  } else {
    as.character(TRADER_backend_version())
  }
  TRADER@mcmc_samples$residual_refresh_count <- if (!is.null(samples$residual_refresh_count)) {
    as.integer(samples$residual_refresh_count)
  } else {
    NA_integer_
  }
  TRADER@mcmc_samples$max_residual_refresh_error <- if (!is.null(samples$max_residual_refresh_error)) {
    as.numeric(samples$max_residual_refresh_error)
  } else {
    NA_real_
  }
  TRADER@mcmc_samples$setup_time_sec <- trader_setup_time_sec
  TRADER@mcmc_samples$cpp_time_sec <- trader_cpp_time_sec

  beta_samples <- TRADER@mcmc_samples$beta
  TRADER@mcmc_mean$beta <- rowMeans(beta_samples)
  TRADER@mcmc_quantile$beta <- t(apply(beta_samples, 1, stats::quantile, probs = c(0.025, 0.975)))

  if (isTRUE(compute_psrf) && ncol(beta_samples) >= 4L) {
    nd <- ncol(beta_samples)
    n2 <- floor(nd / 2)
    first <- seq_len(n2)
    second <- seq.int(nd - n2 + 1L, nd)
    psrf <- vapply(seq_len(p), function(j) {
      s1 <- beta_samples[j, first]
      s2 <- beta_samples[j, second]
      B <- n2 * (mean(s1) - mean(s2))^2 / 2
      W <- (stats::var(s1) + stats::var(s2)) / 2
      sqrt(max(0, 1 - 1 / n2 + B / (W + 5e-324) / n2))
    }, numeric(1))
  } else {
    psrf <- rep(NA_real_, p)
  }
  TRADER@mcmc_samples$PSRF_beta <- psrf
  TRADER
}
