## =========================================================
## Borrow functions: alpha-heterogeneous common-bad CVB, no damping, no freeze
## Compatible with dgp_functions_multisource_mixed_block_priorfix.R
##
## Design principles:
##   - Do NOT change the CVB update formula.
##   - No damping.
##   - No gamma freeze.
##   - Support alpha_method = old / fixed / oracle.
##   - Keep target-only ACVB, VBSS, and classical baselines for comparison.
## =========================================================

suppressPackageStartupMessages({
  library(glmnet)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

safe_log <- function(x, eps = 1e-12) log(pmax(x, eps))
expit <- function(x) plogis(x)
clamp <- function(x, lo, hi) pmin(pmax(x, lo), hi)

solve_reg <- function(A, b) {
  ridge_seq <- c(1e-10, 1e-8, 1e-6, 1e-4, 1e-2)
  for (r in ridge_seq) {
    out <- tryCatch(
      solve(A + diag(r, nrow(A)), b),
      error = function(e) NULL
    )
    if (!is.null(out)) return(out)
  }
  qr.solve(A + diag(1e-2, nrow(A)), b)
}

solve_spd <- function(A, B) {
  ridge_seq <- c(0, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2)
  for (r in ridge_seq) {
    out <- tryCatch({
      R <- chol(A + diag(r, nrow(A)))
      backsolve(R, forwardsolve(t(R), B))
    }, error = function(e) NULL)
    if (!is.null(out)) return(out)
  }
  solve_reg(A, B)
}

get_ns_eff <- function(dat) {
  ## Current theoretical update uses one effective source sample size.
  ## If dat$n_vec exists, use its mean to keep the same scalar update.
  if (!is.null(dat$n_vec)) return(mean(dat$n_vec))
  dat$ns %||% dat$n
}

## =========================================================
## Metrics and prediction
## =========================================================

prediction_mse <- function(beta_hat, dat) {
  if (!is.null(dat$X_test) && !is.null(dat$beta_true)) {
    return(mean((as.numeric(dat$X_test %*% beta_hat) -
                   as.numeric(dat$X_test %*% dat$beta_true))^2))
  }
  mean((as.numeric(dat$X %*% beta_hat) - as.numeric(dat$X %*% dat$beta_true))^2)
}

calc_selection_metrics <- function(sel_hat, beta_true) {
  sel_hat <- as.integer(sel_hat != 0)
  sel_true <- as.integer(beta_true != 0)

  tp <- sum(sel_hat == 1 & sel_true == 1)
  fp <- sum(sel_hat == 1 & sel_true == 0)
  fn <- sum(sel_hat == 0 & sel_true == 1)
  tn <- sum(sel_hat == 0 & sel_true == 0)

  precision <- ifelse(tp + fp > 0, tp / (tp + fp), 0)
  recall <- ifelse(tp + fn > 0, tp / (tp + fn), 0)
  f1 <- ifelse(precision + recall > 0,
               2 * precision * recall / (precision + recall), 0)
  fdr <- ifelse(tp + fp > 0, fp / (tp + fp), 0)
  fpr <- ifelse(fp + tn > 0, fp / (fp + tn), 0)

  list(
    tp = tp, fp = fp, fn = fn, tn = tn,
    precision = precision,
    recall = recall,
    tpr = recall,
    fpr = fpr,
    fdr = fdr,
    f1 = f1,
    acc = mean(sel_hat == sel_true),
    model_size = sum(sel_hat)
  )
}

get_selection_metrics_omega <- function(omega, beta_true, thr = 0.5) {
  calc_selection_metrics(as.integer(omega > thr), beta_true)
}

get_selection_metrics_topk <- function(score, beta_true, k = sum(beta_true != 0)) {
  p <- length(score)
  sel_hat <- rep(0, p)
  k <- max(0, min(as.integer(k), p))
  if (k > 0) {
    ord <- order(score, decreasing = TRUE)
    sel_hat[ord[seq_len(k)]] <- 1
  }
  calc_selection_metrics(sel_hat, beta_true)
}

get_selection_metrics_beta <- function(beta_hat, beta_true, thr = 1e-8) {
  calc_selection_metrics(as.integer(abs(beta_hat) > thr), beta_true)
}

make_selection_table <- function(method, scenario, omega, beta_true) {
  sel_05 <- get_selection_metrics_omega(omega, beta_true, thr = 0.5)
  sel_01 <- get_selection_metrics_omega(omega, beta_true, thr = 0.1)
  sel_top <- get_selection_metrics_topk(omega, beta_true)

  pack_row <- function(rule, sel) {
    data.frame(
      scenario = scenario,
      method = method,
      rule = rule,
      tp = sel$tp,
      fp = sel$fp,
      fn = sel$fn,
      tn = sel$tn,
      model_size = sel$model_size,
      precision = sel$precision,
      recall = sel$recall,
      f1 = sel$f1,
      fdr = sel$fdr,
      fpr = sel$fpr,
      acc = sel$acc,
      stringsAsFactors = FALSE
    )
  }

  rbind(
    pack_row("omega>0.5", sel_05),
    pack_row("omega>0.1", sel_01),
    pack_row("topk", sel_top)
  )
}

## =========================================================
## Common objects for borrowing update
## =========================================================

build_common_objects <- function(omega, gamma, alpha, dat, S0bar_inv) {
  XtX <- dat$XtX
  X <- dat$X
  B <- as.matrix(dat$B)
  n <- dat$n
  ns_eff <- get_ns_eff(dat)
  x_norm2 <- dat$x_norm2

  Omega <- omega %o% omega
  diag(Omega) <- omega

  S1bar_inv <- (ns_eff * x_norm2) / n
  S_inv_all <- S0bar_inv + (S1bar_inv - S0bar_inv) * gamma

  A0 <- XtX * Omega
  diag(A0) <- diag(A0) + S_inv_all
  A0 <- A0 + diag(1e-10, length(omega))

  alpha <- as.numeric(alpha)
  bbar <- as.numeric(B %*% alpha)

  theta <- omega * gamma * bbar
  mu <- as.numeric(X %*% theta)
  r <- dat$y - mu

  list(
    Omega = Omega,
    S1bar_inv = S1bar_inv,
    S_inv_all = S_inv_all,
    A0 = A0,
    bbar = bbar,
    theta = theta,
    mu = mu,
    r = r
  )
}

## =========================================================
## Alpha update: old eta update
## =========================================================

update_alpha_old <- function(alpha, omega, gamma, dat,
                             A_y, B_y, S0bar_inv,
                             eps_eta = 1e-6) {
  X <- dat$X
  B <- as.matrix(dat$B)

  obj <- build_common_objects(omega, gamma, alpha, dat, S0bar_inv)

  wr <- omega * as.numeric(crossprod(X, obj$r))
  tmp <- solve_spd(obj$A0, wr)
  Hr <- obj$r - as.numeric(X %*% (omega * tmp))

  denom <- B_y + 0.5 * sum(obj$r * Hr)
  denom <- max(denom, 1e-12)

  qvec <- gamma * omega * as.numeric(crossprod(X, -Hr))
  numer <- as.numeric(crossprod(B, qvec))

  eta_raw <- 1 - alpha * (A_y + dat$n / 2) * numer / denom
  eta <- pmax(eps_eta, eta_raw)
  alpha_new <- eta / sum(eta)

  list(alpha = alpha_new)
}

## =========================================================
## Omega update
## =========================================================

update_omega_j_old <- function(j, omega, gamma, bbar, dat,
                               S0bar_inv, A_y, B_y, E_logit_rho_z) {
  idx <- setdiff(seq_along(omega), j)
  w <- omega[idx]

  XtX <- dat$XtX
  Xty <- dat$Xty
  x_norm2 <- dat$x_norm2
  yTy <- dat$yTy
  n <- dat$n
  ns_eff <- get_ns_eff(dat)

  S1bar_inv <- (ns_eff * x_norm2) / n
  S_inv_all <- S0bar_inv + (S1bar_inv - S0bar_inv) * gamma
  S_inv_j <- S_inv_all[j]
  S_inv_m <- S_inv_all[idx]

  v_all <- (S_inv_all * gamma) * bbar
  v_j <- v_all[j]
  v_m <- v_all[idx]

  XtX_mm <- XtX[idx, idx, drop = FALSE]
  M_minus <- XtX_mm * (w %o% w)
  diag(M_minus) <- diag(XtX_mm) * w
  M_minus <- M_minus + diag(S_inv_m, length(idx)) + diag(1e-10, length(idx))

  u <- w * Xty[idx]
  ccol <- w * XtX[idx, j]

  SOL <- solve_spd(M_minus, cbind(u, v_m, ccol))
  sol_u <- SOL[, 1]
  sol_v <- SOL[, 2]
  sol_c <- SOL[, 3]

  D_j <- S_inv_j + x_norm2[j] - sum(ccol * sol_c)
  D_j <- pmax(D_j, 1e-10)

  A_j <- Xty[j] - sum(ccol * sol_u)
  B_j <- v_j - sum(ccol * sol_v)

  alpha_term <- sum((gamma * bbar)^2 * S_inv_all)

  T0 <- yTy -
    sum(u * sol_u) -
    2 * sum(u * sol_v) -
    sum(v_m * sol_v) -
    (v_j^2) / pmax(S_inv_j, 1e-12) +
    alpha_term

  Tdiff <- -(A_j + B_j)^2 / pmax(D_j, 1e-12) +
    (v_j^2) / pmax(S_inv_j, 1e-12)

  T1 <- T0 + Tdiff
  T0 <- pmax(T0, 1e-10)
  T1 <- pmax(T1, 1e-10)

  lin <- E_logit_rho_z -
    0.5 * safe_log(D_j) +
    0.5 * safe_log(S_inv_j) -
    (A_y + n / 2) *
    (safe_log(B_y + 0.5 * T1) - safe_log(B_y + 0.5 * T0))

  lin <- clamp(lin, -30, 30)
  expit(lin)
}

## =========================================================
## Gamma update
## =========================================================

update_gamma_j <- function(j, omega, gamma, bbar, dat,
                           S0bar_inv, A_y, B_y, E_logit_rho_c) {
  idx <- setdiff(seq_along(omega), j)
  w <- omega[idx]
  wj <- omega[j]

  XtX <- dat$XtX
  Xty <- dat$Xty
  x_norm2 <- dat$x_norm2
  yTy <- dat$yTy
  n <- dat$n
  ns_eff <- get_ns_eff(dat)

  S1bar_inv <- (ns_eff * x_norm2) / n
  S_inv_all <- S0bar_inv + (S1bar_inv - S0bar_inv) * gamma
  S_inv_m <- S_inv_all[idx]

  v_all <- (S_inv_all * gamma) * bbar
  v_m <- v_all[idx]

  XtX_mm <- XtX[idx, idx, drop = FALSE]
  M_minus <- XtX_mm * (w %o% w)
  diag(M_minus) <- diag(XtX_mm) * w
  M_minus <- M_minus + diag(S_inv_m, length(idx)) + diag(1e-10, length(idx))

  u <- w * Xty[idx]
  r_minus <- u + v_m

  ccol <- w * XtX[idx, j]
  crow <- XtX[j, idx] * w

  SOL <- solve_spd(M_minus, cbind(r_minus, ccol))
  sol_r <- SOL[, 1]
  sol_c <- SOL[, 2]

  cTc <- sum(ccol * sol_c)

  D0j <- S0bar_inv + wj * (x_norm2[j] - cTc)
  D1j <- S1bar_inv[j] + wj * (x_norm2[j] - cTc)
  D0j <- pmax(D0j, 1e-10)
  D1j <- pmax(D1j, 1e-10)

  R0 <- wj * (Xty[j] - sum(crow * sol_r))
  v1j <- S1bar_inv[j] * bbar[j]
  R1 <- R0 + v1j

  alpha_minus <- sum((gamma[idx] * bbar[idx])^2 * S_inv_m)

  Tc0 <- yTy - sum(r_minus * sol_r) - (R0^2) / D0j + alpha_minus
  Tc1 <- yTy - sum(r_minus * sol_r) - (R1^2) / D1j + alpha_minus +
    (S1bar_inv[j] * bbar[j]^2)

  Tc0 <- pmax(Tc0, 1e-10)
  Tc1 <- pmax(Tc1, 1e-10)

  Delta <- -0.5 * safe_log(1 + (S1bar_inv[j] - S0bar_inv) / D0j) -
    (A_y + n / 2) *
    (safe_log(B_y + 0.5 * Tc1) - safe_log(B_y + 0.5 * Tc0))

  lin_g <- E_logit_rho_c + Delta
  lin_g <- clamp(lin_g, -30, 30)

  expit(lin_g)
}

## =========================================================
## Beta estimator: keep original borrowing posterior mean formula
## =========================================================

get_beta_hat <- function(alpha, omega, gamma, dat, S0bar_inv) {
  obj <- build_common_objects(omega, gamma, alpha, dat, S0bar_inv)
  rhs <- omega * dat$Xty + (obj$S_inv_all * gamma) * obj$bbar
  as.numeric(solve_spd(obj$A0, rhs))
}

## =========================================================
## Screening and initialization
## =========================================================

fit_ebic_glmnet_lasso <- function(X, y, gamma_ebic = 0.5) {
  n <- nrow(X)
  p <- ncol(X)

  fit <- glmnet(
    x = X,
    y = y,
    alpha = 1,
    intercept = FALSE,
    standardize = FALSE
  )

  beta_mat <- as.matrix(coef(fit))[-1, , drop = FALSE]
  rss <- colSums((matrix(y, n, ncol(beta_mat)) - X %*% beta_mat)^2)
  dof <- colSums(abs(beta_mat) > 1e-10)

  ebic <- log(pmax(rss / n, 1e-12)) +
    (dof / n) * (log(n) + 2 * gamma_ebic * log(pmax(p, 2)))

  best <- which.min(ebic)
  list(beta = as.numeric(beta_mat[, best]),
       ebic = ebic[best],
       model_size = dof[best])
}

fit_ebic_ncvreg <- function(X, y, penalty = c("MCP", "SCAD"), gamma_ebic = 0.5) {
  penalty <- match.arg(penalty)

  if (!requireNamespace("ncvreg", quietly = TRUE)) {
    return(fit_ebic_glmnet_lasso(X, y, gamma_ebic = gamma_ebic))
  }

  n <- nrow(X)
  p <- ncol(X)

  fit <- ncvreg::ncvreg(
    X = X,
    y = y,
    penalty = penalty,
    intercept = FALSE,
    standardize = FALSE
  )

  beta_mat <- as.matrix(fit$beta)
  if (nrow(beta_mat) == p + 1) beta_mat <- beta_mat[-1, , drop = FALSE]

  if (nrow(beta_mat) != p) {
    stop("Unexpected coefficient dimension returned by ncvreg.")
  }

  rss <- colSums((matrix(y, n, ncol(beta_mat)) - X %*% beta_mat)^2)
  dof <- colSums(abs(beta_mat) > 1e-10)

  ebic <- log(pmax(rss / n, 1e-12)) +
    (dof / n) * (log(n) + 2 * gamma_ebic * log(pmax(p, 2)))

  best <- which.min(ebic)
  list(beta = as.numeric(beta_mat[, best]),
       ebic = ebic[best],
       model_size = dof[best])
}

fit_penalized_baseline <- function(dat,
                                   penalty = c("lasso", "mcp", "scad"),
                                   gamma_ebic = 0.5) {
  penalty <- match.arg(penalty)

  fit <- switch(
    penalty,
    lasso = fit_ebic_glmnet_lasso(dat$X, dat$y, gamma_ebic = gamma_ebic),
    mcp = fit_ebic_ncvreg(dat$X, dat$y, penalty = "MCP", gamma_ebic = gamma_ebic),
    scad = fit_ebic_ncvreg(dat$X, dat$y, penalty = "SCAD", gamma_ebic = gamma_ebic)
  )

  beta_hat <- fit$beta
  sel_beta <- get_selection_metrics_beta(beta_hat, dat$beta_true, thr = 1e-8)

  list(
    method = paste0("no_borrow_", toupper(penalty), "_EBIC"),
    beta_hat = beta_hat,
    pred_mse = prediction_mse(beta_hat, dat),
    sel_beta = sel_beta,
    tp_beta = sel_beta$tp,
    fp_beta = sel_beta$fp,
    fn_beta = sel_beta$fn,
    acc_beta = sel_beta$acc,
    ebic = fit$ebic,
    model_size = fit$model_size
  )
}

screening_set <- function(dat,
                          method = c("lasso", "mcp", "scad", "mix"),
                          gamma_ebic = 0.5) {
  method <- match.arg(method)

  if (method == "lasso") {
    beta <- fit_ebic_glmnet_lasso(dat$X, dat$y, gamma_ebic = gamma_ebic)$beta
    return(which(abs(beta) > 1e-10))
  }

  if (method == "mcp") {
    beta <- fit_ebic_ncvreg(dat$X, dat$y, penalty = "MCP", gamma_ebic = gamma_ebic)$beta
    return(which(abs(beta) > 1e-10))
  }

  if (method == "scad") {
    beta <- fit_ebic_ncvreg(dat$X, dat$y, penalty = "SCAD", gamma_ebic = gamma_ebic)$beta
    return(which(abs(beta) > 1e-10))
  }

  sort(unique(c(
    screening_set(dat, "lasso", gamma_ebic),
    screening_set(dat, "mcp", gamma_ebic),
    screening_set(dat, "scad", gamma_ebic)
  )))
}

make_init_omega <- function(dat,
                            init_method = c("mild", "full", "null", "lasso", "mcp", "scad",
                                            "mix", "source", "custom"),
                            init_omega = NULL,
                            low = 0.10,
                            high = 0.90,
                            gamma_ebic = 0.5) {
  init_method <- match.arg(init_method)
  p <- dat$p

  if (init_method == "custom") {
    if (is.null(init_omega) || length(init_omega) != p) {
      stop("For init_method='custom', init_omega must be a length-p vector.")
    }
    return(clamp(init_omega, 1e-6, 1 - 1e-6))
  }

  if (init_method == "mild") return(rep(0.5, p))
  if (init_method == "full") return(rep(high, p))
  if (init_method == "null") return(rep(low, p))

  sel <- integer(0)

  if (init_method %in% c("lasso", "mcp", "scad", "mix")) {
    sel <- screening_set(dat, method = init_method, gamma_ebic = gamma_ebic)
  }

  if (init_method == "source") {
    if (is.null(dat$B)) stop("dat$B is needed for source initialization.")
    B <- as.matrix(dat$B)

    ## Use the largest absolute source signal as a non-oracle source screening.
    ## This can include spurious source coordinates, which is intentional in
    ## mixed-source diagnostics.
    score <- apply(abs(B), 1, max)

    k <- max(1, min(sum(dat$beta_true != 0), p))
    sel <- order(score, decreasing = TRUE)[seq_len(k)]
  }

  omega <- rep(low, p)
  if (length(sel) > 0) omega[sel] <- high
  clamp(omega, 1e-6, 1 - 1e-6)
}

make_init_alpha <- function(dat,
                            init_alpha = c("uniform", "source1", "random", "random_simplex", "oracle"),
                            alpha_fixed = NULL,
                            alpha_method = c("old", "fixed", "oracle")) {
  init_alpha <- match.arg(init_alpha)
  alpha_method <- match.arg(alpha_method)
  K <- dat$K

  if (alpha_method == "oracle") {
    if (is.null(dat$alpha_true)) stop("dat$alpha_true is needed for alpha_method='oracle'.")
    alpha <- as.numeric(dat$alpha_true)
  } else if (alpha_method == "fixed") {
    if (is.null(alpha_fixed)) alpha_fixed <- rep(1 / K, K)
    alpha <- as.numeric(alpha_fixed)
  } else if (init_alpha == "source1") {
    alpha <- c(1, rep(0, K - 1))
  } else if (init_alpha %in% c("random", "random_simplex")) {
    alpha <- stats::rexp(K, rate = 1)
  } else if (init_alpha == "oracle") {
    if (is.null(dat$alpha_true)) stop("dat$alpha_true is needed for init_alpha='oracle'.")
    alpha <- as.numeric(dat$alpha_true)
  } else {
    alpha <- rep(1 / K, K)
  }

  if (sum(alpha) <= 0 || any(!is.finite(alpha))) {
    alpha <- rep(1 / K, K)
  }
  alpha / sum(alpha)
}

make_init_gamma_borrow <- function(dat,
                                   init_gamma = 0.8,
                                   alpha = rep(1 / dat$K, dat$K),
                                   gamma_ebic = 0.5) {
  p <- dat$p

  if (is.numeric(init_gamma)) {
    if (length(init_gamma) == 1) {
      return(rep(clamp(init_gamma, 1e-4, 1 - 1e-4), p))
    }
    if (length(init_gamma) == p) {
      return(clamp(as.numeric(init_gamma), 1e-4, 1 - 1e-4))
    }
    stop("Numeric init_gamma must have length 1 or dat$p.")
  }

  init_gamma <- match.arg(as.character(init_gamma),
                          choices = c("half", "clean", "agreement"))

  if (init_gamma == "half") return(rep(0.5, p))
  if (init_gamma == "clean") return(rep(0.8, p))

  beta_target <- fit_ebic_glmnet_lasso(dat$X, dat$y, gamma_ebic = gamma_ebic)$beta
  bbar <- as.numeric(as.matrix(dat$B) %*% alpha)
  diff <- abs(bbar - beta_target)
  scale <- stats::quantile(diff, probs = 0.90, names = FALSE, na.rm = TRUE)
  if (!is.finite(scale) || scale <= 0) scale <- max(diff, na.rm = TRUE)
  if (!is.finite(scale) || scale <= 0) return(rep(0.5, p))

  gamma <- 0.90 - 0.80 * pmin(diff / scale, 1)
  clamp(gamma, 0.10, 0.90)
}

## =========================================================
## Borrowing fit: no damping, no freeze
## =========================================================

fit_borrow_model_fixedhyper <- function(dat,
                                        prior_par_borrow,
                                        alpha_method = c("old", "fixed", "oracle"),
                                        alpha_fixed = NULL,
                                        init_alpha = c("uniform", "source1", "random", "random_simplex", "oracle"),
                                        init_method = c("source", "mild", "full", "null", "mcp", "lasso",
                                                        "scad", "mix", "custom"),
                                        init_omega = NULL,
                                        init_gamma = "half",
                                        sigma_beta_sq = 100,
                                        max_iter = 100,
                                        tol = 1e-6,
                                        gamma_ebic = 0.5,
                                        verbose = FALSE) {
  alpha_method <- match.arg(alpha_method)
  init_alpha <- match.arg(init_alpha)
  init_method <- match.arg(init_method)
  A_y <- 1
  B_y <- 1

  A_c <- prior_par_borrow$A_c
  B_c <- prior_par_borrow$B_c
  A_z <- prior_par_borrow$A_z
  B_z <- prior_par_borrow$B_z

  S0bar <- sigma_beta_sq
  S0bar_inv <- 1 / S0bar

  p <- dat$p
  K <- dat$K

  omega <- make_init_omega(
    dat = dat,
    init_method = init_method,
    init_omega = init_omega,
    gamma_ebic = gamma_ebic
  )
  omega_init_vec <- omega

  alpha <- make_init_alpha(
    dat = dat,
    init_alpha = init_alpha,
    alpha_fixed = alpha_fixed,
    alpha_method = alpha_method
  )
  alpha_init_vec <- alpha

  ## Initial value is separated from gamma prior.
  ## The prior enters through E_logit_rho_c; init_gamma only sets the starting point.
  gamma <- make_init_gamma_borrow(
    dat = dat,
    init_gamma = init_gamma,
    alpha = alpha,
    gamma_ebic = gamma_ebic
  )
  gamma_init_vec <- gamma

  E_logit_rho_z <- digamma(A_z) - digamma(B_z)
  E_logit_rho_c <- digamma(A_c) - digamma(B_c)

  n_iter <- max_iter

  for (t in seq_len(max_iter)) {
    omega_old <- omega
    gamma_old <- gamma
    alpha_old <- alpha

    bbar <- as.numeric(as.matrix(dat$B) %*% alpha)

    omega_prop <- omega
    for (j in seq_len(p)) {
      omega_prop[j] <- update_omega_j_old(
        j = j,
        omega = omega,
        gamma = gamma,
        bbar = bbar,
        dat = dat,
        S0bar_inv = S0bar_inv,
        A_y = A_y,
        B_y = B_y,
        E_logit_rho_z = E_logit_rho_z
      )
    }
    omega <- clamp(omega_prop, 1e-6, 1 - 1e-6)

    gamma_prop <- gamma
    for (j in seq_len(p)) {
      gamma_prop[j] <- update_gamma_j(
        j = j,
        omega = omega,
        gamma = gamma,
        bbar = bbar,
        dat = dat,
        S0bar_inv = S0bar_inv,
        A_y = A_y,
        B_y = B_y,
        E_logit_rho_c = E_logit_rho_c
      )
    }
    gamma <- clamp(gamma_prop, 1e-6, 1 - 1e-6)

    if (alpha_method == "old") {
      alpha <- update_alpha_old(alpha, omega, gamma, dat, A_y, B_y, S0bar_inv)$alpha
    } else if (alpha_method == "fixed") {
      alpha <- as.numeric(alpha_fixed) / sum(alpha_fixed)
    } else if (alpha_method == "oracle") {
      alpha <- as.numeric(dat$alpha_true) / sum(dat$alpha_true)
    }

    d_omega <- max(abs(omega - omega_old))
    d_gamma <- max(abs(gamma - gamma_old))
    d_alpha <- max(abs(alpha - alpha_old))

    if (verbose) {
      cat(sprintf("[%s/%s] iter %3d | d_omega=%.3e d_gamma=%.3e d_alpha=%.3e | size05=%d\n",
                  alpha_method, init_method, t, d_omega, d_gamma, d_alpha,
                  sum(omega > 0.5)))
    }

    if (t > 5 && max(d_omega, d_gamma, d_alpha) < tol) {
      n_iter <- t
      break
    }
  }

  beta_hat <- get_beta_hat(alpha, omega, gamma, dat, S0bar_inv)

  sel_05 <- get_selection_metrics_omega(omega, dat$beta_true, thr = 0.5)
  sel_01 <- get_selection_metrics_omega(omega, dat$beta_true, thr = 0.1)
  sel_top <- get_selection_metrics_topk(omega, dat$beta_true)

  method_name <- if (alpha_method == "oracle") {
    "borrow_oracle_alpha"
  } else if (alpha_method == "fixed") {
    "borrow_fixed_alpha"
  } else {
    "borrow_old"
  }

  list(
    method = paste0(method_name, "_init_", init_method),
    alpha = alpha,
    omega = omega,
    gamma = gamma,
    beta_hat = beta_hat,
    pred_mse = prediction_mse(beta_hat, dat),
    sel_05 = sel_05,
    sel_01 = sel_01,
    sel_top = sel_top,
    mean_omega_signal = mean(omega[dat$tag_true == 1]),
    mean_omega_null = mean(omega[dat$tag_true == 0]),
    mean_gamma_signal = mean(gamma[dat$tag_true == 1]),
    mean_gamma_null = mean(gamma[dat$tag_true == 0]),
    fixed_E_logit_rho_z = E_logit_rho_z,
    fixed_E_logit_rho_c = E_logit_rho_c,
    init_method = init_method,
    init_alpha = init_alpha,
    init_gamma = init_gamma,
    omega_init = omega_init_vec,
    alpha_init = alpha_init_vec,
    gamma_init = gamma_init_vec,
    mean_gamma_init_signal = mean(gamma_init_vec[dat$tag_true == 1]),
    mean_gamma_init_null = mean(gamma_init_vec[dat$tag_true == 0]),
    sigma_beta_sq = sigma_beta_sq,
    n_iter = n_iter
  )
}

## =========================================================
## No-borrow VBSS baseline
## =========================================================

fit_noborrow_vbss_fixedhyper <- function(dat,
                                         prior_par_vbss,
                                         tau1_sq = 1,
                                         A_y = 1,
                                         B_y = 1,
                                         init_method = c("prior", "null", "full", "mcp", "lasso",
                                                         "scad", "mix", "mild", "custom"),
                                         init_omega = NULL,
                                         max_iter = 150,
                                         tol = 1e-6,
                                         gamma_ebic = 0.5,
                                         verbose = FALSE) {
  init_method <- match.arg(init_method)

  X <- dat$X
  y <- dat$y
  n <- dat$n
  p <- dat$p
  x2 <- colSums(X^2)

  A_z <- prior_par_vbss$A_z
  B_z <- prior_par_vbss$B_z
  E_logit_rho_z <- digamma(A_z) - digamma(B_z)

  if (init_method == "prior") {
    omega <- rep(clamp(prior_par_vbss$rho_z_mean, 1e-4, 1 - 1e-4), p)
  } else {
    omega <- make_init_omega(
      dat = dat,
      init_method = init_method,
      init_omega = init_omega,
      gamma_ebic = gamma_ebic
    )
  }

  mu <- rep(0, p)
  s2 <- rep(max(var(y), 1e-8) * 0.1, p)
  sigma2 <- max(var(y), 1e-8)

  beta_mean <- omega * mu
  n_iter <- max_iter

  for (t in seq_len(max_iter)) {
    omega_old <- omega
    mu_old <- mu
    sigma2_old <- sigma2

    omega_prop <- omega
    mu_prop <- mu

    for (j in seq_len(p)) {
      r_j <- y - as.numeric(X %*% beta_mean) + X[, j] * beta_mean[j]

      s2[j] <- sigma2 / (x2[j] + 1 / tau1_sq)
      mu_prop[j] <- s2[j] * sum(X[, j] * r_j) / sigma2

      logit_omega_j <- E_logit_rho_z +
        0.5 * log(pmax(s2[j] / (tau1_sq * sigma2), 1e-12)) +
        0.5 * (mu_prop[j]^2 / pmax(s2[j], 1e-12))

      logit_omega_j <- clamp(logit_omega_j, -30, 30)
      omega_prop[j] <- expit(logit_omega_j)
    }

    omega <- clamp(omega_prop, 1e-6, 1 - 1e-6)
    mu <- mu_prop
    beta_mean <- omega * mu

    resid2 <- sum((y - as.numeric(X %*% beta_mean))^2)
    var_part <- sum(x2 * omega * s2)
    prior_part <- sum(omega * (mu^2 + s2) / tau1_sq)

    sigma2 <- (2 * B_y + resid2 + var_part + prior_part) /
      (n + 2 * A_y + sum(omega) + 2)

    d_omega <- max(abs(omega - omega_old))
    d_mu <- max(abs(mu - mu_old))
    d_sig <- abs(sigma2 - sigma2_old)

    if (verbose) {
      cat(sprintf("[VBSS/%s] iter %3d | d_omega=%.3e d_mu=%.3e d_sigma2=%.3e | size05=%d\n",
                  init_method, t, d_omega, d_mu, d_sig, sum(omega > 0.5)))
    }

    if (!is.finite(sigma2)) break
    if (t > 5 && max(d_omega, d_mu, d_sig) < tol) {
      n_iter <- t
      break
    }
  }

  beta_hat <- beta_mean

  sel_05 <- get_selection_metrics_omega(omega, dat$beta_true, thr = 0.5)
  sel_01 <- get_selection_metrics_omega(omega, dat$beta_true, thr = 0.1)
  sel_top <- get_selection_metrics_topk(omega, dat$beta_true)

  list(
    method = paste0("no_borrow_VBSS_init_", init_method),
    omega = omega,
    beta_hat = beta_hat,
    pred_mse = prediction_mse(beta_hat, dat),
    sel_05 = sel_05,
    sel_01 = sel_01,
    sel_top = sel_top,
    mean_omega_signal = mean(omega[dat$tag_true == 1]),
    mean_omega_null = mean(omega[dat$tag_true == 0]),
    fixed_E_logit_rho_z = E_logit_rho_z,
    init_method = init_method,
    n_iter = n_iter
  )
}

## =========================================================
## No-borrow ACVB baseline
## =========================================================

fit_noborrow_acvb <- function(dat,
                              init_method = c("null", "full", "mcp", "lasso", "scad",
                                              "mix", "mild", "custom"),
                              init_omega = NULL,
                              sigma_beta_sq = 1e5,
                              A_y = 1,
                              B_y = 1,
                              pi_prior = NULL,
                              max_iter = 100,
                              tol = 1e-6,
                              gamma_ebic = 0.5,
                              verbose = FALSE) {
  init_method <- match.arg(init_method)

  X <- dat$X
  y <- dat$y
  n <- dat$n
  p <- dat$p
  XtX <- dat$XtX
  Xty <- dat$Xty
  yTy <- dat$yTy

  if (is.null(pi_prior)) pi_prior <- sum(dat$beta_true != 0) / p
  pi_prior <- clamp(pi_prior, 1e-8, 1 - 1e-8)
  logit_pi <- qlogis(pi_prior)

  omega <- make_init_omega(
    dat = dat,
    init_method = init_method,
    init_omega = init_omega,
    gamma_ebic = gamma_ebic
  )

  acvb_beta_given_w <- function(w) {
    Omega <- w %o% w
    diag(Omega) <- w
    G <- XtX * Omega
    diag(G) <- diag(G) + 1 / sigma_beta_sq
    as.numeric(w * solve_spd(G, w * Xty))
  }

  acvb_T <- function(w, j, value) {
    wj <- w
    wj[j] <- value
    beta_j <- acvb_beta_given_w(wj)

    scale_term <- B_y + 0.5 * yTy - 0.5 * sum(Xty * beta_j)
    scale_term <- max(scale_term, 1e-12)

    out <- -(A_y + n / 2) * log(scale_term)

    if (value == 1) {
      out <- out + logit_pi - 0.5 * log(n) - 0.5 * log(sigma_beta_sq)
    }

    out
  }

  n_iter <- max_iter

  for (t in seq_len(max_iter)) {
    omega_old <- omega
    omega_prop <- omega

    for (j in seq_len(p)) {
      T0 <- acvb_T(omega, j, 0)
      T1 <- acvb_T(omega, j, 1)
      omega_prop[j] <- expit(clamp(T1 - T0, -30, 30))
    }

    omega <- clamp(omega_prop, 1e-6, 1 - 1e-6)

    d_omega <- max(abs(omega - omega_old))
    if (verbose) {
      cat(sprintf("[ACVB/%s] iter %3d | d_omega=%.3e | size05=%d\n",
                  init_method, t, d_omega, sum(omega > 0.5)))
    }

    if (t > 2 && d_omega < tol) {
      n_iter <- t
      break
    }
  }

  beta_raw <- acvb_beta_given_w(omega)
  sel <- which(omega > 0.5)
  beta_hat <- rep(0, p)

  if (length(sel) > 0) {
    Xs <- X[, sel, drop = FALSE]
    beta_hat[sel] <- as.numeric(solve_reg(crossprod(Xs), crossprod(Xs, y)))
  }

  sel_05 <- get_selection_metrics_omega(omega, dat$beta_true, thr = 0.5)
  sel_01 <- get_selection_metrics_omega(omega, dat$beta_true, thr = 0.1)
  sel_top <- get_selection_metrics_topk(omega, dat$beta_true)

  list(
    method = paste0("no_borrow_ACVB_init_", init_method),
    omega = omega,
    beta_raw = beta_raw,
    beta_hat = beta_hat,
    pred_mse = prediction_mse(beta_hat, dat),
    sel_05 = sel_05,
    sel_01 = sel_01,
    sel_top = sel_top,
    mean_omega_signal = mean(omega[dat$tag_true == 1]),
    mean_omega_null = mean(omega[dat$tag_true == 0]),
    pi_prior = pi_prior,
    init_method = init_method,
    n_iter = n_iter
  )
}

## =========================================================
## Classical target-only baseline
## =========================================================

fit_noborrow_classical <- function(dat,
                                   method = c("auto", "ols", "lasso", "mcp", "scad"),
                                   gamma_ebic = 0.5) {
  method <- match.arg(method)

  if (method == "auto") {
    method <- ifelse(dat$n > dat$p, "ols", "lasso")
  }

  if (method == "ols") {
    if (dat$n <= dat$p) {
      warning("n <= p: using ridge-stabilized least squares via solve_reg().")
    }

    beta_hat <- as.numeric(solve_reg(dat$XtX, dat$Xty))
    sel_beta <- get_selection_metrics_beta(beta_hat, dat$beta_true, thr = 0.1)

    return(list(
      method = "no_borrow_OLS",
      beta_hat = beta_hat,
      pred_mse = prediction_mse(beta_hat, dat),
      sel_beta = sel_beta,
      tp_beta = sel_beta$tp,
      fp_beta = sel_beta$fp,
      fn_beta = sel_beta$fn,
      acc_beta = sel_beta$acc
    ))
  }

  fit_penalized_baseline(dat, penalty = method, gamma_ebic = gamma_ebic)
}
