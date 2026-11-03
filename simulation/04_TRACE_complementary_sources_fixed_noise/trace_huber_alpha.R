## TRACE Huber initialization: draft_revised.tex and supplement_information.tex
## eqs. target-center-cvb, preliminary-alpha-qp-cvb, huber-threshold-cvb.
## Default tau=1/sqrt(raw ns), weight=omega; diverging is legacy sensitivity only.

.trace_project_simplex <- function(x) {
  x <- as.numeric(x)
  if (length(x) == 1L) return(1)
  u <- sort(x, decreasing = TRUE)
  cssv <- cumsum(u) - 1
  rho <- max(which(u - cssv / seq_along(u) > 0))
  theta <- cssv[rho] / rho
  pmax(x - theta, 0)
}

.trace_huber_loss <- function(residual, tau) {
  a <- abs(residual)
  ifelse(a <= tau, 0.5 * residual^2, tau * a - 0.5 * tau^2)
}

.resolve_trace_huber_tau <- function(dat,
                                     omega,
                                     huber_tau = NULL,
                                     huber_tau_rule = c("fixed", "diverging"),
                                     huber_tau_multiplier = 1) {
  huber_tau_rule <- match.arg(huber_tau_rule)
  if (!is.null(huber_tau)) {
    huber_tau <- as.numeric(huber_tau)[1]
    if (!is.finite(huber_tau) || huber_tau <= 0) {
      stop("huber_tau must be a finite positive number.")
    }
    return(list(tau = huber_tau, rule = "user", s_reference = NA_real_))
  }

  multiplier <- as.numeric(huber_tau_multiplier)[1]
  if (!is.finite(multiplier) || multiplier <= 0) {
    stop("huber_tau_multiplier must be a finite positive number.")
  }
  ns_eff <- max(as.numeric(if (!is.null(dat$n_vec)) mean(dat$n_vec) else dat$ns)[1], 1)
  ## Estimate s_n from the first target-only inclusion probabilities.  This
  ## avoids the oracle dat$nzc while retaining the support-size scale required
  ## by the theorem.  The automatic rule is the geometric-mean rate between
  ## sqrt(log(s_hat) / n_s) and sqrt(log(n_s) / n_s):
  ##
  ##   tau_n = {log(s_hat) log(n_s)}^(1/4) / sqrt(n_s).
  ##
  ## If log(s_hat) is asymptotically equivalent to log(s_n) and
  ## log(s_n) = o(log(n_s)), this satisfies all three stated rate conditions.
  omega_numeric <- as.numeric(omega)
  estimated_support <- sum(is.finite(omega_numeric) & omega_numeric > 0.5)
  estimated_support <- max(estimated_support, 2)
  log_ns <- log(max(ns_eff, exp(1)))
  rate <- if (identical(huber_tau_rule, "fixed")) {
    ## Latest manuscript uses tau proportional to raw ns^(-1/2).
    ## The diverging branch below is retained as an optional sensitivity rule.
    ns_eff^(-0.5)
  } else {
    (log(estimated_support) * log_ns)^(1 / 4) / sqrt(ns_eff)
  }
  list(tau = multiplier * rate, rule = huber_tau_rule,
       ## Retain the existing output-field name for CSV compatibility.
       s_reference = estimated_support)
}

.trace_huber_objective_gradient <- function(alpha, Bmat, beta_target,
                                            weight, tau) {
  residual <- beta_target - as.numeric(Bmat %*% alpha)
  score <- pmax(pmin(residual, tau), -tau)
  list(
    objective = sum(weight * .trace_huber_loss(residual, tau)),
    gradient = -as.numeric(crossprod(Bmat, weight * score)),
    residual = residual
  )
}

.solve_trace_huber_simplex <- function(alpha_start, Bmat, beta_target,
                                       weight, tau, tol = 1e-10,
                                       max_iter = 10000L) {
  alpha <- .trace_project_simplex(alpha_start)
  curvature_bound <- crossprod(Bmat, Bmat * weight)
  lipschitz <- tryCatch(
    max(eigen(curvature_bound, symmetric = TRUE, only.values = TRUE)$values),
    error = function(e) NA_real_
  )
  if (!is.finite(lipschitz) || lipschitz <= .Machine$double.eps) {
    val <- .trace_huber_objective_gradient(
      alpha, Bmat, beta_target, weight, tau
    )
    return(list(alpha = alpha, objective = val$objective, iterations = 0L,
                converged = TRUE, kkt_residual = 0,
                active_set = which(alpha > sqrt(tol)),
                solver = "constant_objective"))
  }

  step <- 1 / (1.001 * lipschitz)
  y <- alpha
  momentum <- 1
  converged <- FALSE
  iterations <- max_iter

  for (iter in seq_len(max_iter)) {
    current <- .trace_huber_objective_gradient(
      alpha, Bmat, beta_target, weight, tau
    )
    at_y <- .trace_huber_objective_gradient(y, Bmat, beta_target, weight, tau)
    candidate <- .trace_project_simplex(y - step * at_y$gradient)
    candidate_eval <- .trace_huber_objective_gradient(
      candidate, Bmat, beta_target, weight, tau
    )
    if (candidate_eval$objective > current$objective + 1e-14) {
      y <- alpha
      momentum <- 1
      candidate <- .trace_project_simplex(alpha - step * current$gradient)
      candidate_eval <- .trace_huber_objective_gradient(
        candidate, Bmat, beta_target, weight, tau
      )
    }

    delta <- max(abs(candidate - alpha))
    momentum_new <- 0.5 * (1 + sqrt(1 + 4 * momentum^2))
    y <- candidate + ((momentum - 1) / momentum_new) * (candidate - alpha)
    alpha <- candidate
    momentum <- momentum_new
    if (delta <= tol * (1 + max(abs(alpha)))) {
      converged <- TRUE
      iterations <- iter
      break
    }
  }

  final <- .trace_huber_objective_gradient(
    alpha, Bmat, beta_target, weight, tau
  )
  projected <- .trace_project_simplex(alpha - step * final$gradient)
  list(
    alpha = alpha,
    objective = final$objective,
    iterations = iterations,
    converged = converged,
    kkt_residual = max(abs(projected - alpha)),
    active_set = which(alpha > sqrt(tol)),
    solver = "monotone_projected_gradient_huber"
  )
}

update_alpha_preliminary_huber <- function(alpha_old, omega, dat,
                                           sigma_beta_sq = 100,
                                           huber_tau = NULL,
                                           huber_tau_rule = c("fixed", "diverging"),
                                           huber_tau_multiplier = 1,
                                           tol = 1e-10,
                                           max_iter = 10000L) {
  alpha_old <- .normalize_simplex(alpha_old, K = dat$K)
  terms <- .compute_alpha_preliminary_qp_terms(
    omega = omega, dat = dat, sigma_beta_sq = sigma_beta_sq
  )
  tau_info <- .resolve_trace_huber_tau(
    dat = dat, omega = omega, huber_tau = huber_tau,
    huber_tau_rule = huber_tau_rule,
    huber_tau_multiplier = huber_tau_multiplier
  )
  ns_eff <- max(as.numeric(if (!is.null(dat$n_vec)) mean(dat$n_vec) else dat$ns)[1], 1)
  # CVB_induction_fixed_alpha_revised eq. preliminary-alpha-qp-cvb.
  weight <- as.numeric(omega)
  Bmat <- as.matrix(dat$B)
  old_eval <- .trace_huber_objective_gradient(
    alpha_old, Bmat, terms$beta_target, weight, tau_info$tau
  )
  sol <- .solve_trace_huber_simplex(
    alpha_start = alpha_old, Bmat = Bmat,
    beta_target = terms$beta_target, weight = weight,
    tau = tau_info$tau, tol = tol, max_iter = max_iter
  )

  list(
    alpha = sol$alpha,
    accepted = TRUE,
    objective_old = old_eval$objective,
    objective_new = sol$objective,
    ## Compatibility aliases for the existing diagnostic table schema.
    qp_old = old_eval$objective,
    qp_new = sol$objective,
    beta_target = terms$beta_target,
    weight = weight,
    huber_tau = tau_info$tau,
    huber_tau_rule = tau_info$rule,
    huber_s_reference = tau_info$s_reference,
    huber_iterations = sol$iterations,
    huber_converged = sol$converged,
    huber_kkt_residual = sol$kkt_residual,
    active_set = sol$active_set,
    removed_set = setdiff(seq_len(dat$K), sol$active_set),
    solver = sol$solver,
    strategy = "preliminary_huber"
  )
}
