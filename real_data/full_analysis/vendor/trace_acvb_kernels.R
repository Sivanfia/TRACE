# Numerical TRACE / ACVB kernels extracted from the delivered simulation package.
# Only truth-based output metrics and oracle initializers were removed.
suppressPackageStartupMessages(library(glmnet))

`.acvb_solve_woodbury` <- function (omega, rhs, dat, sigma_beta_sq, ridge_n = 1e-10) 
{
    omega <- clamp(as.numeric(omega), 0, 1)
    rhs <- as.numeric(rhs)
    lambda <- 1/sigma_beta_sq
    d <- lambda + dat$x_norm2 * omega * (1 - omega)
    d <- pmax(d, 1e-12)
    d_inv <- 1/d
    U <- sweep(dat$X, 2, omega, `*`)
    U_scaled <- sweep(U, 2, sqrt(d_inv), `*`)
    K <- diag(dat$n) + tcrossprod(U_scaled)
    K <- K + diag(ridge_n, dat$n)
    d_rhs <- d_inv * rhs
    middle_rhs <- as.numeric(U %*% d_rhs)
    middle <- tryCatch(solve_spd(K, middle_rhs), error = function(e) solve(K, middle_rhs))
    as.numeric(d_rhs - d_inv * as.numeric(crossprod(U, middle)))
}

`.acvb_update_cached_inverse_legacy` <- function (omega, dat, sigma_beta_sq, A_y, B_y, logit_pi) 
{
    p <- dat$p
    n <- dat$n
    lambda <- 1/sigma_beta_sq
    G <- .build_acvb_G_cached(omega, dat, sigma_beta_sq)
    Ginv <- solve_spd(G, diag(p))
    omega_prop <- omega
    for (j in seq_len(p)) {
        idx <- setdiff(seq_len(p), j)
        w_m <- omega[idx]
        r_m <- w_m * dat$Xty[idx]
        c_m <- w_m * dat$XtX[idx, j]
        SOL <- cbind(.solve_leave_one_cached(Ginv, j, r_m), .solve_leave_one_cached(Ginv, j, c_m))
        sol_r <- SOL[, 1]
        sol_c <- SOL[, 2]
        quad0 <- sum(r_m * sol_r)
        D1 <- pmax(dat$x_norm2[j] + lambda - sum(c_m * sol_c), 1e-10)
        R1 <- dat$Xty[j] - sum(c_m * sol_r)
        quad1 <- quad0 + R1^2/D1
        scale0 <- pmax(B_y + 0.5 * dat$yTy - 0.5 * quad0, 1e-12)
        scale1 <- pmax(B_y + 0.5 * dat$yTy - 0.5 * quad1, 1e-12)
        T0 <- -(A_y + n/2) * safe_log(scale0)
        T1 <- -(A_y + n/2) * safe_log(scale1) + logit_pi - 0.5 * safe_log(n) - 0.5 * safe_log(sigma_beta_sq)
        omega_prop[j] <- expit(T1 - T0)
    }
    list(omega = clamp(omega_prop, .Machine$double.eps, 1 - .Machine$double.eps), min_downdate_denom = NA_real_, 
        n_denom_clamped = NA_integer_, n_system_dimension = p)
}

`.acvb_update_woodbury` <- function (omega, dat, sigma_beta_sq, A_y, B_y, logit_pi, denom_floor = 1e-10, ridge_n = 1e-10) 
{
    X <- dat$X
    y <- dat$y
    n <- dat$n
    p <- dat$p
    lambda <- 1/sigma_beta_sq
    omega <- clamp(as.numeric(omega), 0, 1)
    d <- lambda + dat$x_norm2 * omega * (1 - omega)
    d <- pmax(d, 1e-12)
    scale <- omega/sqrt(d)
    V <- sweep(X, 2, scale, `*`)
    K <- diag(n) + tcrossprod(V)
    K <- K + diag(ridge_n, n)
    rhs <- cbind(X, y)
    Ksol <- tryCatch(solve_spd(K, rhs), error = function(e) solve(K, rhs))
    Kinv_X <- Ksol[, seq_len(p), drop = FALSE]
    Kinv_y <- Ksol[, p + 1L]
    xKx <- colSums(X * Kinv_X)
    xKy <- as.numeric(crossprod(X, Kinv_y))
    yKy <- as.numeric(sum(y * Kinv_y))
    leverage <- scale^2 * xKx
    denom_raw <- 1 - leverage
    denom <- pmax(denom_raw, denom_floor)
    xA <- scale * xKx
    yA <- scale * xKy
    y_Kminus_y <- yKy + yA^2/denom
    x_Kminus_x <- xKx + xA^2/denom
    x_Kminus_y <- xKy + xA * yA/denom
    quad0 <- dat$yTy - y_Kminus_y
    quad0 <- pmin(pmax(quad0, 0), dat$yTy)
    D1 <- pmax(lambda + x_Kminus_x, 1e-10)
    R1 <- x_Kminus_y
    quad1 <- quad0 + R1^2/D1
    scale0 <- pmax(B_y + 0.5 * dat$yTy - 0.5 * quad0, 1e-12)
    scale1 <- pmax(B_y + 0.5 * dat$yTy - 0.5 * quad1, 1e-12)
    T0 <- -(A_y + n/2) * safe_log(scale0)
    T1 <- -(A_y + n/2) * safe_log(scale1) + logit_pi - 0.5 * safe_log(n) - 0.5 * safe_log(sigma_beta_sq)
    list(omega = clamp(expit(T1 - T0), .Machine$double.eps, 1 - .Machine$double.eps), min_downdate_denom = min(denom_raw), 
        n_denom_clamped = sum(denom_raw < denom_floor), n_system_dimension = n)
}

`.alpha_elbo_value` <- function (alpha, omega, gamma, dat, prior_par_borrow, sigma_beta_sq, A_y = 1, B_y = 1, include_source_normalization = TRUE) 
{
    borrow_elbo_surrogate(alpha = alpha, omega = omega, gamma = gamma, dat = dat, prior_par_borrow = prior_par_borrow, 
        sigma_beta_sq = sigma_beta_sq, A_y = A_y, B_y = B_y, include_source_normalization = include_source_normalization)$total
}

`.alpha_qp_objective` <- function (alpha, H, g) 
{
    alpha <- as.numeric(alpha)
    as.numeric(0.5 * sum(alpha * as.numeric(H %*% alpha)) - sum(g * alpha))
}

`.bernoulli_entropy_sum` <- function (q, eps = 1e-12) 
{
    q <- pmin(pmax(as.numeric(q), eps), 1 - eps)
    -sum(q * log(q) + (1 - q) * log(1 - q))
}

`.build_A0_cached` <- function (omega, gamma, dat, S0bar_inv) 
{
    ns_eff <- get_ns_eff(dat)
    S1bar_inv <- (ns_eff * dat$x_norm2)/dat$n
    S_inv_all <- S0bar_inv + (S1bar_inv - S0bar_inv) * gamma
    A0 <- dat$XtX * (omega %o% omega)
    diag(A0) <- dat$x_norm2 * omega + S_inv_all
    A0 + diag(1e-10, length(omega))
}

`.build_acvb_G_cached` <- function (omega, dat, sigma_beta_sq) 
{
    lambda <- 1/sigma_beta_sq
    G <- dat$XtX * (omega %o% omega)
    diag(G) <- dat$x_norm2 * omega + lambda
    G + diag(1e-10, length(omega))
}

`.compute_alpha_em_qp_terms` <- function (omega, gamma, dat, sigma_beta_sq = 100, ridge_M = 1e-10, ridge_H = 1e-10) 
{
    K <- dat$K
    p <- dat$p
    omega <- clamp(as.numeric(omega), 0, 1)
    gamma <- clamp(as.numeric(gamma), 0, 1)
    S0bar_inv <- 1/sigma_beta_sq
    ns_eff <- get_ns_eff(dat)
    S1bar_inv <- (ns_eff * dat$x_norm2)/dat$n
    S_inv <- S0bar_inv + (S1bar_inv - S0bar_inv) * gamma
    S_inv <- pmax(S_inv, 1e-12)
    M <- dat$XtX * (omega %o% omega)
    diag(M) <- dat$x_norm2 * omega + S_inv
    M <- M + diag(ridge_M, p)
    Bmat <- as.matrix(dat$B)
    G1B <- Bmat * as.numeric(gamma * S1bar_inv)
    target_rhs <- as.numeric(omega * dat$Xty)
    factor <- .trace_woodbury_factor(z = omega, diag_add = S_inv + dat$x_norm2 * omega * (1 - omega), dat = dat, 
        ridge = ridge_M)
    solved <- .trace_woodbury_apply(factor, cbind(G1B, target_rhs))
    Minv_G1B <- solved[, seq_len(K), drop = FALSE]
    Minv_target <- solved[, K + 1L]
    H <- crossprod(Bmat, G1B) - crossprod(G1B, Minv_G1B)
    g <- as.numeric(crossprod(G1B, Minv_target))
    yPy <- as.numeric(dat$yTy - sum(target_rhs * Minv_target))
    if (!is.finite(yPy)) 
        yPy <- NA_real_
    H <- (H + t(H))/2
    H[!is.finite(H)] <- 0
    g[!is.finite(g)] <- 0
    ev <- tryCatch(eigen(H, symmetric = TRUE, only.values = TRUE)$values, error = function(e) NA_real_)
    min_ev <- suppressWarnings(min(ev, na.rm = TRUE))
    if (!is.finite(min_ev)) 
        min_ev <- 0
    if (min_ev < ridge_H) {
        H <- H + diag(abs(min_ev) + ridge_H, K)
    }
    else {
        H <- H + diag(ridge_H, K)
    }
    list(H = H, g = g, yPy = yPy, M = M, S_inv = S_inv)
}

`.compute_alpha_preliminary_qp_terms` <- function (omega, dat, sigma_beta_sq = 100, ridge_M = 1e-10, ridge_H = 1e-10) 
{
    omega <- clamp(as.numeric(omega), 0, 1)
    p <- dat$p
    K <- dat$K
    if (length(omega) != p) 
        stop("omega has wrong length in preliminary alpha QP.")
    S0bar_inv <- 1/sigma_beta_sq
    ns_eff <- get_ns_eff(dat)
    S1bar_inv <- (ns_eff * dat$x_norm2)/dat$n
    S1bar_inv <- pmax(S1bar_inv, 1e-12)
    M0 <- dat$XtX * (omega %o% omega)
    diag(M0) <- dat$x_norm2 * omega + S0bar_inv
    M0 <- M0 + diag(ridge_M, p)
    target_rhs <- omega * dat$Xty
    factor <- .trace_woodbury_factor(z = omega, diag_add = rep(S0bar_inv, p) + dat$x_norm2 * omega * (1 - omega), 
        dat = dat, ridge = ridge_M)
    beta_target <- .trace_woodbury_apply(factor, target_rhs)
    Bmat <- as.matrix(dat$B)
    weight <- omega * S1bar_inv
    H <- crossprod(Bmat, Bmat * weight)
    g <- as.numeric(crossprod(Bmat, weight * beta_target))
    H <- (H + t(H))/2
    H[!is.finite(H)] <- 0
    g[!is.finite(g)] <- 0
    ev <- tryCatch(eigen(H, symmetric = TRUE, only.values = TRUE)$values, error = function(e) NA_real_)
    min_ev <- suppressWarnings(min(ev, na.rm = TRUE))
    if (!is.finite(min_ev)) 
        min_ev <- 0
    if (min_ev < ridge_H) {
        H <- H + diag(abs(min_ev) + ridge_H, K)
    }
    else {
        H <- H + diag(ridge_H, K)
    }
    list(H = H, g = g, M0 = M0, beta_target = beta_target, S1bar_inv = S1bar_inv, weight = weight)
}

`.inverse_rank1_coordinate_update` <- function (Minv, j, delta_diag, denom_floor = 1e-10) 
{
    if (!is.finite(delta_diag) || abs(delta_diag) < 1e-15) {
        return(list(inverse = Minv, ok = TRUE, min_abs_denom = Inf))
    }
    mj <- Minv[, j]
    denom <- 1 + delta_diag * Minv[j, j]
    if (!is.finite(denom) || abs(denom) < denom_floor) {
        return(list(inverse = Minv, ok = FALSE, min_abs_denom = abs(denom)))
    }
    out <- Minv - (delta_diag/denom) * tcrossprod(mj)
    out <- 0.5 * (out + t(out))
    list(inverse = out, ok = all(is.finite(out)), min_abs_denom = abs(denom))
}

`.inverse_rank2_omega_update` <- function (Minv, j, delta_omega, omega_before, dat, denom_floor = 1e-10) 
{
    if (!is.finite(delta_omega) || abs(delta_omega) < 1e-15) {
        return(list(inverse = Minv, ok = TRUE, min_abs_denom = Inf))
    }
    p <- length(omega_before)
    ej <- numeric(p)
    ej[j] <- 1
    qj <- omega_before * dat$XtX[, j]
    qj[j] <- 0
    U <- cbind(ej, qj)
    Cmat <- matrix(c(delta_omega * dat$x_norm2[j], delta_omega, delta_omega, 0), nrow = 2L, byrow = TRUE)
    MU <- Minv %*% U
    small_system <- diag(2L) + Cmat %*% crossprod(U, MU)
    small_det <- as.numeric(det(small_system))
    if (!is.finite(small_det) || abs(small_det) < denom_floor) {
        return(list(inverse = Minv, ok = FALSE, min_abs_denom = abs(small_det)))
    }
    correction <- tryCatch(solve(small_system, Cmat), error = function(e) NULL)
    if (is.null(correction) || any(!is.finite(correction))) {
        return(list(inverse = Minv, ok = FALSE, min_abs_denom = abs(small_det)))
    }
    out <- Minv - MU %*% correction %*% t(MU)
    out <- 0.5 * (out + t(out))
    list(inverse = out, ok = all(is.finite(out)), min_abs_denom = abs(small_det))
}

`.make_elbo_row` <- function (iter, stage, alpha_method, init_method, alpha, omega, gamma, dat, prior_par_borrow, sigma_beta_sq, A_y, 
    B_y, d_omega = NA_real_, d_gamma = NA_real_, d_alpha = NA_real_, elbo_prev = NA_real_, include_source_normalization = TRUE) 
{
    e <- borrow_elbo_surrogate(alpha = alpha, omega = omega, gamma = gamma, dat = dat, prior_par_borrow = prior_par_borrow, 
        sigma_beta_sq = sigma_beta_sq, A_y = A_y, B_y = B_y, include_source_normalization = include_source_normalization)
    delta <- if (is.finite(elbo_prev)) 
        e$total - elbo_prev
    else NA_real_
    rel <- if (is.finite(delta)) 
        abs(delta)/(1 + abs(elbo_prev))
    else NA_real_
    param_change_values <- c(d_omega, d_gamma, d_alpha)
    param_change_values <- param_change_values[is.finite(param_change_values)]
    max_param_change_value <- if (length(param_change_values) > 0L) {
        max(param_change_values)
    }
    else {
        NA_real_
    }
    data.frame(iter = iter, stage = stage, objective_name = e$monitor_name, alpha_method = alpha_method, init_method = init_method, 
        elbo = e$total, delta_elbo = delta, rel_delta_elbo = rel, collapsed = e$collapsed, source_norm = e$source_norm, 
        indicator_prior = e$indicator_prior, entropy = e$entropy, logdet = e$logdet, Tval = e$Tval, d_omega = d_omega, 
        d_gamma = d_gamma, d_alpha = d_alpha, max_param_change = max_param_change_value, size05 = sum(omega > 0.5), 
        alpha_1 = alpha[1], alpha_2 = if (length(alpha) >= 2) 
            alpha[2]
        else NA_real_, alpha_3 = if (length(alpha) >= 3) 
            alpha[3]
        else NA_real_, stringsAsFactors = FALSE)
}

`.normalize_simplex` <- function (alpha, K = NULL) 
{
    alpha <- as.numeric(alpha)
    if (!is.null(K) && length(alpha) != K) 
        stop("alpha has wrong length.")
    alpha[!is.finite(alpha)] <- 0
    alpha <- pmax(alpha, 0)
    if (sum(alpha) <= 0) 
        alpha <- rep(1/length(alpha), length(alpha))
    alpha/sum(alpha)
}

`.resolve_trace_huber_tau` <- function (dat, omega, huber_tau = NULL, huber_tau_rule = c("fixed", "diverging"), huber_tau_multiplier = 1) 
{
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
    ns_eff <- max(as.numeric(if (!is.null(dat[["n_vec", exact = TRUE]])) mean(dat[["n_vec", exact = TRUE]]) else dat$ns)[1], 1)
    omega_numeric <- as.numeric(omega)
    estimated_support <- sum(is.finite(omega_numeric) & omega_numeric > 0.5)
    estimated_support <- max(estimated_support, 2)
    log_ns <- log(max(ns_eff, exp(1)))
    rate <- if (identical(huber_tau_rule, "fixed")) {
        ns_eff^(-0.5)
    }
    else {
        (log(estimated_support) * log_ns)^(1/4)/sqrt(ns_eff)
    }
    list(tau = multiplier * rate, rule = huber_tau_rule, s_reference = estimated_support)
}

`.solve_leave_one_cached` <- function (Minv, j, a) 
{
    p <- nrow(Minv)
    if (p <= 1) 
        return(numeric(0))
    idx <- setdiff(seq_len(p), j)
    g <- Minv[j, j]
    if (!is.finite(g) || abs(g) < 1e-12) {
        stop("Cached inverse is numerically singular at coordinate ", j)
    }
    f <- Minv[idx, j]
    as.numeric(Minv[idx, idx, drop = FALSE] %*% a - f * (sum(f * a)/g))
}

`.solve_simplex_qp_active_set` <- function (H, g, tol = 1e-10, ridge = 1e-10, max_iter = NULL) 
{
    H <- as.matrix(H)
    g <- as.numeric(g)
    K <- length(g)
    if (is.null(max_iter)) 
        max_iter <- K + 5L
    H <- (H + t(H))/2
    H <- H + diag(ridge, K)
    active <- seq_len(K)
    removed <- integer(0)
    lambda <- NA_real_
    alpha_active <- rep(NA_real_, K)
    for (iter in seq_len(max_iter)) {
        m <- length(active)
        if (m == 0L) {
            active <- which.max(g)
            alpha_active <- 1
            lambda <- as.numeric(g[active] - H[active, active])
            break
        }
        if (m == 1L) {
            alpha_active <- 1
            lambda <- as.numeric(g[active] - H[active, active])
        }
        else {
            Haa <- H[active, active, drop = FALSE]
            ga <- g[active]
            one <- rep(1, m)
            sol_g <- tryCatch(solve(Haa, ga), error = function(e) solve(Haa + diag(ridge, m), ga))
            sol_1 <- tryCatch(solve(Haa, one), error = function(e) solve(Haa + diag(ridge, m), one))
            denom <- as.numeric(sum(sol_1))
            if (!is.finite(denom) || abs(denom) < 1e-12) 
                denom <- sign(denom + 1e-12) * 1e-12
            lambda <- as.numeric((sum(sol_g) - 1)/denom)
            alpha_active <- as.numeric(sol_g - lambda * sol_1)
        }
        if (all(alpha_active > -tol) || length(active) == 1L) {
            break
        }
        worst_local <- which.min(alpha_active)
        removed <- c(removed, active[worst_local])
        active <- active[-worst_local]
    }
    alpha <- rep(0, K)
    alpha[active] <- pmax(alpha_active, 0)
    alpha <- .normalize_simplex(alpha, K = K)
    grad <- as.numeric(H %*% alpha - g)
    mu <- grad + lambda
    list(alpha = alpha, active_set = active, removed_set = removed, lambda = lambda, mu = mu, objective = .alpha_qp_objective(alpha, 
        H, g), kkt_min_mu_inactive = if (length(setdiff(seq_len(K), active)) > 0) {
        min(mu[setdiff(seq_len(K), active)])
    } else {
        NA_real_
    })
}

`.solve_simplex_qp_global` <- function (H, g, tol = 1e-10, ridge = 1e-10, max_exact_K = 12L) 
{
    H <- as.matrix(H)
    g <- as.numeric(g)
    K <- length(g)
    H <- (H + t(H))/2
    if (K > max_exact_K) {
        out <- .solve_simplex_qp_active_set(H, g, tol = tol, ridge = ridge)
        out$solver <- "active_set_fallback"
        return(out)
    }
    best <- NULL
    all_index <- seq_len(K)
    for (m in seq_len(K)) {
        faces <- combn(all_index, m, simplify = FALSE)
        for (active in faces) {
            Haa <- H[active, active, drop = FALSE]
            ga <- g[active]
            one <- rep(1, m)
            if (m == 1L) {
                alpha_active <- 1
            }
            else {
                sol_g <- tryCatch(solve(Haa, ga), error = function(e) solve(Haa + diag(ridge, m), ga))
                sol_1 <- tryCatch(solve(Haa, one), error = function(e) solve(Haa + diag(ridge, m), one))
                denom <- as.numeric(sum(sol_1))
                if (!is.finite(denom) || abs(denom) < 1e-12) 
                  next
                lambda_face <- as.numeric((sum(sol_g) - 1)/denom)
                alpha_active <- as.numeric(sol_g - lambda_face * sol_1)
            }
            if (any(!is.finite(alpha_active)) || any(alpha_active < -tol)) 
                next
            alpha <- rep(0, K)
            alpha[active] <- pmax(alpha_active, 0)
            if (sum(alpha) <= 0) 
                next
            alpha <- alpha/sum(alpha)
            grad <- as.numeric(H %*% alpha - g)
            lambda <- -mean(grad[active])
            inactive <- setdiff(all_index, active)
            mu <- grad + lambda
            if (length(inactive) > 0 && any(mu[inactive] < -sqrt(tol))) 
                next
            if (max(abs(grad[active] + lambda)) > 10 * sqrt(tol)) 
                next
            objective <- .alpha_qp_objective(alpha, H, g)
            if (is.null(best) || objective < best$objective) {
                best <- list(alpha = alpha, active_set = active, removed_set = inactive, lambda = lambda, mu = mu, 
                  objective = objective, kkt_min_mu_inactive = if (length(inactive) > 0) min(mu[inactive]) else NA_real_, 
                  solver = "exact_face_enumeration")
            }
        }
    }
    if (is.null(best)) {
        best <- .solve_simplex_qp_active_set(H, g, tol = tol, ridge = ridge)
        best$solver <- "active_set_emergency_fallback"
    }
    best
}

`.solve_trace_huber_simplex` <- function (alpha_start, Bmat, beta_target, weight, tau, tol = 1e-10, max_iter = 10000L) 
{
    alpha <- .trace_project_simplex(alpha_start)
    curvature_bound <- crossprod(Bmat, Bmat * weight)
    lipschitz <- tryCatch(max(eigen(curvature_bound, symmetric = TRUE, only.values = TRUE)$values), error = function(e) NA_real_)
    if (!is.finite(lipschitz) || lipschitz <= .Machine$double.eps) {
        val <- .trace_huber_objective_gradient(alpha, Bmat, beta_target, weight, tau)
        return(list(alpha = alpha, objective = val$objective, iterations = 0L, converged = TRUE, kkt_residual = 0, 
            active_set = which(alpha > sqrt(tol)), solver = "constant_objective"))
    }
    step <- 1/(1.001 * lipschitz)
    y <- alpha
    momentum <- 1
    converged <- FALSE
    iterations <- max_iter
    for (iter in seq_len(max_iter)) {
        current <- .trace_huber_objective_gradient(alpha, Bmat, beta_target, weight, tau)
        at_y <- .trace_huber_objective_gradient(y, Bmat, beta_target, weight, tau)
        candidate <- .trace_project_simplex(y - step * at_y$gradient)
        candidate_eval <- .trace_huber_objective_gradient(candidate, Bmat, beta_target, weight, tau)
        if (candidate_eval$objective > current$objective + 1e-14) {
            y <- alpha
            momentum <- 1
            candidate <- .trace_project_simplex(alpha - step * current$gradient)
            candidate_eval <- .trace_huber_objective_gradient(candidate, Bmat, beta_target, weight, tau)
        }
        delta <- max(abs(candidate - alpha))
        momentum_new <- 0.5 * (1 + sqrt(1 + 4 * momentum^2))
        y <- candidate + ((momentum - 1)/momentum_new) * (candidate - alpha)
        alpha <- candidate
        momentum <- momentum_new
        if (delta <= tol * (1 + max(abs(alpha)))) {
            converged <- TRUE
            iterations <- iter
            break
        }
    }
    final <- .trace_huber_objective_gradient(alpha, Bmat, beta_target, weight, tau)
    projected <- .trace_project_simplex(alpha - step * final$gradient)
    list(alpha = alpha, objective = final$objective, iterations = iterations, converged = converged, kkt_residual = max(abs(projected - 
        alpha)), active_set = which(alpha > sqrt(tol)), solver = "monotone_projected_gradient_huber")
}

`.trace_get_ns_eff_core` <- function (dat) 
{
    if (!is.null(dat[["n_vec", exact = TRUE]])) 
        return(mean(dat[["n_vec", exact = TRUE]]))
    dat$ns %||% dat$n
}

`.trace_huber_loss` <- function (residual, tau) 
{
    a <- abs(residual)
    ifelse(a <= tau, 0.5 * residual^2, tau * a - 0.5 * tau^2)
}

`.trace_huber_objective_gradient` <- function (alpha, Bmat, beta_target, weight, tau) 
{
    residual <- beta_target - as.numeric(Bmat %*% alpha)
    score <- pmax(pmin(residual, tau), -tau)
    list(objective = sum(weight * .trace_huber_loss(residual, tau)), gradient = -as.numeric(crossprod(Bmat, weight * 
        score)), residual = residual)
}

`.trace_leave_one_quad` <- function (full_quad, Ha, Hb, h_diag) 
{
    as.numeric(full_quad - Ha * Hb/h_diag)
}

`.trace_project_simplex` <- function (x) 
{
    x <- as.numeric(x)
    if (length(x) == 1L) 
        return(1)
    u <- sort(x, decreasing = TRUE)
    cssv <- cumsum(u) - 1
    rho <- max(which(u - cssv/seq_along(u) > 0))
    theta <- cssv[rho]/rho
    pmax(x - theta, 0)
}

`.trace_woodbury_apply` <- function (factor, rhs) 
{
    rhs_was_vector <- is.null(dim(rhs))
    rhs <- as.matrix(rhs)
    d_rhs <- sweep(rhs, 1L, factor$d_inv, `*`)
    middle <- if (!is.null(factor$Kinv_U)) {
        factor$Kinv_U %*% d_rhs
    }
    else {
        factor$solve_K(factor$U %*% d_rhs)
    }
    correction <- crossprod(factor$U, middle)
    out <- d_rhs - sweep(correction, 1L, factor$d_inv, `*`)
    if (rhs_was_vector) 
        as.numeric(out)
    else out
}

`.trace_woodbury_factor` <- function (z, diag_add, dat, ridge = 1e-10, compute_inverse_diagonal = FALSE) 
{
    z <- as.numeric(z)
    diag_add <- as.numeric(diag_add)
    if (length(z) != dat$p || length(diag_add) != dat$p) {
        stop("z and diag_add must both have length dat$p.")
    }
    d <- pmax(diag_add + ridge, 1e-14)
    d_inv <- 1/d
    U <- sweep(dat$X, 2L, z, `*`)
    U_scaled <- sweep(U, 2L, sqrt(d_inv), `*`)
    K <- diag(dat$n) + tcrossprod(U_scaled)
    chol_K <- tryCatch(chol(K), error = function(e) NULL)
    solve_K <- function(rhs) {
        rhs <- as.matrix(rhs)
        if (!is.null(chol_K)) {
            backsolve(chol_K, forwardsolve(t(chol_K), rhs))
        }
        else {
            solve_spd(K, rhs)
        }
    }
    Kinv_U <- NULL
    h_diag <- NULL
    if (isTRUE(compute_inverse_diagonal)) {
        Kinv_U <- solve_K(U)
        h_diag <- d_inv - d_inv^2 * colSums(U * Kinv_U)
        h_diag <- pmax(h_diag, 1e-14)
    }
    list(U = U, d = d, d_inv = d_inv, Kinv_U = Kinv_U, h_diag = h_diag, solve_K = solve_K, logdet = sum(log(d)) + 
        if (!is.null(chol_K)) {
            2 * sum(log(diag(chol_K)))
        } else {
            as.numeric(determinant(K, logarithm = TRUE)$modulus)
        }, n_system_dimension = dat$n)
}

`.update_gamma_all_cached` <- function (omega, gamma, bbar, dat, S0bar_inv, A_y, B_y, E_logit_rho_c) 
{
    p <- length(gamma)
    omega <- clamp(as.numeric(omega), 0, 1)
    gamma <- clamp(as.numeric(gamma), 0, 1)
    ns_eff <- get_ns_eff(dat)
    S1bar_inv <- (ns_eff * dat$x_norm2)/dat$n
    S_inv <- S0bar_inv + (S1bar_inv - S0bar_inv) * gamma
    v <- gamma * S1bar_inv * bbar
    r <- omega * dat$Xty + v
    Cmat <- sweep(dat$XtX, 1L, omega, `*`)
    factor <- .trace_woodbury_factor(z = omega, diag_add = S_inv + dat$x_norm2 * omega * (1 - omega), dat = dat, 
        compute_inverse_diagonal = TRUE)
    solved <- .trace_woodbury_apply(factor, cbind(r, Cmat))
    Hr <- solved[, 1L]
    HC <- solved[, 1L + seq_len(p), drop = FALSE]
    h <- factor$h_diag
    HC_diag <- diag(HC)
    rr_minus <- .trace_leave_one_quad(sum(r * Hr), Hr, Hr, h)
    cc_minus <- colSums(Cmat * HC) - HC_diag^2/h
    cr_minus <- as.numeric(crossprod(Cmat, Hr)) - HC_diag * Hr/h
    D0 <- pmax(S0bar_inv + omega * dat$x_norm2 - omega^2 * cc_minus, 1e-10)
    D1 <- pmax(S1bar_inv + omega * dat$x_norm2 - omega^2 * cc_minus, 1e-10)
    R0 <- omega * (dat$Xty - cr_minus)
    R1 <- R0 + S1bar_inv * bbar
    alpha_minus <- sum(gamma * S1bar_inv * bbar^2) - gamma * S1bar_inv * bbar^2
    Tc0 <- pmax(dat$yTy - rr_minus - R0^2/D0 + alpha_minus, 1e-10)
    Tc1 <- pmax(dat$yTy - rr_minus - R1^2/D1 + alpha_minus + S1bar_inv * bbar^2, 1e-10)
    Delta <- 0.5 * safe_log(S1bar_inv/S0bar_inv) - 0.5 * safe_log(1 + (S1bar_inv - S0bar_inv)/D0) - (A_y + dat$n/2) * 
        (safe_log(B_y + 0.5 * Tc1) - safe_log(B_y + 0.5 * Tc0))
    clamp(expit(E_logit_rho_c + Delta), .Machine$double.eps, 1 - .Machine$double.eps)
}

`.update_gamma_all_gauss_seidel` <- function (omega, gamma, bbar, dat, S0bar_inv, A_y, B_y, E_logit_rho_c) 
{
    p <- length(gamma)
    if (p <= 1L) 
        return(clamp(gamma, .Machine$double.eps, 1 - .Machine$double.eps))
    gamma <- clamp(as.numeric(gamma), .Machine$double.eps, 1 - .Machine$double.eps)
    ns_eff <- get_ns_eff(dat)
    S1bar_inv <- (ns_eff * dat$x_norm2)/dat$n
    S_inv_all <- S0bar_inv + (S1bar_inv - S0bar_inv) * gamma
    v_all <- gamma * S1bar_inv * bbar
    rebuild_inverse <- function() {
        solve_spd(.build_A0_cached(omega, gamma, dat, S0bar_inv), diag(p))
    }
    Minv <- rebuild_inverse()
    for (j in seq_len(p)) {
        idx <- seq_len(p) != j
        w_m <- omega[idx]
        wj <- omega[j]
        S_inv_m <- S_inv_all[idx]
        v_m <- v_all[idx]
        r_minus <- w_m * dat$Xty[idx] + v_m
        ccol <- w_m * dat$XtX[idx, j]
        crow <- dat$XtX[j, idx] * w_m
        SOL <- tryCatch(cbind(.solve_leave_one_cached(Minv, j, r_minus), .solve_leave_one_cached(Minv, j, ccol)), 
            error = function(e) NULL)
        if (is.null(SOL) || any(!is.finite(SOL))) {
            Minv <- rebuild_inverse()
            SOL <- cbind(.solve_leave_one_cached(Minv, j, r_minus), .solve_leave_one_cached(Minv, j, ccol))
        }
        sol_r <- SOL[, 1]
        sol_c <- SOL[, 2]
        cTc <- sum(ccol * sol_c)
        D0j <- pmax(S0bar_inv + wj * dat$x_norm2[j] - wj^2 * cTc, 1e-10)
        D1j <- pmax(S1bar_inv[j] + wj * dat$x_norm2[j] - wj^2 * cTc, 1e-10)
        R0 <- wj * (dat$Xty[j] - sum(crow * sol_r))
        v1j <- S1bar_inv[j] * bbar[j]
        R1 <- R0 + v1j
        alpha_minus <- sum(gamma[idx] * S1bar_inv[idx] * bbar[idx]^2)
        Tc0 <- pmax(dat$yTy - sum(r_minus * sol_r) - (R0^2)/D0j + alpha_minus, 1e-10)
        Tc1 <- pmax(dat$yTy - sum(r_minus * sol_r) - (R1^2)/D1j + alpha_minus + S1bar_inv[j] * bbar[j]^2, 1e-10)
        Delta <- 0.5 * safe_log(S1bar_inv[j]/S0bar_inv) - 0.5 * safe_log(1 + (S1bar_inv[j] - S0bar_inv)/D0j) - 
            (A_y + dat$n/2) * (safe_log(B_y + 0.5 * Tc1) - safe_log(B_y + 0.5 * Tc0))
        proposal <- clamp(expit(E_logit_rho_c + Delta), 1e-06, 1 - 1e-06)
        gamma_old_j <- gamma[j]
        S_old_j <- S_inv_all[j]
        gamma[j] <- proposal
        S_inv_all[j] <- S0bar_inv + (S1bar_inv[j] - S0bar_inv) * proposal
        v_all[j] <- S1bar_inv[j] * proposal * bbar[j]
        inv_update <- .inverse_rank1_coordinate_update(Minv = Minv, j = j, delta_diag = S_inv_all[j] - S_old_j)
        Minv <- if (isTRUE(inv_update$ok)) 
            inv_update$inverse
        else rebuild_inverse()
    }
    clamp(gamma, .Machine$double.eps, 1 - .Machine$double.eps)
}

`.update_omega_all_cached` <- function (omega, gamma, bbar, dat, S0bar_inv, A_y, B_y, E_logit_rho_z) 
{
    p <- length(omega)
    omega <- clamp(as.numeric(omega), 0, 1)
    gamma <- clamp(as.numeric(gamma), 0, 1)
    ns_eff <- get_ns_eff(dat)
    S1bar_inv <- (ns_eff * dat$x_norm2)/dat$n
    S_inv <- S0bar_inv + (S1bar_inv - S0bar_inv) * gamma
    u <- omega * dat$Xty
    v <- gamma * S1bar_inv * bbar
    Cmat <- sweep(dat$XtX, 1L, omega, `*`)
    factor <- .trace_woodbury_factor(z = omega, diag_add = S_inv + dat$x_norm2 * omega * (1 - omega), dat = dat, 
        compute_inverse_diagonal = TRUE)
    solved <- .trace_woodbury_apply(factor, cbind(u, v, Cmat))
    Hu <- solved[, 1L]
    Hv <- solved[, 2L]
    HC <- solved[, 2L + seq_len(p), drop = FALSE]
    h <- factor$h_diag
    uu_minus <- .trace_leave_one_quad(sum(u * Hu), Hu, Hu, h)
    uv_minus <- .trace_leave_one_quad(sum(u * Hv), Hu, Hv, h)
    vv_minus <- .trace_leave_one_quad(sum(v * Hv), Hv, Hv, h)
    HC_diag <- diag(HC)
    cc_minus <- colSums(Cmat * HC) - HC_diag^2/h
    cu_minus <- as.numeric(crossprod(Cmat, Hu)) - HC_diag * Hu/h
    cv_minus <- as.numeric(crossprod(Cmat, Hv)) - HC_diag * Hv/h
    D <- pmax(S_inv + dat$x_norm2 - cc_minus, 1e-10)
    A <- dat$Xty - cu_minus
    B <- v - cv_minus
    alpha_term <- sum(gamma * S1bar_inv * bbar^2)
    T0 <- dat$yTy - uu_minus - 2 * uv_minus - vv_minus - (gamma * S1bar_inv * bbar^2) + alpha_term
    T1 <- T0 - (A + B)^2/D + (gamma * S1bar_inv * bbar^2)
    T0 <- pmax(T0, 1e-10)
    T1 <- pmax(T1, 1e-10)
    lin <- E_logit_rho_z - 0.5 * safe_log(D) + 0.5 * safe_log(S_inv) - (A_y + dat$n/2) * (safe_log(B_y + 0.5 * 
        T1) - safe_log(B_y + 0.5 * T0))
    clamp(expit(lin), .Machine$double.eps, 1 - .Machine$double.eps)
}

`.update_omega_all_gauss_seidel` <- function (omega, gamma, bbar, dat, S0bar_inv, A_y, B_y, E_logit_rho_z) 
{
    p <- length(omega)
    if (p <= 1L) 
        return(clamp(omega, .Machine$double.eps, 1 - .Machine$double.eps))
    omega <- clamp(as.numeric(omega), .Machine$double.eps, 1 - .Machine$double.eps)
    ns_eff <- get_ns_eff(dat)
    S1bar_inv <- (ns_eff * dat$x_norm2)/dat$n
    S_inv_all <- S0bar_inv + (S1bar_inv - S0bar_inv) * gamma
    v_all <- gamma * S1bar_inv * bbar
    alpha_term <- sum(gamma * S1bar_inv * bbar^2)
    rebuild_inverse <- function() {
        solve_spd(.build_A0_cached(omega, gamma, dat, S0bar_inv), diag(p))
    }
    Minv <- rebuild_inverse()
    for (j in seq_len(p)) {
        idx <- seq_len(p) != j
        w_m <- omega[idx]
        S_inv_j <- S_inv_all[j]
        S_inv_m <- S_inv_all[idx]
        v_j <- v_all[j]
        v_m <- v_all[idx]
        u <- w_m * dat$Xty[idx]
        ccol <- w_m * dat$XtX[idx, j]
        SOL <- tryCatch(cbind(.solve_leave_one_cached(Minv, j, u), .solve_leave_one_cached(Minv, j, v_m), .solve_leave_one_cached(Minv, 
            j, ccol)), error = function(e) NULL)
        if (is.null(SOL) || any(!is.finite(SOL))) {
            Minv <- rebuild_inverse()
            SOL <- cbind(.solve_leave_one_cached(Minv, j, u), .solve_leave_one_cached(Minv, j, v_m), .solve_leave_one_cached(Minv, 
                j, ccol))
        }
        sol_u <- SOL[, 1]
        sol_v <- SOL[, 2]
        sol_c <- SOL[, 3]
        D_j <- pmax(S_inv_j + dat$x_norm2[j] - sum(ccol * sol_c), 1e-10)
        A_j <- dat$Xty[j] - sum(ccol * sol_u)
        B_j <- v_j - sum(ccol * sol_v)
        T0 <- dat$yTy - sum(u * sol_u) - 2 * sum(u * sol_v) - sum(v_m * sol_v) - (gamma[j] * S1bar_inv[j] * bbar[j]^2) + 
            alpha_term
        T1 <- T0 - (A_j + B_j)^2/D_j + (gamma[j] * S1bar_inv[j] * bbar[j]^2)
        T0 <- pmax(T0, 1e-10)
        T1 <- pmax(T1, 1e-10)
        lin <- E_logit_rho_z - 0.5 * safe_log(D_j) + 0.5 * safe_log(S_inv_j) - (A_y + dat$n/2) * (safe_log(B_y + 
            0.5 * T1) - safe_log(B_y + 0.5 * T0))
        proposal <- clamp(expit(lin), .Machine$double.eps, 1 - .Machine$double.eps)
        omega_before <- omega
        delta_omega <- proposal - omega_before[j]
        omega[j] <- proposal
        inv_update <- .inverse_rank2_omega_update(Minv = Minv, j = j, delta_omega = delta_omega, omega_before = omega_before, 
            dat = dat)
        Minv <- if (isTRUE(inv_update$ok)) 
            inv_update$inverse
        else rebuild_inverse()
    }
    clamp(omega, .Machine$double.eps, 1 - .Machine$double.eps)
}

`.update_omega_all_woodbury` <- function (omega, gamma, bbar, dat, S0bar_inv, A_y, B_y, E_logit_rho_z) 
{
    p <- length(omega)
    omega <- clamp(as.numeric(omega), 0, 1)
    gamma <- clamp(as.numeric(gamma), 0, 1)
    ns_eff <- get_ns_eff(dat)
    S1bar_inv <- (ns_eff * dat$x_norm2)/dat$n
    S_inv <- S0bar_inv + (S1bar_inv - S0bar_inv) * gamma
    u <- omega * dat$Xty
    v <- gamma * S1bar_inv * bbar
    Cmat <- sweep(dat$XtX, 1L, omega, `*`)
    factor <- .trace_woodbury_factor(z = omega, diag_add = S_inv + dat$x_norm2 * omega * (1 - omega), dat = dat, 
        compute_inverse_diagonal = TRUE)
    solved <- .trace_woodbury_apply(factor, cbind(u, v, Cmat))
    Hu <- solved[, 1L]
    Hv <- solved[, 2L]
    HC <- solved[, 2L + seq_len(p), drop = FALSE]
    h <- factor$h_diag
    uu_minus <- .trace_leave_one_quad(sum(u * Hu), Hu, Hu, h)
    uv_minus <- .trace_leave_one_quad(sum(u * Hv), Hu, Hv, h)
    vv_minus <- .trace_leave_one_quad(sum(v * Hv), Hv, Hv, h)
    HC_diag <- diag(HC)
    cc_minus <- colSums(Cmat * HC) - HC_diag^2/h
    cu_minus <- as.numeric(crossprod(Cmat, Hu)) - HC_diag * Hu/h
    cv_minus <- as.numeric(crossprod(Cmat, Hv)) - HC_diag * Hv/h
    D <- pmax(S_inv + dat$x_norm2 - cc_minus, 1e-10)
    A <- dat$Xty - cu_minus
    B <- v - cv_minus
    alpha_term <- sum(gamma * S1bar_inv * bbar^2)
    T0 <- dat$yTy - uu_minus - 2 * uv_minus - vv_minus - (gamma * S1bar_inv * bbar^2) + alpha_term
    T1 <- T0 - (A + B)^2/D + (gamma * S1bar_inv * bbar^2)
    T0 <- pmax(T0, 1e-10)
    T1 <- pmax(T1, 1e-10)
    lin <- E_logit_rho_z - 0.5 * safe_log(D) + 0.5 * safe_log(S_inv) - (A_y + dat$n/2) * (safe_log(B_y + 0.5 * 
        T1) - safe_log(B_y + 0.5 * T0))
    clamp(expit(lin), .Machine$double.eps, 1 - .Machine$double.eps)
}

`.update_omega_gamma_interleaved_gauss_seidel` <- function (omega, gamma, bbar, dat, S0bar_inv, A_y, B_y, E_logit_rho_z, E_logit_rho_c) 
{
    p <- length(omega)
    omega <- clamp(as.numeric(omega), .Machine$double.eps, 1 - .Machine$double.eps)
    gamma <- clamp(as.numeric(gamma), .Machine$double.eps, 1 - .Machine$double.eps)
    if (p <= 1L) 
        return(list(omega = omega, gamma = gamma))
    ns_eff <- get_ns_eff(dat)
    S1bar_inv <- (ns_eff * dat$x_norm2)/dat$n
    S_inv_all <- S0bar_inv + (S1bar_inv - S0bar_inv) * gamma
    v_all <- gamma * S1bar_inv * bbar
    rebuild_inverse <- function() {
        solve_spd(.build_A0_cached(omega, gamma, dat, S0bar_inv), diag(p))
    }
    Minv <- rebuild_inverse()
    for (j in seq_len(p)) {
        idx <- seq_len(p) != j
        w_m <- omega[idx]
        S_inv_j <- S_inv_all[j]
        S_inv_m <- S_inv_all[idx]
        v_j <- v_all[j]
        v_m <- v_all[idx]
        u <- w_m * dat$Xty[idx]
        ccol <- w_m * dat$XtX[idx, j]
        SOL <- cbind(.solve_leave_one_cached(Minv, j, u), .solve_leave_one_cached(Minv, j, v_m), .solve_leave_one_cached(Minv, 
            j, ccol))
        sol_u <- SOL[, 1]
        sol_v <- SOL[, 2]
        sol_c <- SOL[, 3]
        D_j <- pmax(S_inv_j + dat$x_norm2[j] - sum(ccol * sol_c), 1e-10)
        A_j <- dat$Xty[j] - sum(ccol * sol_u)
        B_j <- v_j - sum(ccol * sol_v)
        alpha_term <- sum(gamma * S1bar_inv * bbar^2)
        T0 <- dat$yTy - sum(u * sol_u) - 2 * sum(u * sol_v) - sum(v_m * sol_v) - (gamma[j] * S1bar_inv[j] * bbar[j]^2) + 
            alpha_term
        T1 <- T0 - (A_j + B_j)^2/D_j + (gamma[j] * S1bar_inv[j] * bbar[j]^2)
        lin_omega <- E_logit_rho_z - 0.5 * safe_log(D_j) + 0.5 * safe_log(S_inv_j) - (A_y + dat$n/2) * (safe_log(B_y + 
            0.5 * pmax(T1, 1e-10)) - safe_log(B_y + 0.5 * pmax(T0, 1e-10)))
        omega_proposal <- clamp(expit(clamp(lin_omega, -30, 30)), .Machine$double.eps, 1 - .Machine$double.eps)
        omega_before <- omega
        delta_omega <- omega_proposal - omega_before[j]
        omega[j] <- omega_proposal
        inv_update <- .inverse_rank2_omega_update(Minv, j, delta_omega, omega_before, dat)
        Minv <- if (isTRUE(inv_update$ok)) 
            inv_update$inverse
        else rebuild_inverse()
        w_m <- omega[idx]
        wj <- omega[j]
        S_inv_m <- S_inv_all[idx]
        v_m <- v_all[idx]
        r_minus <- w_m * dat$Xty[idx] + v_m
        ccol <- w_m * dat$XtX[idx, j]
        crow <- dat$XtX[j, idx] * w_m
        SOL <- cbind(.solve_leave_one_cached(Minv, j, r_minus), .solve_leave_one_cached(Minv, j, ccol))
        sol_r <- SOL[, 1]
        sol_c <- SOL[, 2]
        cTc <- sum(ccol * sol_c)
        D0j <- pmax(S0bar_inv + wj * dat$x_norm2[j] - wj^2 * cTc, 1e-10)
        D1j <- pmax(S1bar_inv[j] + wj * dat$x_norm2[j] - wj^2 * cTc, 1e-10)
        R0 <- wj * (dat$Xty[j] - sum(crow * sol_r))
        v1j <- S1bar_inv[j] * bbar[j]
        R1 <- R0 + v1j
        alpha_minus <- sum(gamma[idx] * S1bar_inv[idx] * bbar[idx]^2)
        Tc0 <- pmax(dat$yTy - sum(r_minus * sol_r) - (R0^2)/D0j + alpha_minus, 1e-10)
        Tc1 <- pmax(dat$yTy - sum(r_minus * sol_r) - (R1^2)/D1j + alpha_minus + S1bar_inv[j] * bbar[j]^2, 1e-10)
        Delta <- 0.5 * safe_log(S1bar_inv[j]/S0bar_inv) - 0.5 * safe_log(1 + (S1bar_inv[j] - S0bar_inv)/D0j) - 
            (A_y + dat$n/2) * (safe_log(B_y + 0.5 * Tc1) - safe_log(B_y + 0.5 * Tc0))
        gamma_proposal <- clamp(expit(E_logit_rho_c + Delta), .Machine$double.eps, 1 - .Machine$double.eps)
        S_old_j <- S_inv_all[j]
        gamma[j] <- gamma_proposal
        S_inv_all[j] <- S0bar_inv + (S1bar_inv[j] - S0bar_inv) * gamma_proposal
        v_all[j] <- S1bar_inv[j] * gamma_proposal * bbar[j]
        inv_update <- .inverse_rank1_coordinate_update(Minv, j, S_inv_all[j] - S_old_j)
        Minv <- if (isTRUE(inv_update$ok)) 
            inv_update$inverse
        else rebuild_inverse()
    }
    list(omega = clamp(omega, .Machine$double.eps, 1 - .Machine$double.eps), gamma = clamp(gamma, .Machine$double.eps, 
        1 - .Machine$double.eps))
}

`%||%` <- function (x, y) 
if (is.null(x)) y else x

`borrow_elbo_surrogate` <- function (alpha, omega, gamma, dat, prior_par_borrow, sigma_beta_sq = 100, A_y = 1, B_y = 1, include_entropy = TRUE, 
    include_indicator_prior = TRUE, include_source_normalization = TRUE, ridge = 1e-10) 
{
    p <- length(omega)
    if (length(gamma) != p) 
        stop("omega and gamma must have the same length.")
    alpha <- .normalize_simplex(alpha, K = dat$K)
    omega <- clamp(as.numeric(omega), 0, 1)
    gamma <- clamp(as.numeric(gamma), 0, 1)
    S0bar_inv <- 1/sigma_beta_sq
    ns_eff <- get_ns_eff(dat)
    S1bar_inv <- (ns_eff * dat$x_norm2)/dat$n
    S_inv <- S0bar_inv + (S1bar_inv - S0bar_inv) * gamma
    S_inv <- pmax(S_inv, 1e-12)
    bbar <- as.numeric(as.matrix(dat$B) %*% alpha)
    source_rhs <- gamma * S1bar_inv * bbar
    rhs <- omega * dat$Xty + source_rhs
    factor <- .trace_woodbury_factor(z = omega, diag_add = S_inv + dat$x_norm2 * omega * (1 - omega), dat = dat, 
        ridge = ridge)
    sol <- .trace_woodbury_apply(factor, rhs)
    logdet <- factor$logdet
    if (!is.finite(logdet)) 
        logdet <- NA_real_
    quad <- sum(rhs * sol)
    source_quad <- sum(gamma * S1bar_inv * bbar^2)
    Tval <- dat$yTy - quad + source_quad
    Tval <- max(Tval, 1e-10)
    collapsed <- -0.5 * logdet - (A_y + dat$n/2) * safe_log(B_y + 0.5 * Tval)
    source_norm <- if (include_source_normalization) {
        0.5 * sum(gamma * safe_log(S1bar_inv) + (1 - gamma) * safe_log(S0bar_inv))
    }
    else {
        0
    }
    rho_z <- clamp(prior_par_borrow$rho_z_mean, .Machine$double.eps, 1 - .Machine$double.eps)
    rho_c <- clamp(prior_par_borrow$rho_c_mean, .Machine$double.eps, 1 - .Machine$double.eps)
    indicator_prior <- if (include_indicator_prior) {
        sum(omega * log(rho_z) + (1 - omega) * log1p(-rho_z)) + sum(gamma * log(rho_c) + (1 - gamma) * log1p(-rho_c))
    }
    else {
        0
    }
    entropy <- if (include_entropy) {
        .bernoulli_entropy_sum(omega) + .bernoulli_entropy_sum(gamma)
    }
    else {
        0
    }
    total <- collapsed + source_norm + indicator_prior + entropy
    list(total = as.numeric(total), collapsed = as.numeric(collapsed), source_norm = as.numeric(source_norm), indicator_prior = as.numeric(indicator_prior), 
        entropy = as.numeric(entropy), logdet = as.numeric(logdet), Tval = as.numeric(Tval), quad = as.numeric(quad), 
        source_quad = as.numeric(source_quad), monitor_name = "paper_moment_collapsed_variational_objective", alpha = alpha)
}

`build_common_objects` <- function (omega, gamma, alpha, dat, S0bar_inv) 
{
    XtX <- dat$XtX
    X <- dat$X
    B <- as.matrix(dat$B)
    n <- dat$n
    ns_eff <- get_ns_eff(dat)
    x_norm2 <- dat$x_norm2
    Omega <- omega %o% omega
    diag(Omega) <- omega
    S1bar_inv <- (ns_eff * x_norm2)/n
    S_inv_all <- S0bar_inv + (S1bar_inv - S0bar_inv) * gamma
    A0 <- XtX * Omega
    diag(A0) <- diag(A0) + S_inv_all
    A0 <- A0 + diag(1e-10, length(omega))
    alpha <- as.numeric(alpha)
    bbar <- as.numeric(B %*% alpha)
    theta <- omega * gamma * bbar
    mu <- as.numeric(X %*% theta)
    r <- dat$y - mu
    list(Omega = Omega, S1bar_inv = S1bar_inv, S_inv_all = S_inv_all, A0 = A0, bbar = bbar, theta = theta, mu = mu, 
        r = r)
}

`clamp` <- function (x, lo, hi) 
pmin(pmax(x, lo), hi)

`configure_trace_source_precision` <- function (dat, mode = c("high_dimensional_logp", "raw_ns", "target_capped"), cap_ratio = 1) 
{
    mode <- match.arg(mode)
    cap_ratio <- as.numeric(cap_ratio)
    if (length(cap_ratio) != 1L || !is.finite(cap_ratio) || cap_ratio <= 0) {
        stop("source_precision_cap_ratio must be one finite positive number.")
    }
    raw_ns <- trace_get_raw_ns_eff(dat)
    precision_ns <- if (mode == "high_dimensional_logp") {
        if (dat$p < 2) 
            stop("High dimensional precision requires p >= 2.")
        raw_ns/log(dat$p)
    }
    else if (mode == "target_capped") {
        min(raw_ns, cap_ratio * dat$n)
    }
    else {
        raw_ns
    }
    dat$trace_source_precision_mode <- mode
    dat$trace_source_precision_cap_ratio <- cap_ratio
    dat$trace_source_precision_ns_raw <- raw_ns
    dat$trace_source_precision_ns_eff <- precision_ns
    dat$trace_source_precision_power <- precision_ns/raw_ns
    dat
}

`expit` <- function (x) 
plogis(x)

`fit_borrow_model_fixedhyper` <- function (dat, prior_par_borrow, alpha_method = c("old", "fixed", "oracle"), alpha_fixed = NULL, init_alpha = c("uniform", 
    "source1", "random", "oracle"), init_method = c("full", "mild", "source", "null", "mcp", "lasso", "scad", "mix", 
    "custom"), init_omega = NULL, init_gamma = "theory_zero", sigma_beta_sq = 100, max_iter = 100, tol = 1e-06, 
    tol_elbo = 1e-07, min_iter = 5, gamma_ebic = 0.5, verbose = FALSE, fast_cached_update = TRUE, borrow_update_schedule = c("paper_block", 
        "gauss_seidel", "jacobi"), alpha_update_every = 1, alpha_update_strategy = c("em_qp", "elbo_accept", "old"), 
    alpha_freeze_iter = 0, alpha_grid_step = 0.02, alpha_accept_tol = 1e-08, adaptive_damping = TRUE, update_damping = 0.5, 
    cycle_detection_tol = 1e-06, cycle_detection_relative_tol = 0.001, cycle_min_amplitude = NULL, theory_initialization = TRUE, 
    target_init_until_converged = FALSE, target_init_steps = 5L, target_init_max_iter = 1000L, target_init_tol = NULL, 
    preliminary_alpha_update = c("huber", "target_matching_qp", "none"), huber_tau = NULL, huber_tau_rule = c("fixed", 
        "diverging"), huber_tau_multiplier = 1, monitor_elbo = TRUE, elbo_print_every = 1, include_source_normalization = TRUE, 
    stop_by_elbo = FALSE, warn_elbo_drop = TRUE, elbo_drop_tol = 1e-06, stage_label = NULL) 
{
    alpha_method <- match.arg(alpha_method)
    init_alpha <- match.arg(init_alpha)
    init_method <- match.arg(init_method)
    alpha_update_strategy <- match.arg(alpha_update_strategy)
    borrow_update_schedule <- match.arg(borrow_update_schedule)
    preliminary_alpha_update <- match.arg(preliminary_alpha_update)
    huber_tau_rule <- match.arg(huber_tau_rule)
    alpha_update_every <- max(1, as.integer(alpha_update_every))
    alpha_freeze_iter <- max(0L, as.integer(alpha_freeze_iter))
    alpha_grid_step <- as.numeric(alpha_grid_step)
    adaptive_damping <- isTRUE(adaptive_damping)
    update_damping <- as.numeric(update_damping)
    if (length(update_damping) != 1L || !is.finite(update_damping) || update_damping <= 0 || update_damping > 1) {
        stop("update_damping must be a finite scalar in (0, 1].")
    }
    cycle_detection_tol <- as.numeric(cycle_detection_tol)
    cycle_detection_relative_tol <- as.numeric(cycle_detection_relative_tol)
    if (length(cycle_detection_tol) != 1L || !is.finite(cycle_detection_tol) || cycle_detection_tol <= 0) {
        stop("cycle_detection_tol must be a positive finite scalar.")
    }
    if (length(cycle_detection_relative_tol) != 1L || !is.finite(cycle_detection_relative_tol) || cycle_detection_relative_tol < 
        0) {
        stop("cycle_detection_relative_tol must be a nonnegative finite scalar.")
    }
    if (is.null(cycle_min_amplitude)) {
        cycle_min_amplitude <- max(100 * tol, 1e-04)
    }
    else {
        cycle_min_amplitude <- as.numeric(cycle_min_amplitude)
        if (length(cycle_min_amplitude) != 1L || !is.finite(cycle_min_amplitude) || cycle_min_amplitude <= 0) {
            stop("cycle_min_amplitude must be a positive finite scalar.")
        }
    }
    min_iter <- max(1, as.integer(min_iter))
    elbo_print_every <- max(1, as.integer(elbo_print_every))
    A_y <- 1
    B_y <- 1
    A_c <- prior_par_borrow$A_c
    B_c <- prior_par_borrow$B_c
    A_z <- prior_par_borrow$A_z
    B_z <- prior_par_borrow$B_z
    S0bar <- sigma_beta_sq
    S0bar_inv <- 1/S0bar
    p <- dat$p
    if (is.null(stage_label)) 
        stage_label <- alpha_method
    omega <- make_init_omega(dat = dat, init_method = init_method, init_omega = init_omega, gamma_ebic = gamma_ebic)
    alpha <- make_init_alpha(dat = dat, init_alpha = init_alpha, alpha_fixed = alpha_fixed, alpha_method = alpha_method)
    gamma_requested <- init_gamma
    if (isTRUE(theory_initialization)) {
        gamma <- rep(0, p)
        init_gamma_effective <- "theory_zero"
    }
    else {
        gamma <- make_init_gamma_borrow(dat = dat, init_gamma = init_gamma, alpha = alpha, gamma_ebic = gamma_ebic)
        init_gamma_effective <- init_gamma
    }
    rho_z_algorithm <- prior_par_borrow$rho_z_mean %||% (A_z/(A_z + B_z))
    rho_c_algorithm <- prior_par_borrow$rho_c_mean %||% (A_c/(A_c + B_c))
    E_logit_rho_z <- qlogis(clamp(rho_z_algorithm, .Machine$double.eps, 1 - .Machine$double.eps))
    E_logit_rho_c <- qlogis(clamp(rho_c_algorithm, .Machine$double.eps, 1 - .Machine$double.eps))
    elbo_trace <- list()
    elbo_prev <- NA_real_
    theory_init_trace <- list()
    add_theory_trace <- function(step, iter, omega_now, gamma_now, alpha_now) {
        theory_init_trace[[length(theory_init_trace) + 1L]] <<- data.frame(iter = iter, step = step, size05 = sum(omega_now > 
            0.5), gamma_min = min(gamma_now), gamma_max = max(gamma_now), mean_gamma = mean(gamma_now), alpha_1 = alpha_now[1], 
            alpha_2 = if (length(alpha_now) >= 2) 
                alpha_now[2]
            else NA_real_, alpha_3 = if (length(alpha_now) >= 3) 
                alpha_now[3]
            else NA_real_, stringsAsFactors = FALSE)
    }
    add_theory_trace("initial", 0L, omega, gamma, alpha)
    if (isTRUE(monitor_elbo)) {
        row0 <- .make_elbo_row(iter = 0, stage = stage_label, alpha_method = alpha_method, init_method = init_method, 
            alpha = alpha, omega = omega, gamma = gamma, dat = dat, prior_par_borrow = prior_par_borrow, sigma_beta_sq = sigma_beta_sq, 
            A_y = A_y, B_y = B_y, include_source_normalization = include_source_normalization)
        row0$raw_max_param_change <- NA_real_
        row0$cycle_distance_raw <- NA_real_
        row0$adaptive_damping_active <- !adaptive_damping && update_damping < 1
        row0$damping_applied <- FALSE
        row0$damping_factor <- 1
        elbo_trace[[length(elbo_trace) + 1]] <- row0
        elbo_prev <- row0$elbo
        if (verbose) 
            print_elbo_trace_row(row0)
    }
    n_iter <- max_iter
    converged <- FALSE
    convergence_reason <- "max_iter_reached"
    d_omega <- NA_real_
    d_gamma <- NA_real_
    d_alpha <- NA_real_
    max_param_change <- NA_real_
    elbo_drop_count <- 0L
    alpha_update_accept_count <- 0L
    alpha_update_reject_count <- 0L
    last_alpha_step_size <- NA_real_
    last_alpha_update_elbo_old <- NA_real_
    last_alpha_update_elbo_new <- NA_real_
    last_alpha_update_qp_old <- NA_real_
    last_alpha_update_qp_new <- NA_real_
    last_alpha_active_set_size <- NA_integer_
    last_alpha_active_set <- NA_character_
    alpha_eta <- NULL
    last_alpha_update_raw <- NULL
    lag2_state <- NULL
    damping_active <- !adaptive_damping && update_damping < 1
    damping_triggered <- FALSE
    damping_trigger_iter <- NA_integer_
    damping_iterations <- 0L
    cycle_distance_at_trigger <- NA_real_
    raw_max_param_change <- NA_real_
    cycle_distance_raw <- NA_real_
    damping_applied_t <- FALSE
    target_initialization <- NULL
    if (isTRUE(target_init_until_converged) && !isTRUE(theory_initialization)) {
        stop("Target-only convergence initialization requires theory_initialization=TRUE (gamma=0).")
    }
    preliminary_alpha_done <- FALSE
    preliminary_alpha_qp_old <- NA_real_
    preliminary_alpha_qp_new <- NA_real_
    preliminary_alpha_active_set <- NA_character_
    preliminary_alpha_solver <- NA_character_
    preliminary_alpha_beta_target <- NULL
    preliminary_alpha_huber_tau <- NA_real_
    preliminary_alpha_huber_tau_rule <- NA_character_
    preliminary_alpha_huber_s_reference <- NA_real_
    preliminary_alpha_huber_iterations <- NA_integer_
    preliminary_alpha_huber_converged <- NA
    preliminary_alpha_huber_kkt_residual <- NA_real_
    update_omega_once <- function(omega_in, gamma_in, alpha_in) {
        bbar_in <- as.numeric(as.matrix(dat$B) %*% alpha_in)
        if (fast_cached_update) {
            update_fun <- if (borrow_update_schedule == "gauss_seidel") {
                .update_omega_all_gauss_seidel
            }
            else {
                .update_omega_all_cached
            }
            return(update_fun(omega = omega_in, gamma = gamma_in, bbar = bbar_in, dat = dat, S0bar_inv = S0bar_inv, 
                A_y = A_y, B_y = B_y, E_logit_rho_z = E_logit_rho_z))
        }
        omega_prop <- omega_in
        for (j in seq_len(p)) {
            omega_prop[j] <- update_omega_j_old(j = j, omega = omega_in, gamma = gamma_in, bbar = bbar_in, dat = dat, 
                S0bar_inv = S0bar_inv, A_y = A_y, B_y = B_y, E_logit_rho_z = E_logit_rho_z)
        }
        clamp(omega_prop, .Machine$double.eps, 1 - .Machine$double.eps)
    }
    update_gamma_once <- function(omega_in, gamma_in, alpha_in) {
        bbar_in <- as.numeric(as.matrix(dat$B) %*% alpha_in)
        if (fast_cached_update) {
            update_fun <- if (borrow_update_schedule == "gauss_seidel") {
                .update_gamma_all_gauss_seidel
            }
            else {
                .update_gamma_all_cached
            }
            return(update_fun(omega = omega_in, gamma = gamma_in, bbar = bbar_in, dat = dat, S0bar_inv = S0bar_inv, 
                A_y = A_y, B_y = B_y, E_logit_rho_c = E_logit_rho_c))
        }
        gamma_prop <- gamma_in
        for (j in seq_len(p)) {
            gamma_prop[j] <- update_gamma_j(j = j, omega = omega_in, gamma = gamma_in, bbar = bbar_in, dat = dat, 
                S0bar_inv = S0bar_inv, A_y = A_y, B_y = B_y, E_logit_rho_c = E_logit_rho_c)
        }
        clamp(gamma_prop, .Machine$double.eps, 1 - .Machine$double.eps)
    }
    do_regular_alpha_update <- function(t, alpha_in, omega_in, gamma_in) {
        out <- list(alpha = alpha_in, updated = FALSE)
        if (alpha_method != "old") 
            return(out)
        if (t <= alpha_freeze_iter) 
            return(out)
        should_update_alpha <- (t == alpha_freeze_iter + 1L) || ((t - alpha_freeze_iter)%%alpha_update_every == 
            0L) || (t == max_iter)
        if (!should_update_alpha) 
            return(out)
        if (alpha_update_strategy == "old") {
            au <- update_alpha_old(alpha_in, omega_in, gamma_in, dat, A_y, B_y, S0bar_inv)
            out$alpha <- .normalize_simplex(au$alpha, K = dat$K)
            out$raw <- au
            out$eta <- .extract_alpha_eta(au)
            out$accepted <- TRUE
            out$step_size <- NA_real_
        }
        else if (alpha_update_strategy == "elbo_accept") {
            au <- update_alpha_with_elbo_accept(alpha_old = alpha_in, omega = omega_in, gamma = gamma_in, dat = dat, 
                prior_par_borrow = prior_par_borrow, sigma_beta_sq = sigma_beta_sq, A_y = A_y, B_y = B_y, S0bar_inv = S0bar_inv, 
                alpha_eta_old = alpha_eta, include_source_normalization = include_source_normalization, accept_tol = alpha_accept_tol, 
                accept_objective = "point_elbo")
            out$alpha <- au$alpha
            out$eta <- au$eta
            out$raw <- au$raw_update
            out$accepted <- isTRUE(au$accepted)
            out$step_size <- au$step_size
            out$elbo_old <- au$elbo_old
            out$elbo_new <- au$elbo_new
        }
        else {
            au <- update_alpha_em_qp(alpha_old = alpha_in, omega = omega_in, gamma = gamma_in, dat = dat, sigma_beta_sq = sigma_beta_sq)
            out$alpha <- au$alpha
            out$accepted <- TRUE
            out$step_size <- au$step_size
            out$qp_old <- au$qp_old
            out$qp_new <- au$qp_new
            out$active_set <- au$active_set
        }
        out$updated <- TRUE
        out
    }
    for (t in seq_len(max_iter)) {
        omega_old <- omega
        gamma_old <- gamma
        alpha_old <- alpha
        is_theory_first_cycle <- isTRUE(theory_initialization) && t == 1L
        if (is_theory_first_cycle) {
            if (isTRUE(target_init_until_converged)) {
                target_initialization <- trace_target_only_initialize(dat = dat, omega = omega, S0bar_inv = S0bar_inv, 
                  rho_z = rho_z_algorithm, A_y = A_y, B_y = B_y, max_iter = target_init_max_iter, tol = if (is.null(target_init_tol)) 
                    tol
                  else target_init_tol, verbose = verbose)
                if (!isTRUE(target_initialization$converged)) {
                  stop(sprintf("Target-only omega did not converge after %d iterations (change %.4e); alpha/gamma updates were not started.", 
                    target_initialization$n_iter, target_initialization$max_change))
                }
                omega <- target_initialization$omega
                add_theory_trace("omega_target_only_converged", 1L, omega, gamma, alpha)
            }
            else {
                target_initialization <- trace_target_only_fixed_steps(dat = dat, omega = omega, S0bar_inv = S0bar_inv, 
                  rho_z = rho_z_algorithm, A_y = A_y, B_y = B_y, n_steps = target_init_steps, tol = if (is.null(target_init_tol)) 
                    tol
                  else target_init_tol, verbose = verbose)
                omega <- target_initialization$omega
                add_theory_trace("omega_target_only_fixed_steps", 1L, omega, gamma, alpha)
            }
            if (alpha_method == "old" && preliminary_alpha_update != "none") {
                au0 <- if (preliminary_alpha_update == "huber") {
                  update_alpha_preliminary_huber(alpha_old = alpha, omega = omega, dat = dat, sigma_beta_sq = sigma_beta_sq, 
                    huber_tau = huber_tau, huber_tau_rule = huber_tau_rule, huber_tau_multiplier = huber_tau_multiplier)
                }
                else {
                  update_alpha_preliminary_qp(alpha_old = alpha, omega = omega, dat = dat, sigma_beta_sq = sigma_beta_sq)
                }
                alpha <- au0$alpha
                preliminary_alpha_done <- TRUE
                preliminary_alpha_qp_old <- au0$qp_old
                preliminary_alpha_qp_new <- au0$qp_new
                preliminary_alpha_active_set <- paste(au0$active_set, collapse = ";")
                preliminary_alpha_solver <- au0$solver
                preliminary_alpha_beta_target <- au0$beta_target
                preliminary_alpha_huber_tau <- au0$huber_tau %||% NA_real_
                preliminary_alpha_huber_tau_rule <- au0$huber_tau_rule %||% NA_character_
                preliminary_alpha_huber_s_reference <- au0$huber_s_reference %||% NA_real_
                preliminary_alpha_huber_iterations <- au0$huber_iterations %||% NA_integer_
                preliminary_alpha_huber_converged <- au0$huber_converged %||% NA
                preliminary_alpha_huber_kkt_residual <- au0$huber_kkt_residual %||% NA_real_
            }
            else if (alpha_method == "fixed") {
                alpha <- .normalize_simplex(alpha_fixed, K = dat$K)
            }
            else if (alpha_method == "oracle") {
                alpha <- .normalize_simplex(dat$alpha_true, K = dat$K)
            }
            add_theory_trace("alpha_preliminary", 1L, omega, gamma, alpha)
            gamma <- update_gamma_once(omega, gamma, alpha)
            add_theory_trace("gamma_first_formal", 1L, omega, gamma, alpha)
        }
        else {
            if (borrow_update_schedule == "gauss_seidel") {
                if (!isTRUE(fast_cached_update)) {
                  stop("borrow_update_schedule = 'gauss_seidel' requires ", "fast_cached_update = TRUE.")
                }
                joint_update <- .update_omega_gamma_interleaved_gauss_seidel(omega = omega, gamma = gamma, bbar = as.numeric(as.matrix(dat$B) %*% 
                  alpha), dat = dat, S0bar_inv = S0bar_inv, A_y = A_y, B_y = B_y, E_logit_rho_z = E_logit_rho_z, 
                  E_logit_rho_c = E_logit_rho_c)
                omega <- joint_update$omega
                gamma <- joint_update$gamma
            }
            else {
                omega <- update_omega_once(omega, gamma, alpha)
                gamma <- update_gamma_once(omega, gamma, alpha)
            }
            if (alpha_method == "fixed") {
                alpha <- .normalize_simplex(alpha_fixed, K = dat$K)
            }
            else if (alpha_method == "oracle") {
                alpha <- .normalize_simplex(dat$alpha_true, K = dat$K)
            }
            else {
                au <- do_regular_alpha_update(t, alpha, omega, gamma)
                if (isTRUE(au$updated)) {
                  alpha <- au$alpha
                  last_alpha_update_raw <- au$raw %||% last_alpha_update_raw
                  if (!is.null(au$eta)) 
                    alpha_eta <- au$eta
                  if (isTRUE(au$accepted)) {
                    alpha_update_accept_count <- alpha_update_accept_count + 1L
                  }
                  else {
                    alpha_update_reject_count <- alpha_update_reject_count + 1L
                  }
                  last_alpha_step_size <- au$step_size %||% NA_real_
                  last_alpha_update_elbo_old <- au$elbo_old %||% NA_real_
                  last_alpha_update_elbo_new <- au$elbo_new %||% NA_real_
                  last_alpha_update_qp_old <- au$qp_old %||% NA_real_
                  last_alpha_update_qp_new <- au$qp_new %||% NA_real_
                  if (!is.null(au$active_set)) {
                    last_alpha_active_set_size <- length(au$active_set)
                    last_alpha_active_set <- paste(au$active_set, collapse = ";")
                  }
                }
            }
        }
        omega_raw <- omega
        gamma_raw <- gamma
        alpha_raw <- alpha
        raw_max_param_change <- max(max(abs(omega_raw - omega_old)), max(abs(gamma_raw - gamma_old)), max(abs(alpha_raw - 
            alpha_old)))
        cycle_distance_raw <- if (!is.null(lag2_state)) {
            max(max(abs(omega_raw - lag2_state$omega)), max(abs(gamma_raw - lag2_state$gamma)), max(abs(alpha_raw - 
                lag2_state$alpha)))
        }
        else {
            NA_real_
        }
        cycle_threshold_t <- max(cycle_detection_tol, cycle_detection_relative_tol * raw_max_param_change)
        cycle_match <- adaptive_damping && !damping_active && !is_theory_first_cycle && t >= 3L && is.finite(cycle_distance_raw) && 
            raw_max_param_change > cycle_min_amplitude && cycle_distance_raw <= cycle_threshold_t
        if (cycle_match) {
            damping_active <- TRUE
            damping_triggered <- TRUE
            damping_trigger_iter <- t
            cycle_distance_at_trigger <- cycle_distance_raw
            if (isTRUE(verbose)) {
                cat(sprintf("[TRACE adaptive damping] %s iter %d: two-cycle distance=%.3e, raw change=%.3e; eta=%.3f activated.\n", 
                  stage_label, t, cycle_distance_raw, raw_max_param_change, update_damping))
            }
        }
        damping_applied_t <- damping_active && !is_theory_first_cycle && update_damping < 1
        if (damping_applied_t) {
            omega <- clamp((1 - update_damping) * omega_old + update_damping * omega_raw, .Machine$double.eps, 
                1 - .Machine$double.eps)
            gamma <- clamp((1 - update_damping) * gamma_old + update_damping * gamma_raw, .Machine$double.eps, 
                1 - .Machine$double.eps)
            alpha <- .normalize_simplex((1 - update_damping) * alpha_old + update_damping * alpha_raw, K = dat$K)
            damping_iterations <- damping_iterations + 1L
        }
        lag2_state <- list(omega = omega_old, gamma = gamma_old, alpha = alpha_old)
        d_omega <- max(abs(omega - omega_old))
        d_gamma <- max(abs(gamma - gamma_old))
        d_alpha <- max(abs(alpha - alpha_old))
        max_param_change <- max(d_omega, d_gamma, d_alpha)
        rel_delta_elbo <- NA_real_
        delta_elbo <- NA_real_
        if (isTRUE(monitor_elbo)) {
            rowt <- .make_elbo_row(iter = t, stage = if (is_theory_first_cycle) 
                paste0(stage_label, "_theory_first_cycle")
            else stage_label, alpha_method = alpha_method, init_method = init_method, alpha = alpha, omega = omega, 
                gamma = gamma, dat = dat, prior_par_borrow = prior_par_borrow, sigma_beta_sq = sigma_beta_sq, A_y = A_y, 
                B_y = B_y, d_omega = d_omega, d_gamma = d_gamma, d_alpha = d_alpha, elbo_prev = elbo_prev, include_source_normalization = include_source_normalization)
            rowt$raw_max_param_change <- raw_max_param_change
            rowt$cycle_distance_raw <- cycle_distance_raw
            rowt$adaptive_damping_active <- damping_active
            rowt$damping_applied <- damping_applied_t
            rowt$damping_factor <- if (damping_applied_t) 
                update_damping
            else 1
            elbo_trace[[length(elbo_trace) + 1]] <- rowt
            delta_elbo <- rowt$delta_elbo
            rel_delta_elbo <- rowt$rel_delta_elbo
            if (is.finite(delta_elbo) && delta_elbo < -abs(elbo_drop_tol)) {
                elbo_drop_count <- elbo_drop_count + 1L
                if (isTRUE(warn_elbo_drop)) {
                  cat(sprintf("[ELBO drop warning] %s iter %d: dELBO = %.6e. The first theory cycle is not a standard coordinate-ascent sweep, and the monitored surrogate need not be monotone.\n", 
                    stage_label, t, delta_elbo))
                }
            }
            if (verbose && (t%%elbo_print_every == 0 || t == 1)) {
                print_elbo_trace_row(rowt)
            }
            elbo_prev <- rowt$elbo
        }
        else if (verbose) {
            cat(sprintf("[%s/%s fast=%s] iter %3d | d_omega=%.3e d_gamma=%.3e d_alpha=%.3e | size05=%d\n", alpha_method, 
                init_method, fast_cached_update, t, d_omega, d_gamma, d_alpha, sum(omega > 0.5)))
        }
        param_ok <- max_param_change < tol && raw_max_param_change < tol
        elbo_ok <- if (isTRUE(monitor_elbo) && isTRUE(stop_by_elbo)) {
            is.finite(rel_delta_elbo) && rel_delta_elbo < tol_elbo
        }
        else {
            TRUE
        }
        earliest_stop <- if (isTRUE(theory_initialization)) 
            3L
        else 1L
        if (t >= max(min_iter, alpha_freeze_iter + 2L, earliest_stop) && param_ok && elbo_ok) {
            n_iter <- t
            converged <- TRUE
            convergence_reason <- if (isTRUE(stop_by_elbo) && isTRUE(monitor_elbo)) {
                "param_and_elbo_tol"
            }
            else {
                "param_tol"
            }
            break
        }
    }
    beta_hat <- get_beta_hat(alpha, omega, gamma, dat, S0bar_inv)
    method_name <- if (alpha_method == "oracle") {
        "borrow_oracle_alpha"
    }
    else if (alpha_method == "fixed") {
        "borrow_fixed_alpha"
    }
    else {
        "borrow_old"
    }
    elbo_trace_df <- if (length(elbo_trace) > 0) 
        do.call(rbind, elbo_trace)
    else data.frame()
    theory_init_trace_df <- if (length(theory_init_trace) > 0) 
        do.call(rbind, theory_init_trace)
    else data.frame()
    final_elbo <- if (nrow(elbo_trace_df) > 0) 
        tail(elbo_trace_df$elbo, 1)
    else NA_real_
    list(method = paste0(method_name, "_init_", init_method, if (isTRUE(fast_cached_update) && borrow_update_schedule != 
        "gauss_seidel") {
        "_woodbury_elbo_theoryinit"
    } else {
        "_cached_elbo_theoryinit"
    }), alpha = alpha, alpha_eta = alpha_eta, last_alpha_update_raw = last_alpha_update_raw, omega = omega, gamma = gamma, 
        beta_hat = beta_hat, fixed_E_logit_rho_z = E_logit_rho_z, fixed_E_logit_rho_c = E_logit_rho_c, init_method = init_method, 
        init_alpha = init_alpha, init_gamma = init_gamma_effective, init_gamma_requested = gamma_requested, theory_initialization = theory_initialization, 
        target_init_until_converged = target_init_until_converged, target_init_steps = target_init_steps, target_init_max_iter = target_init_max_iter, 
        target_init_tol = target_init_tol, preliminary_alpha_update = preliminary_alpha_update, preliminary_alpha_done = preliminary_alpha_done, 
        preliminary_alpha_qp_old = preliminary_alpha_qp_old, preliminary_alpha_qp_new = preliminary_alpha_qp_new, 
        preliminary_alpha_active_set = preliminary_alpha_active_set, preliminary_alpha_solver = preliminary_alpha_solver, 
        preliminary_alpha_beta_target = preliminary_alpha_beta_target, preliminary_alpha_huber_tau = preliminary_alpha_huber_tau, 
        preliminary_alpha_huber_tau_rule = preliminary_alpha_huber_tau_rule, preliminary_alpha_huber_s_reference = preliminary_alpha_huber_s_reference, 
        preliminary_alpha_huber_iterations = preliminary_alpha_huber_iterations, preliminary_alpha_huber_converged = preliminary_alpha_huber_converged, 
        preliminary_alpha_huber_kkt_residual = preliminary_alpha_huber_kkt_residual, theory_init_trace = theory_init_trace_df, 
        target_initialization = target_initialization, target_init_converged = if (is.null(target_initialization)) NA else target_initialization$converged, 
        target_init_n_iter = if (is.null(target_initialization)) 1L else target_initialization$n_iter, target_init_max_change = if (is.null(target_initialization)) NA_real_ else target_initialization$max_change, 
        target_init_size = if (is.null(target_initialization)) NA_integer_ else sum(target_initialization$omega > 
            0.5), sigma_beta_sq = sigma_beta_sq, n_iter = n_iter, converged = converged, convergence_reason = convergence_reason, 
        d_omega_final = d_omega, d_gamma_final = d_gamma, d_alpha_final = d_alpha, max_param_change_final = max_param_change, 
        stop_by_elbo = stop_by_elbo, final_elbo = final_elbo, elbo_monitor_name = "paper_moment_collapsed_variational_objective", 
        elbo_trace = elbo_trace_df, elbo_drop_count = elbo_drop_count, fast_cached_update = fast_cached_update, 
        trace_linear_algebra_backend = if (isTRUE(fast_cached_update) && borrow_update_schedule != "gauss_seidel") {
            "exact_woodbury_n_by_n"
        } else if (isTRUE(fast_cached_update)) {
            "cached_rank_updates_p_by_p"
        } else {
            "dense_leave_one"
        }, borrow_update_schedule = borrow_update_schedule, alpha_update_every = alpha_update_every, alpha_update_strategy = alpha_update_strategy, 
        alpha_freeze_iter = alpha_freeze_iter, alpha_grid_step = alpha_grid_step, adaptive_damping = adaptive_damping, 
        update_damping = update_damping, damping_triggered = damping_triggered, damping_trigger_iter = damping_trigger_iter, 
        damping_iterations = damping_iterations, cycle_distance_at_trigger = cycle_distance_at_trigger, cycle_detection_tol = cycle_detection_tol, 
        cycle_detection_relative_tol = cycle_detection_relative_tol, cycle_min_amplitude = cycle_min_amplitude, 
        raw_max_param_change_final = raw_max_param_change, alpha_update_accept_count = alpha_update_accept_count, 
        alpha_update_reject_count = alpha_update_reject_count, last_alpha_step_size = last_alpha_step_size, last_alpha_update_elbo_old = last_alpha_update_elbo_old, 
        last_alpha_update_elbo_new = last_alpha_update_elbo_new, last_alpha_update_qp_old = last_alpha_update_qp_old, 
        last_alpha_update_qp_new = last_alpha_update_qp_new, last_alpha_active_set_size = last_alpha_active_set_size, 
        last_alpha_active_set = last_alpha_active_set, tol = tol, tol_elbo = tol_elbo, include_source_normalization = include_source_normalization)
}

`fit_lasso_lambda` <- function (X, y, lambda) 
{
    fit <- glmnet::glmnet(x = X, y = y, alpha = 1, lambda = lambda, intercept = FALSE, standardize = FALSE)
    as.numeric(coef(fit, s = lambda))[-1]
}

`fit_noborrow_acvb` <- function (dat, init_method = c("null", "full", "mcp", "lasso", "scad", "mix", "mild", "custom"), init_omega = NULL, 
    sigma_beta_sq = 0.1, A_y = 1, B_y = 1, pi_prior = NULL, max_iter = 100, tol = 1e-06, gamma_ebic = 0.5, verbose = FALSE, 
    fast_cached_update = TRUE, acceleration = c("woodbury", "cached_inverse", "slow")) 
{
    init_method <- match.arg(init_method)
    acceleration <- match.arg(acceleration)
    if (!isTRUE(fast_cached_update) && acceleration != "slow") 
        acceleration <- "slow"
    X <- dat$X
    y <- dat$y
    n <- dat$n
    p <- dat$p
    Xty <- dat$Xty
    yTy <- dat$yTy
    if (is.null(pi_prior)) 
        pi_prior <- min(0.1, max(1/n, 1/sqrt(p)))
    pi_prior <- clamp(pi_prior, .Machine$double.eps, 1 - .Machine$double.eps)
    logit_pi <- qlogis(pi_prior)
    omega <- make_init_omega(dat = dat, init_method = init_method, init_omega = init_omega, gamma_ebic = gamma_ebic)
    acvb_beta_given_w <- function(w) {
        rhs <- w * Xty
        sol <- if (acceleration == "woodbury") {
            .acvb_solve_woodbury(w, rhs, dat, sigma_beta_sq)
        }
        else {
            G <- .build_acvb_G_cached(w, dat, sigma_beta_sq)
            solve_spd(G, rhs)
        }
        as.numeric(w * sol)
    }
    n_iter <- max_iter
    converged <- FALSE
    convergence_reason <- "max_iter_reached"
    iter_time_sec <- numeric(max_iter)
    min_downdate_denom <- NA_real_
    n_denom_clamped <- 0L
    for (t in seq_len(max_iter)) {
        omega_old <- omega
        omega_prop <- omega
        tm_iter <- system.time({
            if (acceleration == "woodbury" && p > 1L) {
                upd <- .acvb_update_woodbury(omega = omega, dat = dat, sigma_beta_sq = sigma_beta_sq, A_y = A_y, 
                  B_y = B_y, logit_pi = logit_pi)
                omega_prop <- upd$omega
                min_downdate_denom <- upd$min_downdate_denom
                n_denom_clamped <- n_denom_clamped + upd$n_denom_clamped
            }
            else if (acceleration == "cached_inverse" && p > 1L) {
                omega_prop <- .acvb_update_cached_inverse_legacy(omega = omega, dat = dat, sigma_beta_sq = sigma_beta_sq, 
                  A_y = A_y, B_y = B_y, logit_pi = logit_pi)$omega
            }
            else {
                acvb_T <- function(w, j, value) {
                  wj <- w
                  wj[j] <- value
                  beta_j <- acvb_beta_given_w(wj)
                  scale_term <- B_y + 0.5 * yTy - 0.5 * sum(Xty * beta_j)
                  scale_term <- max(scale_term, 1e-12)
                  out <- -(A_y + n/2) * log(scale_term)
                  if (value == 1) {
                    out <- out + logit_pi - 0.5 * log(n) - 0.5 * log(sigma_beta_sq)
                  }
                  out
                }
                for (j in seq_len(p)) {
                  T0 <- acvb_T(omega, j, 0)
                  T1 <- acvb_T(omega, j, 1)
                  omega_prop[j] <- expit(T1 - T0)
                }
            }
            omega <- clamp(omega_prop, .Machine$double.eps, 1 - .Machine$double.eps)
        })
        iter_time_sec[t] <- unname(tm_iter["elapsed"])
        d_omega <- max(abs(omega - omega_old))
        if (verbose) {
            cat(sprintf("[ACVB/%s/%s] iter %3d | d_omega=%.3e | size05=%d | iter_time=%.3fs\n", init_method, acceleration, 
                t, d_omega, sum(omega > 0.5), iter_time_sec[t]))
        }
        if (t > 2L && d_omega < tol) {
            n_iter <- t
            converged <- TRUE
            convergence_reason <- "omega_tol"
            break
        }
    }
    iter_time_sec <- iter_time_sec[seq_len(n_iter)]
    beta_raw <- acvb_beta_given_w(omega)
    sel <- which(omega > 0.5)
    beta_hat <- rep(0, p)
    if (length(sel) > 0) {
        Xs <- X[, sel, drop = FALSE]
        beta_hat[sel] <- as.numeric(solve_reg(crossprod(Xs), crossprod(Xs, y)))
    }
    list(method = paste0("no_borrow_ACVB_init_", init_method, "_", acceleration), omega = omega, beta_raw = beta_raw, 
        beta_hat = beta_hat, pi_prior = pi_prior, sigma_beta_sq = sigma_beta_sq, init_method = init_method, n_iter = n_iter, 
        converged = converged, convergence_reason = convergence_reason, d_omega_final = d_omega, fast_cached_update = fast_cached_update, 
        acvb_acceleration = acceleration, acvb_mean_iter_time_sec = mean(iter_time_sec), acvb_max_iter_time_sec = max(iter_time_sec), 
        acvb_min_downdate_denom = min_downdate_denom, acvb_n_denom_clamped = n_denom_clamped)
}

`fit_oracle_trans_lasso_simple` <- function (dat, A0 = integer(0), method_name = NULL) 
{
    p <- dat$p
    n0 <- dat$n
    if (length(A0) == 0) {
        lambda0 <- sqrt(2 * log(p)/n0)
        beta_hat <- fit_lasso_lambda(dat$X, dat$y, lambda = lambda0)
        if (is.null(method_name)) 
            method_name <- "Oracle_Trans_Lasso_emptyA_target_lasso"
        return(list(method = method_name, beta_hat = beta_hat, A0 = A0))
    }
    X_pool <- do.call(rbind, c(list(dat$X), dat$X_source[A0]))
    y_pool <- c(dat$y, unlist(dat$y_source[A0], use.names = FALSE))
    n_pool <- nrow(X_pool)
    lambda_w <- sqrt(2 * log(p)/n_pool)
    lambda_delta <- sqrt(2 * log(p)/n0)
    w_hat <- fit_lasso_lambda(X_pool, y_pool, lambda = lambda_w)
    delta_hat <- fit_lasso_lambda(dat$X, dat$y - as.numeric(dat$X %*% w_hat), lambda = lambda_delta)
    beta_hat <- w_hat + delta_hat
    if (is.null(method_name)) 
        method_name <- paste0("Oracle_Trans_Lasso_simple_A_", paste(A0, collapse = "_"))
    list(method = method_name, beta_hat = beta_hat, A0 = A0)
}

`fit_trans_lasso_ranked_fallback` <- function (dat) 
{
    n0 <- dat$n
    I_train <- seq_len(floor(n0/2))
    I_valid <- setdiff(seq_len(n0), I_train)
    dat_train <- make_target_subset_for_trans(dat, I_train)
    rank_info <- rank_sources_by_marginal_contrast(dat_train)
    ord <- rank_info$rank
    candidates <- vector("list", dat$K + 1)
    candidates[[1]] <- integer(0)
    for (l in seq_len(dat$K)) candidates[[l + 1]] <- ord[seq_len(l)]
    beta_list <- lapply(seq_along(candidates), function(i) {
        fit_oracle_trans_lasso_simple(dat_train, A0 = candidates[[i]], method_name = paste0("fallback_candidate_", 
            i - 1))$beta_hat
    })
    val_mse <- sapply(beta_list, function(b) {
        mean((dat$y[I_valid] - as.numeric(dat$X[I_valid, , drop = FALSE] %*% b))^2)
    })
    best <- which.min(val_mse)
    list(method = "Trans_Lasso_ranked_fallback", beta_hat = beta_list[[best]], chosen_A = candidates[[best]], rank_info = rank_info, 
        validation_mse = val_mse, note = "Official Trans.lasso was not found; used ranked validation fallback.")
}

`get_beta_hat` <- function (alpha, omega, gamma, dat, S0bar_inv) 
{
    omega <- clamp(as.numeric(omega), 0, 1)
    gamma <- clamp(as.numeric(gamma), 0, 1)
    ns_eff <- get_ns_eff(dat)
    S1bar_inv <- (ns_eff * dat$x_norm2)/dat$n
    S_inv <- S0bar_inv + (S1bar_inv - S0bar_inv) * gamma
    bbar <- as.numeric(as.matrix(dat$B) %*% alpha)
    rhs <- omega * dat$Xty + gamma * S1bar_inv * bbar
    factor <- .trace_woodbury_factor(z = omega, diag_add = S_inv + dat$x_norm2 * omega * (1 - omega), dat = dat, 
        ridge = 1e-10)
    omega * .trace_woodbury_apply(factor, rhs)
}

`get_ns_eff` <- function (dat) 
{
    configured <- dat$trace_source_precision_ns_eff
    if (!is.null(configured) && length(configured) == 1L && is.finite(configured) && configured > 0) {
        return(as.numeric(configured))
    }
    trace_get_raw_ns_eff(dat)
}

`make_init_alpha` <- function (dat, init_alpha = "uniform", alpha_fixed = NULL, alpha_method = "old") 
{
    if (init_alpha != "uniform" || alpha_method != "old") 
        stop("Only learned, uniform-initialized alpha is supported.")
    rep(1/dat$K, dat$K)
}

`make_init_gamma_borrow` <- function (dat, ...) 
rep(0, dat$p)

`make_init_omega` <- function (dat, init_method = "full", init_omega = NULL, ...) 
{
    if (init_method == "full") 
        return(rep(1, dat$p))
    if (init_method == "custom" && length(init_omega) == dat$p) 
        return(init_omega)
    stop("Real-data fits require full or explicit custom initialization.")
}

`make_target_subset_for_trans` <- function (dat, idx) 
{
    out <- dat
    out$X <- dat$X[idx, , drop = FALSE]
    out$y <- dat$y[idx]
    out$n <- length(idx)
    out$X_all <- do.call(rbind, c(list(out$X), dat$X_source))
    out$y_all <- c(out$y, unlist(dat$y_source, use.names = FALSE))
    out$n_vec_trans <- c(out$n, dat$n_source_raw)
    out
}

`print_elbo_trace_row` <- function (row) 
{
    cat(sprintf("[%s/%s] iter %3d | ELBO=% .6e | dELBO=% .3e | rel=% .3e | dpar=% .3e | size05=%3d | alpha=(%.3f, %.3f, %.3f)\n", 
        row$alpha_method[1], row$stage[1], row$iter[1], row$elbo[1], row$delta_elbo[1], row$rel_delta_elbo[1], 
        row$max_param_change[1], row$size05[1], row$alpha_1[1], row$alpha_2[1], row$alpha_3[1]))
}

`rank_sources_by_marginal_contrast` <- function (dat, t_star = NULL) 
{
    K <- dat$K
    p <- dat$p
    if (is.null(t_star)) 
        t_star <- max(1, floor(min(dat$n_vec_trans)^(3/4)))
    t_star <- min(t_star, p)
    m0 <- as.numeric(crossprod(dat$X, dat$y))/dat$n
    R <- rep(NA_real_, K)
    for (k in seq_len(K)) {
        Xk <- dat$X_source[[k]]
        yk <- dat$y_source[[k]]
        mk <- as.numeric(crossprod(Xk, yk))/nrow(Xk)
        Delta <- mk - m0
        idx <- order(abs(Delta), decreasing = TRUE)[seq_len(t_star)]
        R[k] <- sum(Delta[idx]^2)
    }
    list(R = R, rank = order(R), t_star = t_star)
}

`safe_log` <- function (x, eps = 1e-12) 
log(pmax(x, eps))

`solve_reg` <- function (A, b) 
{
    ridge_seq <- c(1e-10, 1e-08, 1e-06, 1e-04, 0.01)
    for (r in ridge_seq) {
        out <- tryCatch(solve(A + diag(r, nrow(A)), b), error = function(e) NULL)
        if (!is.null(out)) 
            return(out)
    }
    qr.solve(A + diag(0.01, nrow(A)), b)
}

`solve_spd` <- function (A, B) 
{
    ridge_seq <- c(0, 1e-10, 1e-08, 1e-06, 1e-04, 0.01)
    for (r in ridge_seq) {
        out <- tryCatch({
            R <- chol(A + diag(r, nrow(A)))
            backsolve(R, forwardsolve(t(R), B))
        }, error = function(e) NULL)
        if (!is.null(out)) 
            return(out)
    }
    solve_reg(A, B)
}

`trace_get_raw_ns_eff` <- function (dat) 
{
    as.numeric(.trace_get_ns_eff_core(dat))
}

`trace_target_only_fixed_steps` <- function (dat, omega, S0bar_inv, rho_z, A_y = 1, B_y = 1, n_steps = 5L, tol = 1e-06, verbose = FALSE) 
{
    stopifnot(length(omega) == dat$p, length(n_steps) == 1L, is.finite(n_steps), n_steps >= 1, n_steps == as.integer(n_steps), 
        length(tol) == 1L, is.finite(tol), tol > 0, is.finite(rho_z), rho_z > 0, rho_z < 1)
    n_steps <- as.integer(n_steps)
    initial_omega <- omega
    rows <- vector("list", n_steps)
    target_map <- function(w) .update_omega_all_woodbury(omega = w, gamma = rep(0, dat$p), bbar = rep(0, dat$p), 
        dat = dat, S0bar_inv = S0bar_inv, A_y = A_y, B_y = B_y, E_logit_rho_z = qlogis(rho_z))
    for (iter in seq_len(n_steps)) {
        proposed <- target_map(omega)
        if (!all(is.finite(proposed))) 
            stop("Nonfinite fixed-step target update.")
        delta <- max(abs(proposed - omega))
        omega <- proposed
        rows[[iter]] <- data.frame(iter = iter, max_omega_change = delta, raw_fixed_point_residual_before = delta, 
            endpoint_fixed_point_residual = NA_real_, damping_active = FALSE, damping_factor = 1, size05 = sum(omega > 
                0.5), sum_omega = sum(omega), gamma_max = 0, alpha_updated = FALSE)
        if (verbose) 
            cat(sprintf("[target fixed] sweep=%d/%d | change=%.4e | size05=%d | gamma=0\n", iter, n_steps, delta, 
                sum(omega > 0.5)))
    }
    endpoint_proposed <- target_map(omega)
    if (!all(is.finite(endpoint_proposed))) 
        stop("Nonfinite fixed-step target endpoint residual.")
    residual <- max(abs(endpoint_proposed - omega))
    rows[[n_steps]]$endpoint_fixed_point_residual <- residual
    list(omega = omega, initial_omega = initial_omega, converged = residual < tol, phase_complete = TRUE, mode = "fixed_steps", 
        n_iter = n_steps, n_steps = n_steps, map_evaluations = n_steps + 1L, max_change = residual, tol = tol, 
        damping_trigger_iter = NA_integer_, damping_active = FALSE, trace = do.call(rbind, rows))
}

`trace_target_only_initialize` <- function (dat, omega, S0bar_inv, rho_z, A_y = 1, B_y = 1, max_iter = 1000L, tol = 1e-06, verbose = FALSE) 
{
    stopifnot(length(omega) == dat$p, length(max_iter) == 1L, is.finite(max_iter), max_iter >= 1, length(tol) == 
        1L, is.finite(tol), tol > 0, is.finite(rho_z), rho_z > 0, rho_z < 1)
    max_iter <- as.integer(max_iter)
    rows <- vector("list", max_iter)
    converged <- FALSE
    initial_omega <- omega
    lag2 <- NULL
    damping_active <- FALSE
    damping_trigger_iter <- NA_integer_
    for (iter in seq_len(max_iter)) {
        proposed <- .update_omega_all_woodbury(omega = omega, gamma = rep(0, dat$p), bbar = rep(0, dat$p), dat = dat, 
            S0bar_inv = S0bar_inv, A_y = A_y, B_y = B_y, E_logit_rho_z = qlogis(rho_z))
        if (!all(is.finite(proposed))) 
            stop("Nonfinite target-only update.")
        residual <- max(abs(proposed - omega))
        cycle_distance <- if (is.null(lag2)) 
            NA_real_
        else max(abs(proposed - lag2))
        if (!damping_active && iter >= 3L && residual > max(100 * tol, 1e-04) && is.finite(cycle_distance) && cycle_distance <= 
            max(1e-06, 0.001 * residual)) {
            damping_active <- TRUE
            damping_trigger_iter <- iter
            if (verbose) 
                cat(sprintf("[target only] detected two-cycle at iter=%d; activating eta=0.5\n", iter))
        }
        eta <- if (damping_active) 
            0.5
        else 1
        converged <- residual < tol
        omega_next <- if (converged) 
            omega
        else omega + eta * (proposed - omega)
        step_delta <- max(abs(omega_next - omega))
        rows[[iter]] <- data.frame(iter = iter, max_omega_change = step_delta, raw_fixed_point_residual = residual, 
            cycle_distance = cycle_distance, damping_active = damping_active, damping_factor = eta, size05 = sum(omega_next > 
                0.5), sum_omega = sum(omega_next), gamma_max = 0, alpha_updated = FALSE)
        lag2 <- omega
        omega <- omega_next
        if (verbose && (iter == 1L || iter%%20L == 0L || converged)) 
            cat(sprintf("[target only] iter=%d | raw residual=%.4e | size05=%d | eta=%.1f | gamma=0\n", iter, residual, 
                sum(omega > 0.5), eta))
        if (converged) 
            break
    }
    list(omega = omega, initial_omega = initial_omega, converged = converged, n_iter = iter, max_change = residual, 
        tol = tol, damping_trigger_iter = damping_trigger_iter, damping_active = damping_active, trace = do.call(rbind, 
            rows[seq_len(iter)]))
}

`update_alpha_em_qp` <- function (alpha_old, omega, gamma, dat, sigma_beta_sq = 100, active_tol = 1e-10, ridge_M = 1e-10, ridge_H = 1e-10) 
{
    alpha_old <- .normalize_simplex(alpha_old, K = dat$K)
    terms <- .compute_alpha_em_qp_terms(omega = omega, gamma = gamma, dat = dat, sigma_beta_sq = sigma_beta_sq, 
        ridge_M = ridge_M, ridge_H = ridge_H)
    sol <- .solve_simplex_qp_global(H = terms$H, g = terms$g, tol = active_tol, ridge = ridge_H)
    list(alpha = sol$alpha, accepted = TRUE, step_size = NA_real_, qp_old = .alpha_qp_objective(alpha_old, terms$H, 
        terms$g), qp_new = .alpha_qp_objective(sol$alpha, terms$H, terms$g), H = terms$H, g = terms$g, active_set = sol$active_set, 
        removed_set = sol$removed_set, lambda = sol$lambda, mu = sol$mu, kkt_min_mu_inactive = sol$kkt_min_mu_inactive, 
        solver = sol$solver, strategy = "em_qp")
}

`update_alpha_old` <- function (alpha, omega, gamma, dat, A_y, B_y, S0bar_inv, eps_eta = 1e-06) 
{
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
    eta_raw <- 1 - alpha * (A_y + dat$n/2) * numer/denom
    eta <- pmax(eps_eta, eta_raw)
    alpha_new <- eta/sum(eta)
    list(alpha = alpha_new)
}

`update_alpha_preliminary_huber` <- function (alpha_old, omega, dat, sigma_beta_sq = 100, huber_tau = NULL, huber_tau_rule = c("fixed", "diverging"), 
    huber_tau_multiplier = 1, tol = 1e-10, max_iter = 10000L) 
{
    alpha_old <- .normalize_simplex(alpha_old, K = dat$K)
    terms <- .compute_alpha_preliminary_qp_terms(omega = omega, dat = dat, sigma_beta_sq = sigma_beta_sq)
    tau_info <- .resolve_trace_huber_tau(dat = dat, omega = omega, huber_tau = huber_tau, huber_tau_rule = huber_tau_rule, 
        huber_tau_multiplier = huber_tau_multiplier)
    ns_eff <- max(as.numeric(if (!is.null(dat[["n_vec", exact = TRUE]])) mean(dat[["n_vec", exact = TRUE]]) else dat$ns)[1], 1)
    weight <- as.numeric(omega)
    Bmat <- as.matrix(dat$B)
    old_eval <- .trace_huber_objective_gradient(alpha_old, Bmat, terms$beta_target, weight, tau_info$tau)
    sol <- .solve_trace_huber_simplex(alpha_start = alpha_old, Bmat = Bmat, beta_target = terms$beta_target, weight = weight, 
        tau = tau_info$tau, tol = tol, max_iter = max_iter)
    list(alpha = sol$alpha, accepted = TRUE, objective_old = old_eval$objective, objective_new = sol$objective, 
        qp_old = old_eval$objective, qp_new = sol$objective, beta_target = terms$beta_target, weight = weight, 
        huber_tau = tau_info$tau, huber_tau_rule = tau_info$rule, huber_s_reference = tau_info$s_reference, huber_iterations = sol$iterations, 
        huber_converged = sol$converged, huber_kkt_residual = sol$kkt_residual, active_set = sol$active_set, removed_set = setdiff(seq_len(dat$K), 
            sol$active_set), solver = sol$solver, strategy = "preliminary_huber")
}

`update_alpha_preliminary_qp` <- function (alpha_old, omega, dat, sigma_beta_sq = 100, active_tol = 1e-10, ridge_M = 1e-10, ridge_H = 1e-10) 
{
    alpha_old <- .normalize_simplex(alpha_old, K = dat$K)
    terms <- .compute_alpha_preliminary_qp_terms(omega = omega, dat = dat, sigma_beta_sq = sigma_beta_sq, ridge_M = ridge_M, 
        ridge_H = ridge_H)
    sol <- .solve_simplex_qp_global(H = terms$H, g = terms$g, tol = active_tol, ridge = ridge_H)
    list(alpha = sol$alpha, accepted = TRUE, qp_old = .alpha_qp_objective(alpha_old, terms$H, terms$g), qp_new = .alpha_qp_objective(sol$alpha, 
        terms$H, terms$g), H = terms$H, g = terms$g, beta_target = terms$beta_target, active_set = sol$active_set, 
        removed_set = sol$removed_set, lambda = sol$lambda, mu = sol$mu, kkt_min_mu_inactive = sol$kkt_min_mu_inactive, 
        solver = sol$solver, strategy = "preliminary_target_matching_qp")
}

`update_alpha_with_elbo_accept` <- function (alpha_old, omega, gamma, dat, prior_par_borrow, sigma_beta_sq = 100, A_y = 1, B_y = 1, S0bar_inv = NULL, 
    alpha_eta_old = NULL, include_source_normalization = TRUE, accept_tol = 1e-08, accept_objective = c("old_dirichlet_surrogate", 
        "point_elbo"), line_search_steps = c(1)) 
{
    accept_objective <- match.arg(accept_objective)
    alpha_old <- .normalize_simplex(alpha_old, K = dat$K)
    eval_alpha_objective <- function(alpha_eval, eta_eval = NULL) {
        if (identical(accept_objective, "old_dirichlet_surrogate")) {
            tmp_fit <- list(alpha = .normalize_simplex(alpha_eval, K = dat$K), omega = omega, gamma = gamma, alpha_eta = eta_eval, 
                sigma_beta_sq = sigma_beta_sq, alpha_update_strategy = "elbo_accept")
            compute_old_alpha_dirichlet_criterion(fit = tmp_fit, dat = dat, prior_par_borrow = prior_par_borrow, 
                sigma_beta_sq = sigma_beta_sq, A_y = A_y, B_y = B_y, include_source_normalization = include_source_normalization)$old_alpha_elbo
        }
        else {
            .alpha_elbo_value(alpha = alpha_eval, omega = omega, gamma = gamma, dat = dat, prior_par_borrow = prior_par_borrow, 
                sigma_beta_sq = sigma_beta_sq, A_y = A_y, B_y = B_y, include_source_normalization = include_source_normalization)
        }
    }
    elbo_old <- eval_alpha_objective(alpha_old, alpha_eta_old)
    if (is.null(S0bar_inv)) 
        S0bar_inv <- 1/sigma_beta_sq
    raw_update <- update_alpha_old(alpha_old, omega, gamma, dat, A_y, B_y, S0bar_inv)
    alpha_cand <- .normalize_simplex(raw_update$alpha, K = dat$K)
    eta_cand <- .extract_alpha_eta(raw_update)
    best_alpha <- alpha_old
    best_elbo <- elbo_old
    best_step <- 0
    best_eta <- alpha_eta_old
    for (lambda in line_search_steps) {
        alpha_try <- (1 - lambda) * alpha_old + lambda * alpha_cand
        alpha_try <- .normalize_simplex(alpha_try, K = dat$K)
        eta_try <- if (abs(lambda - 1) < 1e-12) 
            eta_cand
        else alpha_eta_old
        elbo_try <- eval_alpha_objective(alpha_try, eta_try)
        if (is.finite(elbo_try) && elbo_try > best_elbo) {
            best_alpha <- alpha_try
            best_elbo <- elbo_try
            best_step <- lambda
            best_eta <- eta_try
        }
        if (is.finite(elbo_try) && elbo_try >= elbo_old - abs(accept_tol)) {
            return(list(alpha = alpha_try, eta = eta_try, raw_update = raw_update, accepted = TRUE, step_size = lambda, 
                elbo_old = elbo_old, elbo_new = elbo_try, alpha_proposal = alpha_cand, eta_proposal = eta_cand, 
                accept_objective = accept_objective, strategy = "elbo_accept"))
        }
    }
    list(alpha = alpha_old, eta = alpha_eta_old, raw_update = raw_update, accepted = FALSE, step_size = 0, elbo_old = elbo_old, 
        elbo_new = elbo_old, best_rejected_elbo = best_elbo, best_rejected_alpha = best_alpha, best_rejected_eta = best_eta, 
        best_rejected_step = best_step, alpha_proposal = alpha_cand, eta_proposal = eta_cand, accept_objective = accept_objective, 
        strategy = "elbo_accept")
}

`update_gamma_j` <- function (j, omega, gamma, bbar, dat, S0bar_inv, A_y, B_y, E_logit_rho_c) 
{
    p <- length(omega)
    XtX <- dat$XtX
    Xty <- dat$Xty
    x_norm2 <- dat$x_norm2
    yTy <- dat$yTy
    ns_eff <- get_ns_eff(dat)
    S1bar_inv <- (ns_eff * x_norm2)/dat$n
    S_inv_all <- S0bar_inv + (S1bar_inv - S0bar_inv) * gamma
    idx <- setdiff(seq_len(p), j)
    wj <- omega[j]
    if (length(idx) == 0) {
        cTc <- 0
        base_quad <- 0
        alpha_minus <- 0
        R0 <- wj * Xty[j]
    }
    else {
        w <- omega[idx]
        S_inv_m <- S_inv_all[idx]
        v_m <- gamma[idx] * S1bar_inv[idx] * bbar[idx]
        r_minus <- w * Xty[idx] + v_m
        ccol <- w * XtX[idx, j]
        crow <- XtX[j, idx] * w
        XtX_mm <- XtX[idx, idx, drop = FALSE]
        M_minus <- XtX_mm * (w %o% w)
        diag(M_minus) <- diag(XtX_mm) * w
        M_minus <- M_minus + diag(S_inv_m, length(idx)) + diag(1e-10, length(idx))
        SOL <- solve_spd(M_minus, cbind(r_minus, ccol))
        sol_r <- SOL[, 1]
        sol_c <- SOL[, 2]
        cTc <- sum(ccol * sol_c)
        base_quad <- sum(r_minus * sol_r)
        alpha_minus <- sum(gamma[idx] * S1bar_inv[idx] * bbar[idx]^2)
        R0 <- wj * (Xty[j] - sum(crow * sol_r))
    }
    D0j <- S0bar_inv + wj * x_norm2[j] - wj^2 * cTc
    D1j <- S1bar_inv[j] + wj * x_norm2[j] - wj^2 * cTc
    D0j <- pmax(D0j, 1e-10)
    D1j <- pmax(D1j, 1e-10)
    v1j <- S1bar_inv[j] * bbar[j]
    R1 <- R0 + v1j
    Tc0 <- yTy - base_quad - (R0^2)/D0j + alpha_minus
    Tc1 <- yTy - base_quad - (R1^2)/D1j + alpha_minus + (S1bar_inv[j] * bbar[j]^2)
    Tc0 <- pmax(Tc0, 1e-10)
    Tc1 <- pmax(Tc1, 1e-10)
    Delta <- 0.5 * safe_log(S1bar_inv[j]/S0bar_inv) - 0.5 * safe_log(1 + (S1bar_inv[j] - S0bar_inv)/D0j) - (A_y + 
        dat$n/2) * (safe_log(B_y + 0.5 * Tc1) - safe_log(B_y + 0.5 * Tc0))
    expit(E_logit_rho_c + Delta)
}

`update_omega_j_old` <- function (j, omega, gamma, bbar, dat, S0bar_inv, A_y, B_y, E_logit_rho_z) 
{
    idx <- setdiff(seq_along(omega), j)
    w <- omega[idx]
    ns_eff <- get_ns_eff(dat)
    S1bar_inv <- (ns_eff * dat$x_norm2)/dat$n
    S_inv_all <- S0bar_inv + (S1bar_inv - S0bar_inv) * gamma
    S_inv_j <- S_inv_all[j]
    S_inv_m <- S_inv_all[idx]
    v_all <- gamma * S1bar_inv * bbar
    v_j <- v_all[j]
    v_m <- v_all[idx]
    XtX_mm <- dat$XtX[idx, idx, drop = FALSE]
    M_minus <- XtX_mm * (w %o% w)
    diag(M_minus) <- diag(XtX_mm) * w
    M_minus <- M_minus + diag(S_inv_m, length(idx)) + diag(1e-10, length(idx))
    u <- w * dat$Xty[idx]
    ccol <- w * dat$XtX[idx, j]
    SOL <- solve_spd(M_minus, cbind(u, v_m, ccol))
    sol_u <- SOL[, 1]
    sol_v <- SOL[, 2]
    sol_c <- SOL[, 3]
    D_j <- pmax(S_inv_j + dat$x_norm2[j] - sum(ccol * sol_c), 1e-10)
    A_j <- dat$Xty[j] - sum(ccol * sol_u)
    B_j <- v_j - sum(ccol * sol_v)
    source_quad <- sum(gamma * S1bar_inv * bbar^2)
    T0 <- dat$yTy - sum(u * sol_u) - 2 * sum(u * sol_v) - sum(v_m * sol_v) - (gamma[j] * S1bar_inv[j] * bbar[j]^2) + source_quad
    T1 <- T0 - (A_j + B_j)^2/D_j + (gamma[j] * S1bar_inv[j] * bbar[j]^2)
    lin <- E_logit_rho_z - 0.5 * safe_log(D_j) + 0.5 * safe_log(S_inv_j) - (A_y + dat$n/2) * (safe_log(B_y + 0.5 * 
        pmax(T1, 1e-10)) - safe_log(B_y + 0.5 * pmax(T0, 1e-10)))
    clamp(expit(lin), .Machine$double.eps, 1 - .Machine$double.eps)
}
