
## =========================================================
## TRACE Simulation Setting 1: coordinate-level contamination
## Function file: methods, metrics, and replicate runner
## =========================================================

.get_script_dir <- function() {
  ofiles <- Filter(function(x) is.character(x) && length(x) == 1L && nzchar(x),
                   lapply(sys.frames(), function(x) x$ofile))
  ofile <- if (length(ofiles)) tail(ofiles, 1L)[[1L]] else NULL
  if (!is.null(ofile)) return(dirname(normalizePath(ofile)))
  getwd()
}

.this_dir <- .get_script_dir()
source(file.path(.this_dir, "dgp.R"))

.parent_directories <- function(path, max_up = 10L) {
  path <- normalizePath(path, mustWork = FALSE)
  out <- character(0)

  for (i in seq_len(max_up + 1L)) {
    out <- c(out, path)
    parent <- dirname(path)
    if (identical(parent, path)) break
    path <- parent
  }

  unique(out)
}

.source_borrow_core <- function(this_dir = .this_dir) {
  ## The 01_setting folder is often nested as
  ## borrowcvb_systematic_simulation/borrow_trader_two_groups/01_setting/
  ## 01_current_setting_with_trader.  The project-level _common directory is
  ## therefore three levels above this file.  Search every ancestor rather
  ## than assuming a fixed one- or two-level layout.
  core_env <- Sys.getenv("BORROW_CORE_PATH", unset = "")
  ancestors <- .parent_directories(this_dir, max_up = 10L)

  candidates <- c(
    if (nzchar(core_env)) core_env else character(0),
    unlist(lapply(
      ancestors,
      function(x) c(
        file.path(x, "_common", "borrow_functions_core.R"),
        file.path(x, "borrow_functions_core.R")
      )
    ), use.names = FALSE)
  )

  candidates <- unique(normalizePath(candidates, mustWork = FALSE))
  hit <- candidates[file.exists(candidates)]

  if (length(hit) == 0L) {
    expected_project_path <- file.path(
      dirname(dirname(dirname(this_dir))),
      "_common",
      "borrow_functions_core.R"
    )

    stop(
      paste0(
        "Cannot find borrow_functions_core.R.\n\n",
        "This v9.1 file searched the current folder and up to 10 parent folders.\n",
        "For your current layout, the most likely expected path is:\n",
        normalizePath(expected_project_path, mustWork = FALSE),
        "\n\n",
        "Either place the original core file at that location, or set its exact path before sourcing:\n",
        "Sys.setenv(BORROW_CORE_PATH = '/absolute/path/to/borrow_functions_core.R')\n",
        "source('functions.R')\n\n",
        "Paths searched:\n",
        paste(candidates, collapse = "\n")
      )
    )
  }

  core_file <- hit[1]
  source(core_file, local = .GlobalEnv)

  required_core_functions <- c(
    "make_init_gamma_borrow",
    "get_ns_eff",
    "clamp",
    "solve_spd",
    "make_init_omega",
    "make_init_alpha"
  )
  missing_core_functions <- required_core_functions[
    !vapply(
      required_core_functions,
      exists,
      logical(1),
      mode = "function",
      inherits = TRUE
    )
  ]

  if (length(missing_core_functions) > 0L) {
    stop(
      "The core file was found but is not the expected Borrow-CVB core. Missing functions: ",
      paste(missing_core_functions, collapse = ", "),
      "\nCore file used: ", core_file
    )
  }

  message("[TRACE core] ", core_file)
  invisible(core_file)
}

.core_file_used <- .source_borrow_core(.this_dir)

.make_init_omega_shared <- make_init_omega
make_init_omega <- function(dat, init_method = "full", ...) {
  if (identical(init_method, "full")) return(rep(1, dat$p))
  .make_init_omega_shared(dat = dat, init_method = init_method, ...)
}


## Override the shared-core prediction helper so every method uses the same
## observed-response test MSE.
source(file.path(.this_dir, "prediction_metrics_v9.R"), local = .GlobalEnv)

suppressPackageStartupMessages({
  library(glmnet)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

## Preserve the source sample size reported by the DGP, but allow TRACE's
## finite-sample moment surrogate to use a tempered precision.  This targets
## the soft-gamma failure mode: a moderate gamma multiplied by ns >> n can
## otherwise make a wrong source center dominate the target likelihood.
.trace_get_ns_eff_core <- get_ns_eff

trace_get_raw_ns_eff <- function(dat) {
  as.numeric(.trace_get_ns_eff_core(dat))
}

configure_trace_source_precision <- function(
    dat,
    mode = c("high_dimensional_logp", "raw_ns", "target_capped"),
    cap_ratio = 1) {
  mode <- match.arg(mode)
  cap_ratio <- as.numeric(cap_ratio)
  if (length(cap_ratio) != 1L || !is.finite(cap_ratio) || cap_ratio <= 0) {
    stop("source_precision_cap_ratio must be one finite positive number.")
  }

  raw_ns <- trace_get_raw_ns_eff(dat)
  precision_ns <- if (mode == "high_dimensional_logp") {
    if (dat$p < 2) stop("High dimensional precision requires p >= 2.")
    raw_ns / log(dat$p)
  } else if (mode == "target_capped") {
    min(raw_ns, cap_ratio * dat$n)
  } else {
    raw_ns
  }

  dat$trace_source_precision_mode <- mode
  dat$trace_source_precision_cap_ratio <- cap_ratio
  dat$trace_source_precision_ns_raw <- raw_ns
  dat$trace_source_precision_ns_eff <- precision_ns
  dat$trace_source_precision_power <- precision_ns / raw_ns
  dat
}

## All implemented S1 precision calculations call get_ns_eff().  Once a data
## object has been configured above, they consistently use the tempered value.
## Raw ns is retained separately. The revised default rho_c is constant 0.5;
## Huber tau uses raw ns, never the effective precision coefficient.
get_ns_eff <- function(dat) {
  configured <- dat$trace_source_precision_ns_eff
  if (!is.null(configured) && length(configured) == 1L &&
      is.finite(configured) && configured > 0) {
    return(as.numeric(configured))
  }
  trace_get_raw_ns_eff(dat)
}

.trace_apply_public_method_labels <- function(performance) {
  if (is.null(performance$method_key)) performance$method_key <- performance$method
  missing_key <- is.na(performance$method_key) | !nzchar(as.character(performance$method_key))
  performance$method_key[missing_key] <- performance$method[missing_key]
  labels <- c(
    borrow_cvb_learn_alpha = "TRACE-EM",
    borrow_cvb_em_alpha = "TRACE-EM",
    borrow_cvb_oracle_alpha = "Oracle TRACE",
    borrow_cvb_old_alpha_elboacc = "Legacy alpha diagnostic",
    borrow_cvb_fixed_uniform = "TRACE fixed-uniform diagnostic",
    no_borrow_acvb = "ACVB",
    trans_lasso = "Trans-Lasso",
    trader = "TRADER"
  )
  hit <- labels[as.character(performance$method_key)]
  performance$method <- ifelse(is.na(hit), as.character(performance$method_key), unname(hit))
  performance
}

## Keep the original initializer from borrow_functions_core.R and add the
## theory-aligned option gamma^(0) = 0.  The zero vector is returned exactly,
## without clamping, because the first omega sweep must contain no source term.
## The older ns_quarter option is retained only for backward compatibility.
.make_init_gamma_borrow_core <- make_init_gamma_borrow
make_init_gamma_borrow <- function(dat,
                                   init_gamma = "clean",
                                   alpha = NULL,
                                   gamma_ebic = 0.5) {
  if (is.character(init_gamma) && length(init_gamma) == 1 &&
      init_gamma %in% c("zero", "theory_zero", "target_only_zero")) {
    return(rep(0, dat$p))
  }

  if (is.character(init_gamma) && length(init_gamma) == 1 &&
      init_gamma %in% c("ns_quarter", "ns_quarter_odds", "ns_1_4_odds")) {
    ns_eff <- get_ns_eff(dat)
    odds <- ns_eff^(1 / 4)
    gamma0 <- odds / (1 + odds)
    return(rep(clamp(gamma0, .Machine$double.eps, 1 - .Machine$double.eps), dat$p))
  }

  .make_init_gamma_borrow_core(
    dat = dat,
    init_gamma = init_gamma,
    alpha = alpha,
    gamma_ebic = gamma_ebic
  )
}


## =========================================================
## Fast Borrow-CVB overrides
## =========================================================
## Main speedups:
##   1. Use an exact n-dimensional Woodbury system for paper-block sweeps.
##   2. Optional candidate screening: run Borrow-CVB only on a moderate active pool.
##   3. Update alpha only every few iterations.
## The old p-dimensional cached inverse is retained as a validation fallback.

.solve_leave_one_cached <- function(Minv, j, a) {
  p <- nrow(Minv)
  if (p <= 1) return(numeric(0))
  idx <- setdiff(seq_len(p), j)
  g <- Minv[j, j]
  if (!is.finite(g) || abs(g) < 1e-12) {
    stop("Cached inverse is numerically singular at coordinate ", j)
  }
  f <- Minv[idx, j]
  as.numeric(Minv[idx, idx, drop = FALSE] %*% a - f * (sum(f * a) / g))
}

.build_A0_cached <- function(omega, gamma, dat, S0bar_inv) {
  ns_eff <- get_ns_eff(dat)
  S1bar_inv <- (ns_eff * dat$x_norm2) / dat$n
  S_inv_all <- S0bar_inv + (S1bar_inv - S0bar_inv) * gamma
  A0 <- dat$XtX * (omega %o% omega)
  diag(A0) <- dat$x_norm2 * omega + S_inv_all
  A0 + diag(1e-10, length(omega))
}

.update_omega_all_cached_inverse <- function(omega, gamma, bbar, dat,
                                             S0bar_inv, A_y, B_y, E_logit_rho_z) {
  p <- length(omega)
  if (p <= 1) return(omega)

  XtX <- dat$XtX
  Xty <- dat$Xty
  x_norm2 <- dat$x_norm2
  yTy <- dat$yTy
  n <- dat$n
  ns_eff <- get_ns_eff(dat)

  S1bar_inv <- (ns_eff * x_norm2) / n
  S_inv_all <- S0bar_inv + (S1bar_inv - S0bar_inv) * gamma
  ## Paper identities:
  ## E(S^{-1} C) = E(C S^{-1} C) = Gamma S1^{-1}.
  v_all <- (gamma * S1bar_inv) * bbar
  alpha_term <- sum(gamma * S1bar_inv * bbar^2)

  A0 <- .build_A0_cached(omega, gamma, dat, S0bar_inv)
  Minv <- solve_spd(A0, diag(p))

  omega_prop <- omega
  for (j in seq_len(p)) {
    idx <- setdiff(seq_len(p), j)
    w <- omega[idx]

    S_inv_j <- S_inv_all[j]
    S_inv_m <- S_inv_all[idx]
    v_j <- v_all[j]
    v_m <- v_all[idx]

    u <- w * Xty[idx]
    ccol <- w * XtX[idx, j]

    ok <- TRUE
    SOL <- tryCatch({
      cbind(
        .solve_leave_one_cached(Minv, j, u),
        .solve_leave_one_cached(Minv, j, v_m),
        .solve_leave_one_cached(Minv, j, ccol)
      )
    }, error = function(e) {
      ok <<- FALSE
      NULL
    })

    if (!ok) {
      ## Robust fallback for this coordinate only.
      XtX_mm <- XtX[idx, idx, drop = FALSE]
      M_minus <- XtX_mm * (w %o% w)
      diag(M_minus) <- diag(XtX_mm) * w
      M_minus <- M_minus + diag(S_inv_m, length(idx)) + diag(1e-10, length(idx))
      SOL <- solve_spd(M_minus, cbind(u, v_m, ccol))
    }

    sol_u <- SOL[, 1]
    sol_v <- SOL[, 2]
    sol_c <- SOL[, 3]

    D_j <- S_inv_j + x_norm2[j] - sum(ccol * sol_c)
    D_j <- pmax(D_j, 1e-10)

    A_j <- Xty[j] - sum(ccol * sol_u)
    B_j <- v_j - sum(ccol * sol_v)

    T0 <- yTy -
      sum(u * sol_u) -
      2 * sum(u * sol_v) -
      sum(v_m * sol_v) -
      (gamma[j] * S1bar_inv[j] * bbar[j]^2) +
      alpha_term

    Tdiff <- -(A_j + B_j)^2 / pmax(D_j, 1e-12) +
      (gamma[j] * S1bar_inv[j] * bbar[j]^2)

    T1 <- T0 + Tdiff
    T0 <- pmax(T0, 1e-10)
    T1 <- pmax(T1, 1e-10)

    lin <- E_logit_rho_z -
      0.5 * safe_log(D_j) +
      0.5 * safe_log(S_inv_j) -
      (A_y + n / 2) *
      (safe_log(B_y + 0.5 * T1) - safe_log(B_y + 0.5 * T0))

    lin <- lin
    omega_prop[j] <- expit(lin)
  }

  clamp(omega_prop, .Machine$double.eps, 1 - .Machine$double.eps)
}

.update_gamma_all_cached_inverse <- function(omega, gamma, bbar, dat,
                                             S0bar_inv, A_y, B_y, E_logit_rho_c) {
  p <- length(omega)
  if (p <= 1) return(gamma)

  XtX <- dat$XtX
  Xty <- dat$Xty
  x_norm2 <- dat$x_norm2
  yTy <- dat$yTy
  n <- dat$n
  ns_eff <- get_ns_eff(dat)

  S1bar_inv <- (ns_eff * x_norm2) / n
  S_inv_all <- S0bar_inv + (S1bar_inv - S0bar_inv) * gamma
  v_all <- gamma * S1bar_inv * bbar

  A0 <- .build_A0_cached(omega, gamma, dat, S0bar_inv)
  Minv <- solve_spd(A0, diag(p))

  gamma_prop <- gamma
  for (j in seq_len(p)) {
    idx <- setdiff(seq_len(p), j)
    w <- omega[idx]
    wj <- omega[j]

    S_inv_m <- S_inv_all[idx]
    v_m <- v_all[idx]

    u <- w * Xty[idx]
    r_minus <- u + v_m
    ccol <- w * XtX[idx, j]
    crow <- XtX[j, idx] * w

    ok <- TRUE
    SOL <- tryCatch({
      cbind(
        .solve_leave_one_cached(Minv, j, r_minus),
        .solve_leave_one_cached(Minv, j, ccol)
      )
    }, error = function(e) {
      ok <<- FALSE
      NULL
    })

    if (!ok) {
      XtX_mm <- XtX[idx, idx, drop = FALSE]
      M_minus <- XtX_mm * (w %o% w)
      diag(M_minus) <- diag(XtX_mm) * w
      M_minus <- M_minus + diag(S_inv_m, length(idx)) + diag(1e-10, length(idx))
      SOL <- solve_spd(M_minus, cbind(r_minus, ccol))
    }

    sol_r <- SOL[, 1]
    sol_c <- SOL[, 2]
    cTc <- sum(ccol * sol_c)

    D0j <- S0bar_inv + wj * x_norm2[j] - wj^2 * cTc
    D1j <- S1bar_inv[j] + wj * x_norm2[j] - wj^2 * cTc
    D0j <- pmax(D0j, 1e-10)
    D1j <- pmax(D1j, 1e-10)

    R0 <- wj * (Xty[j] - sum(crow * sol_r))
    v1j <- S1bar_inv[j] * bbar[j]
    R1 <- R0 + v1j

    alpha_minus <- sum(gamma[idx] * S1bar_inv[idx] * bbar[idx]^2)

    Tc0 <- yTy - sum(r_minus * sol_r) - (R0^2) / D0j + alpha_minus
    Tc1 <- yTy - sum(r_minus * sol_r) - (R1^2) / D1j + alpha_minus +
      (S1bar_inv[j] * bbar[j]^2)

    Tc0 <- pmax(Tc0, 1e-10)
    Tc1 <- pmax(Tc1, 1e-10)

    ## Full collapsed marginal likelihood determinant contribution:
    ##   -0.5 log |Sbar| - 0.5 log |Sbar^{-1} + Z X'X Z|.
    ## The first term contributes +0.5 log(S1bar_inv[j] / S0bar_inv)
    ## when c_j changes from 0 to 1.
    Delta <- 0.5 * safe_log(S1bar_inv[j] / S0bar_inv) -
      0.5 * safe_log(1 + (S1bar_inv[j] - S0bar_inv) / D0j) -
      (A_y + n / 2) *
      (safe_log(B_y + 0.5 * Tc1) - safe_log(B_y + 0.5 * Tc0))

    lin_g <- E_logit_rho_c + Delta
    lin_g <- lin_g
    gamma_prop[j] <- expit(lin_g)
  }

  clamp(gamma_prop, .Machine$double.eps, 1 - .Machine$double.eps)
}

## =========================================================
## Exact Woodbury backend for TRACE paper-block updates
## =========================================================
## For U = X diag(z), every TRACE precision matrix used below has the form
##
##   M = diag(d) + U' U.
##
## The factorization and all solves are therefore performed through the exact
## Woodbury identity in target-sample dimension n.  The diagonal ridge is put
## in d, so this represents the same regularized M as the dense implementation.
.trace_woodbury_factor <- function(z,
                                   diag_add,
                                   dat,
                                   ridge = 1e-10,
                                   compute_inverse_diagonal = FALSE) {
  z <- as.numeric(z)
  diag_add <- as.numeric(diag_add)
  if (length(z) != dat$p || length(diag_add) != dat$p) {
    stop("z and diag_add must both have length dat$p.")
  }

  d <- pmax(diag_add + ridge, 1e-14)
  d_inv <- 1 / d
  U <- sweep(dat$X, 2L, z, `*`)
  U_scaled <- sweep(U, 2L, sqrt(d_inv), `*`)
  K <- diag(dat$n) + tcrossprod(U_scaled)

  chol_K <- tryCatch(chol(K), error = function(e) NULL)
  solve_K <- function(rhs) {
    rhs <- as.matrix(rhs)
    if (!is.null(chol_K)) {
      backsolve(chol_K, forwardsolve(t(chol_K), rhs))
    } else {
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

  list(
    U = U,
    d = d,
    d_inv = d_inv,
    Kinv_U = Kinv_U,
    h_diag = h_diag,
    solve_K = solve_K,
    logdet = sum(log(d)) + if (!is.null(chol_K)) {
      2 * sum(log(diag(chol_K)))
    } else {
      as.numeric(determinant(K, logarithm = TRUE)$modulus)
    },
    n_system_dimension = dat$n
  )
}

.trace_woodbury_apply <- function(factor, rhs) {
  rhs_was_vector <- is.null(dim(rhs))
  rhs <- as.matrix(rhs)
  d_rhs <- sweep(rhs, 1L, factor$d_inv, `*`)
  middle <- if (!is.null(factor$Kinv_U)) {
    factor$Kinv_U %*% d_rhs
  } else {
    factor$solve_K(factor$U %*% d_rhs)
  }
  correction <- crossprod(factor$U, middle)
  out <- d_rhs - sweep(correction, 1L, factor$d_inv, `*`)
  if (rhs_was_vector) as.numeric(out) else out
}

.trace_leave_one_quad <- function(full_quad, Ha, Hb, h_diag) {
  as.numeric(full_quad - Ha * Hb / h_diag)
}

.update_omega_all_woodbury <- function(omega, gamma, bbar, dat,
                                       S0bar_inv, A_y, B_y, E_logit_rho_z) {
  p <- length(omega)

  omega <- clamp(as.numeric(omega), 0, 1)
  gamma <- clamp(as.numeric(gamma), 0, 1)
  ns_eff <- get_ns_eff(dat)
  S1bar_inv <- (ns_eff * dat$x_norm2) / dat$n
  S_inv <- S0bar_inv + (S1bar_inv - S0bar_inv) * gamma
  u <- omega * dat$Xty
  v <- gamma * S1bar_inv * bbar
  Cmat <- sweep(dat$XtX, 1L, omega, `*`)

  factor <- .trace_woodbury_factor(
    z = omega,
    diag_add = S_inv + dat$x_norm2 * omega * (1 - omega),
    dat = dat,
    compute_inverse_diagonal = TRUE
  )
  solved <- .trace_woodbury_apply(factor, cbind(u, v, Cmat))
  Hu <- solved[, 1L]
  Hv <- solved[, 2L]
  HC <- solved[, 2L + seq_len(p), drop = FALSE]
  h <- factor$h_diag

  uu_minus <- .trace_leave_one_quad(sum(u * Hu), Hu, Hu, h)
  uv_minus <- .trace_leave_one_quad(sum(u * Hv), Hu, Hv, h)
  vv_minus <- .trace_leave_one_quad(sum(v * Hv), Hv, Hv, h)
  HC_diag <- diag(HC)
  cc_minus <- colSums(Cmat * HC) - HC_diag^2 / h
  cu_minus <- as.numeric(crossprod(Cmat, Hu)) - HC_diag * Hu / h
  cv_minus <- as.numeric(crossprod(Cmat, Hv)) - HC_diag * Hv / h

  D <- pmax(S_inv + dat$x_norm2 - cc_minus, 1e-10)
  A <- dat$Xty - cu_minus
  B <- v - cv_minus
  alpha_term <- sum(gamma * S1bar_inv * bbar^2)
  T0 <- dat$yTy - uu_minus - 2 * uv_minus - vv_minus -
    (gamma * S1bar_inv * bbar^2) + alpha_term
  T1 <- T0 - (A + B)^2 / D + (gamma * S1bar_inv * bbar^2)
  T0 <- pmax(T0, 1e-10)
  T1 <- pmax(T1, 1e-10)

  lin <- E_logit_rho_z - 0.5 * safe_log(D) + 0.5 * safe_log(S_inv) -
    (A_y + dat$n / 2) *
    (safe_log(B_y + 0.5 * T1) - safe_log(B_y + 0.5 * T0))
  clamp(expit(lin), .Machine$double.eps, 1 - .Machine$double.eps)
}

.update_gamma_all_woodbury <- function(omega, gamma, bbar, dat,
                                       S0bar_inv, A_y, B_y, E_logit_rho_c) {
  p <- length(gamma)

  omega <- clamp(as.numeric(omega), 0, 1)
  gamma <- clamp(as.numeric(gamma), 0, 1)
  ns_eff <- get_ns_eff(dat)
  S1bar_inv <- (ns_eff * dat$x_norm2) / dat$n
  S_inv <- S0bar_inv + (S1bar_inv - S0bar_inv) * gamma
  v <- gamma * S1bar_inv * bbar
  r <- omega * dat$Xty + v
  Cmat <- sweep(dat$XtX, 1L, omega, `*`)

  factor <- .trace_woodbury_factor(
    z = omega,
    diag_add = S_inv + dat$x_norm2 * omega * (1 - omega),
    dat = dat,
    compute_inverse_diagonal = TRUE
  )
  solved <- .trace_woodbury_apply(factor, cbind(r, Cmat))
  Hr <- solved[, 1L]
  HC <- solved[, 1L + seq_len(p), drop = FALSE]
  h <- factor$h_diag
  HC_diag <- diag(HC)

  rr_minus <- .trace_leave_one_quad(sum(r * Hr), Hr, Hr, h)
  cc_minus <- colSums(Cmat * HC) - HC_diag^2 / h
  cr_minus <- as.numeric(crossprod(Cmat, Hr)) - HC_diag * Hr / h

  D0 <- pmax(S0bar_inv + omega * dat$x_norm2 - omega^2 * cc_minus, 1e-10)
  D1 <- pmax(S1bar_inv + omega * dat$x_norm2 - omega^2 * cc_minus, 1e-10)
  R0 <- omega * (dat$Xty - cr_minus)
  R1 <- R0 + S1bar_inv * bbar
  alpha_minus <- sum(gamma * S1bar_inv * bbar^2) -
    gamma * S1bar_inv * bbar^2
  Tc0 <- pmax(dat$yTy - rr_minus - R0^2 / D0 + alpha_minus, 1e-10)
  Tc1 <- pmax(
    dat$yTy - rr_minus - R1^2 / D1 + alpha_minus + S1bar_inv * bbar^2,
    1e-10
  )

  Delta <- 0.5 * safe_log(S1bar_inv / S0bar_inv) -
    0.5 * safe_log(1 + (S1bar_inv - S0bar_inv) / D0) -
    (A_y + dat$n / 2) *
    (safe_log(B_y + 0.5 * Tc1) - safe_log(B_y + 0.5 * Tc0))
  clamp(expit(E_logit_rho_c + Delta), .Machine$double.eps, 1 - .Machine$double.eps)
}

## Public fast paper-block kernels now use the exact Woodbury implementation.
## The previous p-dimensional cached inverse remains available above for
## numerical equivalence tests and emergency fallback.
.update_omega_all_cached <- .update_omega_all_woodbury
.update_gamma_all_cached <- .update_gamma_all_woodbury

## =========================================================
## Cyclic Gauss--Seidel updates for the current CAVI fixed-point maps
## =========================================================
## Unlike the cached Jacobi routines above, each coordinate proposal is
## written back immediately.  The inverse of the current posterior precision
## matrix is then updated exactly before the next coordinate is evaluated.
## Omega changes produce a symmetric rank-two matrix update; gamma changes
## produce a diagonal rank-one update.  A direct rebuild is retained as a
## numerical fallback, and every new sweep starts from a freshly factorized
## matrix so rank-update roundoff cannot accumulate across iterations.

.inverse_rank1_coordinate_update <- function(Minv,
                                             j,
                                             delta_diag,
                                             denom_floor = 1e-10) {
  if (!is.finite(delta_diag) || abs(delta_diag) < 1e-15) {
    return(list(inverse = Minv, ok = TRUE, min_abs_denom = Inf))
  }

  mj <- Minv[, j]
  denom <- 1 + delta_diag * Minv[j, j]
  if (!is.finite(denom) || abs(denom) < denom_floor) {
    return(list(inverse = Minv, ok = FALSE, min_abs_denom = abs(denom)))
  }

  out <- Minv - (delta_diag / denom) * tcrossprod(mj)
  out <- 0.5 * (out + t(out))
  list(
    inverse = out,
    ok = all(is.finite(out)),
    min_abs_denom = abs(denom)
  )
}

.inverse_rank2_omega_update <- function(Minv,
                                        j,
                                        delta_omega,
                                        omega_before,
                                        dat,
                                        denom_floor = 1e-10) {
  if (!is.finite(delta_omega) || abs(delta_omega) < 1e-15) {
    return(list(inverse = Minv, ok = TRUE, min_abs_denom = Inf))
  }

  p <- length(omega_before)
  ej <- numeric(p)
  ej[j] <- 1
  qj <- omega_before * dat$XtX[, j]
  qj[j] <- 0
  U <- cbind(ej, qj)
  Cmat <- matrix(
    c(delta_omega * dat$x_norm2[j], delta_omega,
      delta_omega, 0),
    nrow = 2L,
    byrow = TRUE
  )

  MU <- Minv %*% U
  small_system <- diag(2L) + Cmat %*% crossprod(U, MU)
  small_det <- as.numeric(det(small_system))
  if (!is.finite(small_det) || abs(small_det) < denom_floor) {
    return(list(inverse = Minv, ok = FALSE, min_abs_denom = abs(small_det)))
  }

  correction <- tryCatch(
    solve(small_system, Cmat),
    error = function(e) NULL
  )
  if (is.null(correction) || any(!is.finite(correction))) {
    return(list(inverse = Minv, ok = FALSE, min_abs_denom = abs(small_det)))
  }

  out <- Minv - MU %*% correction %*% t(MU)
  out <- 0.5 * (out + t(out))
  list(
    inverse = out,
    ok = all(is.finite(out)),
    min_abs_denom = abs(small_det)
  )
}

.update_omega_all_gauss_seidel <- function(omega, gamma, bbar, dat,
                                            S0bar_inv, A_y, B_y,
                                            E_logit_rho_z) {
  p <- length(omega)
  if (p <= 1L) return(clamp(omega, .Machine$double.eps, 1 - .Machine$double.eps))

  omega <- clamp(as.numeric(omega), .Machine$double.eps, 1 - .Machine$double.eps)
  ns_eff <- get_ns_eff(dat)
  S1bar_inv <- (ns_eff * dat$x_norm2) / dat$n
  S_inv_all <- S0bar_inv + (S1bar_inv - S0bar_inv) * gamma
  v_all <- gamma * S1bar_inv * bbar
  alpha_term <- sum(gamma * S1bar_inv * bbar^2)

  rebuild_inverse <- function() {
    solve_spd(
      .build_A0_cached(omega, gamma, dat, S0bar_inv),
      diag(p)
    )
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

    SOL <- tryCatch(
      cbind(
        .solve_leave_one_cached(Minv, j, u),
        .solve_leave_one_cached(Minv, j, v_m),
        .solve_leave_one_cached(Minv, j, ccol)
      ),
      error = function(e) NULL
    )
    if (is.null(SOL) || any(!is.finite(SOL))) {
      Minv <- rebuild_inverse()
      SOL <- cbind(
        .solve_leave_one_cached(Minv, j, u),
        .solve_leave_one_cached(Minv, j, v_m),
        .solve_leave_one_cached(Minv, j, ccol)
      )
    }

    sol_u <- SOL[, 1]
    sol_v <- SOL[, 2]
    sol_c <- SOL[, 3]
    D_j <- pmax(S_inv_j + dat$x_norm2[j] - sum(ccol * sol_c), 1e-10)
    A_j <- dat$Xty[j] - sum(ccol * sol_u)
    B_j <- v_j - sum(ccol * sol_v)
    T0 <- dat$yTy -
      sum(u * sol_u) -
      2 * sum(u * sol_v) -
      sum(v_m * sol_v) -
      (gamma[j] * S1bar_inv[j] * bbar[j]^2) +
      alpha_term
    T1 <- T0 - (A_j + B_j)^2 / D_j +
      (gamma[j] * S1bar_inv[j] * bbar[j]^2)
    T0 <- pmax(T0, 1e-10)
    T1 <- pmax(T1, 1e-10)

    lin <- E_logit_rho_z -
      0.5 * safe_log(D_j) +
      0.5 * safe_log(S_inv_j) -
      (A_y + dat$n / 2) *
      (safe_log(B_y + 0.5 * T1) - safe_log(B_y + 0.5 * T0))
    proposal <- clamp(expit(lin), .Machine$double.eps, 1 - .Machine$double.eps)

    omega_before <- omega
    delta_omega <- proposal - omega_before[j]
    omega[j] <- proposal
    inv_update <- .inverse_rank2_omega_update(
      Minv = Minv,
      j = j,
      delta_omega = delta_omega,
      omega_before = omega_before,
      dat = dat
    )
    Minv <- if (isTRUE(inv_update$ok)) inv_update$inverse else rebuild_inverse()
  }

  clamp(omega, .Machine$double.eps, 1 - .Machine$double.eps)
}

.update_gamma_all_gauss_seidel <- function(omega, gamma, bbar, dat,
                                            S0bar_inv, A_y, B_y,
                                            E_logit_rho_c) {
  p <- length(gamma)
  if (p <= 1L) return(clamp(gamma, .Machine$double.eps, 1 - .Machine$double.eps))

  gamma <- clamp(as.numeric(gamma), .Machine$double.eps, 1 - .Machine$double.eps)
  ns_eff <- get_ns_eff(dat)
  S1bar_inv <- (ns_eff * dat$x_norm2) / dat$n
  S_inv_all <- S0bar_inv + (S1bar_inv - S0bar_inv) * gamma
  v_all <- gamma * S1bar_inv * bbar

  rebuild_inverse <- function() {
    solve_spd(
      .build_A0_cached(omega, gamma, dat, S0bar_inv),
      diag(p)
    )
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

    SOL <- tryCatch(
      cbind(
        .solve_leave_one_cached(Minv, j, r_minus),
        .solve_leave_one_cached(Minv, j, ccol)
      ),
      error = function(e) NULL
    )
    if (is.null(SOL) || any(!is.finite(SOL))) {
      Minv <- rebuild_inverse()
      SOL <- cbind(
        .solve_leave_one_cached(Minv, j, r_minus),
        .solve_leave_one_cached(Minv, j, ccol)
      )
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
    Tc0 <- pmax(
      dat$yTy - sum(r_minus * sol_r) - (R0^2) / D0j + alpha_minus,
      1e-10
    )
    Tc1 <- pmax(
      dat$yTy - sum(r_minus * sol_r) - (R1^2) / D1j + alpha_minus +
        S1bar_inv[j] * bbar[j]^2,
      1e-10
    )

    Delta <- 0.5 * safe_log(S1bar_inv[j] / S0bar_inv) -
      0.5 * safe_log(1 + (S1bar_inv[j] - S0bar_inv) / D0j) -
      (A_y + dat$n / 2) *
      (safe_log(B_y + 0.5 * Tc1) - safe_log(B_y + 0.5 * Tc0))
    proposal <- clamp(
      expit(E_logit_rho_c + Delta),
      1e-6,
      1 - 1e-6
    )

    gamma_old_j <- gamma[j]
    S_old_j <- S_inv_all[j]
    gamma[j] <- proposal
    S_inv_all[j] <- S0bar_inv + (S1bar_inv[j] - S0bar_inv) * proposal
    v_all[j] <- S1bar_inv[j] * proposal * bbar[j]
    inv_update <- .inverse_rank1_coordinate_update(
      Minv = Minv,
      j = j,
      delta_diag = S_inv_all[j] - S_old_j
    )
    Minv <- if (isTRUE(inv_update$ok)) inv_update$inverse else rebuild_inverse()
  }

  clamp(gamma, .Machine$double.eps, 1 - .Machine$double.eps)
}

.update_omega_gamma_interleaved_gauss_seidel <- function(
    omega, gamma, bbar, dat, S0bar_inv, A_y, B_y,
    E_logit_rho_z, E_logit_rho_c) {
  p <- length(omega)
  omega <- clamp(as.numeric(omega), .Machine$double.eps, 1 - .Machine$double.eps)
  gamma <- clamp(as.numeric(gamma), .Machine$double.eps, 1 - .Machine$double.eps)
  if (p <= 1L) return(list(omega = omega, gamma = gamma))

  ns_eff <- get_ns_eff(dat)
  S1bar_inv <- (ns_eff * dat$x_norm2) / dat$n
  S_inv_all <- S0bar_inv + (S1bar_inv - S0bar_inv) * gamma
  v_all <- gamma * S1bar_inv * bbar

  rebuild_inverse <- function() {
    solve_spd(
      .build_A0_cached(omega, gamma, dat, S0bar_inv),
      diag(p)
    )
  }
  Minv <- rebuild_inverse()

  for (j in seq_len(p)) {
    idx <- seq_len(p) != j

    ## omega_j conditional update using the latest omega_{<j}, gamma_{<j}.
    w_m <- omega[idx]
    S_inv_j <- S_inv_all[j]
    S_inv_m <- S_inv_all[idx]
    v_j <- v_all[j]
    v_m <- v_all[idx]
    u <- w_m * dat$Xty[idx]
    ccol <- w_m * dat$XtX[idx, j]
    SOL <- cbind(
      .solve_leave_one_cached(Minv, j, u),
      .solve_leave_one_cached(Minv, j, v_m),
      .solve_leave_one_cached(Minv, j, ccol)
    )
    sol_u <- SOL[, 1]
    sol_v <- SOL[, 2]
    sol_c <- SOL[, 3]
    D_j <- pmax(S_inv_j + dat$x_norm2[j] - sum(ccol * sol_c), 1e-10)
    A_j <- dat$Xty[j] - sum(ccol * sol_u)
    B_j <- v_j - sum(ccol * sol_v)
    alpha_term <- sum(gamma * S1bar_inv * bbar^2)
    T0 <- dat$yTy - sum(u * sol_u) - 2 * sum(u * sol_v) -
      sum(v_m * sol_v) - (gamma[j] * S1bar_inv[j] * bbar[j]^2) + alpha_term
    T1 <- T0 - (A_j + B_j)^2 / D_j +
      (gamma[j] * S1bar_inv[j] * bbar[j]^2)
    lin_omega <- E_logit_rho_z - 0.5 * safe_log(D_j) +
      0.5 * safe_log(S_inv_j) -
      (A_y + dat$n / 2) *
      (safe_log(B_y + 0.5 * pmax(T1, 1e-10)) -
         safe_log(B_y + 0.5 * pmax(T0, 1e-10)))
    omega_proposal <- clamp(
      expit(clamp(lin_omega, -30, 30)), .Machine$double.eps, 1 - .Machine$double.eps
    )
    omega_before <- omega
    delta_omega <- omega_proposal - omega_before[j]
    omega[j] <- omega_proposal
    inv_update <- .inverse_rank2_omega_update(
      Minv, j, delta_omega, omega_before, dat
    )
    Minv <- if (isTRUE(inv_update$ok)) inv_update$inverse else rebuild_inverse()

    ## gamma_j conditional update immediately after omega_j.  Both updates are
    ## visible to coordinate j+1, giving a cyclic 2p-coordinate schedule.
    w_m <- omega[idx]
    wj <- omega[j]
    S_inv_m <- S_inv_all[idx]
    v_m <- v_all[idx]
    r_minus <- w_m * dat$Xty[idx] + v_m
    ccol <- w_m * dat$XtX[idx, j]
    crow <- dat$XtX[j, idx] * w_m
    SOL <- cbind(
      .solve_leave_one_cached(Minv, j, r_minus),
      .solve_leave_one_cached(Minv, j, ccol)
    )
    sol_r <- SOL[, 1]
    sol_c <- SOL[, 2]
    cTc <- sum(ccol * sol_c)
    D0j <- pmax(S0bar_inv + wj * dat$x_norm2[j] - wj^2 * cTc, 1e-10)
    D1j <- pmax(S1bar_inv[j] + wj * dat$x_norm2[j] - wj^2 * cTc, 1e-10)
    R0 <- wj * (dat$Xty[j] - sum(crow * sol_r))
    v1j <- S1bar_inv[j] * bbar[j]
    R1 <- R0 + v1j
    alpha_minus <- sum(gamma[idx] * S1bar_inv[idx] * bbar[idx]^2)
    Tc0 <- pmax(
      dat$yTy - sum(r_minus * sol_r) - (R0^2) / D0j + alpha_minus,
      1e-10
    )
    Tc1 <- pmax(
      dat$yTy - sum(r_minus * sol_r) - (R1^2) / D1j + alpha_minus +
        S1bar_inv[j] * bbar[j]^2,
      1e-10
    )
    Delta <- 0.5 * safe_log(S1bar_inv[j] / S0bar_inv) -
      0.5 * safe_log(1 + (S1bar_inv[j] - S0bar_inv) / D0j) -
      (A_y + dat$n / 2) *
      (safe_log(B_y + 0.5 * Tc1) - safe_log(B_y + 0.5 * Tc0))
    gamma_proposal <- clamp(
      expit(E_logit_rho_c + Delta), .Machine$double.eps, 1 - .Machine$double.eps
    )
    S_old_j <- S_inv_all[j]
    gamma[j] <- gamma_proposal
    S_inv_all[j] <- S0bar_inv +
      (S1bar_inv[j] - S0bar_inv) * gamma_proposal
    v_all[j] <- S1bar_inv[j] * gamma_proposal * bbar[j]
    inv_update <- .inverse_rank1_coordinate_update(
      Minv, j, S_inv_all[j] - S_old_j
    )
    Minv <- if (isTRUE(inv_update$ok)) inv_update$inverse else rebuild_inverse()
  }

  list(
    omega = clamp(omega, .Machine$double.eps, 1 - .Machine$double.eps),
    gamma = clamp(gamma, .Machine$double.eps, 1 - .Machine$double.eps)
  )
}


## Slow paper-aligned fallback for omega_j.  This overrides the common-core
## implementation so the source moments match Gamma S1^{-1}.
update_omega_j_old <- function(j, omega, gamma, bbar, dat,
                               S0bar_inv, A_y, B_y, E_logit_rho_z) {
  idx <- setdiff(seq_along(omega), j)
  w <- omega[idx]
  ns_eff <- get_ns_eff(dat)
  S1bar_inv <- (ns_eff * dat$x_norm2) / dat$n
  S_inv_all <- S0bar_inv + (S1bar_inv - S0bar_inv) * gamma
  S_inv_j <- S_inv_all[j]
  S_inv_m <- S_inv_all[idx]
  v_all <- gamma * S1bar_inv * bbar
  v_j <- v_all[j]
  v_m <- v_all[idx]
  XtX_mm <- dat$XtX[idx, idx, drop = FALSE]
  M_minus <- XtX_mm * (w %o% w)
  diag(M_minus) <- diag(XtX_mm) * w
  M_minus <- M_minus + diag(S_inv_m, length(idx)) +
    diag(1e-10, length(idx))
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
  T0 <- dat$yTy - sum(u * sol_u) - 2 * sum(u * sol_v) -
    sum(v_m * sol_v) - (gamma[j] * S1bar_inv[j] * bbar[j]^2) +
    source_quad
  T1 <- T0 - (A_j + B_j)^2 / D_j +
    (gamma[j] * S1bar_inv[j] * bbar[j]^2)
  lin <- E_logit_rho_z - 0.5 * safe_log(D_j) +
    0.5 * safe_log(S_inv_j) -
    (A_y + dat$n / 2) *
    (safe_log(B_y + 0.5 * pmax(T1, 1e-10)) -
       safe_log(B_y + 0.5 * pmax(T0, 1e-10)))
  clamp(expit(lin), .Machine$double.eps, 1 - .Machine$double.eps)
}

## Slow fallback for gamma_j with the same full determinant correction as the
## cached update above.  This overrides update_gamma_j() from the common core.
update_gamma_j <- function(j, omega, gamma, bbar, dat,
                           S0bar_inv, A_y, B_y, E_logit_rho_c) {
  p <- length(omega)
  XtX <- dat$XtX
  Xty <- dat$Xty
  x_norm2 <- dat$x_norm2
  yTy <- dat$yTy
  ns_eff <- get_ns_eff(dat)

  S1bar_inv <- (ns_eff * x_norm2) / dat$n
  S_inv_all <- S0bar_inv + (S1bar_inv - S0bar_inv) * gamma

  idx <- setdiff(seq_len(p), j)
  wj <- omega[j]

  if (length(idx) == 0) {
    cTc <- 0
    base_quad <- 0
    alpha_minus <- 0
    R0 <- wj * Xty[j]
  } else {
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

  Tc0 <- yTy - base_quad - (R0^2) / D0j + alpha_minus
  Tc1 <- yTy - base_quad - (R1^2) / D1j + alpha_minus +
    (S1bar_inv[j] * bbar[j]^2)

  Tc0 <- pmax(Tc0, 1e-10)
  Tc1 <- pmax(Tc1, 1e-10)

  Delta <- 0.5 * safe_log(S1bar_inv[j] / S0bar_inv) -
    0.5 * safe_log(1 + (S1bar_inv[j] - S0bar_inv) / D0j) -
    (A_y + dat$n / 2) *
    (safe_log(B_y + 0.5 * Tc1) - safe_log(B_y + 0.5 * Tc0))

  expit(E_logit_rho_c + Delta)
}


## =========================================================
## Algorithm-aligned collapsed variational objective for Borrow-CVB
## =========================================================
## This is the deterministic variational objective corresponding to the
## moment substitutions used in the paper's computational algorithm:
##   E(Z X'X Z) = X'X o Omega,
##   E(S^{-1})   = Gamma S1^{-1} + (I-Gamma) S0^{-1},
##   E(S^{-1} C) = E(C S^{-1} C) = Gamma S1^{-1}.
## It is used as a monitor for the paper's approximate collapsed-CAVI map.

.bernoulli_entropy_sum <- function(q, eps = 1e-12) {
  q <- pmin(pmax(as.numeric(q), eps), 1 - eps)
  -sum(q * log(q) + (1 - q) * log(1 - q))
}

.normalize_simplex <- function(alpha, K = NULL) {
  alpha <- as.numeric(alpha)
  if (!is.null(K) && length(alpha) != K) stop("alpha has wrong length.")
  alpha[!is.finite(alpha)] <- 0
  alpha <- pmax(alpha, 0)
  if (sum(alpha) <= 0) alpha <- rep(1 / length(alpha), length(alpha))
  alpha / sum(alpha)
}

borrow_elbo_surrogate <- function(alpha,
                                  omega,
                                  gamma,
                                  dat,
                                  prior_par_borrow,
                                  sigma_beta_sq = 100,
                                  A_y = 1,
                                  B_y = 1,
                                  include_entropy = TRUE,
                                  include_indicator_prior = TRUE,
                                  include_source_normalization = TRUE,
                                  ridge = 1e-10) {
  p <- length(omega)
  if (length(gamma) != p) stop("omega and gamma must have the same length.")

  alpha <- .normalize_simplex(alpha, K = dat$K)
  omega <- clamp(as.numeric(omega), 0, 1)
  gamma <- clamp(as.numeric(gamma), 0, 1)

  S0bar_inv <- 1 / sigma_beta_sq
  ns_eff <- get_ns_eff(dat)
  S1bar_inv <- (ns_eff * dat$x_norm2) / dat$n
  S_inv <- S0bar_inv + (S1bar_inv - S0bar_inv) * gamma
  S_inv <- pmax(S_inv, 1e-12)

  bbar <- as.numeric(as.matrix(dat$B) %*% alpha)
  source_rhs <- gamma * S1bar_inv * bbar
  rhs <- omega * dat$Xty + source_rhs
  factor <- .trace_woodbury_factor(
    z = omega,
    diag_add = S_inv + dat$x_norm2 * omega * (1 - omega),
    dat = dat,
    ridge = ridge
  )
  sol <- .trace_woodbury_apply(factor, rhs)
  logdet <- factor$logdet
  if (!is.finite(logdet)) logdet <- NA_real_

  quad <- sum(rhs * sol)
  source_quad <- sum(gamma * S1bar_inv * bbar^2)
  Tval <- dat$yTy - quad + source_quad
  Tval <- max(Tval, 1e-10)

  collapsed <- -0.5 * logdet -
    (A_y + dat$n / 2) * safe_log(B_y + 0.5 * Tval)

  source_norm <- if (include_source_normalization) {
    0.5 * sum(
      gamma * safe_log(S1bar_inv) +
        (1 - gamma) * safe_log(S0bar_inv)
    )
  } else {
    0
  }

  rho_z <- clamp(prior_par_borrow$rho_z_mean, .Machine$double.eps, 1 - .Machine$double.eps)
  rho_c <- clamp(prior_par_borrow$rho_c_mean, .Machine$double.eps, 1 - .Machine$double.eps)

  indicator_prior <- if (include_indicator_prior) {
    sum(omega * log(rho_z) + (1 - omega) * log1p(-rho_z)) +
      sum(gamma * log(rho_c) + (1 - gamma) * log1p(-rho_c))
  } else {
    0
  }

  entropy <- if (include_entropy) {
    .bernoulli_entropy_sum(omega) + .bernoulli_entropy_sum(gamma)
  } else {
    0
  }

  total <- collapsed + source_norm + indicator_prior + entropy

  list(
    total = as.numeric(total),
    collapsed = as.numeric(collapsed),
    source_norm = as.numeric(source_norm),
    indicator_prior = as.numeric(indicator_prior),
    entropy = as.numeric(entropy),
    logdet = as.numeric(logdet),
    Tval = as.numeric(Tval),
    quad = as.numeric(quad),
    source_quad = as.numeric(source_quad),
    monitor_name = "paper_moment_collapsed_variational_objective",
    alpha = alpha
  )
}

borrow_algorithm_objective <- borrow_elbo_surrogate

## Paper equation (beta-main-estimator).  This local definition intentionally
## overrides the older common-core formula, which used tilde S^{-1} Gamma in
## the source right-hand side instead of Gamma S1^{-1}.
get_beta_hat <- function(alpha, omega, gamma, dat, S0bar_inv) {
  omega <- clamp(as.numeric(omega), 0, 1)
  gamma <- clamp(as.numeric(gamma), 0, 1)
  ns_eff <- get_ns_eff(dat)
  S1bar_inv <- (ns_eff * dat$x_norm2) / dat$n
  S_inv <- S0bar_inv + (S1bar_inv - S0bar_inv) * gamma
  bbar <- as.numeric(as.matrix(dat$B) %*% alpha)
  rhs <- omega * dat$Xty + gamma * S1bar_inv * bbar
  factor <- .trace_woodbury_factor(
    z = omega,
    diag_add = S_inv + dat$x_norm2 * omega * (1 - omega),
    dat = dat,
    ridge = 1e-10
  )
  omega * .trace_woodbury_apply(factor, rhs)
}

.make_elbo_row <- function(iter,
                           stage,
                           alpha_method,
                           init_method,
                           alpha,
                           omega,
                           gamma,
                           dat,
                           prior_par_borrow,
                           sigma_beta_sq,
                           A_y,
                           B_y,
                           d_omega = NA_real_,
                           d_gamma = NA_real_,
                           d_alpha = NA_real_,
                           elbo_prev = NA_real_,
                           include_source_normalization = TRUE) {
  e <- borrow_elbo_surrogate(
    alpha = alpha,
    omega = omega,
    gamma = gamma,
    dat = dat,
    prior_par_borrow = prior_par_borrow,
    sigma_beta_sq = sigma_beta_sq,
    A_y = A_y,
    B_y = B_y,
    include_source_normalization = include_source_normalization
  )

  delta <- if (is.finite(elbo_prev)) e$total - elbo_prev else NA_real_
  rel <- if (is.finite(delta)) abs(delta) / (1 + abs(elbo_prev)) else NA_real_

  param_change_values <- c(d_omega, d_gamma, d_alpha)
  param_change_values <- param_change_values[is.finite(param_change_values)]
  max_param_change_value <- if (length(param_change_values) > 0L) {
    max(param_change_values)
  } else {
    NA_real_
  }

  data.frame(
    iter = iter,
    stage = stage,
    objective_name = e$monitor_name,
    alpha_method = alpha_method,
    init_method = init_method,
    elbo = e$total,
    delta_elbo = delta,
    rel_delta_elbo = rel,
    collapsed = e$collapsed,
    source_norm = e$source_norm,
    indicator_prior = e$indicator_prior,
    entropy = e$entropy,
    logdet = e$logdet,
    Tval = e$Tval,
    d_omega = d_omega,
    d_gamma = d_gamma,
    d_alpha = d_alpha,
    max_param_change = max_param_change_value,
    size05 = sum(omega > 0.5),
    mean_omega_signal = mean(omega[dat$tag_true == 1]),
    mean_omega_null = mean(omega[dat$tag_true == 0]),
    mean_gamma_signal = mean(gamma[dat$tag_true == 1]),
    mean_gamma_null = mean(gamma[dat$tag_true == 0]),
    alpha_1 = alpha[1],
    alpha_2 = if (length(alpha) >= 2) alpha[2] else NA_real_,
    alpha_3 = if (length(alpha) >= 3) alpha[3] else NA_real_,
    stringsAsFactors = FALSE
  )
}

print_elbo_trace_row <- function(row) {
  cat(sprintf(
    "[%s/%s] iter %3d | ELBO=% .6e | dELBO=% .3e | rel=% .3e | dpar=% .3e | size05=%3d | alpha=(%.3f, %.3f, %.3f)\n",
    row$alpha_method[1], row$stage[1], row$iter[1], row$elbo[1],
    row$delta_elbo[1], row$rel_delta_elbo[1], row$max_param_change[1],
    row$size05[1], row$alpha_1[1], row$alpha_2[1], row$alpha_3[1]
  ))
}

profile_alpha_grid_borrow <- function(fit,
                                      dat,
                                      prior_par_borrow,
                                      sigma_beta_sq = 100,
                                      A_y = 1,
                                      B_y = 1,
                                      grid_step = 0.02,
                                      include_source_normalization = TRUE,
                                      top_n = 20) {
  if (dat$K != 3) stop("profile_alpha_grid_borrow is currently written for K = 3.")
  grid_step <- as.numeric(grid_step)
  vals <- seq(0, 1, by = grid_step)
  grid <- do.call(rbind, lapply(vals, function(a1) {
    a2max <- 1 - a1
    a2vals <- vals[vals <= a2max + 1e-12]
    cbind(alpha_1 = a1, alpha_2 = a2vals, alpha_3 = 1 - a1 - a2vals)
  }))
  grid <- as.data.frame(grid)

  omega <- fit$omega
  gamma <- fit$gamma
  out <- lapply(seq_len(nrow(grid)), function(i) {
    a <- as.numeric(grid[i, c("alpha_1", "alpha_2", "alpha_3")])
    e <- borrow_elbo_surrogate(
      alpha = a,
      omega = omega,
      gamma = gamma,
      dat = dat,
      prior_par_borrow = prior_par_borrow,
      sigma_beta_sq = sigma_beta_sq,
      A_y = A_y,
      B_y = B_y,
      include_source_normalization = include_source_normalization
    )
    c(alpha_1 = a[1], alpha_2 = a[2], alpha_3 = a[3], elbo = e$total, Tval = e$Tval)
  })
  prof <- as.data.frame(do.call(rbind, out))
  prof <- prof[order(prof$elbo, decreasing = TRUE), , drop = FALSE]
  row.names(prof) <- NULL

  current <- fit$alpha
  current_elbo <- borrow_elbo_surrogate(
    alpha = current,
    omega = omega,
    gamma = gamma,
    dat = dat,
    prior_par_borrow = prior_par_borrow,
    sigma_beta_sq = sigma_beta_sq,
    A_y = A_y,
    B_y = B_y,
    include_source_normalization = include_source_normalization
  )$total

  prof$is_current_alpha <- FALSE
  attr(prof, "current_alpha") <- current
  attr(prof, "current_elbo") <- current_elbo
  attr(prof, "best_elbo") <- prof$elbo[1]
  attr(prof, "best_alpha") <- as.numeric(prof[1, c("alpha_1", "alpha_2", "alpha_3")])
  head(prof, top_n)
}

## =========================================================
## Alpha update diagnostics and safeguards
## =========================================================
## alpha_update_strategy = "old": original update_alpha_old().
## alpha_update_strategy = "elbo_accept": use update_alpha_old() as a proposal,
##   then accept only if the chosen old-alpha criterion does not decrease.
##   By default this is a pure accept/reject step with no line search, so the
##   accepted alpha still corresponds to the posterior parameters returned by
##   update_alpha_old().
## alpha_update_strategy = "em_qp": update alpha by the EM-style M-step in the
##   alternative optimization argument.  Given omega and gamma, it minimizes
##       0.5 * alpha^T H_t alpha - g_t^T alpha
##   over the simplex by an active-set method.  This is not a simplex grid.

.alpha_elbo_value <- function(alpha,
                              omega,
                              gamma,
                              dat,
                              prior_par_borrow,
                              sigma_beta_sq,
                              A_y = 1,
                              B_y = 1,
                              include_source_normalization = TRUE) {
  borrow_elbo_surrogate(
    alpha = alpha,
    omega = omega,
    gamma = gamma,
    dat = dat,
    prior_par_borrow = prior_par_borrow,
    sigma_beta_sq = sigma_beta_sq,
    A_y = A_y,
    B_y = B_y,
    include_source_normalization = include_source_normalization
  )$total
}

.alpha_qp_objective <- function(alpha, H, g) {
  alpha <- as.numeric(alpha)
  as.numeric(0.5 * sum(alpha * as.numeric(H %*% alpha)) - sum(g * alpha))
}

.compute_alpha_em_qp_terms <- function(omega,
                                       gamma,
                                       dat,
                                       sigma_beta_sq = 100,
                                       ridge_M = 1e-10,
                                       ridge_H = 1e-10) {
  K <- dat$K
  p <- dat$p
  omega <- clamp(as.numeric(omega), 0, 1)
  gamma <- clamp(as.numeric(gamma), 0, 1)

  S0bar_inv <- 1 / sigma_beta_sq
  ns_eff <- get_ns_eff(dat)
  S1bar_inv <- (ns_eff * dat$x_norm2) / dat$n
  S_inv <- S0bar_inv + (S1bar_inv - S0bar_inv) * gamma
  S_inv <- pmax(S_inv, 1e-12)

  ## M_t = tilde{S}^{-1} + X^T X \odot Omega.
  M <- dat$XtX * (omega %o% omega)
  diag(M) <- dat$x_norm2 * omega + S_inv
  M <- M + diag(ridge_M, p)

  ## Paper equations (H-main) and (g-main):
  ##   H = B' [Gamma S1 - Gamma S1 M^{-1} S1 Gamma] B,
  ##   g = B' Gamma S1 M^{-1} W X'y.
  Bmat <- as.matrix(dat$B)
  G1B <- Bmat * as.numeric(gamma * S1bar_inv)
  target_rhs <- as.numeric(omega * dat$Xty)
  factor <- .trace_woodbury_factor(
    z = omega,
    diag_add = S_inv + dat$x_norm2 * omega * (1 - omega),
    dat = dat,
    ridge = ridge_M
  )
  solved <- .trace_woodbury_apply(factor, cbind(G1B, target_rhs))
  Minv_G1B <- solved[, seq_len(K), drop = FALSE]
  Minv_target <- solved[, K + 1L]

  H <- crossprod(Bmat, G1B) - crossprod(G1B, Minv_G1B)
  g <- as.numeric(crossprod(G1B, Minv_target))
  yPy <- as.numeric(dat$yTy - sum(target_rhs * Minv_target))
  if (!is.finite(yPy)) yPy <- NA_real_

  H <- (H + t(H)) / 2
  H[!is.finite(H)] <- 0
  g[!is.finite(g)] <- 0

  ## Numerically enforce positive semidefiniteness of the small K x K matrix.
  ev <- tryCatch(eigen(H, symmetric = TRUE, only.values = TRUE)$values,
                 error = function(e) NA_real_)
  min_ev <- suppressWarnings(min(ev, na.rm = TRUE))
  if (!is.finite(min_ev)) min_ev <- 0
  if (min_ev < ridge_H) {
    H <- H + diag(abs(min_ev) + ridge_H, K)
  } else {
    H <- H + diag(ridge_H, K)
  }

  list(H = H, g = g, yPy = yPy, M = M, S_inv = S_inv)
}

.solve_simplex_qp_active_set <- function(H,
                                         g,
                                         tol = 1e-10,
                                         ridge = 1e-10,
                                         max_iter = NULL) {
  H <- as.matrix(H)
  g <- as.numeric(g)
  K <- length(g)
  if (is.null(max_iter)) max_iter <- K + 5L
  H <- (H + t(H)) / 2
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
    } else {
      Haa <- H[active, active, drop = FALSE]
      ga <- g[active]
      one <- rep(1, m)

      sol_g <- tryCatch(
        solve(Haa, ga),
        error = function(e) solve(Haa + diag(ridge, m), ga)
      )
      sol_1 <- tryCatch(
        solve(Haa, one),
        error = function(e) solve(Haa + diag(ridge, m), one)
      )
      denom <- as.numeric(sum(sol_1))
      if (!is.finite(denom) || abs(denom) < 1e-12) denom <- sign(denom + 1e-12) * 1e-12
      lambda <- as.numeric((sum(sol_g) - 1) / denom)
      alpha_active <- as.numeric(sol_g - lambda * sol_1)
    }

    if (all(alpha_active > -tol) || length(active) == 1L) {
      break
    }

    ## Remove the index with the most negative component, as described in the
    ## active-set M-step.
    worst_local <- which.min(alpha_active)
    removed <- c(removed, active[worst_local])
    active <- active[-worst_local]
  }

  alpha <- rep(0, K)
  alpha[active] <- pmax(alpha_active, 0)
  alpha <- .normalize_simplex(alpha, K = K)

  grad <- as.numeric(H %*% alpha - g)
  mu <- grad + lambda

  list(
    alpha = alpha,
    active_set = active,
    removed_set = removed,
    lambda = lambda,
    mu = mu,
    objective = .alpha_qp_objective(alpha, H, g),
    kkt_min_mu_inactive = if (length(setdiff(seq_len(K), active)) > 0) {
      min(mu[setdiff(seq_len(K), active)])
    } else {
      NA_real_
    }
  )
}


## Exact global simplex-QP solver for the small number of sources used here.
## For K <= max_exact_K it enumerates all nonempty faces of the simplex and
## checks the full KKT conditions.  This guarantees that the preliminary alpha
## update and the EM M-step use the global minimizer of the convex QP, rather
## than merely a feasible active-set candidate.
.solve_simplex_qp_global <- function(H,
                                     g,
                                     tol = 1e-10,
                                     ridge = 1e-10,
                                     max_exact_K = 12L) {
  H <- as.matrix(H)
  g <- as.numeric(g)
  K <- length(g)
  H <- (H + t(H)) / 2

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
      } else {
        sol_g <- tryCatch(
          solve(Haa, ga),
          error = function(e) solve(Haa + diag(ridge, m), ga)
        )
        sol_1 <- tryCatch(
          solve(Haa, one),
          error = function(e) solve(Haa + diag(ridge, m), one)
        )
        denom <- as.numeric(sum(sol_1))
        if (!is.finite(denom) || abs(denom) < 1e-12) next
        lambda_face <- as.numeric((sum(sol_g) - 1) / denom)
        alpha_active <- as.numeric(sol_g - lambda_face * sol_1)
      }

      if (any(!is.finite(alpha_active)) || any(alpha_active < -tol)) next

      alpha <- rep(0, K)
      alpha[active] <- pmax(alpha_active, 0)
      if (sum(alpha) <= 0) next
      alpha <- alpha / sum(alpha)

      grad <- as.numeric(H %*% alpha - g)
      lambda <- -mean(grad[active])
      inactive <- setdiff(all_index, active)
      mu <- grad + lambda

      if (length(inactive) > 0 && any(mu[inactive] < -sqrt(tol))) next
      if (max(abs(grad[active] + lambda)) > 10 * sqrt(tol)) next

      objective <- .alpha_qp_objective(alpha, H, g)
      if (is.null(best) || objective < best$objective) {
        best <- list(
          alpha = alpha,
          active_set = active,
          removed_set = inactive,
          lambda = lambda,
          mu = mu,
          objective = objective,
          kkt_min_mu_inactive = if (length(inactive) > 0) min(mu[inactive]) else NA_real_,
          solver = "exact_face_enumeration"
        )
      }
    }
  }

  if (is.null(best)) {
    best <- .solve_simplex_qp_active_set(H, g, tol = tol, ridge = ridge)
    best$solver <- "active_set_emergency_fallback"
  }
  best
}

update_alpha_em_qp <- function(alpha_old,
                               omega,
                               gamma,
                               dat,
                               sigma_beta_sq = 100,
                               active_tol = 1e-10,
                               ridge_M = 1e-10,
                               ridge_H = 1e-10) {
  alpha_old <- .normalize_simplex(alpha_old, K = dat$K)
  terms <- .compute_alpha_em_qp_terms(
    omega = omega,
    gamma = gamma,
    dat = dat,
    sigma_beta_sq = sigma_beta_sq,
    ridge_M = ridge_M,
    ridge_H = ridge_H
  )

  sol <- .solve_simplex_qp_global(
    H = terms$H,
    g = terms$g,
    tol = active_tol,
    ridge = ridge_H
  )

  list(
    alpha = sol$alpha,
    accepted = TRUE,
    step_size = NA_real_,
    qp_old = .alpha_qp_objective(alpha_old, terms$H, terms$g),
    qp_new = .alpha_qp_objective(sol$alpha, terms$H, terms$g),
    H = terms$H,
    g = terms$g,
    active_set = sol$active_set,
    removed_set = sol$removed_set,
    lambda = sol$lambda,
    mu = sol$mu,
    kkt_min_mu_inactive = sol$kkt_min_mu_inactive,
    solver = sol$solver,
    strategy = "em_qp"
  )
}


## =========================================================
## Preliminary target center and legacy weighted least-squares helper
## The latest manuscript initializes alpha with the omega-weighted Huber loss.
## update_alpha_preliminary_huber() below uses only beta_target from this helper;
## its H/g weighted least-squares output is retained for optional legacy checks.
## The active regular-cycle alpha update remains update_alpha_em_qp().
## =========================================================
.compute_alpha_preliminary_qp_terms <- function(omega,
                                                dat,
                                                sigma_beta_sq = 100,
                                                ridge_M = 1e-10,
                                                ridge_H = 1e-10) {
  omega <- clamp(as.numeric(omega), 0, 1)
  p <- dat$p
  K <- dat$K
  if (length(omega) != p) stop("omega has wrong length in preliminary alpha QP.")

  S0bar_inv <- 1 / sigma_beta_sq
  ns_eff <- get_ns_eff(dat)
  S1bar_inv <- (ns_eff * dat$x_norm2) / dat$n
  S1bar_inv <- pmax(S1bar_inv, 1e-12)

  ## M0 = S0^{-1} + X'X o Omega(omega).
  M0 <- dat$XtX * (omega %o% omega)
  diag(M0) <- dat$x_norm2 * omega + S0bar_inv
  M0 <- M0 + diag(ridge_M, p)

  target_rhs <- omega * dat$Xty
  factor <- .trace_woodbury_factor(
    z = omega,
    diag_add = rep(S0bar_inv, p) + dat$x_norm2 * omega * (1 - omega),
    dat = dat,
    ridge = ridge_M
  )
  beta_target <- .trace_woodbury_apply(factor, target_rhs)

  Bmat <- as.matrix(dat$B)
  weight <- omega * S1bar_inv

  H <- crossprod(Bmat, Bmat * weight)
  g <- as.numeric(crossprod(Bmat, weight * beta_target))

  H <- (H + t(H)) / 2
  H[!is.finite(H)] <- 0
  g[!is.finite(g)] <- 0

  ev <- tryCatch(
    eigen(H, symmetric = TRUE, only.values = TRUE)$values,
    error = function(e) NA_real_
  )
  min_ev <- suppressWarnings(min(ev, na.rm = TRUE))
  if (!is.finite(min_ev)) min_ev <- 0
  if (min_ev < ridge_H) {
    H <- H + diag(abs(min_ev) + ridge_H, K)
  } else {
    H <- H + diag(ridge_H, K)
  }

  list(
    H = H,
    g = g,
    M0 = M0,
    beta_target = beta_target,
    S1bar_inv = S1bar_inv,
    weight = weight
  )
}

update_alpha_preliminary_qp <- function(alpha_old,
                                        omega,
                                        dat,
                                        sigma_beta_sq = 100,
                                        active_tol = 1e-10,
                                        ridge_M = 1e-10,
                                        ridge_H = 1e-10) {
  alpha_old <- .normalize_simplex(alpha_old, K = dat$K)
  terms <- .compute_alpha_preliminary_qp_terms(
    omega = omega,
    dat = dat,
    sigma_beta_sq = sigma_beta_sq,
    ridge_M = ridge_M,
    ridge_H = ridge_H
  )

  sol <- .solve_simplex_qp_global(
    H = terms$H,
    g = terms$g,
    tol = active_tol,
    ridge = ridge_H
  )

  list(
    alpha = sol$alpha,
    accepted = TRUE,
    qp_old = .alpha_qp_objective(alpha_old, terms$H, terms$g),
    qp_new = .alpha_qp_objective(sol$alpha, terms$H, terms$g),
    H = terms$H,
    g = terms$g,
    beta_target = terms$beta_target,
    active_set = sol$active_set,
    removed_set = sol$removed_set,
    lambda = sol$lambda,
    mu = sol$mu,
    kkt_min_mu_inactive = sol$kkt_min_mu_inactive,
    solver = sol$solver,
    strategy = "preliminary_target_matching_qp"
  )
}


## =========================================================
## TRACE robust preliminary alpha update (paper eqs. 11--13)
## =========================================================
source(file.path(.this_dir, "trace_huber_alpha.R"), local = .GlobalEnv)
source(file.path(.this_dir, "trace_target_initialization.R"), local = .GlobalEnv)

update_alpha_with_elbo_accept <- function(alpha_old,
                                          omega,
                                          gamma,
                                          dat,
                                          prior_par_borrow,
                                          sigma_beta_sq = 100,
                                          A_y = 1,
                                          B_y = 1,
                                          S0bar_inv = NULL,
                                          alpha_eta_old = NULL,
                                          include_source_normalization = TRUE,
                                          accept_tol = 1e-8,
                                          accept_objective = c("old_dirichlet_surrogate", "point_elbo"),
                                          line_search_steps = c(1)) {
  accept_objective <- match.arg(accept_objective)
  alpha_old <- .normalize_simplex(alpha_old, K = dat$K)

  eval_alpha_objective <- function(alpha_eval, eta_eval = NULL) {
    if (identical(accept_objective, "old_dirichlet_surrogate")) {
      tmp_fit <- list(
        alpha = .normalize_simplex(alpha_eval, K = dat$K),
        omega = omega,
        gamma = gamma,
        alpha_eta = eta_eval,
        sigma_beta_sq = sigma_beta_sq,
        alpha_update_strategy = "elbo_accept"
      )
      compute_old_alpha_dirichlet_criterion(
        fit = tmp_fit,
        dat = dat,
        prior_par_borrow = prior_par_borrow,
        sigma_beta_sq = sigma_beta_sq,
        A_y = A_y,
        B_y = B_y,
        include_source_normalization = include_source_normalization
      )$old_alpha_elbo
    } else {
      .alpha_elbo_value(
        alpha = alpha_eval,
        omega = omega,
        gamma = gamma,
        dat = dat,
        prior_par_borrow = prior_par_borrow,
        sigma_beta_sq = sigma_beta_sq,
        A_y = A_y,
        B_y = B_y,
        include_source_normalization = include_source_normalization
      )
    }
  }

  elbo_old <- eval_alpha_objective(alpha_old, alpha_eta_old)

  if (is.null(S0bar_inv)) S0bar_inv <- 1 / sigma_beta_sq
  raw_update <- update_alpha_old(
    alpha_old, omega, gamma, dat, A_y, B_y, S0bar_inv
  )
  alpha_cand <- .normalize_simplex(raw_update$alpha, K = dat$K)
  eta_cand <- .extract_alpha_eta(raw_update)

  best_alpha <- alpha_old
  best_elbo <- elbo_old
  best_step <- 0
  best_eta <- alpha_eta_old

  for (lambda in line_search_steps) {
    alpha_try <- (1 - lambda) * alpha_old + lambda * alpha_cand
    alpha_try <- .normalize_simplex(alpha_try, K = dat$K)

    ## For the old Dirichlet-alpha criterion, the posterior parameters returned
    ## by update_alpha_old() are only coherent for the full proposal lambda = 1.
    ## Therefore the default line_search_steps is c(1).  If a user provides a
    ## smaller lambda, we keep eta unchanged unless lambda is exactly 1.
    eta_try <- if (abs(lambda - 1) < 1e-12) eta_cand else alpha_eta_old
    elbo_try <- eval_alpha_objective(alpha_try, eta_try)

    if (is.finite(elbo_try) && elbo_try > best_elbo) {
      best_alpha <- alpha_try
      best_elbo <- elbo_try
      best_step <- lambda
      best_eta <- eta_try
    }
    if (is.finite(elbo_try) && elbo_try >= elbo_old - abs(accept_tol)) {
      return(list(
        alpha = alpha_try,
        eta = eta_try,
        raw_update = raw_update,
        accepted = TRUE,
        step_size = lambda,
        elbo_old = elbo_old,
        elbo_new = elbo_try,
        alpha_proposal = alpha_cand,
        eta_proposal = eta_cand,
        accept_objective = accept_objective,
        strategy = "elbo_accept"
      ))
    }
  }

  ## If no step passes the non-decrease check, keep alpha and eta unchanged.
  list(
    alpha = alpha_old,
    eta = alpha_eta_old,
    raw_update = raw_update,
    accepted = FALSE,
    step_size = 0,
    elbo_old = elbo_old,
    elbo_new = elbo_old,
    best_rejected_elbo = best_elbo,
    best_rejected_alpha = best_alpha,
    best_rejected_eta = best_eta,
    best_rejected_step = best_step,
    alpha_proposal = alpha_cand,
    eta_proposal = eta_cand,
    accept_objective = accept_objective,
    strategy = "elbo_accept"
  )
}


## Override the slow core fitter from _common/borrow_functions_core.R.
fit_borrow_model_fixedhyper <- function(dat,
                                        prior_par_borrow,
                                        alpha_method = c("old", "fixed", "oracle"),
                                        alpha_fixed = NULL,
                                        init_alpha = c("uniform", "source1", "random", "oracle"),
                                        init_method = c("full", "mild", "source", "null", "mcp", "lasso",
                                                        "scad", "mix", "custom"),
                                        init_omega = NULL,
                                        init_gamma = "theory_zero",
                                        sigma_beta_sq = 100,
                                        max_iter = 100,
                                        tol = 1e-6,
                                        tol_elbo = 1e-7,
                                        min_iter = 5,
                                        gamma_ebic = 0.5,
                                        verbose = FALSE,
                                        fast_cached_update = TRUE,
                                        borrow_update_schedule = c("paper_block", "gauss_seidel", "jacobi"),
                                        alpha_update_every = 1,
                                        alpha_update_strategy = c("em_qp", "elbo_accept", "old"),
                                        alpha_freeze_iter = 0,
                                        alpha_grid_step = 0.02,
                                        alpha_accept_tol = 1e-8,
                                        adaptive_damping = TRUE,
                                        update_damping = 0.5,
                                        cycle_detection_tol = 1e-6,
                                        cycle_detection_relative_tol = 1e-3,
                                        cycle_min_amplitude = NULL,
                                        theory_initialization = TRUE,
    target_init_until_converged = FALSE,
    target_init_steps = 5L,
    target_init_max_iter = 1000L,
    target_init_tol = NULL,
                                        preliminary_alpha_update = c("huber", "target_matching_qp", "none"),
                                        huber_tau = NULL,
                                        huber_tau_rule = c("fixed", "diverging"),
                                        huber_tau_multiplier = 1,
                                        monitor_elbo = TRUE,
                                        elbo_print_every = 1,
                                        include_source_normalization = TRUE,
                                        stop_by_elbo = FALSE,
                                        warn_elbo_drop = TRUE,
                                        elbo_drop_tol = 1e-6,
                                        stage_label = NULL) {
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
  if (length(update_damping) != 1L || !is.finite(update_damping) ||
      update_damping <= 0 || update_damping > 1) {
    stop("update_damping must be a finite scalar in (0, 1].")
  }
  cycle_detection_tol <- as.numeric(cycle_detection_tol)
  cycle_detection_relative_tol <- as.numeric(cycle_detection_relative_tol)
  if (length(cycle_detection_tol) != 1L || !is.finite(cycle_detection_tol) ||
      cycle_detection_tol <= 0) {
    stop("cycle_detection_tol must be a positive finite scalar.")
  }
  if (length(cycle_detection_relative_tol) != 1L ||
      !is.finite(cycle_detection_relative_tol) ||
      cycle_detection_relative_tol < 0) {
    stop("cycle_detection_relative_tol must be a nonnegative finite scalar.")
  }
  if (is.null(cycle_min_amplitude)) {
    cycle_min_amplitude <- max(100 * tol, 1e-4)
  } else {
    cycle_min_amplitude <- as.numeric(cycle_min_amplitude)
    if (length(cycle_min_amplitude) != 1L ||
        !is.finite(cycle_min_amplitude) || cycle_min_amplitude <= 0) {
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
  S0bar_inv <- 1 / S0bar

  p <- dat$p
  if (is.null(stage_label)) stage_label <- alpha_method

  omega <- make_init_omega(
    dat = dat,
    init_method = init_method,
    init_omega = init_omega,
    gamma_ebic = gamma_ebic
  )

  alpha <- make_init_alpha(
    dat = dat,
    init_alpha = init_alpha,
    alpha_fixed = alpha_fixed,
    alpha_method = alpha_method
  )

  ## The initialization requires gamma^(0)=0.
  ## this overrides any user-supplied agreement/clean initialization.
  gamma_requested <- init_gamma
  if (isTRUE(theory_initialization)) {
    gamma <- rep(0, p)
    init_gamma_effective <- "theory_zero"
  } else {
    gamma <- make_init_gamma_borrow(
      dat = dat,
      init_gamma = init_gamma,
      alpha = alpha,
      gamma_ebic = gamma_ebic
    )
    init_gamma_effective <- init_gamma
  }

  ## The theory treats rho_z and rho_c as deterministic hyperparameters.
  ## Consequently the coordinate odds contain logit(rho), not the beta-prior
  ## expectation digamma(A)-digamma(B).
  rho_z_algorithm <- prior_par_borrow$rho_z_mean %||% (A_z / (A_z + B_z))
  rho_c_algorithm <- prior_par_borrow$rho_c_mean %||% (A_c / (A_c + B_c))
  E_logit_rho_z <- qlogis(clamp(rho_z_algorithm, .Machine$double.eps, 1 - .Machine$double.eps))
  E_logit_rho_c <- qlogis(clamp(rho_c_algorithm, .Machine$double.eps, 1 - .Machine$double.eps))

  elbo_trace <- list()
  elbo_prev <- NA_real_
  theory_init_trace <- list()

  add_theory_trace <- function(step, iter, omega_now, gamma_now, alpha_now) {
    theory_init_trace[[length(theory_init_trace) + 1L]] <<- data.frame(
      iter = iter,
      step = step,
      size05 = sum(omega_now > 0.5),
      mean_omega_signal = mean(omega_now[dat$tag_true == 1]),
      mean_omega_null = mean(omega_now[dat$tag_true == 0]),
      gamma_min = min(gamma_now),
      gamma_max = max(gamma_now),
      mean_gamma = mean(gamma_now),
      alpha_1 = alpha_now[1],
      alpha_2 = if (length(alpha_now) >= 2) alpha_now[2] else NA_real_,
      alpha_3 = if (length(alpha_now) >= 3) alpha_now[3] else NA_real_,
      stringsAsFactors = FALSE
    )
  }

  add_theory_trace("initial", 0L, omega, gamma, alpha)

  if (isTRUE(monitor_elbo)) {
    row0 <- .make_elbo_row(
      iter = 0,
      stage = stage_label,
      alpha_method = alpha_method,
      init_method = init_method,
      alpha = alpha,
      omega = omega,
      gamma = gamma,
      dat = dat,
      prior_par_borrow = prior_par_borrow,
      sigma_beta_sq = sigma_beta_sq,
      A_y = A_y,
      B_y = B_y,
      include_source_normalization = include_source_normalization
    )
    row0$raw_max_param_change <- NA_real_
    row0$cycle_distance_raw <- NA_real_
    row0$adaptive_damping_active <- !adaptive_damping && update_damping < 1
    row0$damping_applied <- FALSE
    row0$damping_factor <- 1
    elbo_trace[[length(elbo_trace) + 1]] <- row0
    elbo_prev <- row0$elbo
    if (verbose) print_elbo_trace_row(row0)
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
      } else {
        .update_omega_all_cached
      }
      return(update_fun(
        omega = omega_in,
        gamma = gamma_in,
        bbar = bbar_in,
        dat = dat,
        S0bar_inv = S0bar_inv,
        A_y = A_y,
        B_y = B_y,
        E_logit_rho_z = E_logit_rho_z
      ))
    }

    omega_prop <- omega_in
    for (j in seq_len(p)) {
      omega_prop[j] <- update_omega_j_old(
        j = j,
        omega = omega_in,
        gamma = gamma_in,
        bbar = bbar_in,
        dat = dat,
        S0bar_inv = S0bar_inv,
        A_y = A_y,
        B_y = B_y,
        E_logit_rho_z = E_logit_rho_z
      )
    }
    clamp(omega_prop, .Machine$double.eps, 1 - .Machine$double.eps)
  }

  update_gamma_once <- function(omega_in, gamma_in, alpha_in) {
    bbar_in <- as.numeric(as.matrix(dat$B) %*% alpha_in)
    if (fast_cached_update) {
      update_fun <- if (borrow_update_schedule == "gauss_seidel") {
        .update_gamma_all_gauss_seidel
      } else {
        .update_gamma_all_cached
      }
      return(update_fun(
        omega = omega_in,
        gamma = gamma_in,
        bbar = bbar_in,
        dat = dat,
        S0bar_inv = S0bar_inv,
        A_y = A_y,
        B_y = B_y,
        E_logit_rho_c = E_logit_rho_c
      ))
    }

    gamma_prop <- gamma_in
    for (j in seq_len(p)) {
      gamma_prop[j] <- update_gamma_j(
        j = j,
        omega = omega_in,
        gamma = gamma_in,
        bbar = bbar_in,
        dat = dat,
        S0bar_inv = S0bar_inv,
        A_y = A_y,
        B_y = B_y,
        E_logit_rho_c = E_logit_rho_c
      )
    }
    clamp(gamma_prop, .Machine$double.eps, 1 - .Machine$double.eps)
  }

  do_regular_alpha_update <- function(t, alpha_in, omega_in, gamma_in) {
    out <- list(alpha = alpha_in, updated = FALSE)
    if (alpha_method != "old") return(out)

    if (t <= alpha_freeze_iter) return(out)
    should_update_alpha <-
      (t == alpha_freeze_iter + 1L) ||
      ((t - alpha_freeze_iter) %% alpha_update_every == 0L) ||
      (t == max_iter)
    if (!should_update_alpha) return(out)

    if (alpha_update_strategy == "old") {
      au <- update_alpha_old(alpha_in, omega_in, gamma_in, dat, A_y, B_y, S0bar_inv)
      out$alpha <- .normalize_simplex(au$alpha, K = dat$K)
      out$raw <- au
      out$eta <- .extract_alpha_eta(au)
      out$accepted <- TRUE
      out$step_size <- NA_real_

    } else if (alpha_update_strategy == "elbo_accept") {
      au <- update_alpha_with_elbo_accept(
        alpha_old = alpha_in,
        omega = omega_in,
        gamma = gamma_in,
        dat = dat,
        prior_par_borrow = prior_par_borrow,
        sigma_beta_sq = sigma_beta_sq,
        A_y = A_y,
        B_y = B_y,
        S0bar_inv = S0bar_inv,
        alpha_eta_old = alpha_eta,
        include_source_normalization = include_source_normalization,
        accept_tol = alpha_accept_tol,
        accept_objective = "point_elbo"
      )
      out$alpha <- au$alpha
      out$eta <- au$eta
      out$raw <- au$raw_update
      out$accepted <- isTRUE(au$accepted)
      out$step_size <- au$step_size
      out$elbo_old <- au$elbo_old
      out$elbo_new <- au$elbo_new

    } else {
      au <- update_alpha_em_qp(
        alpha_old = alpha_in,
        omega = omega_in,
        gamma = gamma_in,
        dat = dat,
        sigma_beta_sq = sigma_beta_sq
      )
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
      ## The first theory cycle keeps alpha between omega and gamma by design.
      if (isTRUE(target_init_until_converged)) {
        target_initialization <- trace_target_only_initialize(
          dat=dat,omega=omega,S0bar_inv=S0bar_inv,rho_z=rho_z_algorithm,
          A_y=A_y,B_y=B_y,max_iter=target_init_max_iter,
          tol=if(is.null(target_init_tol)) tol else target_init_tol,verbose=verbose)
        if (!isTRUE(target_initialization$converged)) {
          stop(sprintf("Target-only omega did not converge after %d iterations (change %.4e); alpha/gamma updates were not started.",
                       target_initialization$n_iter,target_initialization$max_change))
        }
        omega <- target_initialization$omega
        add_theory_trace("omega_target_only_converged", 1L, omega, gamma, alpha)
      } else {
        target_initialization <- trace_target_only_fixed_steps(
          dat=dat,omega=omega,S0bar_inv=S0bar_inv,rho_z=rho_z_algorithm,
          A_y=A_y,B_y=B_y,n_steps=target_init_steps,
          tol=if(is.null(target_init_tol)) tol else target_init_tol,verbose=verbose)
        omega <- target_initialization$omega
        add_theory_trace("omega_target_only_fixed_steps", 1L, omega, gamma, alpha)
      }

      ## Step 2 of the special first cycle: optimize alpha before gamma.
      if (alpha_method == "old" && preliminary_alpha_update != "none") {
        au0 <- if (preliminary_alpha_update == "huber") {
          update_alpha_preliminary_huber(
            alpha_old = alpha,
            omega = omega,
            dat = dat,
            sigma_beta_sq = sigma_beta_sq,
            huber_tau = huber_tau,
            huber_tau_rule = huber_tau_rule,
            huber_tau_multiplier = huber_tau_multiplier
          )
        } else {
          update_alpha_preliminary_qp(
            alpha_old = alpha,
            omega = omega,
            dat = dat,
            sigma_beta_sq = sigma_beta_sq
          )
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
      } else if (alpha_method == "fixed") {
        alpha <- .normalize_simplex(alpha_fixed, K = dat$K)
      } else if (alpha_method == "oracle") {
        alpha <- .normalize_simplex(dat$alpha_true, K = dat$K)
      }
      add_theory_trace("alpha_preliminary", 1L, omega, gamma, alpha)

      ## Step 3 of the special first cycle: the first formal gamma update uses
      ## omega^(1) and the newly optimized alpha^(1).
      gamma <- update_gamma_once(omega, gamma, alpha)
      add_theory_trace("gamma_first_formal", 1L, omega, gamma, alpha)

    } else {
      if (borrow_update_schedule == "gauss_seidel") {
        if (!isTRUE(fast_cached_update)) {
          stop(
            "borrow_update_schedule = 'gauss_seidel' requires ",
            "fast_cached_update = TRUE."
          )
        }
        joint_update <- .update_omega_gamma_interleaved_gauss_seidel(
          omega = omega,
          gamma = gamma,
          bbar = as.numeric(as.matrix(dat$B) %*% alpha),
          dat = dat,
          S0bar_inv = S0bar_inv,
          A_y = A_y,
          B_y = B_y,
          E_logit_rho_z = E_logit_rho_z,
          E_logit_rho_c = E_logit_rho_c
        )
        omega <- joint_update$omega
        gamma <- joint_update$gamma
      } else {
        omega <- update_omega_once(omega, gamma, alpha)
        gamma <- update_gamma_once(omega, gamma, alpha)
      }

      if (alpha_method == "fixed") {
        alpha <- .normalize_simplex(alpha_fixed, K = dat$K)
      } else if (alpha_method == "oracle") {
        alpha <- .normalize_simplex(dat$alpha_true, K = dat$K)
      } else {
        au <- do_regular_alpha_update(t, alpha, omega, gamma)
        if (isTRUE(au$updated)) {
          alpha <- au$alpha
          last_alpha_update_raw <- au$raw %||% last_alpha_update_raw
          if (!is.null(au$eta)) alpha_eta <- au$eta
          if (isTRUE(au$accepted)) {
            alpha_update_accept_count <- alpha_update_accept_count + 1L
          } else {
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

    ## Detect an actual period-two return in the undamped full TRACE map.
    ## The first theory cycle is always left unchanged.  Once a cycle is
    ## detected, under-relax the complete (omega, gamma, alpha) sweep and keep
    ## damping active for the rest of this fit.  Convex damping preserves the
    ## simplex constraint and the fixed points of the original map.
    omega_raw <- omega
    gamma_raw <- gamma
    alpha_raw <- alpha
    raw_max_param_change <- max(
      max(abs(omega_raw - omega_old)),
      max(abs(gamma_raw - gamma_old)),
      max(abs(alpha_raw - alpha_old))
    )
    cycle_distance_raw <- if (!is.null(lag2_state)) {
      max(
        max(abs(omega_raw - lag2_state$omega)),
        max(abs(gamma_raw - lag2_state$gamma)),
        max(abs(alpha_raw - lag2_state$alpha))
      )
    } else {
      NA_real_
    }

    cycle_threshold_t <- max(
      cycle_detection_tol,
      cycle_detection_relative_tol * raw_max_param_change
    )
    cycle_match <- adaptive_damping && !damping_active &&
      !is_theory_first_cycle && t >= 3L &&
      is.finite(cycle_distance_raw) &&
      raw_max_param_change > cycle_min_amplitude &&
      cycle_distance_raw <= cycle_threshold_t

    if (cycle_match) {
      damping_active <- TRUE
      damping_triggered <- TRUE
      damping_trigger_iter <- t
      cycle_distance_at_trigger <- cycle_distance_raw
      if (isTRUE(verbose)) {
        cat(sprintf(
          "[TRACE adaptive damping] %s iter %d: two-cycle distance=%.3e, raw change=%.3e; eta=%.3f activated.\n",
          stage_label, t, cycle_distance_raw, raw_max_param_change, update_damping
        ))
      }
    }

    damping_applied_t <- damping_active && !is_theory_first_cycle &&
      update_damping < 1
    if (damping_applied_t) {
      omega <- clamp(
        (1 - update_damping) * omega_old + update_damping * omega_raw,
        .Machine$double.eps, 1 - .Machine$double.eps
      )
      gamma <- clamp(
        (1 - update_damping) * gamma_old + update_damping * gamma_raw,
        .Machine$double.eps, 1 - .Machine$double.eps
      )
      alpha <- .normalize_simplex(
        (1 - update_damping) * alpha_old + update_damping * alpha_raw,
        K = dat$K
      )
      damping_iterations <- damping_iterations + 1L
    }

    ## At the next sweep this is theta^(t-2).
    lag2_state <- list(omega = omega_old, gamma = gamma_old, alpha = alpha_old)

    d_omega <- max(abs(omega - omega_old))
    d_gamma <- max(abs(gamma - gamma_old))
    d_alpha <- max(abs(alpha - alpha_old))
    max_param_change <- max(d_omega, d_gamma, d_alpha)

    rel_delta_elbo <- NA_real_
    delta_elbo <- NA_real_

    if (isTRUE(monitor_elbo)) {
      rowt <- .make_elbo_row(
        iter = t,
        stage = if (is_theory_first_cycle) paste0(stage_label, "_theory_first_cycle") else stage_label,
        alpha_method = alpha_method,
        init_method = init_method,
        alpha = alpha,
        omega = omega,
        gamma = gamma,
        dat = dat,
        prior_par_borrow = prior_par_borrow,
        sigma_beta_sq = sigma_beta_sq,
        A_y = A_y,
        B_y = B_y,
        d_omega = d_omega,
        d_gamma = d_gamma,
        d_alpha = d_alpha,
        elbo_prev = elbo_prev,
        include_source_normalization = include_source_normalization
      )
      rowt$raw_max_param_change <- raw_max_param_change
      rowt$cycle_distance_raw <- cycle_distance_raw
      rowt$adaptive_damping_active <- damping_active
      rowt$damping_applied <- damping_applied_t
      rowt$damping_factor <- if (damping_applied_t) update_damping else 1
      elbo_trace[[length(elbo_trace) + 1]] <- rowt
      delta_elbo <- rowt$delta_elbo
      rel_delta_elbo <- rowt$rel_delta_elbo

      if (is.finite(delta_elbo) && delta_elbo < -abs(elbo_drop_tol)) {
        elbo_drop_count <- elbo_drop_count + 1L
        if (isTRUE(warn_elbo_drop)) {
          cat(sprintf(
            "[ELBO drop warning] %s iter %d: dELBO = %.6e. The first theory cycle is not a standard coordinate-ascent sweep, and the monitored surrogate need not be monotone.\n",
            stage_label, t, delta_elbo
          ))
        }
      }

      if (verbose && (t %% elbo_print_every == 0 || t == 1)) {
        print_elbo_trace_row(rowt)
      }
      elbo_prev <- rowt$elbo
    } else if (verbose) {
      cat(sprintf("[%s/%s fast=%s] iter %3d | d_omega=%.3e d_gamma=%.3e d_alpha=%.3e | size05=%d\n",
                  alpha_method, init_method, fast_cached_update, t,
                  d_omega, d_gamma, d_alpha, sum(omega > 0.5)))
    }

    param_ok <- max_param_change < tol && raw_max_param_change < tol
    elbo_ok <- if (isTRUE(monitor_elbo) && isTRUE(stop_by_elbo)) {
      is.finite(rel_delta_elbo) && rel_delta_elbo < tol_elbo
    } else {
      TRUE
    }

    earliest_stop <- if (isTRUE(theory_initialization)) 3L else 1L
    if (t >= max(min_iter, alpha_freeze_iter + 2L, earliest_stop) && param_ok && elbo_ok) {
      n_iter <- t
      converged <- TRUE
      convergence_reason <- if (isTRUE(stop_by_elbo) && isTRUE(monitor_elbo)) {
        "param_and_elbo_tol"
      } else {
        "param_tol"
      }
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

  elbo_trace_df <- if (length(elbo_trace) > 0) do.call(rbind, elbo_trace) else data.frame()
  theory_init_trace_df <- if (length(theory_init_trace) > 0) do.call(rbind, theory_init_trace) else data.frame()
  final_elbo <- if (nrow(elbo_trace_df) > 0) tail(elbo_trace_df$elbo, 1) else NA_real_

  list(
    method = paste0(
      method_name, "_init_", init_method,
      if (isTRUE(fast_cached_update) && borrow_update_schedule != "gauss_seidel") {
        "_woodbury_elbo_theoryinit"
      } else {
        "_cached_elbo_theoryinit"
      }
    ),
    alpha = alpha,
    alpha_eta = alpha_eta,
    last_alpha_update_raw = last_alpha_update_raw,
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
    init_gamma = init_gamma_effective,
    init_gamma_requested = gamma_requested,
    theory_initialization = theory_initialization,
      target_init_until_converged = target_init_until_converged,
      target_init_steps = target_init_steps,
      target_init_max_iter = target_init_max_iter,
      target_init_tol = target_init_tol,
    preliminary_alpha_update = preliminary_alpha_update,
    preliminary_alpha_done = preliminary_alpha_done,
    preliminary_alpha_qp_old = preliminary_alpha_qp_old,
    preliminary_alpha_qp_new = preliminary_alpha_qp_new,
    preliminary_alpha_active_set = preliminary_alpha_active_set,
    preliminary_alpha_solver = preliminary_alpha_solver,
    preliminary_alpha_beta_target = preliminary_alpha_beta_target,
    preliminary_alpha_huber_tau = preliminary_alpha_huber_tau,
    preliminary_alpha_huber_tau_rule = preliminary_alpha_huber_tau_rule,
    preliminary_alpha_huber_s_reference = preliminary_alpha_huber_s_reference,
    preliminary_alpha_huber_iterations = preliminary_alpha_huber_iterations,
    preliminary_alpha_huber_converged = preliminary_alpha_huber_converged,
    preliminary_alpha_huber_kkt_residual = preliminary_alpha_huber_kkt_residual,
    theory_init_trace = theory_init_trace_df,
    target_initialization = target_initialization,
    target_init_converged = if(is.null(target_initialization)) NA else target_initialization$converged,
    target_init_n_iter = if(is.null(target_initialization)) 1L else target_initialization$n_iter,
    target_init_max_change = if(is.null(target_initialization)) NA_real_ else target_initialization$max_change,
    target_init_size = if(is.null(target_initialization)) NA_integer_ else sum(target_initialization$omega>.5),
    sigma_beta_sq = sigma_beta_sq,
    n_iter = n_iter,
    converged = converged,
    convergence_reason = convergence_reason,
    d_omega_final = d_omega,
    d_gamma_final = d_gamma,
    d_alpha_final = d_alpha,
    max_param_change_final = max_param_change,
    stop_by_elbo = stop_by_elbo,
    final_elbo = final_elbo,
    elbo_monitor_name = "paper_moment_collapsed_variational_objective",
    elbo_trace = elbo_trace_df,
    elbo_drop_count = elbo_drop_count,
    fast_cached_update = fast_cached_update,
    trace_linear_algebra_backend = if (isTRUE(fast_cached_update) &&
                                         borrow_update_schedule != "gauss_seidel") {
      "exact_woodbury_n_by_n"
    } else if (isTRUE(fast_cached_update)) {
      "cached_rank_updates_p_by_p"
    } else {
      "dense_leave_one"
    },
    borrow_update_schedule = borrow_update_schedule,
    alpha_update_every = alpha_update_every,
    alpha_update_strategy = alpha_update_strategy,
    alpha_freeze_iter = alpha_freeze_iter,
    alpha_grid_step = alpha_grid_step,
    adaptive_damping = adaptive_damping,
    update_damping = update_damping,
    damping_triggered = damping_triggered,
    damping_trigger_iter = damping_trigger_iter,
    damping_iterations = damping_iterations,
    cycle_distance_at_trigger = cycle_distance_at_trigger,
    cycle_detection_tol = cycle_detection_tol,
    cycle_detection_relative_tol = cycle_detection_relative_tol,
    cycle_min_amplitude = cycle_min_amplitude,
    raw_max_param_change_final = raw_max_param_change,
    alpha_update_accept_count = alpha_update_accept_count,
    alpha_update_reject_count = alpha_update_reject_count,
    last_alpha_step_size = last_alpha_step_size,
    last_alpha_update_elbo_old = last_alpha_update_elbo_old,
    last_alpha_update_elbo_new = last_alpha_update_elbo_new,
    last_alpha_update_qp_old = last_alpha_update_qp_old,
    last_alpha_update_qp_new = last_alpha_update_qp_new,
    last_alpha_active_set_size = last_alpha_active_set_size,
    last_alpha_active_set = last_alpha_active_set,
    tol = tol,
    tol_elbo = tol_elbo,
    include_source_normalization = include_source_normalization
  )
}

.standardize_score <- function(x) {
  x <- as.numeric(x)
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s <= 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

make_borrow_candidate_set <- function(dat,
                                      candidate_size = 180,
                                      gamma_ebic = 0.5,
                                      min_size = NULL) {
  p <- dat$p
  candidate_size <- min(p, max(1, as.integer(candidate_size)))
  if (is.null(min_size)) min_size <- min(candidate_size, 4 * max(10, ceiling(sqrt(p))))
  candidate_size <- max(candidate_size, min_size)
  candidate_size <- min(candidate_size, p)
  if (candidate_size >= p) return(seq_len(p))

  target_score <- abs(dat$Xty) / sqrt(pmax(dat$x_norm2, 1e-12))
  source_avg_score <- abs(as.numeric(as.matrix(dat$B) %*% rep(1 / dat$K, dat$K)))
  source_max_score <- apply(abs(as.matrix(dat$B)), 1, max)

  lasso_sel <- integer(0)
  lasso_beta <- tryCatch(fit_ebic_glmnet_lasso(dat$X, dat$y, gamma_ebic = gamma_ebic)$beta,
                         error = function(e) rep(0, p))
  lasso_sel <- which(abs(lasso_beta) > 1e-8)

  top_target <- order(target_score, decreasing = TRUE)[seq_len(min(p, ceiling(0.55 * candidate_size)))]
  top_source_avg <- order(source_avg_score, decreasing = TRUE)[seq_len(min(p, ceiling(0.35 * candidate_size)))]
  top_source_max <- order(source_max_score, decreasing = TRUE)[seq_len(min(p, ceiling(0.35 * candidate_size)))]

  cand <- sort(unique(c(top_target, top_source_avg, top_source_max, lasso_sel)))

  if (length(cand) > candidate_size) {
    combined <- .standardize_score(target_score) +
      0.70 * .standardize_score(source_avg_score) +
      0.70 * .standardize_score(source_max_score) +
      2.00 * as.numeric(seq_len(p) %in% lasso_sel)
    ord_c <- cand[order(combined[cand], decreasing = TRUE)]
    cand <- sort(ord_c[seq_len(candidate_size)])
  }

  cand
}

subset_dat_for_borrow <- function(dat, cand) {
  cand <- sort(unique(as.integer(cand)))
  out <- dat
  out$p <- length(cand)
  out$X <- dat$X[, cand, drop = FALSE]
  out$X_test <- dat$X_test[, cand, drop = FALSE]
  out$beta_true <- dat$beta_true[cand]
  out$tag_true <- dat$tag_true[cand]
  out$gamma_true <- dat$gamma_true[cand]
  out$B <- as.matrix(dat$B)[cand, , drop = FALSE]
  out$W_source_true <- as.matrix(dat$W_source_true)[cand, , drop = FALSE]
  out$X_source <- lapply(dat$X_source, function(Xk) Xk[, cand, drop = FALSE])
  out$XtX <- dat$XtX[cand, cand, drop = FALSE]
  out$Xty <- dat$Xty[cand]
  out$x_norm2 <- dat$x_norm2[cand]
  out$SigX <- dat$SigX[cand, cand, drop = FALSE]
  out$common_bad <- which(cand %in% dat$common_bad)
  out$regular_signal <- which(cand %in% dat$regular_signal)
  out$spurious_union <- which(cand %in% dat$spurious_union)
  out$spur_list <- lapply(dat$spur_list, function(s) which(cand %in% s))
  out$borrow_candidate_original_index <- cand
  out$scenario_description <- paste0(dat$scenario_description %||% dat$scenario_title,
                                     " Candidate-screened Borrow-CVB with ", length(cand), " variables.")
  out
}

expand_borrow_fit <- function(fit, dat_full, cand) {
  p <- dat_full$p
  omega_full <- rep(1e-6, p)
  gamma_full <- rep(1e-6, p)
  beta_full <- rep(0, p)

  omega_full[cand] <- fit$omega
  gamma_full[cand] <- fit$gamma
  beta_full[cand] <- fit$beta_hat

  fit$omega <- omega_full
  fit$gamma <- gamma_full
  fit$beta_hat <- beta_full
  fit$pred_mse <- prediction_mse(beta_full, dat_full)
  fit$sel_05 <- get_selection_metrics_omega(omega_full, dat_full$beta_true, thr = 0.5)
  fit$sel_01 <- get_selection_metrics_omega(omega_full, dat_full$beta_true, thr = 0.1)
  fit$sel_top <- get_selection_metrics_topk(omega_full, dat_full$beta_true)
  fit$mean_omega_signal <- mean(omega_full[dat_full$tag_true == 1])
  fit$mean_omega_null <- mean(omega_full[dat_full$tag_true == 0])
  fit$mean_gamma_signal <- mean(gamma_full[dat_full$tag_true == 1])
  fit$mean_gamma_null <- mean(gamma_full[dat_full$tag_true == 0])
  fit$p_screened <- length(cand)
  fit$borrow_candidate_original_index <- cand
  fit$method <- paste0(fit$method, "_screen", length(cand))
  fit
}

fit_borrow_with_screening <- function(dat,
                                      ...,
                                      use_candidate_screening = FALSE,
                                      candidate_size = 180,
                                      fast_cached_update = TRUE) {
  args <- list(...)
  args$fast_cached_update <- fast_cached_update
  if (isTRUE(args$theory_initialization %||% TRUE) && isTRUE(use_candidate_screening)) {
    stop("Theory target-only initialization requires use_candidate_screening=FALSE.")
  }

  if (!isTRUE(use_candidate_screening) || is.null(candidate_size) || candidate_size >= dat$p) {
    args$dat <- dat
    return(do.call(fit_borrow_model_fixedhyper, args))
  }

  gamma_ebic <- args$gamma_ebic %||% 0.5
  cand <- make_borrow_candidate_set(dat, candidate_size = candidate_size, gamma_ebic = gamma_ebic)
  dat_small <- subset_dat_for_borrow(dat, cand)

  ## If custom initial omega/gamma are supplied on the full scale, subset them.
  if (!is.null(args$init_omega) && length(args$init_omega) == dat$p) args$init_omega <- args$init_omega[cand]
  if (!is.null(args$init_gamma) && is.numeric(args$init_gamma) && length(args$init_gamma) == dat$p) args$init_gamma <- args$init_gamma[cand]

  args$dat <- dat_small
  fit_small <- do.call(fit_borrow_model_fixedhyper, args)
  expand_borrow_fit(fit_small, dat_full = dat, cand = cand)
}


## =========================================================
## Fast no-borrow ACVB override
## =========================================================
## The default backend uses an exact Woodbury representation of the same ACVB
## Jacobi update.  It replaces a p x p inverse and p repeated submatrix slices
## by one n x n factorization and vectorized Sherman-Morrison downdates.  The
## legacy cached-inverse and original slow backends remain available for
## numerical validation.

.build_acvb_G_cached <- function(omega, dat, sigma_beta_sq) {
  lambda <- 1 / sigma_beta_sq
  G <- dat$XtX * (omega %o% omega)
  diag(G) <- dat$x_norm2 * omega + lambda
  G + diag(1e-10, length(omega))
}

.acvb_solve_woodbury <- function(omega,
                                 rhs,
                                 dat,
                                 sigma_beta_sq,
                                 ridge_n = 1e-10) {
  omega <- clamp(as.numeric(omega), 0, 1)
  rhs <- as.numeric(rhs)
  lambda <- 1 / sigma_beta_sq

  ## G = D + U'U, where U = X diag(omega) and
  ## D_jj = lambda + ||X_j||^2 omega_j (1 - omega_j).
  d <- lambda + dat$x_norm2 * omega * (1 - omega)
  d <- pmax(d, 1e-12)
  d_inv <- 1 / d
  U <- sweep(dat$X, 2, omega, `*`)
  U_scaled <- sweep(U, 2, sqrt(d_inv), `*`)
  K <- diag(dat$n) + tcrossprod(U_scaled)
  K <- K + diag(ridge_n, dat$n)

  d_rhs <- d_inv * rhs
  middle_rhs <- as.numeric(U %*% d_rhs)
  middle <- tryCatch(
    solve_spd(K, middle_rhs),
    error = function(e) solve(K, middle_rhs)
  )

  as.numeric(d_rhs - d_inv * as.numeric(crossprod(U, middle)))
}

.acvb_update_woodbury <- function(omega,
                                   dat,
                                   sigma_beta_sq,
                                   A_y,
                                   B_y,
                                   logit_pi,
                                   denom_floor = 1e-10,
                                   ridge_n = 1e-10) {
  X <- dat$X
  y <- dat$y
  n <- dat$n
  p <- dat$p
  lambda <- 1 / sigma_beta_sq
  omega <- clamp(as.numeric(omega), 0, 1)

  d <- lambda + dat$x_norm2 * omega * (1 - omega)
  d <- pmax(d, 1e-12)
  scale <- omega / sqrt(d)

  ## K = I_n + X diag(omega^2 / d) X'.
  V <- sweep(X, 2, scale, `*`)
  K <- diag(n) + tcrossprod(V)
  K <- K + diag(ridge_n, n)

  ## One factorization solves all p columns and y at once.
  rhs <- cbind(X, y)
  Ksol <- tryCatch(
    solve_spd(K, rhs),
    error = function(e) solve(K, rhs)
  )
  Kinv_X <- Ksol[, seq_len(p), drop = FALSE]
  Kinv_y <- Ksol[, p + 1L]

  xKx <- colSums(X * Kinv_X)
  xKy <- as.numeric(crossprod(X, Kinv_y))
  yKy <- as.numeric(sum(y * Kinv_y))

  ## Vectorized leave-one-coordinate quantities obtained by an exact rank-one
  ## downdate from K to K_{-j}.
  leverage <- scale^2 * xKx
  denom_raw <- 1 - leverage
  denom <- pmax(denom_raw, denom_floor)
  xA <- scale * xKx
  yA <- scale * xKy

  y_Kminus_y <- yKy + yA^2 / denom
  x_Kminus_x <- xKx + xA^2 / denom
  x_Kminus_y <- xKy + xA * yA / denom

  quad0 <- dat$yTy - y_Kminus_y
  quad0 <- pmin(pmax(quad0, 0), dat$yTy)
  D1 <- pmax(lambda + x_Kminus_x, 1e-10)
  R1 <- x_Kminus_y
  quad1 <- quad0 + R1^2 / D1

  scale0 <- pmax(B_y + 0.5 * dat$yTy - 0.5 * quad0, 1e-12)
  scale1 <- pmax(B_y + 0.5 * dat$yTy - 0.5 * quad1, 1e-12)
  T0 <- -(A_y + n / 2) * safe_log(scale0)
  T1 <- -(A_y + n / 2) * safe_log(scale1) +
    logit_pi - 0.5 * safe_log(n) - 0.5 * safe_log(sigma_beta_sq)

  list(
    omega = clamp(expit(T1 - T0), .Machine$double.eps, 1 - .Machine$double.eps),
    min_downdate_denom = min(denom_raw),
    n_denom_clamped = sum(denom_raw < denom_floor),
    n_system_dimension = n
  )
}

.acvb_update_cached_inverse_legacy <- function(omega,
                                               dat,
                                               sigma_beta_sq,
                                               A_y,
                                               B_y,
                                               logit_pi) {
  p <- dat$p
  n <- dat$n
  lambda <- 1 / sigma_beta_sq
  G <- .build_acvb_G_cached(omega, dat, sigma_beta_sq)
  Ginv <- solve_spd(G, diag(p))
  omega_prop <- omega

  for (j in seq_len(p)) {
    idx <- setdiff(seq_len(p), j)
    w_m <- omega[idx]
    r_m <- w_m * dat$Xty[idx]
    c_m <- w_m * dat$XtX[idx, j]
    SOL <- cbind(
      .solve_leave_one_cached(Ginv, j, r_m),
      .solve_leave_one_cached(Ginv, j, c_m)
    )
    sol_r <- SOL[, 1]
    sol_c <- SOL[, 2]
    quad0 <- sum(r_m * sol_r)
    D1 <- pmax(dat$x_norm2[j] + lambda - sum(c_m * sol_c), 1e-10)
    R1 <- dat$Xty[j] - sum(c_m * sol_r)
    quad1 <- quad0 + R1^2 / D1

    scale0 <- pmax(B_y + 0.5 * dat$yTy - 0.5 * quad0, 1e-12)
    scale1 <- pmax(B_y + 0.5 * dat$yTy - 0.5 * quad1, 1e-12)
    T0 <- -(A_y + n / 2) * safe_log(scale0)
    T1 <- -(A_y + n / 2) * safe_log(scale1) +
      logit_pi - 0.5 * safe_log(n) - 0.5 * safe_log(sigma_beta_sq)
    omega_prop[j] <- expit(T1 - T0)
  }

  list(
    omega = clamp(omega_prop, .Machine$double.eps, 1 - .Machine$double.eps),
    min_downdate_denom = NA_real_,
    n_denom_clamped = NA_integer_,
    n_system_dimension = p
  )
}

fit_noborrow_acvb <- function(dat,
                              init_method = c("null", "full", "mcp", "lasso", "scad",
                                              "mix", "mild", "custom"),
                              init_omega = NULL,
                              sigma_beta_sq = 0.1,
                              A_y = 1,
                              B_y = 1,
                              pi_prior = NULL,
                              max_iter = 100,
                              tol = 1e-6,
                              gamma_ebic = 0.5,
                              verbose = FALSE,
                              fast_cached_update = TRUE,
                              acceleration = c("woodbury", "cached_inverse", "slow")) {
  init_method <- match.arg(init_method)
  acceleration <- match.arg(acceleration)
  if (!isTRUE(fast_cached_update) && acceleration != "slow") acceleration <- "slow"

  X <- dat$X
  y <- dat$y
  n <- dat$n
  p <- dat$p
  Xty <- dat$Xty
  yTy <- dat$yTy

  if (is.null(pi_prior)) pi_prior <- min(.1,max(1/n,1/sqrt(p)))
  pi_prior <- clamp(pi_prior, .Machine$double.eps, 1 - .Machine$double.eps)
  logit_pi <- qlogis(pi_prior)
  omega <- make_init_omega(
    dat = dat,
    init_method = init_method,
    init_omega = init_omega,
    gamma_ebic = gamma_ebic
  )

  acvb_beta_given_w <- function(w) {
    rhs <- w * Xty
    sol <- if (acceleration == "woodbury") {
      .acvb_solve_woodbury(w, rhs, dat, sigma_beta_sq)
    } else {
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
        upd <- .acvb_update_woodbury(
          omega = omega,
          dat = dat,
          sigma_beta_sq = sigma_beta_sq,
          A_y = A_y,
          B_y = B_y,
          logit_pi = logit_pi
        )
        omega_prop <- upd$omega
        min_downdate_denom <- upd$min_downdate_denom
        n_denom_clamped <- n_denom_clamped + upd$n_denom_clamped

      } else if (acceleration == "cached_inverse" && p > 1L) {
        omega_prop <- .acvb_update_cached_inverse_legacy(
          omega = omega,
          dat = dat,
          sigma_beta_sq = sigma_beta_sq,
          A_y = A_y,
          B_y = B_y,
          logit_pi = logit_pi
        )$omega

      } else {
        ## Slow fallback, matching the original implementation.
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
      cat(sprintf(
        "[ACVB/%s/%s] iter %3d | d_omega=%.3e | size05=%d | iter_time=%.3fs\n",
        init_method, acceleration, t, d_omega, sum(omega > 0.5), iter_time_sec[t]
      ))
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

  sel_05 <- get_selection_metrics_omega(omega, dat$beta_true, thr = 0.5)
  sel_01 <- get_selection_metrics_omega(omega, dat$beta_true, thr = 0.1)
  sel_top <- get_selection_metrics_topk(omega, dat$beta_true)

  list(
    method = paste0("no_borrow_ACVB_init_", init_method, "_", acceleration),
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
    sigma_beta_sq = sigma_beta_sq,
    init_method = init_method,
    n_iter = n_iter,
    converged = converged,
    convergence_reason = convergence_reason,
    d_omega_final = d_omega,
    fast_cached_update = fast_cached_update,
    acvb_acceleration = acceleration,
    acvb_mean_iter_time_sec = mean(iter_time_sec),
    acvb_max_iter_time_sec = max(iter_time_sec),
    acvb_min_downdate_denom = min_downdate_denom,
    acvb_n_denom_clamped = n_denom_clamped
  )
}

.time_it <- function(expr) {
  tm <- system.time(value <- eval.parent(substitute(expr)))
  list(value = value, elapsed = unname(tm["elapsed"]))
}

# Pinned official implementation; missing/invalid backend stops the run.
source(file.path(.this_dir, "translasso_official.R"), local = .GlobalEnv)

fit_lasso_lambda <- function(X, y, lambda) {
  fit <- glmnet::glmnet(
    x = X, y = y, alpha = 1, lambda = lambda,
    intercept = FALSE, standardize = FALSE
  )
  as.numeric(coef(fit, s = lambda))[-1]
}

make_target_subset_for_trans <- function(dat, idx) {
  out <- dat
  out$X <- dat$X[idx, , drop = FALSE]
  out$y <- dat$y[idx]
  out$n <- length(idx)
  out$X_all <- do.call(rbind, c(list(out$X), dat$X_source))
  out$y_all <- c(out$y, unlist(dat$y_source, use.names = FALSE))
  out$n_vec_trans <- c(out$n, dat$n_source_raw)
  out
}

fit_oracle_trans_lasso_simple <- function(dat, A0 = integer(0), method_name = NULL) {
  p <- dat$p
  n0 <- dat$n
  if (length(A0) == 0) {
    lambda0 <- sqrt(2 * log(p) / n0)
    beta_hat <- fit_lasso_lambda(dat$X, dat$y, lambda = lambda0)
    if (is.null(method_name)) method_name <- "Oracle_Trans_Lasso_emptyA_target_lasso"
    return(list(method = method_name, beta_hat = beta_hat, A0 = A0))
  }
  X_pool <- do.call(rbind, c(list(dat$X), dat$X_source[A0]))
  y_pool <- c(dat$y, unlist(dat$y_source[A0], use.names = FALSE))
  n_pool <- nrow(X_pool)
  lambda_w <- sqrt(2 * log(p) / n_pool)
  lambda_delta <- sqrt(2 * log(p) / n0)
  w_hat <- fit_lasso_lambda(X_pool, y_pool, lambda = lambda_w)
  delta_hat <- fit_lasso_lambda(dat$X, dat$y - as.numeric(dat$X %*% w_hat), lambda = lambda_delta)
  beta_hat <- w_hat + delta_hat
  if (is.null(method_name)) method_name <- paste0("Oracle_Trans_Lasso_simple_A_", paste(A0, collapse = "_"))
  list(method = method_name, beta_hat = beta_hat, A0 = A0)
}

beta_metrics <- function(beta_hat, dat) {
  if (is.null(beta_hat) || anyNA(beta_hat)) {
    return(list(pred_mse = NA_real_))
  }
  list(
    pred_mse = prediction_mse(beta_hat, dat)
  )
}

selection_from_beta <- function(beta_hat, dat, thr = 1e-4) {
  if (is.null(beta_hat) || anyNA(beta_hat)) return(calc_selection_metrics(rep(0, dat$p), dat$beta_true))
  get_selection_metrics_beta(beta_hat, dat$beta_true, thr = thr)
}

source_alignment_metrics <- function(dat, alpha_hat = NULL) {
  S_tr <- setdiff(which(dat$tag_true == 1), dat$common_bad)
  empty_out <- list(
    alpha_l1_error = if (is.null(alpha_hat) || anyNA(alpha_hat)) NA_real_ else sum(abs(alpha_hat - dat$alpha_true)),
    weighted_summary_sse = NA_real_,
    best_single_summary_sse = NA_real_,
    source1_summary_sse = NA_real_,
    source2_summary_sse = NA_real_,
    source3_summary_sse = NA_real_
  )
  if (length(S_tr) == 0) return(empty_out)
  if (is.null(alpha_hat)) alpha_hat <- rep(NA_real_, dat$K)
  if (anyNA(alpha_hat)) {
    weighted_sse <- NA_real_
  } else {
    weighted_sse <- sum((as.numeric(dat$B[S_tr, , drop = FALSE] %*% alpha_hat) - dat$beta_true[S_tr])^2)
  }
  single_sse <- sapply(seq_len(dat$K), function(k) {
    sum((dat$B[S_tr, k] - dat$beta_true[S_tr])^2)
  })
  list(
    alpha_l1_error = if (anyNA(alpha_hat)) NA_real_ else sum(abs(alpha_hat - dat$alpha_true)),
    weighted_summary_sse = weighted_sse,
    best_single_summary_sse = min(single_sse),
    source1_summary_sse = single_sse[1],
    source2_summary_sse = if (dat$K >= 2) single_sse[2] else NA_real_,
    source3_summary_sse = if (dat$K >= 3) single_sse[3] else NA_real_
  )
}



borrow_formula_beta <- function(dat,
                                fit,
                                thr = 0.5,
                                source_prec_ratio = 0.25,
                                source_prec_cap = 0.75,
                                ridge = 1e-6) {
  p <- dat$p
  beta_hat <- rep(0, p)

  if (is.null(fit$omega) || is.null(fit$gamma) || is.null(fit$alpha)) {
    if (!is.null(fit$beta_hat)) return(fit$beta_hat)
    return(beta_hat)
  }

  selected <- which(fit$omega > thr)
  if (length(selected) == 0) return(beta_hat)

  alpha <- as.numeric(fit$alpha)
  if (sum(alpha) <= 0 || any(!is.finite(alpha))) {
    alpha <- rep(1 / dat$K, dat$K)
  } else {
    alpha <- alpha / sum(alpha)
  }

  bbar <- as.numeric(as.matrix(dat$B) %*% alpha)

  gamma_s <- as.numeric(fit$gamma[selected])
  gamma_s <- pmin(pmax(gamma_s, 1e-6), 1 - 1e-6)

  ## Active-set posterior precision formula:
  ## beta_S = (X_S'X_S + Sbar_S^{-1})^{-1}
  ##          (X_S'y + Sbar_S^{-1} (C B alpha)_S).
  ## Here Sbar_S^{-1} is diagonal with entries lambda_s.
  ## gamma_j gates borrowing: common-bad coordinates should have gamma_j near zero.
  lambda0 <- source_prec_ratio * dat$n
  lambda0 <- min(lambda0, source_prec_cap * dat$n)
  lambda_s <- lambda0 * gamma_s

  Xs <- dat$X[, selected, drop = FALSE]
  y <- dat$y
  k <- length(selected)

  A <- crossprod(Xs) + diag(lambda_s, nrow = k, ncol = k) + ridge * diag(k)
  rhs <- as.numeric(crossprod(Xs, y)) + lambda_s * bbar[selected]

  beta_s <- solve(A, rhs)
  beta_hat[selected] <- as.numeric(beta_s)
  beta_hat
}

target_refit_beta <- function(dat,
                              fit,
                              thr = 0.5,
                              ridge = 1e-6,
                              beta_thr = 1e-6) {
  p <- dat$p
  beta_hat <- rep(0, p)
  
  ## For ACVB, use omega > 0.5 as the selected set.
  if (!is.null(fit$omega)) {
    selected <- which(fit$omega > thr)
  } else if (!is.null(fit$beta_hat)) {
    ## For lasso / trans-lasso style methods, use nonzero beta support.
    selected <- which(abs(fit$beta_hat) > beta_thr)
  } else {
    return(beta_hat)
  }
  
  if (length(selected) == 0) {
    return(beta_hat)
  }
  
  Xs <- dat$X[, selected, drop = FALSE]
  y <- dat$y
  k <- length(selected)
  
  A <- crossprod(Xs) + ridge * diag(k)
  rhs <- as.numeric(crossprod(Xs, y))
  
  beta_s <- tryCatch({
    R <- chol(A)
    backsolve(R, forwardsolve(t(R), rhs))
  }, error = function(e) {
    solve_spd(A, rhs)
  })
  
  beta_hat[selected] <- as.numeric(beta_s)
  beta_hat
}



## =========================================================
## Collapsed log marginal likelihood diagnostic
## =========================================================
## This function evaluates the expression
## log p(y | z, c, alpha)
## requested for the performance table.  By default it uses
## z_j = 1(omega_j > threshold) and c_j = 1(gamma_j > threshold),
## so that it matches the displayed formula with Z and C as indicator
## matrices.  Set mode = "variational" only for a plug-in diagnostic
## using z_j = omega_j and c_j = gamma_j.

compute_collapsed_logp_y_zc_alpha <- function(alpha,
                                              omega,
                                              gamma,
                                              dat,
                                              sigma_beta_sq = 100,
                                              A_y = 1,
                                              B_y = 1,
                                              mode = c("selected", "variational"),
                                              threshold = 0.5,
                                              ridge = 1e-10) {
  mode <- match.arg(mode)
  if (is.null(alpha) || is.null(omega) || is.null(gamma)) {
    return(list(
      logp = NA_real_, Tval = NA_real_, logdet_S_inv = NA_real_,
      logdet_M = NA_real_, source_quad = NA_real_, quad = NA_real_,
      mode = mode, threshold = threshold,
      selected_z_size = NA_real_, selected_c_size = NA_real_
    ))
  }

  alpha <- .normalize_simplex(alpha, K = dat$K)
  omega <- as.numeric(omega)
  gamma <- as.numeric(gamma)
  if (length(omega) != dat$p || length(gamma) != dat$p) {
    stop("omega and gamma must have length dat$p in compute_collapsed_logp_y_zc_alpha().")
  }

  if (mode == "selected") {
    z <- as.numeric(omega > threshold)
    cvec <- as.numeric(gamma > threshold)
  } else {
    z <- pmin(pmax(omega, 0), 1)
    cvec <- pmin(pmax(gamma, 0), 1)
  }

  S0bar_inv <- 1 / sigma_beta_sq
  ns_eff <- get_ns_eff(dat)
  S1bar_inv <- (ns_eff * dat$x_norm2) / dat$n
  S_inv <- S0bar_inv + (S1bar_inv - S0bar_inv) * cvec
  S_inv <- pmax(S_inv, 1e-12)

  bbar <- as.numeric(as.matrix(dat$B) %*% alpha)
  ## Binary selected mode and paper-moment variational mode share this form:
  ## E(S^{-1} C) = Gamma S1^{-1}.
  rhs <- z * dat$Xty + cvec * S1bar_inv * bbar
  factor <- .trace_woodbury_factor(
    z = z,
    diag_add = S_inv + dat$x_norm2 * z * (1 - z),
    dat = dat,
    ridge = ridge
  )
  sol <- .trace_woodbury_apply(factor, rhs)
  logdet_M <- factor$logdet
  logdet_S_inv <- if (mode == "selected") {
    sum(safe_log(S_inv))
  } else {
    sum(cvec * safe_log(S1bar_inv) +
          (1 - cvec) * safe_log(S0bar_inv))
  }

  quad <- sum(rhs * sol)
  source_quad <- sum(cvec * S1bar_inv * bbar^2)
  Tval <- dat$yTy + source_quad - quad
  Tval <- max(Tval, 1e-10)

  scale <- B_y + 0.5 * Tval
  logp <- -0.5 * dat$n * log(2 * pi) +
    0.5 * logdet_S_inv -
    0.5 * logdet_M -
    (A_y + dat$n / 2) * safe_log(scale) +
    A_y * safe_log(B_y) -
    lgamma(A_y) +
    lgamma(A_y + dat$n / 2)

  list(
    logp = as.numeric(logp),
    Tval = as.numeric(Tval),
    logdet_S_inv = as.numeric(logdet_S_inv),
    logdet_M = as.numeric(logdet_M),
    source_quad = as.numeric(source_quad),
    quad = as.numeric(quad),
    mode = mode,
    threshold = threshold,
    selected_z_size = sum(z > 0),
    selected_c_size = sum(cvec > 0)
  )
}



## =========================================================
## Method-specific criterion functions
## =========================================================
## All Borrow-CVB methods use one paper-moment collapsed variational objective
## as their primary criterion.  Alpha is an EM parameter, so no q(alpha) or
## Dirichlet entropy belongs to the primary ELBO.  Legacy alpha diagnostics are
## retained in separate columns only.
## * ACVB: target-only collapsed VB objective.
## * Trans-Lasso: no variational ELBO; the criterion ELBO is NA and the logp is
##   the target-only post-selection marginal log likelihood on its selected support.

.extract_alpha_eta <- function(obj) {
  if (is.null(obj) || !is.list(obj)) return(NULL)
  candidates <- c("eta", "eta_new", "alpha_eta", "dirichlet_eta", "eta_alpha", "eta_vec")
  for (nm in candidates) {
    if (!is.null(obj[[nm]])) {
      eta <- as.numeric(obj[[nm]])
      if (all(is.finite(eta)) && length(eta) > 0 && all(eta > 0)) return(eta)
    }
  }
  NULL
}

.dirichlet_entropy <- function(eta) {
  eta <- as.numeric(eta)
  if (length(eta) == 0 || any(!is.finite(eta)) || any(eta <= 0)) return(NA_real_)
  eta0 <- sum(eta)
  K <- length(eta)
  as.numeric(sum(lgamma(eta)) - lgamma(eta0) + (eta0 - K) * digamma(eta0) - sum((eta - 1) * digamma(eta)))
}

compute_alpha_em_objective <- function(alpha,
                                       omega,
                                       gamma,
                                       dat,
                                       sigma_beta_sq = 100,
                                       A_y = 1,
                                       B_y = 1) {
  terms <- .compute_alpha_em_qp_terms(
    omega = omega,
    gamma = gamma,
    dat = dat,
    sigma_beta_sq = sigma_beta_sq
  )
  alpha <- .normalize_simplex(alpha, K = dat$K)
  H <- terms$H
  g <- terms$g
  yPy <- terms$yPy
  qp <- .alpha_qp_objective(alpha, H, g)
  rss_alpha <- as.numeric(yPy - 2 * sum(g * alpha) + sum(alpha * as.numeric(H %*% alpha)))
  rss_alpha <- max(rss_alpha, 1e-10)
  scale <- max(B_y + 0.5 * rss_alpha, 1e-12)
  F <- - (A_y + dat$n / 2) * safe_log(scale)
  list(
    F = as.numeric(F),
    qp = as.numeric(qp),
    rss_alpha = as.numeric(rss_alpha),
    yPy = as.numeric(yPy)
  )
}

compute_old_alpha_dirichlet_criterion <- function(fit,
                                                  dat,
                                                  prior_par_borrow,
                                                  sigma_beta_sq = 100,
                                                  A_y = 1,
                                                  B_y = 1,
                                                  include_source_normalization = TRUE) {
  point <- borrow_elbo_surrogate(
    alpha = fit$alpha,
    omega = fit$omega,
    gamma = fit$gamma,
    dat = dat,
    prior_par_borrow = prior_par_borrow,
    sigma_beta_sq = sigma_beta_sq,
    A_y = A_y,
    B_y = B_y,
    include_source_normalization = include_source_normalization
  )
  eta <- fit$alpha_eta %||% .extract_alpha_eta(fit$last_alpha_update_raw)
  H_alpha <- .dirichlet_entropy(eta)
  old_elbo <- point$total
  if (is.finite(H_alpha)) old_elbo <- old_elbo + H_alpha
  list(
    old_alpha_elbo = as.numeric(old_elbo),
    old_alpha_plugin_elbo = as.numeric(point$total),
    old_alpha_dirichlet_entropy = as.numeric(H_alpha),
    old_alpha_eta_sum = if (!is.null(eta)) sum(eta) else NA_real_
  )
}

compute_target_logp_y_z <- function(z,
                                    dat,
                                    sigma_beta_sq = 100,
                                    A_y = 1,
                                    B_y = 1,
                                    mode = c("selected", "variational"),
                                    threshold = 0.5,
                                    ridge = 1e-10) {
  mode <- match.arg(mode)
  if (is.null(z)) {
    return(list(logp = NA_real_, Tval = NA_real_, logdet_M = NA_real_, selected_z_size = NA_real_))
  }
  z <- as.numeric(z)
  if (length(z) != dat$p) stop("z must have length dat$p in compute_target_logp_y_z().")
  if (mode == "selected") {
    z <- as.numeric(z > threshold)
  } else {
    z <- pmin(pmax(z, 0), 1)
  }
  lambda <- 1 / sigma_beta_sq
  M <- dat$XtX * (z %o% z)
  diag(M) <- diag(M) + lambda
  M <- M + diag(ridge, dat$p)
  rhs <- z * dat$Xty
  sol <- tryCatch(solve_spd(M, rhs), error = function(e) solve(M, rhs))
  logdet_M <- as.numeric(determinant(M, logarithm = TRUE)$modulus)
  logdet_S_inv <- dat$p * safe_log(lambda)
  quad <- sum(rhs * sol)
  Tval <- max(dat$yTy - quad, 1e-10)
  scale <- max(B_y + 0.5 * Tval, 1e-12)
  logp <- -0.5 * dat$n * log(2 * pi) +
    0.5 * logdet_S_inv -
    0.5 * logdet_M -
    (A_y + dat$n / 2) * safe_log(scale) +
    A_y * safe_log(B_y) -
    lgamma(A_y) +
    lgamma(A_y + dat$n / 2)
  list(
    logp = as.numeric(logp),
    Tval = as.numeric(Tval),
    logdet_M = as.numeric(logdet_M),
    selected_z_size = sum(z > 0)
  )
}

compute_acvb_elbo_surrogate <- function(omega,
                                        dat,
                                        sigma_beta_sq = 100,
                                        pi_prior = NULL,
                                        A_y = 1,
                                        B_y = 1,
                                        include_entropy = TRUE,
                                        include_indicator_prior = TRUE,
                                        ridge = 1e-10) {
  if (is.null(omega)) return(list(total = NA_real_))
  omega <- clamp(as.numeric(omega), 0, 1)
  if (length(omega) != dat$p) stop("omega must have length dat$p in compute_acvb_elbo_surrogate().")
  if (is.null(pi_prior)) pi_prior <- mean(omega)
  pi_prior <- clamp(pi_prior, .Machine$double.eps, 1 - .Machine$double.eps)
  lambda <- 1 / sigma_beta_sq
  M <- dat$XtX * (omega %o% omega)
  diag(M) <- dat$x_norm2 * omega + lambda
  M <- M + diag(ridge, dat$p)
  rhs <- omega * dat$Xty
  sol <- tryCatch(solve_spd(M, rhs), error = function(e) solve(M, rhs))
  logdet_M <- as.numeric(determinant(M, logarithm = TRUE)$modulus)
  logdet_S_inv <- dat$p * safe_log(lambda)
  quad <- sum(rhs * sol)
  Tval <- max(dat$yTy - quad, 1e-10)
  collapsed <- 0.5 * logdet_S_inv - 0.5 * logdet_M -
    (A_y + dat$n / 2) * safe_log(B_y + 0.5 * Tval)
  prior <- if (include_indicator_prior) sum(omega) * qlogis(pi_prior) else 0
  entropy <- if (include_entropy) .bernoulli_entropy_sum(omega) else 0
  list(
    total = as.numeric(collapsed + prior + entropy),
    collapsed = as.numeric(collapsed),
    indicator_prior = as.numeric(prior),
    entropy = as.numeric(entropy),
    Tval = as.numeric(Tval),
    logdet_M = as.numeric(logdet_M)
  )
}

compute_method_criteria <- function(name,
                                    fit,
                                    dat,
                                    prior_par_borrow = NULL,
                                    source_prec_ratio = 0.25,
                                    source_prec_cap = 0.75,
                                    refit_ridge = 1e-6,
                                    collapsed_logp_mode = c("selected", "variational"),
                                    collapsed_logp_threshold = 0.5,
                                    include_source_normalization = TRUE) {
  collapsed_logp_mode <- match.arg(collapsed_logp_mode)
  out <- list(
    criterion_name = NA_character_,
    criterion_elbo = NA_real_,
    criterion_logp = NA_real_,
    point_alpha_elbo = NA_real_,
    point_alpha_logp = NA_real_,
    old_alpha_elbo = NA_real_,
    old_alpha_plugin_elbo = NA_real_,
    old_alpha_dirichlet_entropy = NA_real_,
    old_alpha_eta_sum = NA_real_,
    em_alpha_F = NA_real_,
    em_alpha_qp_objective = NA_real_,
    em_alpha_rss = NA_real_,
    acvb_elbo = NA_real_,
    target_logp_y_z = NA_real_,
    criterion_direction = "larger_is_better"
  )

  is_borrow <- grepl("borrow_cvb", name) && !is.null(fit$omega) && !is.null(fit$gamma) && !is.null(fit$alpha)
  if (is_borrow) {
    if (is.null(prior_par_borrow)) {
      out$criterion_name <- "borrow_criterion_not_computed_missing_prior"
      return(out)
    }
    sigsq <- fit$sigma_beta_sq %||% 100
    point <- borrow_elbo_surrogate(
      alpha = fit$alpha,
      omega = fit$omega,
      gamma = fit$gamma,
      dat = dat,
      prior_par_borrow = prior_par_borrow,
      sigma_beta_sq = sigsq,
      include_source_normalization = include_source_normalization
    )
    logp <- compute_collapsed_logp_y_zc_alpha(
      alpha = fit$alpha,
      omega = fit$omega,
      gamma = fit$gamma,
      dat = dat,
      sigma_beta_sq = sigsq,
      mode = collapsed_logp_mode,
      threshold = collapsed_logp_threshold
    )
    em <- tryCatch(compute_alpha_em_objective(
      alpha = fit$alpha,
      omega = fit$omega,
      gamma = fit$gamma,
      dat = dat,
      sigma_beta_sq = sigsq
    ), error = function(e) list(F = NA_real_, qp = NA_real_, rss_alpha = NA_real_))
    old <- if (identical(fit$alpha_update_strategy, "old") || grepl("old_alpha", name)) {
      compute_old_alpha_dirichlet_criterion(
        fit = fit,
        dat = dat,
        prior_par_borrow = prior_par_borrow,
        sigma_beta_sq = sigsq,
        include_source_normalization = include_source_normalization
      )
    } else {
      list(old_alpha_elbo = NA_real_, old_alpha_plugin_elbo = NA_real_, old_alpha_dirichlet_entropy = NA_real_, old_alpha_eta_sum = NA_real_)
    }
    out$point_alpha_elbo <- point$total
    out$point_alpha_logp <- logp$logp
    out$em_alpha_F <- em$F
    out$em_alpha_qp_objective <- em$qp
    out$em_alpha_rss <- em$rss_alpha
    out$old_alpha_elbo <- old$old_alpha_elbo
    out$old_alpha_plugin_elbo <- old$old_alpha_plugin_elbo
    out$old_alpha_dirichlet_entropy <- old$old_alpha_dirichlet_entropy
    out$old_alpha_eta_sum <- old$old_alpha_eta_sum

    out$criterion_name <- "paper_moment_collapsed_variational_objective"
    out$criterion_elbo <- point$total
    out$criterion_logp <- logp$logp
    return(out)
  }

  if (!is.null(fit$omega)) {
    sigsq <- fit$sigma_beta_sq %||% 100
    acvb <- compute_acvb_elbo_surrogate(
      omega = fit$omega,
      dat = dat,
      sigma_beta_sq = sigsq,
      pi_prior = fit$pi_prior %||% NULL
    )
    tlogp <- compute_target_logp_y_z(
      z = fit$omega,
      dat = dat,
      sigma_beta_sq = sigsq,
      mode = collapsed_logp_mode,
      threshold = collapsed_logp_threshold
    )
    out$criterion_name <- "target_acvb_collapsed_elbo"
    out$criterion_elbo <- acvb$total
    out$criterion_logp <- tlogp$logp
    out$acvb_elbo <- acvb$total
    out$target_logp_y_z <- tlogp$logp
    return(out)
  }

  if (!is.null(fit$beta_hat)) {
    z <- as.numeric(abs(fit$beta_hat) > 1e-8)
    tlogp <- compute_target_logp_y_z(
      z = z,
      dat = dat,
      sigma_beta_sq = 100,
      mode = "selected",
      threshold = 0.5
    )
    out$criterion_name <- "post_selection_target_logp_no_elbo"
    out$criterion_elbo <- NA_real_
    out$criterion_logp <- tlogp$logp
    out$target_logp_y_z <- tlogp$logp
    return(out)
  }

  out
}

pack_fit_row <- function(name,
                         fit,
                         dat,
                         elapsed = NA_real_,
                         prior_par_borrow = NULL,
                         include_source_normalization = TRUE,
                         source_prec_ratio = 0.25,
                         source_prec_cap = 0.75,
                         refit_ridge = 1e-6,
                         collapsed_logp_mode = c("selected", "variational"),
                         collapsed_logp_threshold = 0.5) {
  collapsed_logp_mode <- match.arg(collapsed_logp_mode)

  beta_raw <- fit$beta_hat
  if (is.null(beta_raw)) beta_raw <- rep(NA_real_, dat$p)

  is_borrow <- grepl("borrow_cvb", name) &&
    !is.null(fit$omega) && !is.null(fit$gamma) && !is.null(fit$alpha)

  ## Target-only post-selection refit on the selected support.
  ## For ACVB this is the main fair post-selection estimator.
  ## For Borrow-CVB it is reported as a second beta estimator.
  beta_target_refit <- if (!is.null(fit$omega) || !is.null(fit$beta_hat)) {
    target_refit_beta(dat = dat, fit = fit, thr = 0.5, ridge = refit_ridge)
  } else {
    beta_raw
  }

  ## Source-assisted posterior precision formula for Borrow-CVB.
  ## For non-Borrow methods, use target-only refit if available.
  if (is_borrow) {
    beta_formula <- borrow_formula_beta(
      dat = dat,
      fit = fit,
      thr = 0.5,
      source_prec_ratio = source_prec_ratio,
      source_prec_cap = source_prec_cap,
      ridge = refit_ridge
    )
  } else {
    beta_formula <- beta_target_refit
  }

  met_raw <- beta_metrics(beta_raw, dat)
  met_target_refit <- beta_metrics(beta_target_refit, dat)
  met_formula <- beta_metrics(beta_formula, dat)

  ## Keep selection metrics consistent with the variational inclusion probability.
  ## Do not use beta_formula != 0 to count TP/FP/FN.
  if (!is.null(fit$omega)) {
    sel <- get_selection_metrics_omega(fit$omega, dat$beta_true, thr = 0.5)
  } else {
    sel <- fit$sel_05 %||% fit$sel_beta %||% selection_from_beta(beta_target_refit, dat)
  }

  transferable <- setdiff(which(dat$tag_true == 1), dat$common_bad)
  align <- source_alignment_metrics(dat, fit$alpha)

  logp_diag <- if (is_borrow) {
    compute_collapsed_logp_y_zc_alpha(
      alpha = fit$alpha,
      omega = fit$omega,
      gamma = fit$gamma,
      dat = dat,
      sigma_beta_sq = fit$sigma_beta_sq %||% 100,
      A_y = 1,
      B_y = 1,
      mode = collapsed_logp_mode,
      threshold = collapsed_logp_threshold
    )
  } else {
    list(logp = NA_real_)
  }

  criteria <- compute_method_criteria(
    name = name,
    fit = fit,
    dat = dat,
    prior_par_borrow = prior_par_borrow,
    source_prec_ratio = source_prec_ratio,
    source_prec_cap = source_prec_cap,
    refit_ridge = refit_ridge,
    collapsed_logp_mode = collapsed_logp_mode,
    collapsed_logp_threshold = collapsed_logp_threshold,
    include_source_normalization = include_source_normalization
  )

  data.frame(
    scenario_id = dat$scenario_id,
    method = name,
    method_name = fit$method %||% name,
    translasso_backend = fit$translasso_backend %||% NA_character_,
    translasso_commit = fit$translasso_commit %||% NA_character_,
    translasso_source_md5 = fit$translasso_source_md5 %||% NA_character_,
    translasso_adapter_md5 = fit$translasso_adapter_md5 %||% NA_character_,
    translasso_split = fit$translasso_split %||% NA_character_,
    translasso_numerics = fit$translasso_numerics %||% NA_character_,
    translasso_rng_policy = fit$translasso_rng_policy %||% NA_character_,
    estimator_definition = if (is_borrow) "paper_model_averaged_beta" else "method_returned_beta",

    pred_mse_raw = met_raw$pred_mse,

    pred_mse_target_refit = met_target_refit$pred_mse,

    pred_mse_formula = met_formula$pred_mse,

    ## Main columns use the returned estimator (TRACE: paper model average).
    ## Post-selection refits remain explicitly named diagnostic columns.
    pred_mse = met_raw$pred_mse,
    pred_mse_definition = PRED_MSE_DEFINITION_V9,
    pred_mse_target = PRED_MSE_TARGET_V9,

    tp = sel$tp,
    fp = sel$fp,
    fn = sel$fn,
    precision = sel$precision,
    recall = sel$recall,
    f1 = sel$f1,
    model_size = sel$model_size,

    mean_omega_signal = if (!is.null(fit$omega)) mean(fit$omega[dat$tag_true == 1]) else NA_real_,
    mean_omega_null = if (!is.null(fit$omega)) mean(fit$omega[dat$tag_true == 0]) else NA_real_,
    mean_omega_common_bad = if (!is.null(fit$omega) && length(dat$common_bad) > 0) mean(fit$omega[dat$common_bad]) else NA_real_,
    mean_omega_transferable_signal = if (!is.null(fit$omega) && length(transferable) > 0) mean(fit$omega[transferable]) else NA_real_,
    mean_omega_spurious_null = if (!is.null(fit$omega) && length(dat$spurious_union) > 0) mean(fit$omega[dat$spurious_union]) else NA_real_,
    mean_gamma_signal = if (!is.null(fit$gamma)) mean(fit$gamma[dat$tag_true == 1]) else NA_real_,
    mean_gamma_common_bad = if (!is.null(fit$gamma) && length(dat$common_bad) > 0) mean(fit$gamma[dat$common_bad]) else NA_real_,
    mean_gamma_transferable_signal = if (!is.null(fit$gamma) && length(transferable) > 0) mean(fit$gamma[transferable]) else NA_real_,
    mean_gamma_spurious_null = if (!is.null(fit$gamma) && length(dat$spurious_union) > 0) mean(fit$gamma[dat$spurious_union]) else NA_real_,

    alpha_hat = if (!is.null(fit$alpha)) paste(round(fit$alpha, 4), collapse = ";") else NA_character_,
    alpha_l1_error = align$alpha_l1_error,
    weighted_summary_sse = align$weighted_summary_sse,
    best_single_summary_sse = align$best_single_summary_sse,
    source1_summary_sse = align$source1_summary_sse,
    source2_summary_sse = align$source2_summary_sse,
    source3_summary_sse = align$source3_summary_sse,

    collapsed_logp_y_zc_alpha = logp_diag$logp,

    criterion_name = criteria$criterion_name,
    criterion_elbo = criteria$criterion_elbo,
    criterion_logp = criteria$criterion_logp,
    criterion_direction = criteria$criterion_direction,
    point_alpha_elbo = criteria$point_alpha_elbo,
    point_alpha_logp = criteria$point_alpha_logp,
    old_alpha_elbo = criteria$old_alpha_elbo,
    old_alpha_plugin_elbo = criteria$old_alpha_plugin_elbo,
    old_alpha_dirichlet_entropy = criteria$old_alpha_dirichlet_entropy,
    old_alpha_eta_sum = criteria$old_alpha_eta_sum,
    em_alpha_F = criteria$em_alpha_F,
    em_alpha_qp_objective = criteria$em_alpha_qp_objective,
    em_alpha_rss = criteria$em_alpha_rss,
    acvb_elbo = criteria$acvb_elbo,
    target_logp_y_z = criteria$target_logp_y_z,

    source_prec_ratio = source_prec_ratio,
    source_prec_cap = source_prec_cap,
    refit_ridge = refit_ridge,
    time_sec = elapsed,
    n_iter = fit$n_iter %||% NA_integer_,
    stringsAsFactors = FALSE
  )
}


prediction_diagnostics <- function(fit,
                                 dat,
                                 method_key = NULL,
                                 source_prec_ratio = 0.25,
                                 source_prec_cap = 0.75,
                                 refit_ridge = 1e-6) {
  beta_raw <- fit$beta_hat
  if (is.null(beta_raw) || anyNA(beta_raw)) {
    return(c(
      pred_mse_raw = NA_real_,
      pred_mse_target_refit = NA_real_,
      pred_mse_formula = NA_real_,
      pred_mse = NA_real_
    ))
  }

  is_borrow <- !is.null(method_key) && grepl("borrow_cvb", method_key) &&
    !is.null(fit$omega) && !is.null(fit$gamma) && !is.null(fit$alpha)

  beta_target_refit <- if (!is.null(fit$omega) || !is.null(fit$beta_hat)) {
    target_refit_beta(dat = dat, fit = fit, thr = 0.5, ridge = refit_ridge)
  } else {
    beta_raw
  }

  if (is_borrow) {
    beta_formula <- borrow_formula_beta(
      dat = dat,
      fit = fit,
      thr = 0.5,
      source_prec_ratio = source_prec_ratio,
      source_prec_cap = source_prec_cap,
      ridge = refit_ridge
    )
  } else {
    beta_formula <- beta_target_refit
  }

  c(
    pred_mse_raw = unname(prediction_mse(beta_raw, dat)),

    pred_mse_target_refit = unname(prediction_mse(beta_target_refit, dat)),

    pred_mse_formula = unname(prediction_mse(beta_formula, dat)),

    pred_mse = unname(prediction_mse(beta_formula, dat))
  )
}


run_one <- function(
    seed = 3036,
    z_prior = "dimension_mild",
    gamma_prior = "constant_half",
    pi_omega = NULL,

    source_precision_mode = c("high_dimensional_logp", "raw_ns", "target_capped"),
    source_precision_cap_ratio = 1,

    init_method_borrow = "full",
    init_alpha_borrow = "uniform",
    init_gamma_borrow = "theory_zero",
    warm_start_learned_from_fixed = FALSE,
    theory_initialization = TRUE,
    target_init_until_converged = FALSE,
    target_init_steps = 5L,
    target_init_max_iter = 1000L,
    target_init_tol = NULL,
    preliminary_alpha_update = "huber",
    huber_tau = NULL,
    huber_tau_rule = c("fixed", "diverging"),
    huber_tau_multiplier = 1,
    use_candidate_screening = FALSE,
    candidate_size = 180,
    fast_cached_update = TRUE,
    borrow_update_schedule = c("paper_block", "gauss_seidel", "jacobi"),
    warm_start_iter = 25,

    sigma_beta_sq_borrow = 0.1,
    max_iter_borrow = 150,
    tol_borrow = 1e-6,

    sigma_beta_sq_acvb = 100,
    max_iter_acvb = 100,
    tol_acvb = 1e-6,
    acvb_acceleration = c("woodbury", "cached_inverse", "slow"),

    source_prec_ratio = 0.25,
    source_prec_cap = 0.75,
    refit_ridge = 1e-6,
    include_source_normalization = TRUE,

    gamma_ebic = 0.5,
    gamma_ebic_classical = 0.25,

    run_fixed_alpha = TRUE,
    run_oracle_alpha = TRUE,
    run_single_source_alpha = FALSE,
    run_vbss = FALSE,
    run_mcp = TRUE,
    run_scad = FALSE,
    run_oracle_trans_simple = TRUE,
    run_naive_trans_simple = TRUE,
    verbose = FALSE,
    ...
) {
  acvb_acceleration <- match.arg(acvb_acceleration)
  borrow_update_schedule <- match.arg(borrow_update_schedule)
  source_precision_mode <- match.arg(source_precision_mode)
  if (is.null(target_init_tol)) target_init_tol <- tol_borrow
  preliminary_alpha_update <- match.arg(
    preliminary_alpha_update,
    choices = c("huber", "target_matching_qp", "none")
  )
  huber_tau_rule <- match.arg(huber_tau_rule)
  if (isTRUE(theory_initialization) && isTRUE(warm_start_learned_from_fixed)) {
    if (isTRUE(verbose)) {
      cat("[run_one] theory_initialization=TRUE: disabling fixed-alpha warm start so the first cycle starts from gamma^(0)=0.\n")
    }
    warm_start_learned_from_fixed <- FALSE
  }

  ## Only pass DGP arguments to make_dat().
  dots <- list(...)
  dgp_names <- names(formals(make_dat))
  dgp_dots <- dots[names(dots) %in% dgp_names]
  dat <- do.call(make_dat, c(list(seed = seed), dgp_dots))
  dat <- configure_trace_source_precision(
    dat = dat,
    mode = source_precision_mode,
    cap_ratio = source_precision_cap_ratio
  )

  nzc_true <- sum(dat$beta_true != 0)

  prior_par_borrow <- make_prior_params_borrow_commonbad(
    n = dat$n,
    p = dat$p,
    nzc = nzc_true,
    z_prior = z_prior,
    gamma_prior = gamma_prior,
    ns_eff = trace_get_raw_ns_eff(dat),
    min_A_z = 2,
    min_conc_z = 50
  )

  prior_par_vbss <- make_prior_params_vbss_commonbad(
    n = dat$n,
    p = dat$p,
    nzc = nzc_true,
    z_prior = z_prior,
    min_A_z = 2,
    min_conc_z = 50
  )

  fits <- list()
  times <- list()

  fit_borrow_fast_wrapper <- function(...) {
    fit_borrow_with_screening(
      ...,
      use_candidate_screening = use_candidate_screening,
      candidate_size = candidate_size,
      fast_cached_update = fast_cached_update,
      borrow_update_schedule = borrow_update_schedule,
      theory_initialization = theory_initialization,
      target_init_until_converged = target_init_until_converged,
      target_init_steps = target_init_steps,
      target_init_max_iter = target_init_max_iter,
      target_init_tol = target_init_tol,
      preliminary_alpha_update = preliminary_alpha_update,
      huber_tau = huber_tau,
      huber_tau_rule = huber_tau_rule,
      huber_tau_multiplier = huber_tau_multiplier
    )
  }

  ## =========================================================
  ## Borrow-CVB learned alpha, using the old/stable update path.
  ## Optional warm start: fixed uniform alpha first, then learned alpha.
  ## =========================================================
  fit_fixed_warm <- NULL
  time_fixed_warm <- 0

  if (warm_start_learned_from_fixed) {
    tmp_fixed_warm <- .time_it(fit_borrow_fast_wrapper(
      dat = dat,
      prior_par_borrow = prior_par_borrow,
      alpha_method = "fixed",
      alpha_fixed = rep(1 / dat$K, dat$K),
      init_alpha = "uniform",
      init_method = init_method_borrow,
      init_gamma = init_gamma_borrow,
      sigma_beta_sq = sigma_beta_sq_borrow,
      max_iter = min(warm_start_iter, max_iter_borrow),
      tol = tol_borrow,
      gamma_ebic = gamma_ebic,
      verbose = verbose
    ))
    fit_fixed_warm <- tmp_fixed_warm$value
    time_fixed_warm <- tmp_fixed_warm$elapsed

    tmp_learned <- .time_it(fit_borrow_fast_wrapper(
      dat = dat,
      prior_par_borrow = prior_par_borrow,
      alpha_method = "old",
      init_alpha = init_alpha_borrow,
      init_method = "custom",
      init_omega = fit_fixed_warm$omega,
      init_gamma = fit_fixed_warm$gamma,
      sigma_beta_sq = sigma_beta_sq_borrow,
      max_iter = max_iter_borrow,
      tol = tol_borrow,
      gamma_ebic = gamma_ebic,
      verbose = verbose
    ))
    fit_learned <- tmp_learned$value
    fit_learned$method <- paste0(fit_learned$method, "_warm_fixed_alpha")
    fit_learned$warm_start_from <- "fixed_uniform_alpha"
    fit_learned$total_time_with_warm_start <- time_fixed_warm + tmp_learned$elapsed

    fits$borrow_cvb_learn_alpha <- fit_learned
    times$borrow_cvb_learn_alpha <- time_fixed_warm + tmp_learned$elapsed
  } else {
    tmp_learned <- .time_it(fit_borrow_fast_wrapper(
      dat = dat,
      prior_par_borrow = prior_par_borrow,
      alpha_method = "old",
      init_alpha = init_alpha_borrow,
      init_method = init_method_borrow,
      init_gamma = init_gamma_borrow,
      sigma_beta_sq = sigma_beta_sq_borrow,
      max_iter = max_iter_borrow,
      tol = tol_borrow,
      gamma_ebic = gamma_ebic,
      verbose = verbose
    ))
    fits$borrow_cvb_learn_alpha <- tmp_learned$value
    times$borrow_cvb_learn_alpha <- tmp_learned$elapsed
  }

  if (run_fixed_alpha) {
    if (!is.null(fit_fixed_warm)) {
      fits$borrow_cvb_fixed_uniform <- fit_fixed_warm
      times$borrow_cvb_fixed_uniform <- time_fixed_warm
    } else {
      tmp_fixed <- .time_it(fit_borrow_fast_wrapper(
        dat = dat,
        prior_par_borrow = prior_par_borrow,
        alpha_method = "fixed",
        alpha_fixed = rep(1 / dat$K, dat$K),
        init_alpha = "uniform",
        init_method = init_method_borrow,
        init_gamma = init_gamma_borrow,
        sigma_beta_sq = sigma_beta_sq_borrow,
        max_iter = max_iter_borrow,
        tol = tol_borrow,
        gamma_ebic = gamma_ebic,
        verbose = verbose
      ))
      fits$borrow_cvb_fixed_uniform <- tmp_fixed$value
      times$borrow_cvb_fixed_uniform <- tmp_fixed$elapsed
    }
  }

  if (run_oracle_alpha) {
    tmp_oracle <- .time_it(fit_borrow_fast_wrapper(
      dat = dat,
      prior_par_borrow = prior_par_borrow,
      alpha_method = "oracle",
      init_alpha = "oracle",
      init_method = init_method_borrow,
      init_gamma = init_gamma_borrow,
      sigma_beta_sq = sigma_beta_sq_borrow,
      max_iter = max_iter_borrow,
      tol = tol_borrow,
      gamma_ebic = gamma_ebic,
      verbose = verbose
    ))
    fits$borrow_cvb_oracle_alpha <- tmp_oracle$value
    times$borrow_cvb_oracle_alpha <- tmp_oracle$elapsed
  }

  if (run_single_source_alpha) {
    for (k in seq_len(dat$K)) {
      alpha_k <- rep(0, dat$K); alpha_k[k] <- 1
      nm <- paste0("borrow_cvb_source", k, "_only")
      tmp <- .time_it(fit_borrow_fast_wrapper(
        dat = dat,
        prior_par_borrow = prior_par_borrow,
        alpha_method = "fixed",
        alpha_fixed = alpha_k,
        init_alpha = "uniform",
        init_method = init_method_borrow,
        init_gamma = init_gamma_borrow,
        sigma_beta_sq = sigma_beta_sq_borrow,
        max_iter = max_iter_borrow,
        tol = tol_borrow,
        gamma_ebic = gamma_ebic,
        verbose = verbose
      ))
      fits[[nm]] <- tmp$value
      times[[nm]] <- tmp$elapsed
    }
  }

  ## Target-only ACVB.
  pi_prior_acvb <- prior_par_borrow$rho_z_mean
  tmp <- .time_it(fit_noborrow_acvb(
    dat = dat,
    init_method = init_method_borrow,
    pi_prior = pi_prior_acvb,
    max_iter = max_iter_acvb,
    sigma_beta_sq = sigma_beta_sq_acvb,
    tol = tol_acvb,
    gamma_ebic = gamma_ebic,
    verbose = verbose,
    fast_cached_update = TRUE,
    acceleration = acvb_acceleration
  ))
  fits$no_borrow_acvb <- tmp$value
  times$no_borrow_acvb <- tmp$elapsed

  if (run_vbss) {
    tmp <- .time_it(fit_noborrow_vbss_fixedhyper(
      dat = dat,
      prior_par_vbss = prior_par_vbss,
      init_method = "mix",
      max_iter = max_iter_acvb,
      tol = tol_acvb,
      gamma_ebic = gamma_ebic,
      verbose = verbose
    ))
    fits$no_borrow_vbss <- tmp$value
    times$no_borrow_vbss <- tmp$elapsed
  }

  tmp <- .time_it(fit_noborrow_classical(dat = dat, method = "lasso", gamma_ebic = gamma_ebic_classical))
  fits$no_borrow_lasso <- tmp$value
  times$no_borrow_lasso <- tmp$elapsed

  if (run_mcp) {
    tmp <- .time_it(fit_noborrow_classical(dat = dat, method = "mcp", gamma_ebic = gamma_ebic_classical))
    fits$no_borrow_mcp <- tmp$value
    times$no_borrow_mcp <- tmp$elapsed
  }

  if (run_scad) {
    tmp <- .time_it(fit_noborrow_classical(dat = dat, method = "scad", gamma_ebic = gamma_ebic_classical))
    fits$no_borrow_scad <- tmp$value
    times$no_borrow_scad <- tmp$elapsed
  }

  tmp <- .time_it(fit_trans_lasso_main(dat))
  fits$trans_lasso <- tmp$value
  times$trans_lasso <- tmp$elapsed

  if (run_oracle_trans_simple) {
    tmp <- .time_it(fit_oracle_trans_lasso_simple(
      dat,
      A0 = dat$A_trans_oracle,
      method_name = if (length(dat$A_trans_oracle) == 0) {
        "Oracle_Trans_Lasso_simple_emptyA"
      } else {
        paste0("Oracle_Trans_Lasso_simple_A_", paste(dat$A_trans_oracle, collapse = "_"))
      }
    ))
    fits$oracle_trans_lasso_simple <- tmp$value
    times$oracle_trans_lasso_simple <- tmp$elapsed
  }

  if (run_naive_trans_simple) {
    tmp <- .time_it(fit_oracle_trans_lasso_simple(
      dat,
      A0 = seq_len(dat$K),
      method_name = "Naive_Trans_Lasso_simple_all_sources"
    ))
    fits$naive_trans_lasso_simple <- tmp$value
    times$naive_trans_lasso_simple <- tmp$elapsed
  }

  performance <- do.call(rbind, lapply(names(fits), function(nm) {
    pack_fit_row(
      nm,
      fits[[nm]],
      dat,
      elapsed = times[[nm]] %||% NA_real_,
      prior_par_borrow = prior_par_borrow,
      include_source_normalization = include_source_normalization,
      source_prec_ratio = source_prec_ratio,
      source_prec_cap = source_prec_cap,
      refit_ridge = refit_ridge
    )
  }))
  performance <- .trace_apply_public_method_labels(performance)

  decomp <- do.call(rbind, lapply(names(fits), function(nm) {
    prediction_diagnostics(
      fits[[nm]],
      dat = dat,
      method_key = nm,
      source_prec_ratio = source_prec_ratio,
      source_prec_cap = source_prec_cap,
      refit_ridge = refit_ridge
    )
  }))
  decomp <- data.frame(scenario_id = dat$scenario_id, method = names(fits), decomp, row.names = NULL)

  S <- which(dat$beta_true != 0)
  beta_formula_borrow <- borrow_formula_beta(
    dat = dat,
    fit = fits$borrow_cvb_learn_alpha,
    thr = 0.5,
    source_prec_ratio = source_prec_ratio,
    source_prec_cap = source_prec_cap,
    ridge = refit_ridge
  )

  beta_refit_borrow <- target_refit_beta(
    dat = dat,
    fit = fits$borrow_cvb_learn_alpha,
    thr = 0.5,
    ridge = refit_ridge
  )

  beta_refit_acvb <- target_refit_beta(
    dat = dat,
    fit = fits$no_borrow_acvb,
    thr = 0.5,
    ridge = refit_ridge
  )

  signal_table <- data.frame(
    j = S,
    beta_true = dat$beta_true[S],
    gamma_true = dat$gamma_true[S],
    is_common_bad = S %in% dat$common_bad,
    B_source_1 = dat$B[S, 1],
    B_source_2 = if (dat$K >= 2) dat$B[S, 2] else NA_real_,
    B_source_3 = if (dat$K >= 3) dat$B[S, 3] else NA_real_,
    W_source_1 = dat$W_source_true[S, 1],
    W_source_2 = if (dat$K >= 2) dat$W_source_true[S, 2] else NA_real_,
    W_source_3 = if (dat$K >= 3) dat$W_source_true[S, 3] else NA_real_,
    omega_borrow = fits$borrow_cvb_learn_alpha$omega[S],
    gamma_borrow = fits$borrow_cvb_learn_alpha$gamma[S],
    beta_borrow_raw = fits$borrow_cvb_learn_alpha$beta_hat[S],
    beta_borrow_target_refit = beta_refit_borrow[S],
    beta_borrow_formula = beta_formula_borrow[S],
    beta_acvb = beta_refit_acvb[S],
    beta_acvb_raw = fits$no_borrow_acvb$beta_hat[S],
    beta_acvb_refit = beta_refit_acvb[S],
    beta_trans = fits$trans_lasso$beta_hat[S]
  )

  scenario <- data.frame(
    scenario_id = dat$scenario_id,
    scenario_title = dat$scenario_title,
    seed = seed,
    n = dat$n,
    p = dat$p,
    nzc = nzc_true,
    rho = dat$rho,
    snr = dat$snr,
    noise_mode = dat$noise_mode %||% NA_character_,
    sigma_true = dat$sigma_true %||% NA_real_,
    snr_requested = dat$snr_requested %||% NA_real_,
    realized_snr = dat$realized_snr %||% dat$snr,
    signal_profile = dat$signal_profile %||% NA_character_,
    theory_signal_magnitude = dat$theory_signal_magnitude %||% NA_real_,
    common_bad_signal_magnitude = dat$common_bad_signal_magnitude %||% NA_real_,
    transferable_signal_magnitude = dat$transferable_signal_magnitude %||% NA_real_,
    beta_min = dat$beta_min %||% NA_real_,
    beta_max = dat$beta_max %||% NA_real_,
    beta_min_information = dat$beta_min_information %||% NA_real_,
    beta_min_log_benchmark = dat$beta_min_log_benchmark %||% NA_real_,
    beta_min_margin_ratio = dat$beta_min_margin_ratio %||% NA_real_,
    common_bad_gap_population = dat$common_bad_gap_population %||% NA_real_,
    transferable_count = dat$transferable_count %||% length(dat$regular_signal),
    source_weight_identification_eigenvalue =
      dat$source_weight_identification_eigenvalue %||% NA_real_,
    source_weight_identification_pass =
      dat$source_weight_identification_pass %||% NA,
    exact_combination_max_error =
      dat$exact_combination_max_error %||% NA_real_,
    shared_spurious_count =
      dat$shared_spurious_count %||% length(dat$spurious_union),
    shared_spurious_gap_population =
      dat$shared_spurious_gap_population %||% NA_real_,
    common_bad_frac = dat$common_bad_frac %||% NA_real_,
    common_bad_count = dat$common_bad_count %||% length(dat$common_bad),
    K = dat$K,
    ns_summary = dat$ns,
    n_source_raw = paste(dat$n_source_raw, collapse = ";"),
    source_type = dat$source_type,
    source_scale = paste(dat$source_scale, collapse = ";"),
    common_bad = paste(dat$common_bad, collapse = ";"),
    common_bad_scale = paste(dat$common_bad_scale, collapse = ";"),
    n_spur_each = paste(dat$n_spur_each, collapse = ";"),
    spur_size = paste(dat$spur_size, collapse = ";"),
    alpha_true = paste(round(dat$alpha_true, 4), collapse = ";"),
    A_trans_oracle = paste(dat$A_trans_oracle, collapse = ";"),
    z_prior = z_prior,
    z_prior_rule = prior_par_borrow$z_prior_rule %||% z_prior,
    z_prior_uses_true_nzc = prior_par_borrow$z_prior_uses_true_nzc %||% NA,
    gamma_prior = gamma_prior,
    source_precision_mode = dat$trace_source_precision_mode,
    source_precision_cap_ratio = dat$trace_source_precision_cap_ratio,
    source_precision_ns_raw = dat$trace_source_precision_ns_raw,
    source_precision_ns_eff = dat$trace_source_precision_ns_eff,
    source_precision_power = dat$trace_source_precision_power,
    pi_omega = if (is.null(pi_omega)) NA_real_ else pi_omega,
    init_method_borrow = init_method_borrow,
    init_alpha_borrow = init_alpha_borrow,
    init_gamma_borrow = paste(init_gamma_borrow, collapse = ";"),
    warm_start_learned_from_fixed = warm_start_learned_from_fixed,
    theory_initialization = theory_initialization,
      target_init_until_converged = target_init_until_converged,
      target_init_steps = target_init_steps,
      target_init_max_iter = target_init_max_iter,
      target_init_tol = target_init_tol,
    preliminary_alpha_update = preliminary_alpha_update,
    use_candidate_screening = use_candidate_screening,
    candidate_size = if (is.null(candidate_size)) NA_integer_ else candidate_size,
    fast_cached_update = fast_cached_update,
    borrow_update_schedule = borrow_update_schedule,
    warm_start_iter = warm_start_iter,
    sigma_beta_sq_borrow = sigma_beta_sq_borrow,
    max_iter_borrow = max_iter_borrow,
    tol_borrow = tol_borrow,
    sigma_beta_sq_acvb = sigma_beta_sq_acvb,
    max_iter_acvb = max_iter_acvb,
    tol_acvb = tol_acvb,
    source_prec_ratio = source_prec_ratio,
    source_prec_cap = source_prec_cap,
    refit_ridge = refit_ridge,
    gamma_ebic = gamma_ebic,
    gamma_ebic_classical = gamma_ebic_classical,
    stringsAsFactors = FALSE
  )

  list(
    scenario = scenario,
    dat = dat,
    fits = fits,
    performance = performance,
    decomp = decomp,
    signal_table = signal_table
  )
}

run_replicates <- function(nrep = 100,
                           base_seed = 3036,
                           save_dir = "results_simulation4_old_method_revised",
                           ...) {
  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)

  all_perf <- list()
  all_decomp <- list()
  all_scen <- list()

  for (r in seq_len(nrep)) {
    seed_r <- base_seed + r - 1
    cat("\n================ simulation4 replicate ", r, "/", nrep, " seed=", seed_r, " ================\n", sep = "")
    ans <- run_one(seed = seed_r, ...)

    perf_r <- ans$performance
    perf_r$rep <- r
    perf_r$seed <- seed_r

    decomp_r <- ans$decomp
    decomp_r$rep <- r
    decomp_r$seed <- seed_r

    scen_r <- ans$scenario
    scen_r$rep <- r

    all_perf[[r]] <- perf_r
    all_decomp[[r]] <- decomp_r
    all_scen[[r]] <- scen_r

    write.csv(perf_r, file.path(save_dir, paste0("performance_rep", r, ".csv")), row.names = FALSE)
  }

  performance <- do.call(rbind, all_perf)
  decomp <- do.call(rbind, all_decomp)
  scenarios <- do.call(rbind, all_scen)

  write.csv(performance, file.path(save_dir, "performance_all.csv"), row.names = FALSE)
  write.csv(decomp, file.path(save_dir, "beta_decomposition_all.csv"), row.names = FALSE)
  write.csv(scenarios, file.path(save_dir, "scenarios_all.csv"), row.names = FALSE)

  metric_cols <- intersect(
    c("pred_mse",
      "pred_mse_raw",
      "pred_mse_target_refit",
      "pred_mse_formula",
      "tp", "fp", "fn", "precision", "recall", "f1", "model_size", "time_sec",
      "alpha_l1_error", "weighted_summary_sse", "best_single_summary_sse",
      "mean_gamma_common_bad", "mean_gamma_transferable_signal"),
    names(performance)
  )

  summary_mean <- aggregate(performance[, metric_cols], by = list(method = performance$method), FUN = function(x) mean(x, na.rm = TRUE))
  summary_sd <- aggregate(performance[, metric_cols], by = list(method = performance$method), FUN = function(x) sd(x, na.rm = TRUE))

  write.csv(summary_mean, file.path(save_dir, "summary_mean.csv"), row.names = FALSE)
  write.csv(summary_sd, file.path(save_dir, "summary_sd.csv"), row.names = FALSE)

  list(
    performance = performance,
    decomp = decomp,
    scenarios = scenarios,
    summary_mean = summary_mean,
    summary_sd = summary_sd
  )
}

## =========================================================
## Small ELBO diagnostic runner
## =========================================================
## This runner intentionally avoids the large grid.  It keeps n = 300,
## p = 600, and nzc in c(5, 10), and saves one ELBO trace per method.

.add_elbo_to_perf <- function(perf, fit) {
  perf$final_elbo <- fit$final_elbo %||% NA_real_
  perf$elbo_monitor_name <- fit$elbo_monitor_name %||%
    "paper_moment_collapsed_variational_objective"
  if (!is.null(fit$elbo_trace) && nrow(fit$elbo_trace) > 0) {
    perf$initial_elbo <- fit$elbo_trace$elbo[1]
    perf$total_elbo_gain <- tail(fit$elbo_trace$elbo, 1) - fit$elbo_trace$elbo[1]
    finite_drop <- fit$elbo_trace$delta_elbo[is.finite(fit$elbo_trace$delta_elbo)]
    perf$max_elbo_drop <- if (length(finite_drop) > 0) min(finite_drop, na.rm = TRUE) else NA_real_
  } else {
    perf$initial_elbo <- NA_real_
    perf$total_elbo_gain <- NA_real_
    perf$max_elbo_drop <- NA_real_
  }
  perf$elbo_drop_count <- fit$elbo_drop_count %||% NA_integer_
  perf$converged <- fit$converged %||% NA
  perf$convergence_reason <- fit$convergence_reason %||% NA_character_
  perf$borrow_update_schedule <- fit$borrow_update_schedule %||% NA_character_
  perf$trace_linear_algebra_backend <- fit$trace_linear_algebra_backend %||% NA_character_
  perf$alpha_update_every <- fit$alpha_update_every %||% NA_integer_
  perf$alpha_update_strategy <- fit$alpha_update_strategy %||% NA_character_
  perf$alpha_freeze_iter <- fit$alpha_freeze_iter %||% NA_integer_
  perf$alpha_update_accept_count <- fit$alpha_update_accept_count %||% NA_integer_
  perf$alpha_update_reject_count <- fit$alpha_update_reject_count %||% NA_integer_
  perf$last_alpha_step_size <- fit$last_alpha_step_size %||% NA_real_
  perf$last_alpha_update_qp_old <- fit$last_alpha_update_qp_old %||% NA_real_
  perf$last_alpha_update_qp_new <- fit$last_alpha_update_qp_new %||% NA_real_
  perf$last_alpha_active_set_size <- fit$last_alpha_active_set_size %||% NA_integer_
  perf$last_alpha_active_set <- fit$last_alpha_active_set %||% NA_character_
  perf$theory_initialization <- fit$theory_initialization %||% NA
  perf$preliminary_alpha_update <- fit$preliminary_alpha_update %||% NA_character_
  perf$preliminary_alpha_done <- fit$preliminary_alpha_done %||% NA
  perf$preliminary_alpha_qp_old <- fit$preliminary_alpha_qp_old %||% NA_real_
  perf$preliminary_alpha_qp_new <- fit$preliminary_alpha_qp_new %||% NA_real_
  perf$preliminary_alpha_active_set <- fit$preliminary_alpha_active_set %||% NA_character_
  perf$preliminary_alpha_solver <- fit$preliminary_alpha_solver %||% NA_character_
  perf$preliminary_alpha_huber_tau <- fit$preliminary_alpha_huber_tau %||% NA_real_
  perf$preliminary_alpha_huber_tau_rule <- fit$preliminary_alpha_huber_tau_rule %||% NA_character_
  perf$preliminary_alpha_huber_s_reference <- fit$preliminary_alpha_huber_s_reference %||% NA_real_
  perf$preliminary_alpha_huber_iterations <- fit$preliminary_alpha_huber_iterations %||% NA_integer_
  perf$preliminary_alpha_huber_converged <- fit$preliminary_alpha_huber_converged %||% NA
  perf$preliminary_alpha_huber_kkt_residual <- fit$preliminary_alpha_huber_kkt_residual %||% NA_real_
  perf$d_omega_final <- fit$d_omega_final %||% NA_real_
  perf$d_gamma_final <- fit$d_gamma_final %||% NA_real_
  perf$d_alpha_final <- fit$d_alpha_final %||% NA_real_
  perf$max_param_change_final <- fit$max_param_change_final %||% NA_real_
  perf$borrow_adaptive_damping <- fit$adaptive_damping %||% NA
  perf$borrow_update_damping <- fit$update_damping %||% NA_real_
  perf$damping_triggered <- fit$damping_triggered %||% NA
  perf$damping_trigger_iter <- fit$damping_trigger_iter %||% NA_integer_
  perf$damping_iterations <- fit$damping_iterations %||% NA_integer_
  perf$cycle_distance_at_trigger <- fit$cycle_distance_at_trigger %||% NA_real_
  perf$raw_max_param_change_final <- fit$raw_max_param_change_final %||% NA_real_
  perf$stop_by_elbo <- fit$stop_by_elbo %||% NA
  perf$acvb_acceleration <- fit$acvb_acceleration %||% NA_character_
  perf$acvb_mean_iter_time_sec <- fit$acvb_mean_iter_time_sec %||% NA_real_
  perf$acvb_max_iter_time_sec <- fit$acvb_max_iter_time_sec %||% NA_real_
  perf$acvb_min_downdate_denom <- fit$acvb_min_downdate_denom %||% NA_real_
  perf$acvb_n_denom_clamped <- fit$acvb_n_denom_clamped %||% NA_integer_
  perf
}

run_one_elbo <- function(
    seed = 3036,
    n = 300,
    p = 600,
    nzc = 10,
    rho = 0.3,
    snr = 3.0,
    common_bad_frac = 0.2,
    common_bad_count = NULL,

    z_prior = "dimension_mild",
    gamma_prior = "constant_half",
    source_precision_mode = c("high_dimensional_logp", "raw_ns", "target_capped"),
    source_precision_cap_ratio = 1,
    init_method_borrow = "full",
    init_alpha_borrow = "uniform",
    init_gamma_borrow = "theory_zero",
    warm_start_learned_from_fixed = FALSE,
    theory_initialization = TRUE,
    target_init_until_converged = FALSE,
    target_init_steps = 5L,
    target_init_max_iter = 1000L,
    target_init_tol = NULL,
    preliminary_alpha_update = "huber",
    huber_tau = NULL,
    huber_tau_rule = c("fixed", "diverging"),
    huber_tau_multiplier = 1,
    use_candidate_screening = FALSE,
    candidate_size = 220,
    fast_cached_update = TRUE,
    borrow_update_schedule = c("paper_block", "gauss_seidel", "jacobi"),
    warm_start_iter = 20,

    sigma_beta_sq_borrow = 0.1,
    max_iter_borrow = 80,
    tol_borrow = 1e-6,
    tol_elbo = 1e-7,
    stop_by_elbo = FALSE,
    min_iter = 8,
    borrow_adaptive_damping = TRUE,
    borrow_update_damping = 0.5,
    borrow_cycle_detection_tol = 1e-6,
    borrow_cycle_detection_relative_tol = 1e-3,
    borrow_cycle_min_amplitude = NULL,
    alpha_update_every_learned = 1,
    alpha_update_strategy_learned = c("elbo_accept", "em_qp", "old"),
    alpha_update_strategy_old = c("elbo_accept", "old"),
    run_old_alpha = TRUE,
    run_em_alpha = TRUE,
    alpha_freeze_iter_learned = 0,
    alpha_grid_step_learned = 0.02,
    alpha_accept_tol = 1e-8,

    source_prec_ratio = 0.25,
    source_prec_cap = 0.75,
    refit_ridge = 1e-6,
    gamma_ebic = 0.5,

    run_current_alpha_every5 = TRUE,
    run_fixed_alpha = TRUE,
    run_oracle_alpha = TRUE,
    run_acvb = FALSE,
    acvb_acceleration = c("woodbury", "cached_inverse", "slow"),
    run_trans_lasso = FALSE,
    run_oracle_trans_simple = FALSE,
    run_naive_trans_simple = FALSE,
    run_alpha_grid_profile = TRUE,
    alpha_grid_step = 0.02,

    monitor_elbo = TRUE,
    elbo_print_every = 1,
    include_source_normalization = TRUE,
    collapsed_logp_mode = c("selected", "variational"),
    collapsed_logp_threshold = 0.5,
    verbose = TRUE,
    save_dir = NULL,
    ...
) {
  borrow_update_schedule <- match.arg(borrow_update_schedule)
  source_precision_mode <- match.arg(source_precision_mode)
  if (is.null(target_init_tol)) target_init_tol <- tol_borrow
  collapsed_logp_mode <- match.arg(collapsed_logp_mode)
  acvb_acceleration <- match.arg(acvb_acceleration)
  alpha_update_strategy_learned <- match.arg(alpha_update_strategy_learned)
  alpha_update_strategy_old <- match.arg(alpha_update_strategy_old)
  preliminary_alpha_update <- match.arg(
    preliminary_alpha_update,
    choices = c("huber", "target_matching_qp", "none")
  )
  huber_tau_rule <- match.arg(huber_tau_rule)

  if (isTRUE(theory_initialization) && isTRUE(warm_start_learned_from_fixed)) {
    if (isTRUE(verbose)) {
      cat("[run_one_elbo] theory_initialization=TRUE: disabling the old fixed-alpha warm start so gamma^(0)=0 is preserved.\n")
    }
    warm_start_learned_from_fixed <- FALSE
  }

  ## Keep algorithm-only arguments away from make_dat().  This is important
  ## when run_one_elbo() is called by do.call() with a merged list of DGP and
  ## algorithm options.  Only arguments that are formal arguments of make_dat()
  ## are forwarded through dgp_dots.
  dots <- list(...)
  dgp_names <- names(formals(make_dat))
  dgp_dots <- dots[names(dots) %in% dgp_names]
  for (nm in c("seed", "n", "p", "nzc", "rho", "snr",
               "common_bad_frac", "common_bad_count")) {
    dgp_dots[[nm]] <- NULL
  }

  ignored_dots <- setdiff(names(dots), names(dgp_dots))
  ignored_dots <- setdiff(ignored_dots, c("seed", "n", "p", "nzc", "rho", "snr",
                                          "common_bad_frac", "common_bad_count"))
  if (length(ignored_dots) > 0 && isTRUE(verbose)) {
    cat("[run_one_elbo] Ignoring non-DGP arguments in ...: ",
        paste(ignored_dots, collapse = ", "), "\n", sep = "")
  }

  dat <- do.call(
    make_dat,
    c(
      list(
        seed = seed,
        n = n,
        p = p,
        nzc = nzc,
        rho = rho,
        snr = snr,
        common_bad_frac = common_bad_frac,
        common_bad_count = common_bad_count
      ),
      dgp_dots
    )
  )
  dat <- configure_trace_source_precision(
    dat = dat,
    mode = source_precision_mode,
    cap_ratio = source_precision_cap_ratio
  )

  nzc_true <- sum(dat$beta_true != 0)
  prior_par_borrow <- make_prior_params_borrow_commonbad(
    n = dat$n,
    p = dat$p,
    nzc = nzc_true,
    z_prior = z_prior,
    gamma_prior = gamma_prior,
    ns_eff = trace_get_raw_ns_eff(dat),
    min_A_z = 2,
    min_conc_z = 50
  )

  if (!is.null(save_dir) && !dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)

  fit_borrow_diag <- function(..., stage_label) {
    fit_borrow_with_screening(
      ...,
      use_candidate_screening = use_candidate_screening,
      candidate_size = candidate_size,
      fast_cached_update = fast_cached_update,
      borrow_update_schedule = borrow_update_schedule,
      monitor_elbo = monitor_elbo,
      elbo_print_every = elbo_print_every,
      include_source_normalization = include_source_normalization,
      tol_elbo = tol_elbo,
      stop_by_elbo = stop_by_elbo,
      min_iter = min_iter,
      adaptive_damping = borrow_adaptive_damping,
      update_damping = borrow_update_damping,
      cycle_detection_tol = borrow_cycle_detection_tol,
      cycle_detection_relative_tol = borrow_cycle_detection_relative_tol,
      cycle_min_amplitude = borrow_cycle_min_amplitude,
      theory_initialization = theory_initialization,
      target_init_until_converged = target_init_until_converged,
      target_init_steps = target_init_steps,
      target_init_max_iter = target_init_max_iter,
      target_init_tol = target_init_tol,
      preliminary_alpha_update = preliminary_alpha_update,
      huber_tau = huber_tau,
      huber_tau_rule = huber_tau_rule,
      huber_tau_multiplier = huber_tau_multiplier,
      stage_label = stage_label
    )
  }

  fits <- list()
  times <- list()

  cat("\n================ data setting ================\n")
  cat(sprintf("seed=%d | n=%d | p=%d | nzc=%d | rho=%.2f | snr=%.2f | common_bad_count=%d | common_bad_frac=%.3f\n",
              seed, dat$n, dat$p, nzc_true, dat$rho, dat$snr,
              dat$common_bad_count, dat$common_bad_frac))
  cat(sprintf("alpha_true=(%s) | source_scale=(%s)\n",
              paste(round(dat$alpha_true, 4), collapse = ", "),
              paste(round(dat$source_scale, 4), collapse = ", ")))
  cat(sprintf("use_candidate_screening=%s | candidate_size=%s\n",
              use_candidate_screening, candidate_size))
  cat(sprintf(
    "target prior=%s (rho_z=%.8g, uses true nzc=%s) | source precision ns: raw=%.1f, effective=%.1f\n",
    prior_par_borrow$z_prior_rule %||% z_prior,
    prior_par_borrow$rho_z_mean,
    prior_par_borrow$z_prior_uses_true_nzc %||% NA,
    dat$trace_source_precision_ns_raw,
    dat$trace_source_precision_ns_eff
  ))
  cat(sprintf("theory_initialization=%s | init_gamma_effective=%s | preliminary_alpha_update=%s\n",
              theory_initialization,
              if (isTRUE(theory_initialization)) "theory_zero" else paste(init_gamma_borrow, collapse = ";"),
              preliminary_alpha_update))

  ## ---------------------------------------------------------
  ## Fixed uniform alpha.  This is also used as the warm start
  ## for learned-alpha if warm_start_learned_from_fixed = TRUE.
  ## ---------------------------------------------------------
  fit_fixed_warm <- NULL
  time_fixed_warm <- 0

  if (isTRUE(run_fixed_alpha) || isTRUE(warm_start_learned_from_fixed)) {
    cat("\n================ fit: fixed uniform alpha ================\n")
    tmp <- .time_it(fit_borrow_diag(
      dat = dat,
      prior_par_borrow = prior_par_borrow,
      alpha_method = "fixed",
      alpha_fixed = rep(1 / dat$K, dat$K),
      init_alpha = "uniform",
      init_method = init_method_borrow,
      init_gamma = init_gamma_borrow,
      sigma_beta_sq = sigma_beta_sq_borrow,
      max_iter = if (isTRUE(run_fixed_alpha)) max_iter_borrow else max(warm_start_iter, min_iter),
      tol = tol_borrow,
      gamma_ebic = gamma_ebic,
      verbose = verbose,
      alpha_update_every = 1,
      stage_label = "fixed_uniform"
    ))
    fit_fixed_warm <- tmp$value
    time_fixed_warm <- tmp$elapsed
    if (isTRUE(run_fixed_alpha)) {
      fits$borrow_cvb_fixed_uniform <- fit_fixed_warm
      times$borrow_cvb_fixed_uniform <- time_fixed_warm
    }
  }

  ## ---------------------------------------------------------
  ## Learned alpha methods.  We now fit the old Dirichlet-style alpha update
  ## and the EM-QP alpha update as separate methods, because their alpha
  ## objectives are not the same.
  ## ---------------------------------------------------------
  if (isTRUE(warm_start_learned_from_fixed) && !is.null(fit_fixed_warm)) {
    init_method_learned <- "custom"
    init_omega_learned <- fit_fixed_warm$omega
    init_gamma_learned <- fit_fixed_warm$gamma
  } else {
    init_method_learned <- init_method_borrow
    init_omega_learned <- NULL
    init_gamma_learned <- init_gamma_borrow
  }

  if (isTRUE(run_old_alpha)) {
    cat("\n================ fit: old learned alpha with ELBO accept/reject ================\n")
    tmp <- .time_it(fit_borrow_diag(
      dat = dat,
      prior_par_borrow = prior_par_borrow,
      alpha_method = "old",
      init_alpha = init_alpha_borrow,
      init_method = init_method_learned,
      init_omega = init_omega_learned,
      init_gamma = init_gamma_learned,
      sigma_beta_sq = sigma_beta_sq_borrow,
      max_iter = max_iter_borrow,
      tol = tol_borrow,
      gamma_ebic = gamma_ebic,
      verbose = verbose,
      alpha_update_every = alpha_update_every_learned,
      alpha_update_strategy = alpha_update_strategy_old,
      alpha_freeze_iter = alpha_freeze_iter_learned,
      alpha_grid_step = alpha_grid_step_learned,
      alpha_accept_tol = alpha_accept_tol,
      stage_label = paste0("old_alpha_elboacc_every", alpha_update_every_learned)
    ))
    fits$borrow_cvb_old_alpha_elboacc <- tmp$value
    fits$borrow_cvb_old_alpha_elboacc$warm_start_from <- if (warm_start_learned_from_fixed) "fixed_uniform_alpha" else NA_character_
    times$borrow_cvb_old_alpha_elboacc <- tmp$elapsed + if (warm_start_learned_from_fixed) time_fixed_warm else 0
  }

  if (isTRUE(run_em_alpha)) {
    cat("\n================ fit: EM-QP learned alpha ================\n")
    tmp <- .time_it(fit_borrow_diag(
      dat = dat,
      prior_par_borrow = prior_par_borrow,
      alpha_method = "old",
      init_alpha = init_alpha_borrow,
      init_method = init_method_learned,
      init_omega = init_omega_learned,
      init_gamma = init_gamma_learned,
      sigma_beta_sq = sigma_beta_sq_borrow,
      max_iter = max_iter_borrow,
      tol = tol_borrow,
      gamma_ebic = gamma_ebic,
      verbose = verbose,
      alpha_update_every = alpha_update_every_learned,
      alpha_update_strategy = "em_qp",
      alpha_freeze_iter = alpha_freeze_iter_learned,
      alpha_grid_step = alpha_grid_step_learned,
      alpha_accept_tol = alpha_accept_tol,
      stage_label = paste0("em_alpha_every", alpha_update_every_learned)
    ))
    fits$borrow_cvb_em_alpha <- tmp$value
    fits$borrow_cvb_em_alpha$warm_start_from <- if (warm_start_learned_from_fixed) "fixed_uniform_alpha" else NA_character_
    times$borrow_cvb_em_alpha <- tmp$elapsed + if (warm_start_learned_from_fixed) time_fixed_warm else 0
  }

  ## ---------------------------------------------------------
  ## Current update frequency, if you want to compare with the
  ## existing alpha_update_every = 5 behavior.
  ## ---------------------------------------------------------
  if (isTRUE(run_current_alpha_every5) && alpha_update_every_learned != 5) {
    cat("\n================ fit: learned alpha, current every-5 update ================\n")
    tmp <- .time_it(fit_borrow_diag(
      dat = dat,
      prior_par_borrow = prior_par_borrow,
      alpha_method = "old",
      init_alpha = init_alpha_borrow,
      init_method = init_method_borrow,
      init_gamma = init_gamma_borrow,
      sigma_beta_sq = sigma_beta_sq_borrow,
      max_iter = max_iter_borrow,
      tol = tol_borrow,
      gamma_ebic = gamma_ebic,
      verbose = verbose,
      alpha_update_every = 1,
      alpha_update_strategy = "old",
      alpha_freeze_iter = 0,
      alpha_grid_step = alpha_grid_step_learned,
      alpha_accept_tol = alpha_accept_tol,
      stage_label = "learned_alpha_every5_current"
    ))
    fits$borrow_cvb_learn_alpha_every5 <- tmp$value
    times$borrow_cvb_learn_alpha_every5 <- tmp$elapsed
  }

  if (isTRUE(run_oracle_alpha)) {
    cat("\n================ fit: oracle alpha ================\n")
    tmp <- .time_it(fit_borrow_diag(
      dat = dat,
      prior_par_borrow = prior_par_borrow,
      alpha_method = "oracle",
      init_alpha = "oracle",
      init_method = init_method_borrow,
      init_gamma = init_gamma_borrow,
      sigma_beta_sq = sigma_beta_sq_borrow,
      max_iter = max_iter_borrow,
      tol = tol_borrow,
      gamma_ebic = gamma_ebic,
      verbose = verbose,
      alpha_update_every = 1,
      stage_label = "oracle_alpha"
    ))
    fits$borrow_cvb_oracle_alpha <- tmp$value
    times$borrow_cvb_oracle_alpha <- tmp$elapsed
  }

  if (isTRUE(run_acvb)) {
    cat("\n================ fit: no-borrow ACVB ================\n")
    pi_prior_acvb <- prior_par_borrow$rho_z_mean
    tmp <- .time_it(fit_noborrow_acvb(
      dat = dat,
      init_method = init_method_borrow,
      pi_prior = pi_prior_acvb,
      max_iter = max_iter_borrow,
      sigma_beta_sq = sigma_beta_sq_borrow,
      tol = tol_borrow,
      gamma_ebic = gamma_ebic,
      verbose = verbose,
      fast_cached_update = TRUE,
      acceleration = acvb_acceleration
    ))
    fits$no_borrow_acvb <- tmp$value
    times$no_borrow_acvb <- tmp$elapsed
  }

  if (isTRUE(run_trans_lasso)) {
    cat("\n================ fit: Trans-Lasso ================\n")
    tmp <- .time_it(fit_trans_lasso_main(dat))
    fits$trans_lasso <- tmp$value
    times$trans_lasso <- tmp$elapsed
  }

  if (isTRUE(run_oracle_trans_simple)) {
    cat("\n================ fit: oracle simple Trans-Lasso ================\n")
    tmp <- .time_it(fit_oracle_trans_lasso_simple(
      dat,
      A0 = dat$A_trans_oracle,
      method_name = if (length(dat$A_trans_oracle) == 0) {
        "Oracle_Trans_Lasso_simple_emptyA"
      } else {
        paste0("Oracle_Trans_Lasso_simple_A_", paste(dat$A_trans_oracle, collapse = "_"))
      }
    ))
    fits$oracle_trans_lasso_simple <- tmp$value
    times$oracle_trans_lasso_simple <- tmp$elapsed
  }

  if (isTRUE(run_naive_trans_simple)) {
    cat("\n================ fit: naive simple Trans-Lasso ================\n")
    tmp <- .time_it(fit_oracle_trans_lasso_simple(
      dat,
      A0 = seq_len(dat$K),
      method_name = "Naive_Trans_Lasso_simple_all_sources"
    ))
    fits$naive_trans_lasso_simple <- tmp$value
    times$naive_trans_lasso_simple <- tmp$elapsed
  }

  performance <- do.call(rbind, lapply(names(fits), function(nm) {
    row <- pack_fit_row(
      nm,
      fits[[nm]],
      dat,
      elapsed = times[[nm]] %||% NA_real_,
      prior_par_borrow = prior_par_borrow,
      include_source_normalization = include_source_normalization,
      source_prec_ratio = source_prec_ratio,
      source_prec_cap = source_prec_cap,
      refit_ridge = refit_ridge,
      collapsed_logp_mode = collapsed_logp_mode,
      collapsed_logp_threshold = collapsed_logp_threshold
    )
    .add_elbo_to_perf(row, fits[[nm]])
  }))
  performance <- .trace_apply_public_method_labels(performance)

  trace_fit_for_diag <- fits$borrow_cvb_em_alpha %||% fits$borrow_cvb_learn_alpha
  diagnostic_index <- sort(unique(c(which(dat$tag_true == 1), dat$spurious_union)))
  source_center_true <- as.numeric(dat$W_source_true %*% dat$alpha_true)
  coordinate_diagnostics <- data.frame(
    j = diagnostic_index,
    coordinate_type = ifelse(
      diagnostic_index %in% dat$common_bad,
      "active_nontransferable",
      ifelse(diagnostic_index %in% dat$spurious_union,
             "source_only_null", "transferable_active")
    ),
    beta_true = dat$beta_true[diagnostic_index],
    gamma_true = dat$gamma_true[diagnostic_index],
    source_center_true = source_center_true[diagnostic_index],
    source_center_gap = source_center_true[diagnostic_index] - dat$beta_true[diagnostic_index],
    omega_TRACE = if (!is.null(trace_fit_for_diag$omega)) trace_fit_for_diag$omega[diagnostic_index] else NA_real_,
    gamma_TRACE = if (!is.null(trace_fit_for_diag$gamma)) trace_fit_for_diag$gamma[diagnostic_index] else NA_real_,
    stringsAsFactors = FALSE
  )

  scenario <- data.frame(
    scenario_id = dat$scenario_id,
    seed = seed,
    n = dat$n,
    p = dat$p,
    nzc = nzc_true,
    rho = dat$rho,
    snr = dat$snr,
    noise_mode = dat$noise_mode %||% NA_character_,
    sigma_true = dat$sigma_true %||% NA_real_,
    snr_requested = dat$snr_requested %||% NA_real_,
    realized_snr = dat$realized_snr %||% dat$snr,
    signal_profile = dat$signal_profile %||% NA_character_,
    theory_signal_magnitude = dat$theory_signal_magnitude %||% NA_real_,
    common_bad_signal_magnitude = dat$common_bad_signal_magnitude %||% NA_real_,
    transferable_signal_magnitude = dat$transferable_signal_magnitude %||% NA_real_,
    beta_min = dat$beta_min %||% NA_real_,
    beta_max = dat$beta_max %||% NA_real_,
    beta_min_information = dat$beta_min_information %||% NA_real_,
    beta_min_log_benchmark = dat$beta_min_log_benchmark %||% NA_real_,
    beta_min_margin_ratio = dat$beta_min_margin_ratio %||% NA_real_,
    common_bad_gap_population = dat$common_bad_gap_population %||% NA_real_,
    transferable_count = dat$transferable_count %||% length(dat$regular_signal),
    source_weight_identification_eigenvalue =
      dat$source_weight_identification_eigenvalue %||% NA_real_,
    source_weight_identification_pass =
      dat$source_weight_identification_pass %||% NA,
    exact_combination_max_error =
      dat$exact_combination_max_error %||% NA_real_,
    shared_spurious_count =
      dat$shared_spurious_count %||% length(dat$spurious_union),
    shared_spurious_gap_population =
      dat$shared_spurious_gap_population %||% NA_real_,
    common_bad_count = dat$common_bad_count,
    common_bad_frac = dat$common_bad_frac,
    K = dat$K,
    alpha_true = paste(round(dat$alpha_true, 4), collapse = ";"),
    source_scale = paste(round(dat$source_scale, 4), collapse = ";"),
    source_target_population_correlation = paste(
      round(dat$source_target_population_correlation %||% NA_real_, 6),
      collapse = ";"
    ),
    dgp_variant = dat$dgp_variant %||% NA_character_,
    common_bad_mode = dat$common_bad_mode %||% NA_character_,
    independent_contamination_factors =
      dat$independent_contamination_factors %||% NA,
    common_bad_scale = paste(dat$common_bad_scale %||% NA_real_, collapse = ";"),
    spurious_union_count = length(dat$spurious_union),
    common_spurious_count = dat$common_spurious_count %||% 0L,
    support_swap_spur_ratio = dat$support_swap_spur_ratio %||% NA_real_,
    support_swap_size = dat$support_swap_size %||% NA_real_,
    target_model_config_aligned = TRUE,
    target_z_prior = prior_par_borrow$z_prior %||% z_prior,
    target_z_prior_rule = prior_par_borrow$z_prior_rule %||% z_prior,
    target_z_prior_uses_true_nzc = prior_par_borrow$z_prior_uses_true_nzc %||% NA,
    target_pi_prior = prior_par_borrow$rho_z_mean %||% NA_real_,
    target_init_method = init_method_borrow,
    target_candidate_rule = if (isTRUE(use_candidate_screening)) "screened" else "full_p",
    trace_linear_algebra_backend = if (isTRUE(fast_cached_update) &&
                                         borrow_update_schedule != "gauss_seidel") {
      "exact_woodbury_n_by_n"
    } else if (isTRUE(fast_cached_update)) {
      "cached_rank_updates_p_by_p"
    } else {
      "dense_leave_one"
    },
    source_precision_mode = dat$trace_source_precision_mode,
    source_precision_cap_ratio = dat$trace_source_precision_cap_ratio,
    source_precision_ns_raw = dat$trace_source_precision_ns_raw,
    source_precision_ns_eff = dat$trace_source_precision_ns_eff,
    source_precision_power = dat$trace_source_precision_power,
    theory_initialization = theory_initialization,
      target_init_until_converged = target_init_until_converged,
      target_init_steps = target_init_steps,
      target_init_max_iter = target_init_max_iter,
      target_init_tol = target_init_tol,
    preliminary_alpha_update = preliminary_alpha_update,
    huber_tau = if (is.null(huber_tau)) NA_real_ else huber_tau,
    huber_tau_rule = huber_tau_rule,
    huber_tau_multiplier = huber_tau_multiplier,
    init_gamma_borrow = paste(init_gamma_borrow, collapse = ";"),
    borrow_update_schedule = borrow_update_schedule,
    use_candidate_screening = use_candidate_screening,
    candidate_size = candidate_size,
    sigma_beta_sq_borrow = sigma_beta_sq_borrow,
    max_iter_borrow = max_iter_borrow,
    tol_borrow = tol_borrow,
    tol_elbo = tol_elbo,
    stop_by_elbo = stop_by_elbo,
    include_source_normalization = include_source_normalization,
    run_old_alpha = run_old_alpha,
    run_em_alpha = run_em_alpha,
    alpha_update_strategy_learned = alpha_update_strategy_learned,
    alpha_update_strategy_old = alpha_update_strategy_old,
    run_acvb = run_acvb,
    acvb_acceleration = acvb_acceleration,
    run_trans_lasso = run_trans_lasso,
    run_oracle_trans_simple = run_oracle_trans_simple,
    run_naive_trans_simple = run_naive_trans_simple,
    stringsAsFactors = FALSE
  )

  ## Save traces and final tables.
  alpha_profiles <- list()
  if (!is.null(save_dir)) {
    write.csv(scenario, file.path(save_dir, "scenario.csv"), row.names = FALSE)
    write.csv(performance, file.path(save_dir, "performance.csv"), row.names = FALSE)
    write.csv(coordinate_diagnostics, file.path(save_dir, "coordinate_diagnostics.csv"), row.names = FALSE)

    for (nm in names(fits)) {
      if (!is.null(fits[[nm]]$elbo_trace) && nrow(fits[[nm]]$elbo_trace) > 0) {
        tr <- fits[[nm]]$elbo_trace
        tr$method <- nm
        write.csv(tr, file.path(save_dir, paste0("elbo_trace_", nm, ".csv")), row.names = FALSE)
      }
      if (!is.null(fits[[nm]]$theory_init_trace) && nrow(fits[[nm]]$theory_init_trace) > 0) {
        tr0 <- fits[[nm]]$theory_init_trace
        tr0$method <- nm
        write.csv(tr0, file.path(save_dir, paste0("theory_init_trace_", nm, ".csv")), row.names = FALSE)
      }
    }
  }

  if (isTRUE(run_alpha_grid_profile) && "borrow_cvb_em_alpha" %in% names(fits)) {
    cat("\n================ alpha grid profile at learned-alpha final omega/gamma ================\n")
    prof <- profile_alpha_grid_borrow(
      fit = fits$borrow_cvb_em_alpha,
      dat = dat,
      prior_par_borrow = prior_par_borrow,
      sigma_beta_sq = sigma_beta_sq_borrow,
      grid_step = alpha_grid_step,
      include_source_normalization = include_source_normalization,
      top_n = 20
    )
    alpha_profiles$em_alpha_final_state <- prof
    print(prof)
    cat(sprintf("current alpha=(%s), current ELBO=%.6e\n",
                paste(round(attr(prof, "current_alpha"), 4), collapse = ", "),
                attr(prof, "current_elbo")))
    cat(sprintf("best grid alpha=(%s), best grid ELBO=%.6e, gap=%.6e\n",
                paste(round(attr(prof, "best_alpha"), 4), collapse = ", "),
                attr(prof, "best_elbo"),
                attr(prof, "best_elbo") - attr(prof, "current_elbo")))
    if (!is.null(save_dir)) {
      write.csv(prof, file.path(save_dir, "alpha_grid_profile_top.csv"), row.names = FALSE)
      meta <- data.frame(
        current_alpha = paste(round(attr(prof, "current_alpha"), 6), collapse = ";"),
        current_elbo = attr(prof, "current_elbo"),
        best_alpha = paste(round(attr(prof, "best_alpha"), 6), collapse = ";"),
        best_elbo = attr(prof, "best_elbo"),
        elbo_gap_best_minus_current = attr(prof, "best_elbo") - attr(prof, "current_elbo")
      )
      write.csv(meta, file.path(save_dir, "alpha_grid_profile_meta.csv"), row.names = FALSE)
    }
  }

  cat("\n================ performance summary ================\n")
  show_cols <- intersect(
    c("method", "pred_mse", "tp", "fp", "fn", "f1", "model_size",
      "mean_omega_signal", "mean_omega_null", "mean_gamma_common_bad",
      "mean_gamma_transferable_signal", "alpha_hat", "alpha_l1_error",
      "final_elbo", "elbo_drop_count", "converged", "convergence_reason",
      "time_sec", "n_iter", "alpha_update_every", "alpha_update_strategy", "alpha_freeze_iter",
      "alpha_update_accept_count", "alpha_update_reject_count", "last_alpha_step_size",
      "last_alpha_update_qp_old", "last_alpha_update_qp_new",
      "last_alpha_active_set_size", "last_alpha_active_set",
      "theory_initialization", "preliminary_alpha_update", "preliminary_alpha_done",
      "preliminary_alpha_qp_old", "preliminary_alpha_qp_new", "preliminary_alpha_active_set",
      "preliminary_alpha_solver"),
    names(performance)
  )
  print(performance[, show_cols, drop = FALSE])

  list(
    scenario = scenario,
    dat = dat,
    prior_par_borrow = prior_par_borrow,
    fits = fits,
    performance = performance,
    coordinate_diagnostics = coordinate_diagnostics,
    alpha_profiles = alpha_profiles
  )
}

run_small_elbo_experiment <- function(nzc_values = c(5, 10),
                                      base_seed = 3036,
                                      save_root = "results_elbo_small",
                                      n = 300,
                                      p = 600,
                                      ...) {
  ## Keep the diagnostic design fixed, as requested.
  if (!identical(as.integer(n), 300L)) stop("This diagnostic script fixes n = 300.")
  if (!identical(as.integer(p), 600L)) stop("This diagnostic script fixes p = 600.")
  if (!all(as.integer(nzc_values) %in% c(5L, 10L))) {
    stop("This diagnostic script only allows nzc in c(5, 10).")
  }

  ## Arguments forwarded through ... may accidentally contain n, p, nzc,
  ## seed, or save_dir. Remove them here because this wrapper supplies those
  ## arguments explicitly to run_one_elbo(). Otherwise R throws:
  ## formal argument "n" matched by multiple actual arguments.
  dots <- list(...)
  for (nm in c("n", "p", "nzc", "seed", "save_dir")) dots[[nm]] <- NULL

  if (!dir.exists(save_root)) dir.create(save_root, recursive = TRUE)
  all_perf <- list()
  all_scen <- list()
  outputs <- list()

  for (ii in seq_along(nzc_values)) {
    nzc_i <- as.integer(nzc_values[ii])
    seed_i <- base_seed + ii - 1L
    save_dir_i <- file.path(save_root, paste0("nzc", nzc_i, "_seed", seed_i))
    cat("

########################################################
")
    cat(sprintf("### Small ELBO diagnostic: nzc=%d, seed=%d
", nzc_i, seed_i))
    cat("########################################################
")

    ans <- do.call(
      run_one_elbo,
      c(
        list(
          seed = seed_i,
          n = as.integer(n),
          p = as.integer(p),
          nzc = nzc_i,
          save_dir = save_dir_i
        ),
        dots
      )
    )

    perf_i <- ans$performance
    perf_i$nzc_setting <- nzc_i
    perf_i$seed <- seed_i
    scen_i <- ans$scenario
    scen_i$nzc_setting <- nzc_i
    scen_i$seed <- seed_i

    all_perf[[ii]] <- perf_i
    all_scen[[ii]] <- scen_i
    outputs[[paste0("nzc", nzc_i)]] <- ans
  }

  performance_all <- do.call(rbind, all_perf)
  scenarios_all <- do.call(rbind, all_scen)
  write.csv(performance_all, file.path(save_root, "performance_all.csv"), row.names = FALSE)
  write.csv(scenarios_all, file.path(save_root, "scenarios_all.csv"), row.names = FALSE)

  cat("\n================ combined small-experiment summary ================\n")
  show_cols <- intersect(
    c("nzc_setting", "seed", "method", "pred_mse", "f1", "model_size",
      "alpha_hat", "alpha_l1_error", "final_elbo", "elbo_drop_count",
      "converged", "convergence_reason", "n_iter"),
    names(performance_all)
  )
  print(performance_all[, show_cols, drop = FALSE])

  invisible(list(
    performance = performance_all,
    scenarios = scenarios_all,
    outputs = outputs
  ))
}


## =========================================================
## TRADER baseline extension for the current Borrow-CVB DGP
## =========================================================
source(file.path(.this_dir, "trader_bridge.R"))

.run_one_elbo_without_trader <- run_one_elbo

run_one_elbo <- function(...,
                         run_trader = TRUE,
                         trader_source_input = c("B", "raw"),
                         trader_estimator = c("spikeslab", "glmnet"),
                         trader_allow_estimator_fallback = FALSE,
                         trader_corr_threshold = 0.1,
                         trader_max_var = NULL,
                         trader_n_mcmc = 10000,
                         trader_n_burnin = 5000,
                         trader_thin = 2,
                         trader_n_chain = 1,
                         trader_fold = 1.1,
                         trader_sign_threshold = 0.95,
                         trader_keep_samples = FALSE,
                         trader_verbose = TRUE,
                         trader_seed = NULL,
                         save_dir = NULL) {
  trader_source_input <- match.arg(trader_source_input)
  trader_estimator <- match.arg(trader_estimator)
  dots <- list(...)

  ans <- do.call(
    .run_one_elbo_without_trader,
    c(dots, list(save_dir = save_dir))
  )

  if (!isTRUE(run_trader)) {
    ans$scenario$run_trader <- FALSE
    if (!is.null(save_dir)) {
      write.csv(ans$scenario, file.path(save_dir, "scenario.csv"), row.names = FALSE)
    }
    return(ans)
  }

  seed_base <- dots$seed %||% 123L
  if (is.null(trader_seed)) trader_seed <- as.integer(seed_base) + 700000L

  cat("\n================ fit: TRADER baseline ================\n")
  fit_trader <- fit_trader_baseline(
    dat = ans$dat,
    source_input = trader_source_input,
    estimator = trader_estimator,
    profile = "current",
    corr_threshold = trader_corr_threshold,
    max_var = trader_max_var,
    allow_estimator_fallback = trader_allow_estimator_fallback,
    n_mcmc = trader_n_mcmc,
    n_burnin = trader_n_burnin,
    thin = trader_thin,
    n_chain = trader_n_chain,
    fold = trader_fold,
    sign_threshold = trader_sign_threshold,
    seed = trader_seed,
    keep_samples = trader_keep_samples,
    verbose = trader_verbose
  )

  trader_row <- make_trader_performance_row(
    fit = fit_trader,
    dat = ans$dat,
    template = ans$performance
  )
  ans$performance <- .bind_rows_fill(ans$performance, trader_row)
  ans$performance <- .trace_apply_public_method_labels(ans$performance)
  ans$fits$trader <- fit_trader

  ans$scenario$run_trader <- TRUE
  ans$scenario$trader_source_input <- trader_source_input
  ans$scenario$trader_estimator <- trader_estimator
  ans$scenario$trader_allow_estimator_fallback <- trader_allow_estimator_fallback
  ans$scenario$trader_corr_threshold <- trader_corr_threshold
  ans$scenario$trader_max_var <- fit_trader$max_var_used
  ans$scenario$trader_max_var_rule <- fit_trader$max_var_rule
  ans$scenario$trader_n_mcmc <- trader_n_mcmc
  ans$scenario$trader_n_burnin <- trader_n_burnin
  ans$scenario$trader_thin <- trader_thin
  ans$scenario$trader_sign_threshold <- trader_sign_threshold

  if (!is.null(save_dir)) {
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    write.csv(ans$performance, file.path(save_dir, "performance.csv"), row.names = FALSE)
    write.csv(ans$scenario, file.path(save_dir, "scenario.csv"), row.names = FALSE)
    trader_diag <- data.frame(
      kept_sources = paste(fit_trader$kept_sources, collapse = ";"),
      dropped_sources = paste(fit_trader$dropped_sources, collapse = ";"),
      source_correlations = paste(round(fit_trader$source_correlations, 8), collapse = ";"),
      eta_mean = paste(round(fit_trader$eta_mean, 8), collapse = ";"),
      source_input = fit_trader$source_input,
      target_estimator = fit_trader$target_estimator,
      source_estimators = paste(fit_trader$source_estimators, collapse = ";"),
      max_var_used = fit_trader$max_var_used,
      max_var_rule = fit_trader$max_var_rule,
      preprocess_time_sec = fit_trader$preprocess_time_sec,
      mcmc_time_sec = fit_trader$mcmc_time_sec,
      total_time_sec = fit_trader$total_time_sec,
      max_psrf = fit_trader$max_psrf,
      mean_psrf = fit_trader$mean_psrf,
      stringsAsFactors = FALSE
    )
    write.csv(trader_diag, file.path(save_dir, "trader_diagnostics.csv"), row.names = FALSE)
    if (isTRUE(trader_keep_samples) && !is.null(fit_trader$beta_samples)) {
      saveRDS(fit_trader$beta_samples, file.path(save_dir, "trader_beta_samples.rds"))
    }
  }

  ans
}

source(file.path(.this_dir, "trace_theory_diagnostics.R"), local = .GlobalEnv)
