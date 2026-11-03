## Prior constructors used by TRACE and ACVB. These are algorithm helpers,
## independent of the scenario DGP, and must be loaded before run_one_elbo().

.safe_prob <- function(x, low = 1e-6, high = 1 - 1e-6) {
  pmin(pmax(x, low), high)
}

.make_z_prior_params_alpha_commonbad <- function(
    n, p, nzc = 8,
    z_prior = c("dimension_mild", "dimension_mild_dense", "high_dimensional", "truth", "adaptive",
                "anti_collapse", "n_inverse"),
    min_A_z = 2,
    min_conc_z = 50) {
  z_prior <- match.arg(z_prior)
  truth_sparsity <- nzc / p

  if (z_prior == "high_dimensional") {
    rho_z_mean <- 1 / (p * max(p^(3 / 2), sqrt(n)))
  } else if (z_prior == "dimension_mild") {
    rho_z_mean <- min(0.10, max(1 / max(n, 2), 1 / sqrt(max(p, 2))))
  } else if (z_prior == "dimension_mild_dense") {
    rho_z_mean <- min(
      0.10,
      max(1 / max(n, 2), 1.5 / sqrt(max(p, 2)))
    )
  }else if (z_prior == "truth") {
    rho_z_mean <- truth_sparsity
  } else if (z_prior == "n_inverse") {
    rho_z_mean <- 1 / max(n, 2)
  } else if (z_prior == "anti_collapse") {
    rho_z_mean <- min(0.30, max(4.0 * truth_sparsity, 1 / n))
  } else {
    if (n > p) {
      rho_z_mean <- min(0.20, max(truth_sparsity, 1 / n))
    } else {
      rho_z_mean <- min(0.12, max(2.0 * truth_sparsity, 1 / n))
    }
  }

  if (z_prior != "high_dimensional") rho_z_mean <- .safe_prob(rho_z_mean)
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
      "min(0.20, max(1/n, 1.5/sqrt(p)))"
    }else {
      z_prior
    },
    z_prior_uses_true_nzc = z_prior %in% c("truth", "adaptive", "anti_collapse")
  )
}

make_prior_params_borrow_alpha_commonbad <- function(
    n, p, nzc = 8,
    z_prior = c("dimension_mild","dimension_mild_dense",  "high_dimensional", "truth", "adaptive",
                "anti_collapse", "n_inverse"),
    gamma_prior = c("constant_half", "ns_inverse_quarter", "n_inverse",
                    "moderate", "strong", "weak"),
    ns_eff = NULL,
    min_A_z = 2,
    min_conc_z = 50) {
  z_prior <- match.arg(z_prior)
  gamma_prior <- match.arg(gamma_prior)
  zpar <- .make_z_prior_params_alpha_commonbad(
    n = n, p = p, nzc = nzc, z_prior = z_prior,
    min_A_z = min_A_z, min_conc_z = min_conc_z
  )

  if (gamma_prior == "constant_half") {
    rho_c_mean <- 0.5
    conc_c <- 50
  } else if (gamma_prior == "ns_inverse_quarter") {
    if (is.null(ns_eff)) ns_eff <- n
    rho_c_mean <- .safe_prob(as.numeric(ns_eff)^(-1 / 4))
    conc_c <- max(50, 2 / rho_c_mean, 2 / (1 - rho_c_mean))
  } else if (gamma_prior == "n_inverse") {
    rho_c_mean <- .safe_prob(1 - 1 / max(n, 2))
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
  c(zpar, list(
    A_c = A_c,
    B_c = B_c,
    rho_c_mean = rho_c_mean,
    conc_c = conc_c,
    E_logit_rho_c = qlogis(rho_c_mean),
    gamma_prior = gamma_prior
  ))
}

make_prior_params_vbss_alpha_commonbad <- function(
    n, p, nzc = 8,
    z_prior = c("dimension_mild", "high_dimensional", "dimension_mild_dense", "truth", "adaptive",
                "anti_collapse", "n_inverse"),
    min_A_z = 2,
    min_conc_z = 50) {
  z_prior <- match.arg(z_prior)
  .make_z_prior_params_alpha_commonbad(
    n = n, p = p, nzc = nzc, z_prior = z_prior,
    min_A_z = min_A_z, min_conc_z = min_conc_z
  )
}

make_prior_params_borrow_commonbad <- make_prior_params_borrow_alpha_commonbad
make_prior_params_vbss_commonbad <- make_prior_params_vbss_alpha_commonbad
