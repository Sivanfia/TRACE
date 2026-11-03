## =========================================================
## Unified prediction metric for 01_setting
## =========================================================
## All methods use the observed-response test error
##
##   pred_mse = mean((y_test - X_test %*% beta_hat)^2)
##
## This includes the irreducible test noise.  The previous core helper could
## calculate the noiseless signal error, which made Borrow-CVB/ACVB/Trans-Lasso
## incomparable with TRADER.
## =========================================================

PRED_MSE_DEFINITION_V9 <-
  "mean((y_test - X_test %*% beta_hat)^2)"

PRED_MSE_TARGET_V9 <- "observed_test_response"

response_prediction_mse <- function(beta_hat, dat) {
  if (is.null(beta_hat) || anyNA(beta_hat)) {
    return(NA_real_)
  }

  beta_hat <- as.numeric(beta_hat)

  if (is.null(dat$X_test) || is.null(dat$y_test)) {
    stop(
      "Unified pred_mse requires both dat$X_test and dat$y_test. ",
      "The training error is not used as a fallback."
    )
  }

  X_test <- as.matrix(dat$X_test)
  y_test <- as.numeric(dat$y_test)

  if (ncol(X_test) != length(beta_hat)) {
    stop(
      "Dimension mismatch in pred_mse: ncol(X_test)=", ncol(X_test),
      " but length(beta_hat)=", length(beta_hat), "."
    )
  }

  if (nrow(X_test) != length(y_test)) {
    stop(
      "Dimension mismatch in pred_mse: nrow(X_test)=", nrow(X_test),
      " but length(y_test)=", length(y_test), "."
    )
  }

  pred <- as.numeric(X_test %*% beta_hat)
  value <- mean((y_test - pred)^2)

  if (!is.finite(value)) {
    stop("Unified pred_mse returned a non-finite value.")
  }

  unname(value)
}

## Backward-compatible name used throughout functions.R and the shared core.
## R resolves this global binding at evaluation time, so this definition
## overrides the older signal-prediction helper after borrow_functions_core.R
## has been sourced.
prediction_mse <- response_prediction_mse

signal_prediction_mse_v9 <- function(beta_hat, dat) {
  if (is.null(beta_hat) || anyNA(beta_hat)) return(NA_real_)
  beta_hat <- as.numeric(beta_hat)
  if (is.null(dat$X_test) || is.null(dat$beta_true)) return(NA_real_)
  mean(
    (
      as.numeric(dat$X_test %*% beta_hat) -
        as.numeric(dat$X_test %*% dat$beta_true)
    )^2
  )
}
