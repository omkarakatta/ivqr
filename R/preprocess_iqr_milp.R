### Meta -------------------------
###
### Title: Compute IQR estimator by preprocessing MILP
###
### Description: The IQR estimator is computed by solving a MILP.
### This function pre-processes the data and fixes the sign of the residuals of
### outliers to speed up the procedure.
###
### Author: Omkar A. Katta
###

### preprocess_iqr_milp -------------------------
#' Compute IQR estimator by preprocessing MILP
#'
#' Fix the sign of outliers' residuals to solve the IQR MILP more quickly
#'
#' @param Y Dependent variable (vector of length n)
#' @param X Exogenous variable (n by p_X matrix)
#' @param D Endogenous variable (n by p_D matrix)
#' @param Z Instrumental variable (n by p_Z matrix)
#' @param tau Quantile (number strictly between 0 and 1)
#' @param prop_alpha_initial Proportion of residuals that are \emph{not} fixed
#'  at the start;thus, 1 - \code{alpha_initial} is the proportion of
#'  observations whose residuals are indeed fixed at the start;
#'
#'
preprocess_iqr_milp <- function(Y,
                                D,
                                X,
                                Z,
                                tau,
                                prop_alpha_initial = 0.7) {

  # Get dimensions of data
  n <- length(Y)
  n_D <- nrow(D)
  n_X <- nrow(X)
  n_Z <- nrow(Z)
  p_D <- ncol(D)
  p_X <- ncol(X)
  p_Z <- ncol(Z)

  # Determine preliminary residuals
  resid <- quantreg::rq(Y ~ X + D - 1, tau = tau)$residuals

  # Determine initial residual bounds
  alpha_initial <- quantile(abs(resid), prop_alpha_initial)

  # Start the while loop
  alpha <- width_initial
  status <- "TIME_LIMIT"

  # Continue the while loop if the program took too long to solve or if
  # the objective (i.e., absolute value of beta_Z) is not 0. Note that
  # if the program is infeasible, i.e., the objective was NULL, we also
  # continue the while loop (we mechanically set the objective to be nonzero
  # later in the code).
  while (status == "TIME_LIMIT" | obj != 0) {
    # Fix the most negative and most positive residuals
    O_neg <- which(Y < -1 * alpha)
    O_pos <- which(Y > alpha)
  }
}

