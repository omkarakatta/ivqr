### Meta -------------------------
###
### Title: Create simulation according to Chen and Lee (2017)
###
### Description: Chen and Lee (2017) describe a simulation study in Section 4 of
### their paper. This function constructs the same simulation, and also extends
### the simulation from the original case of 3 endogenous variables to the more
### general case of at least 3 endogeneous variables.
###
### Author: Omkar A. Katta
###

### chen_lee -------------------------
#' Simulate data according to Chen and Lee (2017)
#'
#' Given the number of observations and number of endogeneous variables, create
#' an outcome variable defined by a location scale model where the coefficients
#' on the endogenous variables are supposed to be 1.
#'
#' The error term in the location scale model that underpins this simulation
#' design is defined in terms of the endogeneous variables (which is why we
#' call these variables "endogenous"). To properly estimate the coefficients on
#' the endogeneous variables, we require instruments that are uncorrelated with
#' the errors, related to the endogeneous variables, and only related to the
#' outcome variable through their association with these endogeneous variables.
#'
#' This function creates errors, endogeneous variables, instruments, and an
#' outcome variable such that the above terms are satisfied.
#' The errors are drawn independently of the instruments from a multivariate
#' normal distribution. The instruments are drawn from a standard normal
#' normal distribution. The endogeneous variables are multiples of the
#' cumulative distribution function of the shocked instruments.
#' The error in the true model for the outcome variable is defined in terms
#' of the endogeneous variables.
#'
#' The original Chen and Lee simulation design used 3 endogeneous variables.
#' This design allows for an arbitrary number of endogeneous variables.
#' To allow fewer endogeneous variables, say 2 endogeneous variables, we simply
#' omit the third endogeneous variable from the original Chen and Lee
#' simulation before constructing our outcome variable.
#'
#' @param n Number of observations; defaults to 500 (numeric)
#' @param p_D Number of endogeneous variables; defaults to 3
#'  (numeric)
#'
#' @return A named list:
#'  \enumerate{
#'    \item Y: outcome variable (n by 1 matrix)
#'    \item D: endogeneous variable (n by p_D matrix)
#'    \item Z: instruments (n by p_D matrix)
#'    \item X: matrix of 1's (n by 1 matrix)
#'    \item errors: matrix of errors and shocks (n by (p_D + 1) matrix); first
#'      column is the vector of errors on the location scale model; all other
#'      columns are shocks to the instruments when defining D.
#'  }
chen_lee <- function(n = 500, p_D = 3) {

  # If p_D is less than 3, we will first create a Chen and Lee simulation with
  # p_D = 3 and then remove the extra endogeneous variables and extra errors in
  # Steps 3, 4, and 5.
  actual_p <- p_D
  if (p_D < 3) {
    p_D <- 3 # actual_p != p_D iff p_D < 3
    warning("p_D is less than 3")
  }

  msg <- "n must be a positive integer."
  send_note_if(msg, round(n) != n || n <= 0, stop, call. = FALSE)

  # Step 1: Create shocks and errors.
  # Each instrument will be associated with a shock (see second equation of (16)
  # in the Chen and Lee (2017) paper), and there will also be an error in the
  # the scale model (see first equation of (16) in the Chen and Lee (2017)
  # paper).
  # Hence, we need to draw from n vector-valued observations from a
  # (p_D + 1)-dimensional multivariate normal distribution.
  # We begin by computing the var-cov matrix for these shocks and errors.
  V <- diag(p_D + 1)
  tmp <- c(0.4, 0.6, -0.2, seq(-1, 1, length.out = p_D - 3))
  new_values <- sqrt(0.4^2 + 0.6^2 + 0.2^2) * tmp / norm(tmp, "2")
  V[1, 2:ncol(V)] <- new_values
  V[2:ncol(V), 1] <- new_values
  # Note: when p_D = 3, new_values is the same as tmp => V is same as first
  # display on page 557 of Chen and Lee (2017).

  # Step 2: Draw errors from multivariate normal.
  errors <- MASS::mvrnorm(n, mu = rep(0, nrow(V)), Sigma = 0.25 *  V)
  eps <- errors[, 1] # first column is error in the scale model
  nu <- errors[, -1] # all other columns are shocks to instruments to define D

  # Step 3: Define Z (instruments).
  Z <- matrix(stats::rnorm(n * p_D), ncol = p_D) # n by p_D matrix of instruments
  # note: I remove extra instruments in Step 4 after I define D

  # Step 4: Define D (endogeneous variables).
  # D is the product of some coefficient and the CDF of shocked instrument.
  coef_D <- c(1, 2.5, 1.5, seq(1, 2, length.out = p_D - 3))
  scale_D <- matrix(rep(coef_D, n), ncol = p_D, byrow = T)
  D <- stats::pnorm(Z + nu) * scale_D # element-wise product, not matrix multiplication
  D <- D[, seq_len(actual_p), drop = FALSE] # remove extra if p_D < 3
  Z <- Z[, seq_len(actual_p), drop = FALSE] # remove extra if p_D < 3

  # Step 5: Define Y (outcome variable)
  beta_D <- rep(1, actual_p) # coefficients on endogeneous variables!
  beta_D_errors <- c(1, 0.25, 0.15, seq(0.1, 1, length.out = p_D - 3))
  beta_D_errors <- beta_D_errors[seq_len(actual_p)] # remove extra if p_D < 3
  X <- matrix(rep(1, n), ncol = 1)
  Y <- X + D %*% beta_D + (0.5 + D %*% beta_D_errors) * eps

  list("Y" = Y, "D" = D, "Z" = Z, "X" = X, "errors" = errors)
}
