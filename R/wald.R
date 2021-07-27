### Title: Create Wald-type confidence interval
###
### Description: To approximate the confidence interval from line_confint, we
### consider a Wald-type confidence interval (see (3.11) from Chernozhukov and
### Hansen (2006)). This confidence interval may give us a good guess for where
### the bounds of the true confidence interval may lie.
###
### Author: Omkar A. Katta

### wald_univariate -------------------------
#' Compute Wald-type univariate confidence interval
#'
#' See (3.11) from Chernozhukov and Hansen (2006).
#'
#' @param index Index of the coefficient of interest
#'  (numeric between 1 and p_D if \code{endogeneous} is TRUE;
#'  numeric between 1 and p_X if \code{endogeneous} is FALSE)
#' @param endogeneous If TRUE (default), \code{index} refers to the
#'  coefficients on the endogeneous variables; if FALSE, \code{index} refers to
#'  the coefficients on the exogeneous variables (boolean)
#' @param center MILP point estimate that is the center of the univariate
#'  confidence interval (scalar)
#' @param resid MILP residuals associated with point estimate
#' @param alpha Alpha level of the test; defaults to 0.1; used to compute
#'  critical value and the Hall and Sheather bandwidth (for estimating
#'  J_{Psi}(tau) in CH (2006))
#' @param tau Quantile (number between 0 and 1)
#' @param D Matrix of endogeneous variables
#' @param X Matrix of covariates (including intercept)
#' @param Z Matrix of instruments (only relevant for obtaining \code{Phi})
#' @param Phi Transformed matrix of instruments; defaults to linear projection
#'  of D on X and Z
wald_univariate <- function(center,
                            endogeneous,
                            index,
                            resid,
                            alpha = 0.1,
                            tau,
                            D,
                            X,
                            Z,
                            Phi = linear_projection(D, X, Z)) {

  start_clock <- Sys.time()
  out <- list() # initialize list of results to return

  # Get dimensions of data
  n <- nrow(D)
  p_D <- ncol(D)
  p_X <- ncol(X)
  stopifnot(nrow(X) == n)
  stopifnot(nrow(Phi) == n)

  # critical value
  crit_val <- stats::qnorm(1 - alpha / 2)
  out$crit_val <- crit_val

  # Matrices
  B <- cbind(Phi, X)
  C <- cbind(D, X)

  # Hall and Sheather (1988) bandwidth
  tmp_a <- n ^ (1 / 3)
  tmp_b <- stats::qnorm(1 - 0.5 * alpha) ^ (2 / 3)
  tmp_c <- 1.5 * (stats::dnorm(stats::qnorm(tau)) ^ 2)
  tmp_d <- 2 * (stats::qnorm(tau) ^ 2) + 1
  hs <- tmp_a * tmp_b * ((tmp_c / tmp_d) ^ (1 / 3))
  out$hs <- hs
  # Powell
  bw <- hs
  # Note that the 1 / (2 * n * bw) is negated in the formula for B_tilde
  Psi <- diag(as.numeric(abs(resid) < bw), nrow = n, ncol = n)

  # Find S and J in formula (3.11) of CH (2006)
  S <- tau * (1 - tau) / n * t(B) %*% B # note that the 'n' will be negated
  J <- 1 / (2 * n * bw) * t(B) %*% Psi %*% C # note that the 'n' will be negated
  varcov <- 1 / n * solve(J) %*% S %*% t(solve(J)) # variance of estimator
  out$S <- S
  out$J <- J
  out$varcov <- varcov
  stopifnot(nrow(varcov) == p_D + p_X)
  stopifnot(ncol(varcov) == p_D + p_X)

  # Get standard error for the coefficient of interest
  if (!endogeneous) {
    index <- p_D + index
  }
  se <- sqrt(varcov[index, index])
  out$se <- se

  # Confidence interval
  lower <- center - crit_val * se
  upper <- center + crit_val * se
  out$center <- center
  out$lower <- lower
  out$upper <- upper

  # Coda
  end_clock <- Sys.time()
  elapsed <- difftime(end_clock, start_clock, units = "mins")
  out$time_mins <- elapsed
  out
}
