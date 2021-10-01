### Title: Compute warm-starting solutions
###
### Description: Warm-starting the IQR MILP helps us find an incumbent solution
### from the start. These functions help us find these warm-starting solutions.
### For an exploration, see:
### ../scrap/education-autor/iqr_milp/warmstart/warmstart.R
###
### Author: Omkar A. Katta

### compute_warmstart -------------------------
#' Compute warm-starting solutions
#'
#' The result of this function should be fed into the \code{start} argument of
#' \code{iqr_milp}.
#'
#' If \code{method} is NULL, we won't warm-start the program at all, i.e., we
#' return NULL.
#'
#' If \code{method} contains the phrase "naive3", we construct the
#' warm-starting solution with two quantile regressions.
#' First, we run a preliminary quantile regression of Y on D and X to obtain
#' \code{beta_D}. Then, we concentrate out D from Y according to \code{beta_D}.
#' After regressing the concentrated-out Y on Phi and X, we obtain
#' \code{beta_Phi} and \code{beta_X}. Using this information, we can construct
#' our decision vector.
#'
#' If \code{method} contains the phrase "full", then we not only specify the
#' coefficients in the decision vector, we also specify the dual variables,
#' residuals, and the binary variables.
#'
#'
#' @param Y Dependent variable (vector of length n)
#' @param X Exogenous variable (including constant vector) (n by p_X matrix)
#' @param D Endogenous variable (n by p_D matrix)
#' @param Z Instrumental variable (n by p_Z matrix)
#' @param Phi Transformation of X and Z to be used in the program;
#'  defaults to the linear projection of D on X and Z (n by p_D matrix)
#' @param tau Quantile (numeric between 0 and 1)
#' @param method Determines how we compute the warm-starting solution;
#'  defaults to NULL, see 'Details'
#'
#' @return Either NULL or a vector of values to warm-start the program; the
#'  output is meant to be fed into the \code{start} parameter of \code{iqr_milp}
compute_warmstart <- function(Y,
                              X,
                              D,
                              Z,
                              Phi = linear_projection(D, X, Z),
                              tau,
                              method = NULL) {
  n <- nrow(Y)
  if (is.null(method)) {
    NULL
  } else if (grepl("naive3", method)) {
    prelim <- quantreg::rq(Y ~ D + X - 1, tau = tau)
    prelim_coef <- coef(prelim)
    beta_D <- prelim_coef[seq_len(ncol(D))]
    Y_tilde <- Y - D %*% beta_D
    qr <- quantreg::rq(Y_tilde ~ Phi + X - 1, tau = tau)
    qr_coef <- coef(qr)
    beta_Phi <- qr_coef[seq_len(ncol(Phi))]
    beta_Phi_plus <- sapply(beta_Phi, function(i) max(i, 0))
    beta_Phi_minus <- sapply(beta_Phi, function(i) max(-i, 0))
    beta_X <- qr_coef[seq_len(ncol(X)) + ncol(Phi)]
    start <- c(beta_X, beta_Phi_plus, beta_Phi_minus, beta_D, rep(NA, 5 * n))
    if (grepl("full", method)) {
      dual <- qr$dual
      resid <- qr$residuals
      resid_plus <- sapply(resid, function(i) {
        i <- ifelse(isTRUE(all.equal(as.numeric(i), 0)), 0, i)
        max(i, 0)
      })
      resid_minus <- sapply(resid, function(i) {
        i <- ifelse(isTRUE(all.equal(as.numeric(i), 0)), 0, i)
        max(-i, 0)
      })
      k <- rep(NA, length = "n")
      l <- rep(NA, length = "n")
      k[resid_plus > 0] <- 1
      l[resid_plus > 0] <- 0
      k[resid_minus > 0] <- 0
      l[resid_minus > 0] <- 1
      k[resid_plus == 0 & resid_minus == 0] <- 0 # resid = 0 => a \in (0, 1)
      l[resid_plus == 0 & resid_minus == 0] <- 0 # resid = 0 => a \in (0, 1)
      start <- c(beta_X, beta_Phi_plus, beta_Phi_minus, beta_D,
                 resid_plus, resid_minus, dual, k, l)
    }
    start
  } else {
    stop("Unknown argument for `method`")
  }
}
