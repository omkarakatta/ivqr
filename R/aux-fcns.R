### Meta -------------------------
###
### Title: Auxiliary Functions
###
### Description: Functions that are shared by main functions
###
### Author: Omkar A. Katta
###


### manipulate_qr_residuals -------------------------
#' Manipulate residuals from a quantile regression
#'
#' Apply a function to the residuals of a quantile regression
#'
#' The purpose of this function is to define the big M constant in
#' \code{iqr_milp}, \code{miqcp_proj}, and potentially other
#' functions.
#'
#' The formula relies on variables defined outside the scope of the function.
#'
#' @param string_formula A formula to be passed to \code{quantreg::rq} in the
#'  the form of a string (string)
#' @param tau Quantile of interest (number between 0 and 1)
#' @param factor Number to be multiplied to the output of \code{FUN};
#'  defaults to 1 (numeric)
#' @param FUN A function to manipulate the residuals;
#'  defaults to \code{stats::sd}
#' @param ... Arguments passed to \code{FUN}
#'
#' @return Result of \code{factor} multiplied by the output of \code{FUN}
#'  evaluated at the quantile regression specified by \code{string_formula} and
#'  \code{tau}
manipulate_qr_residuals <- function(string_formula,
                                    tau,
                                    factor = 1,
                                    FUN = stats::sd,
                                    ...) {
  qr <- quantreg::rq(stats::as.formula(string_formula), tau)
  resid <- qr$resid
  fun_resid <- FUN(resid, ...)
  factor * fun_resid
}

### linear_projection -------------------------
#' Linearly project a vector onto the space spanned by other vectors
#'
#' Find the fitted values from a regression of each column of \code{Y} on
#' regressors, which are given by matrices in \code{...}.
#'
#' Ensure that the rows of \code{Y} and \code{...} are the same!
#'
#' @param Y Vector or matrix that will be projected (n by p_Y matrix)
#' @param ... "Regressors" whose span is the destination of the projection
#'  (matrix of n rows)
#'
#' @return Linear projection of \code{Y} on \code{...}
#'  (matrix of dimension nrow(Y) by ncol(Y))
linear_projection <- function(Y, ...) {
  X <- cbind(...)
  coef <- solve(t(X) %*% X) %*% t(X) %*% Y
  X %*% coef
}
