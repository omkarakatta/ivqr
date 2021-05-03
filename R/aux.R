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
#' \code{\link{iqr_milp}}, \code{\link{\miqcp_proj}}, and potentially other
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
manipulate_qr_residuals <- function(string_formula,
                                    tau,
                                    factor = 1,
                                    FUN = stats::sd,
                                    ...) {
  qr <- quantreg::rq(as.formula(string_formula), tau)
  resid <- qr$resid
  fun_resid <- FUN(resid, ...)
  factor * fun_resid
}

