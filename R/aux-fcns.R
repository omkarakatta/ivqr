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
#'
#' @export
linear_projection <- function(Y, ...) {
  X <- cbind(...)
  coef <- solve(t(X) %*% X) %*% t(X) %*% Y
  X %*% coef
}

### round_to_magnitude -------------------------
#' Round up the nearest order of magnitude
#'
#' For example, if \code{x} is 12, this function returns 100.
#' If \code{x} is 0.12, this function returns 0.1.
#'
#' @param x Number to be rounded (numeric)
#'
#' @return Nearest order of magnitude larger than \code{x}
round_to_magnitude <- function(x) {
  10 ^ (ceiling(log10(x)))
}

### p_val_interpolation -------------------------
#' Perform P-value-weighted Linear Interpolation
#'
#' In line searches, we should interpolate between the two closest values that
#' bound the true boundary of the confidence interval according to the p-value.
#'
#' When performing a line search in a given direction, the final two beta
#' values will surround the true beta value that is the bound of the confidence
#' interval. Interpolating between these two beta values will get us closer to
#' the true beta value. The weights of this interpolation is given by the
#' p-values of the final two beta values from the line search as well as the
#' p-value of the true beta value, which is the alpha-level.
#'
#' @param old_p_val,new_p_val Values of the p-value associated with
#'  \code{old_beta} and \code{new_beta} respectively
#' @param old_beta,new_beta Values of the coefficient that are inside and
#'  outside of the confidence interval (both can't be inside the CI nor can both
#'  be outside)
#' @param alpha Level of the hypothesis test
#'
#' @return Named list:
#'  \enumerate{
#'    \item beta_border: New beta value that is the result of interpolation
#'    \item pi: p-value-based weights
#'  }
#'
#' @seealso \code{\link{line_confint}}, \code{\link{line_confint_interpolation}}
p_val_interpolation <- function(old_p_val,
                                new_p_val,
                                old_beta,
                                new_beta,
                                alpha) {
  pair_p_val <- c(old_p_val, new_p_val)

  msg <- "`old_beta` and `new_beta` are on same side of the confidence interval." #nolint
  both_reject <- old_p_val < alpha & new_p_val < alpha
  both_accept <- old_p_val > alpha & new_p_val > alpha
  send_note_if(msg, both_reject | both_accept, warning)

  pair_beta <- c(old_beta, new_beta)
  ordered <- order(pair_p_val) # ordered[1] = index of smaller p-value

  # construct p-value-based weights
  pi <- (alpha - pair_p_val[ordered[1]]) /
    (pair_p_val[ordered[2]] - pair_p_val[ordered[1]])

  # compute weighted sum of beta values
  beta_border <- (1 - pi) * pair_beta[ordered[1]] + pi * pair_beta[ordered[2]]

  # return new beta value and weight
  list(beta_border = beta_border, pi = pi)
}

### parse_single_log -------------------------
#' Parse information from a Gurobi log file
#'
#' Obtain information on incumbent solutions, Gurobi version, etc. from a
#' Gurobi log file.
#'
#' When \code{information} is "incumbent", this function will retrieve the
#' times at which Gurobi has found an incumbent solution either through
#' heuristics or through solving the problem at the node.
#'
#' @param log_path Path to a Gurobi log file
#' @param information Vector with values:
#'  \enumerate{
#'    \item "incumbent"
#'  }
#'
#' @export
parse_single_log <- function(log_path, information = "incumbent") {
  if (!file.exists(log_path)) {
    stop(paste("Does not exist:", log_path))
  } else {
    # TODO: check this is an actual log file
    log <- readLines(log_path)
  }
  if ("incumbent" %in% information) {
    incumbent_index <- grepl(pattern = "^(H|\\*)", log)
    incumbent_log <- log[incumbent_index]
    # remove everything before and includeing last space
    incumbent_time <- gsub(".*? ", "", incumbent_log)
    # remove all non-numeric characters
    incumbent_time_numeric <- as.numeric(gsub("[^0-9]+", "", incumbent_time))
    result <- incumbent_time_numeric
  }
  result
}

### parse_mult_logs -------------------------
#' Parse multiple log files
#'
#' Use \code{parse_single_log} across multiple log files.
#' Only valid for \code{information = "incumbent"}.
#'
#' @param log_dir Directory of log files
#' @param expr Expression to limit files in \code{log_dir}; defaults to
#'  searching for "log" extension
#' @param information Set to "incumbent"
#' @param value How should the information be returned? defaults to "list",
#'  also can be "data.frame"
parse_mult_logs <- function(log_dir,
                            expr = "log$",
                            information = "incumbent",
                            value = "list") {
  if (!dir.exists(log_dir)) {
    stop(paste("Directory not found:", log_dir))
  }
  files <- list.files(log_dir, pattern = expr)
  info <- lapply(paste0(log_dir, "/", files),
                 function(f){
                   parse_single_log(f, information = "incumbent")
                 })
  names(info) <- files
  result <- info
  if (value == "data.frame") {
    max_soln <- max(sapply(info, length))
    # extend each entry in list to have same number of entries
    extended_info <- lapply(info,
      function(i) {
        l <- length(i)
        new <- max_soln - l
        c(i, rep(NA, new))
      }
    )
    result <- do.call(cbind, extended_info)
  }
  result
}
