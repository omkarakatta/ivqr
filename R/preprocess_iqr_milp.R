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
#' Obtain the residuals from a quantile regression of \code{Y} on \code{D} and
#' \code{X}. Use these residuals to guide which ones are the outliers, and fix
#' their dual variables. In particular, fix the sign of the dual variables
#' associated with residuals whose magnitude is larger than the
#' \eqn{\alpha}-percentile of the absolute residuals.
#' (Note that the initial value of \eqn{\alpha} is given by
#' \code{prop_alpha_initial}.)
#' Then, solve the preprocessed MILP.
#'
#' If the minimized objective of the preprocessed MILP is not 0,
#' then too many dual variables have been fixed.
#' So, relax the number of fixed dual variables by multiplying the bandwidth
#' \eqn{alpha} by \code{r}. Note that \code{r} should be at least 1 to ensure
#' we are relaxing our preprocessing constraints. We then solve the relaxed
#' preprocessed MILP.
#'
#' We repeatedly relax our preprocessing until the minimized objective is
#' finally 0.
#'
#' Note: no warm-starting with \code{start} argument is allowed.
#'
#' @param Y Dependent variable (vector of length n)
#' @param X Exogenous variable (including constant vector) (n by p_X matrix)
#' @param D Endogenous variable (n by p_D matrix)
#' @param Z Instrumental variable (n by p_Z matrix)
#' @param Phi Transformation of X and Z to be used in the program;
#'  defaults to the linear projection of D on X and Z (matrix with n rows)
#' @param tau Quantile (number strictly between 0 and 1)
#' @param M A large number that bounds the absolute value of the residuals
#'  (a positive number); defaults to 2 times the largest absolute residual from
#'  quantile regression of Y on X and D
#' @param prop_alpha_initial Initial value of the bandwidth \eqn{alpha};
#' @param TimeLimit Maximum time (in seconds) spent on a linear program;
#'  defaults to heuristic (numeric)
#' @param globalTimeLimit Maximum time (in seconds) spent on the entire
#' preprocessing function; defaults to Inf (i.e., no time limit) (numeric)
#'  thus, 1 - \code{prop_alpha_initial} is the proportion of observations whose
#'  dual variables (i.e., sign of the residuals) are fixed at the start
#'  (number between 0 and 1)
#' @param r Rate at which we increase the bandwidth (\eqn{alpha}) at every
#'  iteration (number greater than 1)
#' @param show_iterations If TRUE, print the iteration number to the console;
#'  defaults to FALSE (boolean)
#' @param LogFileExt Extension of Gurobi log file; If \code{LogFileName} is
#'  empty, then Gurobi log won't be saved and this argument will be ignored;
#'  defaults to "log" (string)
#' @param ... Arguments that will be passed to \code{\link{iqr_milp}}
#'
#' @return A named list of
#'  \enumerate{
#'    \item final_fit: output of \code{\link{iqr_milp}} from the final
#'      iteration
#'    \item time: time elapsed for this function to run
#'    \item iteration: number of iterations before objective is 0
#'    \item O_neg: indices for residuals that are fixed to be negative
#'      by the final iteration
#'    \item O_pos: indices for residuals that are fixed to be positive
#'      by the final iteration
#'  }
#'
#' @export
preprocess_iqr_milp <- function(Y,
                                D,
                                X,
                                Z,
                                Phi = linear_projection(D, X, Z),
                                tau,
                                M = NULL,
                                TimeLimit = NULL,
                                globalTimeLimit = Inf,
                                prop_alpha_initial = 0.7,
                                r = 1.25,
                                show_iterations = FALSE,
                                LogFileExt = ".log",
                                start = NULL,
                                ...) {

  out <- list() # Initialize list of results to return

  # Start the clock
  clock_start <- Sys.time()

  # Get dimensions of data
  n <- length(Y)
  n_D <- nrow(D)
  n_X <- nrow(X)
  n_Z <- nrow(Z)
  p_D <- ncol(D)
  p_X <- ncol(X)
  p_Z <- ncol(Z)

  if (!is.null(start)) {
    msg <- "No warm-starting allowed with preprocessing"
    send_note_if(msg, TRUE, warning)
    start <- NULL
  }

  # If there are no endogeneous variables, return quantile regression results:
  if (p_D == 0) {
    msg <- paste("p_D is 0 -- running QR instead of IQR MILP...")
    send_note_if(msg, TRUE, warning)
    qr <- quantreg::rq(Y ~ X - 1, tau = tau)
    out$final_fit <- qr
    return(out)
  }

  # Determine preliminary residuals
  if (p_X == 0) {
    resid <- quantreg::rq(Y ~ D - 1, tau = tau)$residuals
  } else if (p_D == 0) {
    resid <- quantreg::rq(Y ~ X - 1, tau = tau)$residuals
  } else {
    resid <- quantreg::rq(Y ~ X + D - 1, tau = tau)$residuals
  }

  # Determine initial residual bounds
  alpha_initial <- stats::quantile(abs(resid), prop_alpha_initial)

  # Determine M before iqr_milp to avoid rerunning quantreg::rq
  if (is.null(M)) {
    # by default, M = 2 * max(resid from QR of Y on X and D)
    # TODO: update heuristic for choosing M
    # TODO: update documentation with default M
    max_qr <- max(abs(resid))
    M <- 2 * max_qr
  }

  # Start the while loop
  alphawidth <- alpha_initial
  status <- "TIME_LIMIT"
  num_fixed_vars_per_iteration <- c()
  time_limit_per_iteration <- c()
  time_elapsed_per_iteration <- c()
  final_objective_per_iteration <- c()

  # Continue the while loop if the program took too long to solve or if
  # the objective (i.e., absolute value of beta_Z) is not 0. Note that
  # if the program is infeasible, i.e., the objective was NULL, we also
  # continue the while loop (we mechanically set the objective to be nonzero
  # later in the code).
  counter <- 0
  while (status == "TIME_LIMIT" || obj != 0) {
    counter <- counter + 1
    while_start_time <- Sys.time()
    # TODO: Is it possible for two iterations of the while loop to have the same number of fixed variables?
    # Fix the most negative and most positive residuals
    O_neg <- which(resid < -1 * alphawidth)
    O_pos <- which(resid > alphawidth)
    O <- c(O_neg, O_pos)
    send_note_if(paste("Alpha:", alphawidth), show_iterations, message)
    # TODO: are we fixing the "dual" variables? I think we can improve the message below
    send_note_if(paste("Number of Fixed Dual Variables:", length(O)), show_iterations, message)
    num_fixed_vars_per_iteration <- c(num_fixed_vars_per_iteration, length(O))
    # Heuristic for time limit
    if (length(O) == 0) {
      TT <- Inf
    } else if (is.null(TimeLimit)) {
      num_free <- n - length(O)
      TT <- exp(num_free / 200 + p_D / 5 + num_free * p_D / 1000) * 4
      if (TT > globalTimeLimit) {
        TT <- globalTimeLimit
      }
    } else {
      TT <- TimeLimit
    }
    send_note_if(paste("TT:", TT), show_iterations, message)
    time_limit_per_iteration <- c(time_limit_per_iteration, TT)
    # IQR
    fit <- iqr_milp(Y = Y,
                    X = X,
                    D = D,
                    Z = Z,
                    Phi = Phi,
                    tau = tau,
                    O_neg = O_neg,
                    O_pos = O_pos,
                    TimeLimit = TT,
                    M = M, # by default, M = 2 * max(resid from QR of Y on X and D)
                    LogFileExt = paste0("_", counter, LogFileExt),
                    ...)
    final_objective_per_iteration <- c(final_objective_per_iteration,
                                       ifelse(is.null(fit$objval),
                                              "NULL",
                                              as.character(round(fit$objval, 3))))
    if (is.null(fit$objval)) {
      obj <- 0.5
    } else {
      obj <- fit$objval
    }
    status <- fit$status
    alphawidth <- alphawidth * r
    if (show_iterations) {
      print(paste("Iteration", counter, "complete"))
    }
    if (TT == Inf & obj != 0) {
      warning("Nonzero Coefficients on Instruments")
      break # exit while loop
    }
    current <- Sys.time()
    while_elapsed <- difftime(current, while_start_time)
    time_elapsed_per_iteration <- c(time_elapsed_per_iteration, while_elapsed)

    elapsed_time <- difftime(current, clock_start, units = "secs")
    if (as.numeric(elapsed_time) > globalTimeLimit) {
      warning(paste("Global Time Limit of", globalTimeLimit, "reached."))
      break # exit while loop
    }
  }

  # Stop the clock
  clock_end <- Sys.time()
  elapsed_time <- difftime(clock_end, clock_start, units = "mins")

  # Return results
  out$status <- fit$status
  out$final_fit <- fit
  out$time <- elapsed_time # mins
  out$O_neg <- O_neg
  out$O_pos <- O_pos
  out$iterations <- counter
  out$num_fixed_vars_per_iteration <- num_fixed_vars_per_iteration
  out$time_limit_per_iteration <- time_limit_per_iteration
  out$time_elapsed_per_iteration <- time_elapsed_per_iteration
  out$final_objective_per_iteration <- final_objective_per_iteration

  return(out)
}

