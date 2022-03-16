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

### h_to_beta -------------------------

#' Find the coefficients given the active basis
#'
#' Given the active basis, find the coefficients that solve the IQR problem.
#'
#' @param h Active basis in terms of the data provided
#' @param Y Dependent variable (vector of length n)
#' @param X Exogenous variable (including constant vector) (n by p_X matrix)
#' @param D Endogenous variable (n by p_D matrix)
#' @param Z Instrumental variable (n by p_Z matrix)
#' @param Phi Transformation of X and Z to be used in the program;
#'  defaults to the linear projection of D on X and Z (matrix with n rows)
#'
#' @return A named list of coefficients
#' \enumerate{
#'   \item \code{beta_D}: coefficients on the endogeneous variables
#'   \item \code{beta_X}: coefficients on the exogenous variables
#'   \item \code{beta_Phi}: coefficients on the transformed instruments
#' }
#'
#' @family mcmc_subsampling
h_to_beta <- function(h, Y, X, D, Z, Phi = linear_projection(D, X, Z)) {
  # dimensions
  p_Phi <- ncol(Phi)
  p_D <- ncol(D)
  stopifnot(p_Phi == p_D)
  p_X <- ncol(X)
  p <- p_Phi + p_X

  design <- cbind(X, Phi)
  designh <- design[h, ]
  a <- solve(designh, Y[h])
  B <- solve(designh, D[h, , drop = FALSE])

  # Solve system from the lower rows, corresponding to beta_Phi = 0
  a_Phi <- a[(p_X + 1):p]
  B_Phi <- B[(p_X + 1):p, ]

  beta_D <- solve(B_Phi, a_Phi)
  beta_X_Phi <- a - B %*% beta_D
  beta_X <- beta_X_Phi[1:p_X]
  beta_Phi <- beta_X_Phi[(p_X + 1):p]

  list("beta_D" = beta_D,
       "beta_X" = beta_X,
       "beta_Phi" = beta_Phi)
}

### run_concentrated_qr

#' @param beta_D Endogenous coefficients (vector of length p_D)
#' @param Y Dependent variable (vector of length n)
#' @param X Exogenous variable (including constant vector) (n by p_X matrix)
#' @param D Endogenous variable (n by p_D matrix)
#' @param Z Instrumental variable (n by p_Z matrix)
#' @param Phi Transformation of X and Z to be used in the program;
#'  defaults to the linear projection of D on X and Z (matrix with n rows)
#' @param tau Quantile (numeric)
#'
#' @return Results of QR(Y - D %*% beta_D ~ X + Phi) for a given \code{tau}
run_concentrated_qr <- function(beta_D, Y, X, D, Z, Phi = linear_projection(D, X, Z), tau) {
  quantreg::rq(Y - D %*% beta_D ~ X + Phi - 1, tau = tau)
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

# compute_foc_conditions -------------------------------------------------------

#' @param h Indices in active basis
#' @param Y,X,D,Phi Data
#' @param tau Quantile
#'
#' @return A p by n matrix
compute_foc_conditions <- function(
  h,
  beta_D = NULL,
  beta_X = NULL,
  Y, X, D, Phi,
  tau
) {
  if (is.null(beta_D)) {
    coef_full <- h_to_beta(h, Y = Y, X = X, D = D, Phi = Phi)
    beta_D <- coef_full$beta_D
    beta_X <- coef_full$beta_X
  }
  residuals <- Y - D %*% beta_D - X %*% beta_X
  design <- cbind(X, Phi)
  designh_inv <- solve(design[h, , drop = FALSE])

  xi_i <- vector("list", length = nrow(Y)) # create `n` matrices of dim 1 by `p`
  for (i in seq_len(nrow(Y))) {
    xi_i[[i]] <-
      as.numeric(!is.element(i, h)) * # if index is in active basis, set xi to 0
      (tau - as.numeric(residuals[i] < 0)) *
      design[i, ] %*%
      designh_inv
  }
  t(do.call(rbind, xi_i)) # p by n matrix; note: don't cbind row matrices
}

### ot -------------------------

#' Compute OT between identical uniform distributions
#'
#' Compute transport map between U(1,...,n-p) and U(1,...,n-p) where C(i, j) =
#' norm difference of pre[, i] and post[, j] where i and j are between 1 and
#' n-p.
#'
#' We use Gurobi to solve the OT problem if \code{method} is "gurobi".
#' We use the transport package to solve the OT problem if \code{method} is "transport".
ot <- function(pre, post, params = list(OutputFlag = 0), method = "gurobi") {
  n_minus_p <- ncol(pre)
  stopifnot(n_minus_p == ncol(post))

  # prelims
  num_decision_vars <- n_minus_p^2

  # compute cost matrix
  c_ij <- matrix(0, nrow = n_minus_p, ncol = n_minus_p)
  for (col in seq_len(n_minus_p)) {
    # get diagonals
    c_ij[col, col] <- sum(abs(pre[, col] - post[, col]))
    for (i in seq_len(n_minus_p - col)) {
      row <- i + col
      # create a lower-triangular matrix with the cost
      c_ij[row, col] <- sum(abs(pre[, row] - post[, col]))
    }
  }
  # fill out the cost matrix
  c_ij[upper.tri(c_ij)] <- c_ij[lower.tri(c_ij)]

  if (tolower(method) == "gurobi") {
    # create constraints
    const_pre <- vector("list", length = n_minus_p)
    const_post <- vector("list", length = n_minus_p)
    for (i in seq_len(n_minus_p)) {
      zeros_left <- matrix(0, nrow = n_minus_p, ncol = i - 1)
      ones <- matrix(1, nrow = n_minus_p, ncol = 1)
      zeros_right <- matrix(0, nrow = n_minus_p, ncol = n_minus_p - i)
      a_mat <- cbind(zeros_left, ones, zeros_right)
      const_pre[[i]] <- c(a_mat)
      const_post[[i]] <- c(t(a_mat))
    }
    a_mat_pre <- do.call(rbind, const_pre)
    a_mat_post <- do.call(rbind, const_post)
    a_mat <- rbind(a_mat_pre, a_mat_post)
    stopifnot(ncol(a_mat) == num_decision_vars)
    stopifnot(nrow(a_mat) == n_minus_p * 2)

    # create program
    model <- list()
    model$obj <- c(c_ij) # turn into vector (go down each column)
    model$A <- a_mat
    model$sense <- rep("=", length = n_minus_p * 2)
    model$rhs <- rep(1, length = n_minus_p * 2)
    model$vtype <- rep("C", num_decision_vars)

    sol <- gurobi::gurobi(model, params)
    status <- sol$status
    t_ij <- matrix(sol$x, nrow = n_minus_p, ncol = n_minus_p)
    map <- apply(t_ij, 1, function(x) which(x == 1))
  } else if (tolower(method) == "transport") {
    unif <- rep(1, length = n_minus_p)
    sol <- transport::transport(unif, unif, costm = c_ij)
    model <- NA
    status <- NA
    map <- sol$to
  }

  list(
    model = model, # gurobi-specific
    status = status, # gurobi-specific
    sol = sol,
    c_ij = c_ij,
    map = map
  )
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
    incumbent_warmstart_bool <- sum(grepl(pattern = "Loaded user MIP start with objective", log)) > 0
    if (incumbent_warmstart_bool) {
      # find start of table, then add three rows
      # first row below start is subheader; second row is empty; third row is warm-start!
      start_index <- grep(pattern = "Unexpl", log) + 3
      incumbent_index[start_index] <- TRUE
    }
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
#'
#' @export
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

### Parse sol files -------------------------
#' Parse a single solution (.sol) file
#'
#' Obtain solution vector from a single .sol file.
#'
#' @param sol_path Path to .sol file
#'
#' @export
parse_single_sol <- function(sol_path) {
  if (!file.exists(sol_path)) {
    stop(paste("Does not exist:", sol_path))
  } else {
    # TODO: check this is an actual sol file
    sol <- readLines(sol_path)
  }
  sol <- sol[3:length(sol)] # remove first two comments
  varnames <- gsub(" .*", "", sol) # get names of decision variables
  values <- as.numeric(gsub(".* ", "", sol)) # get values of each variable
  names(values) <- varnames
  result <- values
  result
}

### Parse multiple sol files -------------------------
#' Parse multiple solution (.sol) files
#'
#' Obtain solution from multiple .sol files.
#'
#' @param sol_dir Path to directory with .sol files
#' @param expr Expression to limit files in \code{sol_dir}; defaults to
#'  searching for "sol" extension
#' @param value How should the information be returned? defaults to "list",
#'  accepts "data.frame" if the number of decision veariables are the same
#'  across solution files
#'
#' @export
parse_mult_sols <- function(sol_dir, expr = "sol$", value = "list") {
  if (!dir.exists(sol_dir)) {
    stop(paste("Directory not found:", sol_dir))
  }
  files <- list.files(sol_dir, pattern = expr)
  info <- lapply(paste0(sol_dir, "/", files), parse_single_sol)
  numbers_and_extension <- gsub(".*_[^0-9]*", "", files)
  numbers <- as.numeric(gsub(".sol$", "", numbers_and_extension))
  names(info) <- numbers
  result <- info
  if (value == "data.frame") {
    if (length(unique(sapply(result, length))) == 1) {
      result <- do.call(cbind, info)[, as.character(sort(numbers))]
    } else {
      warning("Number of decision variables are different; returning list")
    }
  }
  result
}

### Connect sol and log files -------------------------
parse_logsol <- function(log_path, sol_dir, expr = "sol$") {
  print("Parsing log for times to incumbent solutions...")
  time <- parse_single_log(log_path)
  print("Parsing sol for incumbent solutions...")
  sols <- parse_mult_sols(sol_dir, expr = expr, value = "data.frame")
  list(time = time, sols = sols, result = rbind(sols, time))
}
