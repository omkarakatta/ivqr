### Meta -------------------------
###
### Title: Univariate inference under weak identification with a mixed integer
### quadratically constrained program (MIQCP)
###
### Description: After inverting the null on the full weakly identified vector
### of endogenous coefficients, we project the resulting multivariate confidence
### region on the axis corresponding to the coefficient.
### See display (25) - (31) of Jan 25 draft.
###
### Author: Omkar A. Katta
###

### miqcp_proj -------------------------
#' Subvector inference under weak identification with MIQCP
#'
#' Invert the null on the full weakly identified vector of endogeneous
#' coefficients and project the resulting multivariate confidence region on the
#' axis of a specified coefficient.
#'
#' The maximization and minimization of the projection is formulated as a
#' Mixed Integer Quadratically Constrained Program (MIQCP).
#'
#' The index of projection is given by \code{projection_index} and
#' \code{endogeneous}.  For example, if \code{projection_index} is 1 and
#' \code{endogeneous} is TRUE, we project on the axis of \eqn{\beta_{D, 1}}.
#' If \code{projection_index} is 2 and \code{endogeneous} is FALSE, we project
#' on the axis of \eqn{\beta_{X, 2}}.
#'
#' Under weak identification, we invert the null hypothesis on the full weakly
#' identified vector of endogeneous coefficients.  To include exogeneous
#' coefficients in the null, specify the desired indices of the exogeneous
#' variables in \code{beta_X_indices}.
#'
#' @param projection_index Index associated with the coefficient of interest
#'  (numeric between 1 and p_D if \code{endogeneous} is TRUE;
#'  numeric between 1 and p_X if \code{endogeneous} is FALSE)
#' @param endogeneous If TRUE (default), the function will project on the axis
#'  of the corresponding endogeneous variable's coefficient; if FALSE, the
#'  function will project on the axis of the corresponding exogeneous variable's
#'  coefficient (boolean)
#' @param beta_X_indices Indices of endogeneous variables to be included in the
#'  null hypothesis; if NULL (default), none of the coefficients on the
#'  endogeneous variables will be specified in the null (numeric vector)
#' @param alpha Alpha level; defaults to 10% (numeric between 0 and 1)
#' @param sense Maximize or minimize beta_D,j (either "max" or "min")
#' @param orthogonalize_statistic If TRUE, \eqn{\tilde{B}} will be used in
#'  numerator of test statistic; defaults to FALSE; for advanced users only
#' @param homoskedasticity If TRUE, assume density of error at 0 is constant;
#'  defaults to FALSE (boolean)
#' @param kernel Only active if \code{homoskedasticity} is FALSE; either
#'  "Powell" (default) to use the Powell estimator or
#'  "Gaussian" to use a Gaussian kernel; only used when
#'  \code{homoskedasticity} is FALSE
#' @param residuals Residuals from IQR MILP program; if NULL (default), use
#'  naive residuals from quantile regression
#' @param sparse If TRUE (default), use sparse matrices # TODO: incorporate into iqr_milp
#' @inheritParams iqr_milp
#'
#' @importFrom methods as
#'
#' @return A named list of # TODO: update
#'  \enumerate{
#'    \item proj: Gurobi model that was solved
#'    \item params: Gurobi parameters used
#'    \item result: solution to MILP returned by Gurobi
#'    \item status: status of Gurobi's solution
#'    \item beta_X: coefficients on exogenous variables
#'    \item beta_D: coefficients on endogenous variables
#'    \item u: positive part of residuals
#'    \item v: negative part of residuals
#'    \item a: dual variable
#'    \item k: binary variable associated with u
#'    \item l: binary variable associated with v
#'    \item resid: residuals (u - v)
#'    \item objval: value of objective function (beta_D,j)
#'    \item M: big M constant used for complementary slackness conditions
#'    \item Phi_J,Phi_J_minus,X_K,X_K_minus,B,B_tilde,Psi: matrices used in
#'      program
#'    \item projection_index,endogeneous: coefficient onto which the
#'      multivariate confidence region was projected
#'    \item homoskedasticity,kernel: indicates estimator of residual density
#'  }
miqcp_proj <- function(projection_index,
                       endogeneous = TRUE,
                       beta_X_indices = NULL,
                       alpha = 0.1,
                       sense,
                       Y,
                       X,
                       D,
                       Z,
                       Phi = linear_projection(D, X, Z),
                       tau,
                       orthogonalize_statistic = FALSE,
                       homoskedasticity = FALSE,
                       kernel = "Powell",
                       residuals = NULL,
                       O_neg = NULL,
                       O_pos = NULL,
                       M = NULL,
                       TimeLimit = 300,
                       params = list(FeasibilityTol = 1e-6,
                                     OutputFlag = 0),
                       sparse = TRUE, # TODO: document
                       quietly = TRUE,
                       show_progress = TRUE,
                       LogFileName = "",
                       LogFileExt = ".log") {

  out <- list() # Initialize list of results to return

  # Check that alpha is between 0 and 1 (exclusive)
  msg <- "`alpha` must be between 0 and 1."
  send_note_if(msg, alpha <= 0 || alpha >= 1, stop, call. = FALSE)

  # Check that kernel is appropriate
  kernel <- tolower(kernel)
  if (kernel != "powell" & kernel != "gaussian" & !homoskedasticity) {
    stop(paste0("`kernel` should either be 'Powell' or 'Gaussian', not ",
                kernel))
  }

  # Get dimensions of data
  n <- length(Y)
  n_D <- nrow(D)
  n_X <- nrow(X)
  n_Z <- nrow(Z)
  n_Phi <- nrow(Phi)
  p_D <- ncol(D)
  p_X <- ncol(X)

  # Ensure that number of observations is the same
  stopifnot(all.equal(n, n_D))
  stopifnot(all.equal(n, n_X))
  stopifnot(all.equal(n, n_Z))
  stopifnot(all.equal(n, n_Phi))

  # Ensure that there are some RHS variables
  stopifnot(p_D + p_X > 0)

  # Check beta_X_indices and create K
  if (is.null(beta_X_indices)) {
    beta_X_indices <- NULL
  } else {
    stopifnot(sum(beta_X_indices < 0) == 0) # no negative indices
    stopifnot(all.equal(round(beta_X_indices), beta_X_indices)) # indices are integers
    stopifnot(max(beta_X_indices) < p_X) # no indices greater than p_X
    stopifnot(length(beta_X_indices) < p_X) # number of indices specified can't exceed p_X
  }
  K <- rep(FALSE, p_X)
  K[beta_X_indices] <- TRUE
  cardinality_K <- sum(K)

  # Create J
  J <- rep(TRUE, p_D)
  cardinality_J <- sum(J)

  # Check that projection_index is an integer that is between 1 and p_D or p_X (inclusive)
  stopifnot(projection_index > 0) # can't be negative
  stopifnot(length(projection_index) == 1) # we project on axis of only one coefficient
  stopifnot(round(projection_index) == projection_index) # must be integer
  if (endogeneous) {
    stopifnot(projection_index <= p_D)
  } else {
    stopifnot(any(projection_index %in% beta_X_indices))
  }
  out$projection_index <- projection_index
  out$endogeneous <- endogeneous

  # Create basic matrices
  D_J <- D[, J, drop = FALSE]
  D_J_minus <- D[, !J, drop = FALSE]
  X_K <- X[, K, drop = FALSE]
  X_K_minus <- X[, !K, drop = FALSE]
  # Z_J <- Z[, J, drop = FALSE] # not used
  # Z_J_minus <- Z[, !J, drop = FALSE] # not used
  Phi_J <- Phi[, J, drop = FALSE] # by default, Phi = projection of D on X and Z
  Phi_J_minus <- Phi[, !J, drop = FALSE]

  out$Phi_J <- Phi_J
  # out$Phi_J_minus <- Phi_J_minus # don't need to return; this will always be empty since |J| = p_D under weak identification
  out$X_K <- X_K
  out$X_K_minus <- X_K_minus

  # Get residuals
  if (is.null(residuals)) {
    if (p_X == 0) {
      residuals <- quantreg::rq(Y ~ D - 1, tau = tau)$residuals
    } else if (p_D == 0) {
      residuals <- quantreg::rq(Y ~ X - 1, tau = tau)$residuals
    } else {
      residuals <- quantreg::rq(Y ~ X + D - 1, tau = tau)$residuals
    }
  }
  stopifnot(is.numeric(residuals))
  stopifnot(length(residuals) == n)

  if (is.null(M)) {
    # by default, M = 10 * sd(resid from QR of Y on X and D)
    M <- 10 * stats::sd(residuals)
  }
  out$M <- M

  # Decision variables in order from left/top to right/bottom:
  # 1. beta_X
  # 2. beta_D
  # 3. u
  # 4. v
  # 5. a
  # 6. k
  # 7. l

  num_decision_vars <- p_X + p_D + 5 * n

  # Create vector of 1s
  ones <- rep(1, n)
  ones_dv <- rep(1, num_decision_vars) # dv stands for decision variable to remind us that this vector has dimension equal to the number of decision variables

  # Objective: maximize or minimize beta_D,projection_index
  beta_D_obj <- rep(0, p_D)
  beta_X_obj <- rep(0, p_X)
  if (endogeneous) {
    beta_D_obj[projection_index] <- 1
  } else {
    beta_X_obj[projection_index] <- 1
  }
  obj <- c(beta_X_obj,    # beta_X
           beta_D_obj,    # beta_D
           rep(0, n),     # u
           rep(0, n),     # v
           rep(0, n),     # a
           rep(0, n),     # k
           rep(0, n))     # l
  stopifnot(length(obj) == num_decision_vars)

  # Primal Feasibility Constraint (25)
  A_pf <- cbind(X,                  # beta_X
                D,                  # beta_D
                diag(1, nrow = n),  # u
                -diag(1, nrow = n), # v
                diag(0, nrow = n),  # a
                diag(0, nrow = n),  # k
                diag(0, nrow = n))  # l
  b_pf <- Y
  sense_pf <- rep("=", n)

  stopifnot(ncol(A_pf) == num_decision_vars)
  stopifnot(nrow(A_pf) == n)
  stopifnot(length(b_pf) == n)
  stopifnot(length(sense_pf) == n)
  if (sparse) {
    A_pf <- as(A_pf, "sparseMatrix")
  }
  msg <- paste("Primal Feasibility Complete.")
  send_note_if(msg, show_progress, message)

  # Dual Feasibility Constraint (26)
  A_df_X <- cbind(matrix(0, nrow = p_X - cardinality_K, ncol = p_X),  # beta_X
                  matrix(0, nrow = p_X - cardinality_K, ncol = p_D),  # beta_D
                  matrix(0, nrow = p_X - cardinality_K, ncol = n),    # u
                  matrix(0, nrow = p_X - cardinality_K, ncol = n),    # v
                  t(X_K_minus),                                       # a
                  matrix(0, nrow = p_X - cardinality_K, ncol = n),    # k
                  matrix(0, nrow = p_X - cardinality_K, ncol = n))    # l
  b_df_X <- (1 - tau) * t(X) %*% ones
  sense_df_X <- rep("=", p_X)

  stopifnot(ncol(A_df_X) == num_decision_vars)
  stopifnot(nrow(A_df_X) == p_X - cardinality_K)
  stopifnot(length(b_df_X) == p_X - cardinality_K)
  stopifnot(length(sense_df_X) == p_X - cardinality_K)
  if (sparse) {
    A_df_X <- as(A_df_X, "sparseMatrix")
  }
  msg <- paste("Dual Feasibility for X Complete.")
  send_note_if(msg, show_progress, message)

  # Complementary Slackness (28)
  A_cs_uk <- cbind(matrix(0, nrow = n, ncol = p_X),       # beta_X
                   matrix(0, nrow = n, ncol = p_D),       # beta_D
                   diag(1, nrow = n, ncol = n),           # u
                   matrix(0, nrow = n, ncol = n),         # v
                   matrix(0, nrow = n, ncol = n),         # a
                   -M * diag(1, nrow = n, ncol = n),      # k
                   matrix(0, nrow = n, ncol = n))         # l
  b_cs_uk <- rep(0, n)
  sense_cs_uk <- rep("<=", n)

  stopifnot(ncol(A_cs_uk) == num_decision_vars)
  stopifnot(nrow(A_cs_uk) == n)
  stopifnot(length(b_cs_uk) == n)
  stopifnot(length(sense_cs_uk) == n)
  msg <- paste("Complementary Slackness for u and k Complete.")
  if (sparse) {
    A_cs_uk <- as(A_cs_uk, "sparseMatrix")
  }
  send_note_if(msg, show_progress, message)

  A_cs_vl <- cbind(matrix(0, nrow = n, ncol = p_X),       # beta_X
                   matrix(0, nrow = n, ncol = p_D),       # beta_D
                   matrix(0, nrow = n, ncol = n),         # u
                   diag(1, nrow = n, ncol = n),           # v
                   matrix(0, nrow = n, ncol = n),         # a
                   matrix(0, nrow = n, ncol = n),         # k
                   -M * diag(1, nrow = n, ncol = n))      # l
  b_cs_vl <- rep(0, n)
  sense_cs_vl <- rep("<=", n)

  stopifnot(ncol(A_cs_vl) == num_decision_vars)
  stopifnot(nrow(A_cs_vl) == n)
  stopifnot(length(b_cs_vl) == n)
  stopifnot(length(sense_cs_vl) == n)
  if (sparse) {
    A_cs_vl <- as(A_cs_vl, "sparseMatrix")
  }
  msg <- paste("Complementary Slackness for v and l Complete.")
  send_note_if(msg, show_progress, message)

  # Complementary Slackness (29)
  A_cs_ak <- cbind(matrix(0, nrow = n, ncol = p_X),     # beta_X
                   matrix(0, nrow = n, ncol = p_D),     # beta_D
                   matrix(0, nrow = n, ncol = n),       # u
                   matrix(0, nrow = n, ncol = n),       # v
                   diag(1, nrow = n, ncol = n),         # a
                   -diag(1, nrow = n, ncol = n),        # k
                   matrix(0, nrow = n, ncol = n))       # l
  b_cs_ak <- rep(0, n)
  sense_cs_ak <- rep(">=", n)

  stopifnot(ncol(A_cs_ak) == num_decision_vars)
  stopifnot(nrow(A_cs_ak) == n)
  stopifnot(length(b_cs_ak) == n)
  stopifnot(length(sense_cs_ak) == n)
  if (sparse) {
    A_cs_ak <- as(A_cs_ak, "sparseMatrix")
  }
  msg <- paste("Complementary Slackness for a and k Complete.")
  send_note_if(msg, show_progress, message)

  A_cs_al <- cbind(matrix(0, nrow = n, ncol = p_X),   # beta_X
                   matrix(0, nrow = n, ncol = p_D),   # beta_D
                   matrix(0, nrow = n, ncol = n),     # u
                   matrix(0, nrow = n, ncol = n),     # v
                   diag(1, nrow = n, ncol = n),       # a
                   matrix(0, nrow = n, ncol = n),     # k
                   diag(1, nrow = n, ncol = n))       # l
  b_cs_al <- rep(1, n)
  sense_cs_al <- rep("<=", n)

  stopifnot(ncol(A_cs_al) == num_decision_vars)
  stopifnot(nrow(A_cs_al) == n)
  stopifnot(length(b_cs_al) == n)
  stopifnot(length(sense_cs_al) == n)
  if (sparse) {
    A_cs_al <- as(A_cs_al, "sparseMatrix")
  }
  msg <- paste("Complementary Slackness for a and l Complete.")
  send_note_if(msg, show_progress, message)

  # Non-negativity and Boundedness Constraints (30) and (31)
  lb <- c(rep(-Inf, p_X), # beta_X
          rep(-Inf, p_D), # beta_D
          rep(0, n),      # u
          rep(0, n),      # v
          rep(0, n),      # a
          rep(0, n),      # k
          rep(0, n))      # l
  ub <- c(rep(Inf, p_X),  # beta_X
          rep(Inf, p_D),  # beta_D
          rep(Inf, n),    # u
          rep(Inf, n),    # v
          rep(1, n),      # a
          rep(1, n),      # k
          rep(1, n))      # l

  stopifnot(length(lb) == num_decision_vars)
  stopifnot(length(ub) == num_decision_vars)
  msg <- "Non-negativity and Boundedness Constraints Complete."
  send_note_if(msg, show_progress, message)

  # Integrality Constraint (see vtype) (18)
  vtype <- c(rep("C", p_X), # beta_X
             rep("C", p_D), # beta_D
             rep("C", n),   # u
             rep("C", n),   # v
             rep("C", n),   # a
             rep("B", n),   # k
             rep("B", n))   # l

  stopifnot(length(vtype) == num_decision_vars)
  msg <- "Integrality Constraints Complete."
  send_note_if(msg, show_progress, message)

  # Pre-processing: fix residuals of outliers
  O_neg <- sort(O_neg)
  O_pos <- sort(O_pos)
  O <- c(O_neg, O_pos)        # indices of fixed residuals
  if (!is.null(O)) {
    # If a residual is positive, then the associated k must be 1, which means
    # the dual variable, a, must also be 1. Accordingly, the associated l must
    # be 0.
    # If a residual is negative, then the associated l must be 1, which means
    # the dual variable, a, must also be 1. Accordingly, the associated k must
    # be 0.
    fixed <- rep(0, n)
    fixed[O] <- 1
    fixed_mat <- diag(fixed)

    A_pp_a <- cbind(matrix(0, nrow = n, ncol = p_X),    # beta_X
                    matrix(0, nrow = n, ncol = p_D),    # beta_D
                    matrix(0, nrow = n, ncol = n),      # u
                    matrix(0, nrow = n, ncol = n),      # v
                    fixed_mat,                          # a
                    matrix(0, nrow = n, ncol = n),      # k
                    matrix(0, nrow = n, ncol = n))      # l
    b_a_fixed <- rep(0, n)
    b_a_fixed[O_pos] <- 1
    b_a_fixed[O_neg] <- 0
    b_pp_a <- b_a_fixed
    sense_pp_a <- rep("=", n)

    stopifnot(ncol(A_pp_a) == num_decision_vars)
    stopifnot(nrow(A_pp_a) == n)
    stopifnot(length(b_pp_a) == n)
    stopifnot(length(sense_pp_a) == n)
    if (sparse) {
      A_pp_a <- as(A_pp_a, "sparseMatrix")
    }
    msg <- "Pre-processing for a Complete."
    send_note_if(msg, show_progress, message)

    A_pp_k <- cbind(matrix(0, nrow = n, ncol = p_X),    # beta_X
                    matrix(0, nrow = n, ncol = p_D),    # beta_D
                    matrix(0, nrow = n, ncol = n),      # u
                    matrix(0, nrow = n, ncol = n),      # v
                    matrix(0, nrow = n, ncol = n),      # a
                    fixed_mat,                          # k
                    matrix(0, nrow = n, ncol = n))      # l
    b_k_fixed <- rep(0, n)
    b_k_fixed[O_pos] <- 1
    b_k_fixed[O_neg] <- 0
    b_pp_k <- b_k_fixed
    sense_pp_k <- rep("=", n)

    stopifnot(ncol(A_pp_k) == num_decision_vars)
    stopifnot(nrow(A_pp_k) == n)
    stopifnot(length(b_pp_k) == n)
    stopifnot(length(sense_pp_k) == n)
    if (sparse) {
      A_pp_k <- as(A_pp_k, "sparseMatrix")
    }
    msg <- "Pre-processing for k Complete."
    send_note_if(msg, show_progress, message)

    A_pp_l <- cbind(matrix(0, nrow = n, ncol = p_X),  # beta_X
                    matrix(0, nrow = n, ncol = p_D),  # beta_D
                    matrix(0, nrow = n, ncol = n),    # u
                    matrix(0, nrow = n, ncol = n),    # v
                    matrix(0, nrow = n, ncol = n),    # a
                    matrix(0, nrow = n, ncol = n),    # k
                    fixed_mat)                        # l
    b_l_fixed <- rep(0, n)
    b_l_fixed[O_pos] <- 0
    b_l_fixed[O_neg] <- 1
    b_pp_l <- b_l_fixed
    sense_pp_l <- rep("=", n)

    stopifnot(ncol(A_pp_l) == num_decision_vars)
    stopifnot(nrow(A_pp_l) == n)
    stopifnot(length(b_pp_l) == n)
    stopifnot(length(sense_pp_l) == n)
    if (sparse) {
      A_pp_l <- as(A_pp_l, "sparseMatrix")
    }
    msg <- "Pre-processing for l Complete."
    send_note_if(msg, show_progress, message)
  } else {
    A_pp_a <- c()
    b_pp_a <- c()
    sense_pp_a <- c()
    A_pp_k <- c()
    b_pp_k <- c()
    sense_pp_k <- c()
    A_pp_l <- c()
    b_pp_l <- c()
    sense_pp_l <- c()
  }

  # Quadratic Constraint: (27) and Corollary 2
  B <- cbind(Phi_J, X_K)
  stopifnot(nrow(B) == n)
  stopifnot(ncol(B) == cardinality_J + cardinality_K)
  B_minus <- cbind(Phi_J_minus, X_K_minus)
  stopifnot(nrow(B_minus) == n)
  stopifnot(ncol(B_minus) == p_D - cardinality_J + p_X - cardinality_K)
  C_minus <- cbind(D_J_minus, X_K_minus)
  stopifnot(nrow(C_minus) == n)
  stopifnot(ncol(C_minus) == p_D - cardinality_J + p_X - cardinality_K)
  if (homoskedasticity) {
    Psi <- diag(1, nrow = n)
  } else {
    # Hall and Sheather (1988) bandwidth
    tmp_a <- n ^ (1 / 3)
    tmp_b <- stats::qnorm(1 - 0.5 * alpha) ^ (2 / 3)
    tmp_c <- 1.5 * (stats::dnorm(stats::qnorm(tau)) ^ 2)
    tmp_d <- 2 * (stats::qnorm(tau) ^ 2) + 1
    hs <- tmp_a * tmp_b * ((tmp_c / tmp_d) ^ (1 / 3))
    if (kernel == "powell") {
      bw <- hs
      # TODO: double-check whether this is correct
      # Note that the 1 / (2 * n * bw) is negated in the formula for B_tilde
      Psi <- diag(as.numeric(abs(residuals) < bw), nrow = n, ncol = n)
    } else if (kernel == "gaussian") {
      bw <- hs
      # Note that the 1 / (n * bw) is negated in the formula for B_tilde
      Psi <- diag(stats::dnorm(residuals / bw), nrow = n, ncol = n)
    } else {
      stop(
       "Let `homoskedasticity` be TRUE or choose an appropriate `kernel`."
      )
    }
  }
  G_minus <- Psi %*% C_minus
  B_tilde <- B - B_minus %*% solve(t(G_minus) %*% B_minus) %*% t(G_minus) %*% B
  out$B <- B
  out$B_tilde <- B_tilde
  out$Psi <- Psi
  out$homoskedasticity <- homoskedasticity
  out$kernel <- ifelse(homoskedasticity, "homoskedasticity", kernel)

  if (orthogonalize_statistic) {
    tmp <- B_tilde %*% solve(t(B_tilde) %*% B_tilde) %*% t(B_tilde)
  } else {
    tmp <- B %*% solve(t(B_tilde) %*% B_tilde) %*% t(B)
  }
  crit_value <- stats::qchisq(1 - alpha, p_D)

  qc <- matrix(0, nrow = num_decision_vars, ncol = num_decision_vars)
  qc[(p_X + p_D + 2 * n + 1):(p_X + p_D + 3 * n),
     (p_X + p_D + 2 * n + 1):(p_X + p_D + 3 * n)] <- tmp # quadratic
  q <- as.numeric(-2 * (1 - tau) * qc %*% ones_dv) # linear
  qc_rhs_1 <- -1 * (1 - tau)^2 * t(ones_dv) %*% qc %*% ones_dv
  qc_rhs_2 <- tau * (1 - tau) * crit_value
  qc_rhs <- qc_rhs_1 + qc_rhs_2 # constant
  if (sparse) {
    qc <- as(qc, "sparseMatrix")
    q <- as(q, "sparseVector")
  }

  msg <- "Quadratic Constraints Complete."
  send_note_if(msg, show_progress, message)

  # Putting it all together
  proj <- list()
  proj$obj <- obj
  proj$A <- rbind(A_pf,    # Primal Feasibility
                  A_df_X,  # Dual Feasibility - X
                  A_cs_uk, # Complementary Slackness - u and k
                  A_cs_vl, # Complementary Slackness - v and l
                  A_cs_ak, # Complementary Slackness - a and k
                  A_cs_al, # Complementary Slackness - a and l
                  A_pp_a,  # Pre-processing - fixing a
                  A_pp_k,  # Pre-processing - fixing k
                  A_pp_l)  # Pre-processing - fixing l
  # message(paste("A:", nrow(proj$A), ncol(proj$A)))
  # message(paste("A_pf:", nrow(A_pf), ncol(A_pf)))
  # message(paste("A_df_X:", nrow(A_df_X), ncol(A_df_X)))
  # message(paste("A_df_Phi:", nrow(A_df_Phi), ncol(A_df_Phi)))
  # message(paste("A_cs_uk:", nrow(A_cs_uk), ncol(A_cs_uk)))
  # message(paste("A_cs_vl:", nrow(A_cs_vl), ncol(A_cs_vl)))
  # message(paste("A_cs_ak:", nrow(A_cs_ak), ncol(A_cs_ak)))
  # message(paste("A_cs_al:", nrow(A_cs_al), ncol(A_cs_al)))
  # message(paste("A_pp_a:", nrow(A_pp_a), ncol(A_pp_a)))
  # message(paste("A_pp_k:", nrow(A_pp_k), ncol(A_pp_k)))
  # message(paste("A_pp_l:", nrow(A_pp_l), ncol(A_pp_l)))

  # out$b_pf <- b_pf # DEBUG
  # out$b_df_X <- b_df_X # DEBUG
  # out$b_cs_uk <- b_cs_uk # DEBUG
  # out$b_cs_vl <- b_cs_vl # DEBUG
  # out$b_cs_ak <- b_cs_ak # DEBUG
  # out$b_cs_al <- b_cs_al # DEBUG
  # out$b_pp_a <- b_pp_a # DEBUG
  # out$b_pp_k <- b_pp_k # DEBUG
  # out$b_pp_l <- b_pp_l # DEBUG

  proj$rhs <- c(b_pf,    # Primal Feasibility
                b_df_X,  # Dual Feasibility - X
                b_cs_uk, # Complementary Slackness - u and k
                b_cs_vl, # Complementary Slackness - v and l
                b_cs_ak, # Complementary Slackness - a and k
                b_cs_al, # Complementary Slackness - a and l
                b_pp_a,  # Pre-processing - fixing a
                b_pp_k,  # Pre-processing - fixing k
                b_pp_l)  # Pre-processing - fixing l
  # message(paste("b:", length(proj$rhs)))

  proj$sense <- c(sense_pf,   # Primal Feasibility
                 sense_df_X,  # Dual Feasibility - X
                 sense_cs_uk, # Complementary Slackness - u and k
                 sense_cs_vl, # Complementary Slackness - v and l
                 sense_cs_ak, # Complementary Slackness - a and k
                 sense_cs_al, # Complementary Slackness - a and l
                 sense_pp_a,  # Pre-processing - fixing a
                 sense_pp_k,  # Pre-processing - fixing k
                 sense_pp_l)  # Pre-processing - fixing l
  # message(paste("sense:", length(proj$sense)))

  # Quadratic Constraint
  proj$quadcon[[1]]$Qc <- qc
  proj$quadcon[[1]]$q <- q
  proj$quadcon[[1]]$rhs <- qc_rhs

  proj$lb <- lb
  proj$ub <- ub
  proj$vtype <- vtype
  proj$modelsense <- sense
  params$TimeLimit <- TimeLimit
  if (LogFileName != "") {
    params$LogFile <- paste0(LogFileName, LogFileExt)
  }

  # out$proj <- proj # DEBUG
  # return(out) # DEBUG

  result <- gurobi::gurobi(proj, params)

  msg <- paste("Mixed Integer Quadratic Program Complete.")
  send_note_if(msg, show_progress, message)

  # Return results
  msg <- paste("Status of MIQCP Projection program:",
               result$status,
               "| Objective:",
               format(result$objval, scientific = F, digits = 10))
  send_note_if(msg, !quietly, message)  # Print status of program if !quietly

  out$proj <- proj
  out$params <- params
  out$result <- result
  out$status <- result$status
  if (result$status %in% c("OPTIMAL", "SUBOPTIMAL")) {
    answer <- result$x
    out$beta_X <- answer[1:p_X]
    out$beta_D <- answer[(p_X + 1):(p_X + p_D)]
    out$u <- answer[(p_X + p_D + 1):(p_X + p_D + n)]
    out$v <- answer[(p_X + p_D + n + 1):(p_X + p_D + 2 * n)]
    out$a <- answer[(p_X + p_D + 2 * n + 1):(p_X + p_D + 3 * n)]
    out$k <- answer[(p_X + p_D + 3 * n + 1):(p_X + p_D + 4 * n)]
    out$l <- answer[(p_X + p_D + 4 * n + 1):(p_X + p_D + 5 * n)]

    out$resid <- out$u - out$v
    out$objval <- result$objval
  }

  return(out)
}
