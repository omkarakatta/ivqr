### Meta -------------------------
###
### Title: Compute inverse quantile regression estimator
###
### Description: Solve a mixed integer linear program to compute the inverse
### quantile regression estimator.
###
### Author: Omkar A. Katta
###

### iqr_milp -------------------------
#' Compute inverse quantile regression estimator
#'
#' Solve a mixed integer linear program to compute the inverse
#' quantile regression estimator.
#'
#' @param Y Dependent variable (vector of length n)
#' @param X Exogenous variable (including constant vector) (n by p_X matrix)
#' @param D Endogenous variable (n by p_D matrix)
#' @param Z Instrumental variable (n by p_Z matrix)
#' @param Phi Transformation of X and Z to be used in the program;
#'  defaults to the linear projection of D on X and Z (matrix with n rows)
#' @param tau Quantile (number between 0 and 1)
#' @param O_neg,O_pos Indices for residuals whose sign is fixed to be negative
#'  and positive, respectively (vectors)
#' @param M A large number that bounds the absolute value of the residuals
#'  (a positive number); defaults to 2 times the largest absolute residual from
#'  quantile regression of Y on X and D
#' @param TimeLimit Maximum time (in seconds) spent on a linear program;
#'  defaults to 300, will be appended to \code{params}
#' @param VarHintVal_bool If TRUE, instead of fixing the binary variables as
#'  determined by \code{O_neg,O_pos}, we will "hint" these values; defaults to
#'  FALSE; see
#'  \url{https://www.gurobi.com/documentation/9.1/refman/varhintval.html#attr:VarHintVal};
#'  overrides \code{attributes}
#' @param VarHintPri_bool If TRUE, we specify the strength of the "hints"
#'  according to the magnitude of the residual from a naive QR;
#'  only valid if VarHintVal_bool is TRUE; overrides \code{attributes};
#'  defaults to FALSE;
#'  see \url{https://www.gurobi.com/documentation/9.1/refman/varhintpri.html#attr:VarHintPri}
#' @param BranchPriority_bool If TRUE, we specify which variables will be
#'  branched on first according to the magnitude of the residual from a naive QR;
#'  defaults to FALSE;
#'  see \url{https://www.gurobi.com/documentation/9.1/refman/branchpriority.html}
#' @param attributes named list of Gurobi attributes, see
#'  \url{https://www.gurobi.com/documentation/9.1/refman/attributes.html} and
#'  \url{https://www.gurobi.com/documentation/9.1/refman/the_model_argument.html};
#'  defaults to empty list
#' @param params Gurobi parameters, see
#'  \url{https://www.gurobi.com/documentation/9.1/refman/parameter_descriptions.html}
#' @param start Gurobi attribute, see
#'  \url{https://www.gurobi.com/documentation/9.1/examples/mip_starts.html}; If
#'  NULL (default), no starting solution will be provided; see \code{compute_warmstart}
#' @param fix Fix decision variables; If NULL (default), no starting solution will
#'  be provided; if NA, the respective decision variable won't be fixed;
#'  decision variables are:
#'  \enumerate{
#'    \item beta_X
#'    \item beta_Phi_plus
#'    \item beta_Phi_minus ; note: abs(beta_Phi) = beta_Phi_plus + beta_Phi_minus
#'    \item beta_D
#'    \item u
#'    \item v
#'    \item a
#'    \item k
#'    \item l
#'  }
#' @param sparse If TRUE (default), use sparse matrices
#' @param quietly If TRUE (default), supress messages during execution (boolean)
#' @param show_progress If TRUE (default), sends progress messages during
#'  execution (boolean)
#' @param LogFileName Name of Gurobi log file; If string is empty (default),
#'  Gurobi log won't be saved (string)
#' @param LogFileExt Extension of Gurobi log file; If \code{LogFileName} is
#'  empty, then Gurobi log won't be saved and this argument will be ignored;
#'  defaults to "log" (string)
#'
#' @importFrom methods as
#'
#' @return A named list of
#'  \enumerate{
#'    \item iqr: Gurobi model that was solved
#'    \item params: Gurobi parameters used
#'    \item result: solution to MILP returned by Gurobi
#'    \item status: status of Gurobi's solution
#'    \item beta_X: coefficients on exogenous variables
#'    \item beta_Phi_plus: positive part of coefficients on instruments
#'    \item beta_Phi_minus: negative part of coefficients on instruments
#'    \item beta_D: coefficients on endogenous variables
#'    \item u: positive part of residuals
#'    \item v: negative part of residuals
#'    \item a: dual variable
#'    \item k: binary variable associated with u
#'    \item l: binary variable associated with v
#'    \item beta_Phi: coefficients on instruments
#'      (beta_Phi_plus - beta_Phi_minus)
#'    \item resid: residuals (u - v)
#'    \item objval: value of objective function (absolute value of beta_Phi)
#'  }
#'
#' @export
iqr_milp <- function(Y,
                     X,
                     D,
                     Z,
                     Phi = linear_projection(D, X, Z),
                     tau,
                     O_neg = NULL,
                     O_pos = NULL,
                     M = NULL,
                     TimeLimit = 300,
                     VarHintVal_bool = FALSE,
                     VarHintPri_bool = FALSE,
                     BranchPriority_bool = FALSE,
                     sparse = TRUE,
                     attributes = list(),
                     params = list(FeasibilityTol = 1e-6,
                                   LogToConsole = 0),
                     start = compute_warmstart(Y = Y,
                                               X = X,
                                               D = D,
                                               Z = Z,
                                               Phi = Phi,
                                               tau = tau,
                                               method = NULL),
                     fix = NULL,
                     quietly = TRUE,
                     show_progress = TRUE,
                     LogFileName = "",
                     LogFileExt = ".log") {

  out <- list() # Initialize list of results to return

  send_note_if(paste("TimeLimit:", TimeLimit, "secs"), show_progress, message)

  # Get dimensions of data
  n <- length(Y)
  n_D <- nrow(D)
  n_X <- nrow(X)
  n_Z <- nrow(Z)
  n_Phi <- nrow(Phi)
  p_D <- ncol(D)
  p_X <- ncol(X)
  p_Z <- ncol(Z)
  p_Phi <- ncol(Phi)

  # Ensure that number of observations is the same
  stopifnot(all.equal(n, n_D))
  stopifnot(all.equal(n, n_X))
  stopifnot(all.equal(n, n_Z))
  stopifnot(all.equal(n, n_Phi))

  # If there are no endogeneous variables, return quantile regression results:
  if (p_D == 0) {
    msg <- paste("p_D is 0 -- running QR instead of IQR MILP...")
    send_note_if(msg, TRUE, warning)
    qr <- quantreg::rq(Y ~ X - 1, tau = tau)
    out <- qr
    return(out)
  }

  # Create vector of 1s
  ones <- rep(1, n)

  out$Phi <- Phi # by default, Phi = projection of D on X and Z

  if (is.null(M)) {
    # by default, M = 2 * max(resid from QR of Y on X and D)
    # TODO: update heuristic for choosing M
    # TODO: update documentation with default M
    if (p_X == 0) {
      max_qr <- max(abs(quantreg::rq(Y ~ D - 1, tau = tau)$residuals))
    } else if (p_D == 0) {
      max_qr <- max(abs(quantreg::rq(Y ~ X - 1, tau = tau)$residuals))
    } else {
      max_qr <- max(abs(quantreg::rq(Y ~ X + D - 1, tau = tau)$residuals))
    }
    M <- 2 * max_qr
  }
  out$M <- M


  # Decision variables in order from left/top to right/bottom:
  # 1. beta_X
  # 2. beta_Phi_plus
  # 3. beta_Phi_minus ; note: abs(beta_Phi) = beta_Phi_plus + beta_Phi_minus
  # 4. beta_D
  # 5. u
  # 6. v
  # 7. a
  # 8. k
  # 9. l

  num_decision_vars <- p_X + 2 * p_Phi + p_D + 5 * n

  # Objective: minimize absolute value of \beta_Phi
  obj <- c(rep(0, p_X),   # beta_X
           rep(1, p_Phi), # beta_Phi_plus
           rep(1, p_Phi), # beta_Phi_minus
           rep(0, p_D),   # beta_D
           rep(0, n),     # u
           rep(0, n),     # v
           rep(0, n),     # a
           rep(0, n),     # k
           rep(0, n))     # l
  stopifnot(length(obj) == num_decision_vars)

  # Fix decision variables according to `fix`
  if (is.null(fix)) {
    A_fix <- c()
    b_fix <- c()
    sense_fix <- c()
  } else {
    stopifnot(length(fix) - num_decision_vars == 0)
    not_na <- !is.na(fix)
    A_fix <- diag(1, num_decision_vars)[not_na, ]
    b_fix <- fix[not_na]
    sense_fix <- rep('=', sum(not_na))
    stopifnot(ncol(A_fix) == num_decision_vars)
    stopifnot(nrow(A_fix) == sum(not_na))
    stopifnot(length(b_fix) == sum(not_na))
  }

  # Primal Feasibility Constraint (11)
  A_pf <- cbind(X,                  # beta_X
                Phi,                # beta_Phi_plus
                -Phi,               # beta_Phi_minus
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
  msg <- paste("Primal Feasibility Complete.")
  send_note_if(msg, show_progress, message)

  # Dual Feasibility Constraint (13) and (14)
  A_df_X <- cbind(matrix(0, nrow = p_X, ncol = p_X),  # beta_X
                  matrix(0, nrow = p_X, ncol = p_Phi),  # beta_Phi_plus
                  matrix(0, nrow = p_X, ncol = p_Phi),  # beta_Phi_minus
                  matrix(0, nrow = p_X, ncol = p_D),  # beta_D
                  matrix(0, nrow = p_X, ncol = n),  # u
                  matrix(0, nrow = p_X, ncol = n),  # v
                  t(X),                 # a
                  matrix(0, nrow = p_X, ncol = n),  # k
                  matrix(0, nrow = p_X, ncol = n))  # l
  b_df_X <- (1 - tau) * t(X) %*% ones
  sense_df_X <- rep("=", p_X)

  stopifnot(ncol(A_df_X) == num_decision_vars)
  stopifnot(nrow(A_df_X) == p_X)
  stopifnot(length(b_df_X) == p_X)
  stopifnot(length(sense_df_X) == p_X)
  msg <- paste("Dual Feasibility for X Complete.")
  send_note_if(msg, show_progress, message)


  A_df_Phi <- cbind(matrix(0, nrow = p_Phi, ncol = p_X),    # beta_X
                    matrix(0, nrow = p_Phi, ncol = p_Phi),  # beta_Phi_plus
                    matrix(0, nrow = p_Phi, ncol = p_Phi),  # beta_Phi_minus
                    matrix(0, nrow = p_Phi, ncol = p_D),    # beta_D
                    matrix(0, nrow = p_Phi, ncol = n),      # u
                    matrix(0, nrow = p_Phi, ncol = n),      # v
                    t(Phi),                                 # a
                    matrix(0, nrow = p_Phi, ncol = n),      # k
                    matrix(0, nrow = p_Phi, ncol = n))      # l
  b_df_Phi <- (1 - tau) * t(Phi) %*% ones
  sense_df_Phi <- rep("=", p_Phi)

  stopifnot(ncol(A_df_Phi) == num_decision_vars)
  stopifnot(nrow(A_df_Phi) == p_Phi)
  stopifnot(length(b_df_Phi) == p_Phi)
  stopifnot(length(sense_df_Phi) == p_Phi)
  msg <- paste("Dual Feasibility for Phi Complete.")
  send_note_if(msg, show_progress, message)

  # Complementary Slackness (16) and (17)
  A_cs_uk <- cbind(matrix(0, nrow = n, ncol = p_X),       # beta_X
                   matrix(0, nrow = n, ncol = p_Phi),     # beta_Phi_plus
                   matrix(0, nrow = n, ncol = p_Phi),     # beta_Phi_minus
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
  send_note_if(msg, show_progress, message)

  A_cs_vl <- cbind(matrix(0, nrow = n, ncol = p_X),       # beta_X
                   matrix(0, nrow = n, ncol = p_Phi),     # beta_Phi_plus
                   matrix(0, nrow = n, ncol = p_Phi),     # beta_Phi_minus
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
  msg <- paste("Complementary Slackness for v and l Complete.")
  send_note_if(msg, show_progress, message)

  A_cs_ak <- cbind(matrix(0, nrow = n, ncol = p_X),     # beta_X
                   matrix(0, nrow = n, ncol = p_Phi),   # beta_Phi_plus
                   matrix(0, nrow = n, ncol = p_Phi),   # beta_Phi_minus
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
  msg <- paste("Complementary Slackness for a and k Complete.")
  send_note_if(msg, show_progress, message)

  A_cs_al <- cbind(matrix(0, nrow = n, ncol = p_X), # beta_X
                   matrix(0, nrow = n, ncol = p_Phi), # beta_Phi_plus
                   matrix(0, nrow = n, ncol = p_Phi), # beta_Phi_minus
                   matrix(0, nrow = n, ncol = p_D), # beta_D
                   matrix(0, nrow = n, ncol = n), # u
                   matrix(0, nrow = n, ncol = n), # v
                   diag(1, nrow = n, ncol = n), # a
                   matrix(0, nrow = n, ncol = n), # k
                   diag(1, nrow = n, ncol = n)) # l
  b_cs_al <- rep(1, n)
  sense_cs_al <- rep("<=", n)

  stopifnot(ncol(A_cs_al) == num_decision_vars)
  stopifnot(nrow(A_cs_al) == n)
  stopifnot(length(b_cs_al) == n)
  stopifnot(length(sense_cs_al) == n)
  msg <- paste("Complementary Slackness for a and l Complete.")
  send_note_if(msg, show_progress, message)

  # Non-negativity and Boundedness Constraints (12) and (15)
  lb <- c(rep(-Inf, p_X), # beta_X
          rep(0, p_Phi),    # beta_Phi_plus
          rep(0, p_Phi),    # beta_Phi_minus
          rep(-Inf, p_D), # beta_D
          rep(0, n),      # u
          rep(0, n),      # v
          rep(0, n),      # a
          rep(0, n),      # k
          rep(0, n))      # l
  ub <- c(rep(Inf, p_X),  # beta_X
          rep(Inf, p_Phi),  # beta_Phi_plus
          rep(Inf, p_Phi),  # beta_Phi_minus
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
             rep("C", p_Phi), # beta_Phi_plus
             rep("C", p_Phi), # beta_Phi_minus
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
  if (!is.null(O) & !VarHintVal_bool) {
    # If a residual is positive, then the associated k must be 1, which means
    # the dual variable, a, must also be 1. Accordingly, the associated l must
    # be 0.
    # If a residual is negative, then the associated l must be 1, which means
    # the dual variable, a, must also be 0. Accordingly, the associated k must
    # be 0.
    fixed <- rep(0, n)
    fixed[O] <- 1
    fixed_mat <- diag(fixed)

    A_pp_a <- cbind(matrix(0, nrow = n, ncol = p_X),    # beta_X
                    matrix(0, nrow = n, ncol = p_Phi),  # beta_Phi_plus
                    matrix(0, nrow = n, ncol = p_Phi),  # beta_Phi_minus
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
    msg <- "Pre-processing for a Complete."
    send_note_if(msg, show_progress, message)

    A_pp_k <- cbind(matrix(0, nrow = n, ncol = p_X),    # beta_X
                    matrix(0, nrow = n, ncol = p_Phi),  # beta_Phi_plus
                    matrix(0, nrow = n, ncol = p_Phi),  # beta_Phi_minus
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
    msg <- "Pre-processing for k Complete."
    send_note_if(msg, show_progress, message)

    A_pp_l <- cbind(matrix(0, nrow = n, ncol = p_X),  # beta_X
                    matrix(0, nrow = n, ncol = p_Phi),  # beta_Phi_plus
                    matrix(0, nrow = n, ncol = p_Phi),  # beta_Phi_minus
                    matrix(0, nrow = n, ncol = p_D),  # beta_D
                    matrix(0, nrow = n, ncol = n),    # u
                    matrix(0, nrow = n, ncol = n),    # v
                    matrix(0, nrow = n, ncol = n),    # a
                    matrix(0, nrow = n, ncol = n),    # k
                    fixed_mat)        # l
    b_l_fixed <- rep(0, n)
    b_l_fixed[O_pos] <- 0
    b_l_fixed[O_neg] <- 1
    b_pp_l <- b_l_fixed
    sense_pp_l <- rep("=", n)

    stopifnot(ncol(A_pp_l) == num_decision_vars)
    stopifnot(nrow(A_pp_l) == n)
    stopifnot(length(b_pp_l) == n)
    stopifnot(length(sense_pp_l) == n)
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

  # Pre-processing alternative: "hint" or prioritize branching
  if (VarHintPri_bool | BranchPriority_bool) {
    # set hint priorities: larger absolute residuals get largest priorities
    # NOTE: the `resid` is obtained from code used in `preprocess_iqr_milp`
    # TODO: remove possible redundancy of computing QR to get residuals between `iqr_milp` and `preprocess_iqr_milp`
    if (p_X == 0) {
      resid <- quantreg::rq(Y ~ D - 1, tau = tau)$residuals
    } else if (p_D == 0) {
      resid <- quantreg::rq(Y ~ X - 1, tau = tau)$residuals
    } else {
      resid <- quantreg::rq(Y ~ X + D - 1, tau = tau)$residuals
    }
    hints_pri <- order(abs(resid))
  }
  if (!is.null(O) & VarHintVal_bool) {
    # use VarHintVal
    a_hints <- rep(NA, n)
    a_hints[O_pos] <- 1
    a_hints[O_neg] <- 0

    k_hints <- rep(NA, n)
    k_hints[O_pos] <- 1
    k_hints[O_neg] <- 0

    l_hints <- rep(NA, n)
    l_hints[O_pos] <- 0
    l_hints[O_neg] <- 1

    VarHintVal <- c(rep(NA, num_decision_vars - 3 * n), a_hints, k_hints, l_hints)

    if (VarHintPri_bool) {
      VarHintPri <- c(rep(0, num_decision_vars - 3 * n), hints_pri, hints_pri, hints_pri)
    } else {
      VarHintPri <- NULL
    }
  } else {
    VarHintVal <- NULL
    VarHintPri <- NULL
  }
  if (!VarHintVal_bool & VarHintPri_bool) {
    message("`VarHintVal_bool` is FALSE; ignoring `VarHintPri_bool`")
  }
  if (!is.null(O) & BranchPriority_bool) {
    BranchPriority <- c(rep(0, num_decision_vars - 3 * n), hints_pri, hints_pri, hints_pri)
  } else {
    BranchPriority <- NULL
  }

  # Putting it all together
  iqr <- list()
  iqr$obj <- obj
  iqr$A <- rbind(A_fix,   # `fix`
                 A_pf,    # Primal Feasibility
                 A_df_X,  # Dual Feasibility - X
                 A_df_Phi,  # Dual Feasibility - Phi
                 A_cs_uk, # Complementary Slackness - u and k
                 A_cs_vl, # Complementary Slackness - v and l
                 A_cs_ak, # Complementary Slackness - a and k
                 A_cs_al, # Complementary Slackness - a and l
                 A_pp_a,  # Pre-processing - fixing a
                 A_pp_k,  # Pre-processing - fixing k
                 A_pp_l)  # Pre-processing - fixing l
  # message(paste("A:", nrow(iqr$A), ncol(iqr$A)))
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
  if (sparse) {
    iqr$A <- as(iqr$A, "sparseMatrix")
  }

  iqr$rhs <- c(b_fix,   # `fix`
               b_pf,    # Primal Feasibility
               b_df_X,  # Dual Feasibility - X
               b_df_Phi,  # Dual Feasibility - Phi
               b_cs_uk, # Complementary Slackness - u and k
               b_cs_vl, # Complementary Slackness - v and l
               b_cs_ak, # Complementary Slackness - a and k
               b_cs_al, # Complementary Slackness - a and l
               b_pp_a,  # Pre-processing - fixing a
               b_pp_k,  # Pre-processing - fixing k
               b_pp_l)  # Pre-processing - fixing l
  # message(paste("b:", length(iqr$rhs)))

  iqr$sense <- c(sense_fix,   # `fix`
                 sense_pf,    # Primal Feasibility
                 sense_df_X,  # Dual Feasibility - X
                 sense_df_Phi,  # Dual Feasibility - Phi
                 sense_cs_uk, # Complementary Slackness - u and k
                 sense_cs_vl, # Complementary Slackness - v and l
                 sense_cs_ak, # Complementary Slackness - a and k
                 sense_cs_al, # Complementary Slackness - a and l
                 sense_pp_a,  # Pre-processing - fixing a
                 sense_pp_k,  # Pre-processing - fixing k
                 sense_pp_l)  # Pre-processing - fixing l
  # message(paste("sense:", length(iqr$sense)))

  iqr$constrnames <- c(
    paste0("fix", seq_along(sense_fix), recycle0 = TRUE),        # sum(not_na)
    paste0("pf", seq_along(sense_pf), recycle0 = TRUE),          # n
    paste0("df_X", seq_along(sense_df_X), recycle0 = TRUE),      # p_X
    paste0("df_Phi", seq_along(sense_df_Phi), recycle0 = TRUE),  # p_Phi
    paste0("cs_uk", seq_along(sense_cs_uk), recycle0 = TRUE),    # n
    paste0("cs_vl", seq_along(sense_cs_vl), recycle0 = TRUE),    # n
    paste0("cs_ak", seq_along(sense_cs_ak), recycle0 = TRUE),    # n
    paste0("cs_al", seq_along(sense_cs_al), recycle0 = TRUE),    # n
    paste0("pp_a", seq_along(sense_pp_a), recycle0 = TRUE),      # n if !is.null(O)
    paste0("pp_k", seq_along(sense_pp_k), recycle0 = TRUE),      # n if !is.null(O)
    paste0("pp_l", seq_along(sense_pp_l), recycle0 = TRUE)       # n if !is.null(O)
  )

  iqr$varnames <- c(
    paste0("beta_X", seq_len(p_X), recycle0 = TRUE),
    paste0("beta_Phi_plus", seq_len(p_Phi), recycle0 = TRUE),
    paste0("beta_Phi_minus", seq_len(p_Phi), recycle0 = TRUE),
    paste0("beta_D", seq_len(p_D), recycle0 = TRUE),
    paste0("u", seq_len(n), recycle0 = TRUE),
    paste0("v", seq_len(n), recycle0 = TRUE),
    paste0("a", seq_len(n), recycle0 = TRUE),
    paste0("k", seq_len(n), recycle0 = TRUE),
    paste0("l", seq_len(n), recycle0 = TRUE)
  )

  for (i in seq_along(attributes)) {
    # TODO: error check the attribute names?
    # TODO: send warning if attribute already exists?
    att_name <- names(attributes)[[i]]
    att_val <- attributes[[att_name]]
    iqr[[att_name]] <- att_val
  }
  if (!is.null(start)) {
    # TODO: error-check starting solution
    # TODO: figure out api to help people specify warmstart solutions
    iqr$start <- start
    if (ncol(iqr$A) != length(iqr$start)) {
      warning(paste("ncol of A:", ncol(iqr$A)))
      warning(paste("length of start:", length(iqr$start)))
      stop("ncol of A doesn't match length of start.")
    }
  }

  iqr$varhintval <- VarHintVal # TODO: send warning if attribute already exists?
  iqr$varhintpri <- VarHintPri
  iqr$branchpriority <- BranchPriority
  iqr$lb <- lb
  iqr$ub <- ub
  iqr$vtype <- vtype
  iqr$modelsense <- "min"
  params$TimeLimit <- TimeLimit
  if (LogFileName != "") {
    params$LogFile <- paste0(LogFileName, LogFileExt)
  }
  result <- gurobi::gurobi(iqr, params)

  msg <- paste("Mixed Integer Linear Program Complete.")
  send_note_if(msg, show_progress, message)

  # Return results
  msg <- paste("Status of IQR program:",
               result$status,
               "| Objective:",
               format(result$objval, scientific = F, digits = 10))
  send_note_if(msg, !quietly, message)  # Print status of program if !quietly

  out$iqr <- iqr
  out$params <- params
  out$attributes <- attributes
  out$result <- result
  out$status <- result$status
  if (result$status %in% c("OPTIMAL", "SUBOPTIMAL")) {
    answer <- result$x
    if (p_X > 0) {
      out$beta_X <- answer[1:p_X]
    } else {
      out$beta_X <- NA
    }
    if (p_Phi > 0) {
      out$beta_Phi_plus <- answer[(p_X + 1):(p_X + p_Phi)]
      out$beta_Phi_minus <- answer[(p_X + p_Phi + 1):(p_X + 2*p_Phi)]
    } else {
      out$beta_Phi_plus <- NA
      out$beta_Phi_minus <- NA
    }
    if (p_D > 0) {
      out$beta_D <- answer[(p_X + 2*p_Phi + 1):(p_X + 2*p_Phi + p_D)]
    } else {
      out$beta_D <- NA
    }
    out$u <- answer[(p_X + 2*p_Phi + p_D + 1):(p_X + 2*p_Phi + p_D + n)]
    out$v <- answer[(p_X + 2*p_Phi + p_D + n + 1):(p_X + 2*p_Phi + p_D + 2*n)]
    out$a <- answer[(p_X + 2*p_Phi + p_D + 2*n + 1):(p_X + 2*p_Phi + p_D + 3*n)]
    out$k <- answer[(p_X + 2*p_Phi + p_D + 3*n + 1):(p_X + 2*p_Phi + p_D + 4*n)]
    out$l <- answer[(p_X + 2*p_Phi + p_D + 4*n + 1):(p_X + 2*p_Phi + p_D + 5*n)]

    out$beta_Phi <- out$beta_Phi_plus - out$beta_Phi_minus
    out$resid <- out$u - out$v
    out$objval <- result$objval
  }

  return(out)
}
