### Title: Early-stopping
###
### Description: Solving the IQR MILP to optimality may take some time. Can we
### instead stop the program early and get bounds on the IQR estimator? The
### idea is to solve for the largest and smallest values of each beta_D such
### that the L1 norm of the transformed instruments is less than some epsilon.
### These largest and smallest values of beta_D create an interval that
### contains the true IQR estimator.
### However, this new program is very similar to the original one. So, if we
### are okay with solving this new program, we may as well solve the original
### one to optimality. So instead of running this new program, we will solve the
### continuous relaxation of this program. We also want to find good values of
### `M`, which bounds the magnitude of the residuals. To do this, we will
### compare the residuals of the early-stopping solution with an even earlier
### early-stopping solution.
###
### Author: Omkar A. Katta

### iqr_milp_es -------------------------
# Modifications to iqr_milp:
# 1. change objective to max/min beta_D_j:
#   - beta_D_index
#   - modelsense
# 2. make all decision variables continuous
#   - change vtype
# 3. use individual M_l and M_k
#   - introduce M_k and M_l arguments
#   - change A_uk, A_vl constraints # TODO: do this
# 4. require norm of beta_Phi to be less than something
iqr_milp_es <- function(beta_D_index, # which coefficient are we interested in?
                        bound = 0, # bound on norm of beta_Phi
                        modelsense = "min",
                        Y, X, D, Z,
                        Phi = linear_projection(D, X, Z),
                        tau,
                        O_neg = which(M_l == 0),
                        O_pos = which(M_k == 0),
                        M_k = NULL,
                        M_l = NULL,
                        M = NULL, # used to fill in for NULL M_k and M_l
                        TimeLimit = 300,
                        VarHintVal_bool = FALSE,
                        VarHintPri_bool = FALSE,
                        BranchPriority_bool = FALSE,
                        sparse = TRUE,
                        attributes = list(),
                        params = list(FeasibilityTol = 1e-6,
                                      LogToConsole = 0),
                        start_method = NULL,
                        start = NULL,
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
  n_Phi <- nrow(Phi)
  p_D <- ncol(D)
  p_X <- ncol(X)
  p_Phi <- ncol(Phi)

  # NOTE: modification from iqr_milp
  stopifnot(beta_D_index <= p_D)
  stopifnot(length(M_k) == n)
  stopifnot(length(M_l) == n)
  stopifnot(bound >= 0)
  stopifnot(modelsense %in% c("max", "min"))
  out$beta_D_index <- beta_D_index
  out$bound <- bound
  out$modelsense <- modelsense

  # Ensure that number of observations is the same
  stopifnot(all.equal(n, n_D))
  stopifnot(all.equal(n, n_X))
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

  if (is.null(M_k) | is.null(M_l)) {
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
    if (is.null(M_k)) {
      M_k <- M
    }
    if (is.null(M_l)) {
      M_l <- M
    }
  }
  out$M_k <- M_k
  out$M_l <- M_l


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

  # Objective: minimize / maximize beta_D[beta_D_index]
  # NOTE: modification from iqr_milp
  tmp <- rep(0, p_D)
  tmp[beta_D_index] <- 1
  obj <- c(rep(0, p_X),   # beta_X
           rep(0, p_Phi), # beta_Phi_plus
           rep(0, p_Phi), # beta_Phi_minus
           tmp,           # beta_D
           rep(0, n),     # u
           rep(0, n),     # v
           rep(0, n),     # a
           rep(0, n),     # k
           rep(0, n))     # l
  stopifnot(length(obj) == num_decision_vars)

  # Constrain ||beta_Phi|| < bound
  # NOTE: modification from early-stopping
  A_es <- c(rep(0, p_X),   # beta_X
            rep(1, p_Phi), # beta_Phi_plus
            rep(1, p_Phi), # beta_Phi_minus
            rep(0, p_D),   # beta_D
            rep(0, n),     # u
            rep(0, n),     # v
            rep(0, n),     # a
            rep(0, n),     # k
            rep(0, n))     # l
  b_es <- bound
  sense_es <- "<="

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
    sense_fix <- rep("=", sum(not_na))
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
                   - diag(M_k),                           # k
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
                   - diag(M_l))                           # l
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
             rep("C", n),   # k # NOTE: modificaiton from iqr_milp
             rep("C", n))   # l # NOTE: modificaiton from iqr_milp

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
  iqr$A <- rbind(A_es,    # NOTE: modification from iqr_milp
                 A_fix,   # `fix`
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

  iqr$rhs <- c(b_es,    # NOTE: modification from iqr_milp
               b_fix,   # `fix`
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

  iqr$sense <- c(sense_es,    # NOTE: modification from iqr_milp
                 sense_fix,   # `fix`
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
    paste0("es", seq_along(sense_es), recycle0 = TRUE),          # 1
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
  } else if (!is.null(start_method)) {
    iqr$start <- compute_warmstart(Y = Y,
                                   X = X,
                                   D = D,
                                   Z = Z,
                                   Phi = Phi,
                                   tau = tau,
                                   method = start_method)
  }

  iqr$varhintval <- VarHintVal # TODO: send warning if attribute already exists?
  iqr$varhintpri <- VarHintPri
  iqr$branchpriority <- BranchPriority
  iqr$lb <- lb
  iqr$ub <- ub
  iqr$vtype <- vtype
  iqr$modelsense <- modelsense # NOTE: modification from iqr_milp
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


### iqr_milp_es_gp -------------------------
# from ~/BFI/5_IVQR_GP/fromGuillaume/early_stopping/bounds as preprocessing_OK.R
iqr_milp_es_gp <- function(beta_index = 2,
                           bound = 0,
                           modelsense = "min",
                           Y, D, X, Z,
                           tau = 0.5,
                           O_pos = NULL, O_neg = NULL,
                           cuts = NULL, eps = 1e-14,
                           M_k, M_l, M = 10,
                           start = NULL,
                           TimeLimit = 300, FeasibilityTol = 1e-6, MIPFocus=0, heuristics = 0.05, LogToConsole = 0) {
	
	 #TODO: add condition that makes M_k and M_l all equal to M if they are NULL
	
	p_D = dim(D)[2] ; if(is.null(p_D)){p_D<-0}; if(dim(as.matrix(D))[2]==1){p_D<-1}
	n = dim(X)[1]
	p_X = dim(X)[2]
	p_Z = dim(Z)[2]

	ones = rep(1,n)
	
	O = c(O_neg, O_pos)
	
	if(is.null(O)==F){O <- sort(O)} #put index of pos and neg outlier together
	I = setdiff(1:n, O)
	Iu = setdiff(1:n, O_neg) #no pos error variable u if error is negative
	Iv = setdiff(1:n, O_pos) #no neg error variable v if error is positive

	n_O = length(O)
	n_I = length(I) # number of observations that are neither fixed as positive or negative residuals
	n_Iu = length(Iu)
	n_Iv = length(Iv)

	if(is.null(O)==F){
	a_fix = rep(NA, n); a_fix[O_pos] = 1; a_fix[O_neg] = 0; a_O <- a_fix[O]
	}


	#B_Z_plus, B_Z_min, B_D, B_X, u_Iu, v_Iv, a_I, k_I, l_I
	#obj_core = c(rep(1, 2*p_Z), rep(0, p_D), rep(0, p_X), rep(0, n_Iu+n_Iv), rep(0, n_I), rep(0, n_I+n_I)) # old objective
	beta_D_obj = rep(0, p_D); beta_D_obj[beta_index] <- 1 # new objective
	obj_core = c(rep(0, 2*p_Z), beta_D_obj, rep(0, p_X), rep(0, n_Iu+n_Iv), rep(0, n_I), rep(0, n_I+n_I))
	
	#constraint on objective

	# constrain ||beta_Phi|| = 0 (or some epsilon)
	A_obj = c(rep(1, 2*p_Z), rep(0, p_D), rep(0, p_X), rep(0, n_Iu+n_Iv), rep(0, n_I), rep(0, n_I+n_I))
	b_obj = bound # TODO: add this as an argument to the parameter to allow this to be > 0 so that the feasible set of solutions contain the early-stopping solution
	sense_obj = "<="
	
	#primal feasibility
	A_pf = cbind(Z, -Z, D, X, diag(as.numeric(is.element(1:n, Iu)))[,is.element(1:n, Iu)], -diag(as.numeric(is.element(1:n, Iv)))[,is.element(1:n, Iv)], matrix(0, n, n_I+n_I+n_I)) 	
 	b_pf = Y
 	sense_pf = rep("=", n)
 	
 	#dual feasibility
 	A_df_X = cbind(matrix(0, p_X, 2*p_Z+p_D+p_X), matrix(0, p_X, n_Iu+n_Iv), t(X[I,]), matrix(0, p_X, n_I+n_I)) 
 	if(is.null(O)){ b_df_X = (1-tau)*t(X)%*%ones }  
  	if(is.null(O)==F){ b_df_X = (1-tau)*t(X)%*%ones - t(X[O,])%*%a_O }
  	sense_df_X = rep("=", p_X)
 		
 	A_df_Z = cbind(matrix(0, p_Z, 2*p_Z+p_D+p_X), matrix(0, p_Z, n_Iu+n_Iv), t(Z[I,]), matrix(0, p_Z, n_I+n_I)) 
 	if(is.null(O)){ b_df_Z = (1-tau)*t(Z)%*%ones }
 	if(is.null(O)==F){ b_df_Z = (1-tau)*t(Z)%*%ones - t(Z[O,])%*%a_O }
  	sense_df_Z = rep("=", p_Z)
 			
 	
 	
	#complementary slackness .. do I need the eps? 
	
  # note: all of this is for the n_I observations, not for all n observations
	A_vl_I = cbind(matrix(0, n_I, (2*p_Z + p_D + p_X)), matrix(0,n_I,n_Iu), diag(as.numeric(is.element(Iv, I)))[is.element(Iv, I),], matrix(0, n_I, n_I), matrix(0, n_I, n_I), - diag(M_l[I])*diag(n_I)) #v<=lM # Each v,l has a unique M_l (TODO: remove diag(n_I))
	b_vl_I = rep(0 + eps, n_I)
	sense_vl_I = rep("<=", n_I)   
	
	A_uk_I = cbind(matrix(0, n_I, (2*p_Z + p_D + p_X)), diag(as.numeric(is.element(Iu, I)))[is.element(Iu, I),], matrix(0,n_I,n_Iv), matrix(0, n_I, n_I), - diag(M_k[I])*diag(n_I), matrix(0, n_I, n_I)) #u<=kM # Each u,k has a unique M_l (TODO: remove diag(n_I))
	b_uk_I = rep(0 + eps, n_I)
	sense_uk_I = rep("<=", n_I)  
	
	
	A_ak_I = cbind(matrix(0, n_I, (2*p_Z + p_D + p_X + n_Iu + n_Iv)), diag(n_I), -diag(n_I), matrix(0, n_I, n_I)) #a_i >= k_i 
	b_ak_I = rep(0, n_I) 
	sense_ak_I =  rep(">=", n_I) 

	A_al_I = cbind(matrix(0, n_I, (2*p_Z + p_D + p_X + n_Iu + n_Iv)), diag(n_I), matrix(0, n_I, n_I), diag(n_I)) #a_i <= 1 - l_i 
	b_al_I = rep(1, n_I) 
	sense_al_I =  rep("<=", n_I) 

	
	A_core = rbind(A_obj, A_pf, A_df_X, A_df_Z, A_vl_I, A_uk_I, A_ak_I, A_al_I)
	b_core = c(b_obj, b_pf, b_df_X, b_df_Z, b_vl_I, b_uk_I, b_ak_I, b_al_I)
	sense_core = c(sense_obj, sense_pf, sense_df_X, sense_df_Z, sense_vl_I, sense_uk_I, sense_ak_I, sense_al_I)
	
	#B_Z_plus, B_Z_min, B_D, B_X, u_Iu, v_Iv, a_I, k_I, l_I
	lb_core = c(rep(0, 2*p_Z), rep(-Inf, p_D+p_X), rep(0, n_Iu+n_Iv), rep(0, 3*n_I))
	ub_core = c(rep(Inf, 2*p_Z), rep(Inf, p_D+p_X), rep(Inf, n_Iu+n_Iv), rep(1, 3*n_I))
	vtype_core = "C" #c(rep("C", 2*p_Z+p_D+p_X+n_Iu+n_Iv+n_I), rep("B", 2*n_I)) # we solve continuous relaxation
	
	
	model_core = list()
	model_core$obj = obj_core
	model_core$A = A_core
	model_core$rhs = b_core
	model_core$modelsense = modelsense
	model_core$lb = lb_core
	model_core$ub = ub_core
	model_core$vtype = vtype_core
	model_core$sense = sense_core
	model_core$start = start
	
	
	model <- model_core
	
	 
	params <- list(TimeLimit = TimeLimit, FeasibilityTol = FeasibilityTol, MIPFocus = MIPFocus, heuristics = heuristics, LogToConsole = LogToConsole)
  		
  	out = gurobi::gurobi(model, params) 
  	
  	if(out$status != "TIME_LIMIT"){
  	objval = out$objval
  	x = out$x
  	B_Z_plus = x[1:p_Z]
  	B_Z_min = x[(p_Z+1):(2*p_Z)]
  	if(p_D>0){ B_D = x[(2*p_Z+1):(2*p_Z+p_D)] }; if(p_D==0){ B_D = NULL }
  	B_X = x[(2*p_Z+p_D+1):(2*p_Z+p_D+p_X)]
  	u_Iu = x[(2*p_Z+p_D+p_X+1):(2*p_Z+p_D+p_X+n_Iu)] 
  	v_Iv = x[(2*p_Z+p_D+p_X+n_Iu+1):(2*p_Z+p_D+p_X+n_Iu+n_Iv)] 
  	a_I = x[(2*p_Z+p_D+p_X+n_Iu+n_Iv+1):(2*p_Z+p_D+p_X+n_Iu+n_Iv+n_I)] 
  	k_I = x[(2*p_Z+p_D+p_X+n_Iu+n_Iv+n_I+1):(2*p_Z+p_D+p_X+n_Iu+n_Iv+2*n_I)] 
  	l_I = x[(2*p_Z+p_D+p_X+n_Iu+n_Iv+2*n_I+1):(2*p_Z+p_D+p_X+n_Iu+n_Iv+3*n_I)] 
	#put full dual back together
	if(is.null(O)){a = a_I}
	if(is.null(O)==F){ a_fix[I] = a_I; a <- a_fix }
	
	
	regression_out = list("objval" = objval, "B_Z" =  B_Z_plus-B_Z_min, "B_D" = B_D, "B_X" = B_X, "a" = a, "x" = x, "status" = out$status)
	
	return(regression_out)
	}
	
	if(out$status == "TIME_LIMIT"){
		regression_out = list("status" = out$status, "out" = out)
		return(regression_out)
	}

}
