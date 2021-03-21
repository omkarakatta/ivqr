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
#' @param X Exogenous variable (n by p_X matrix)
#' @param D Endogenous variable (n by p_D matrix)
#' @param Z Instrumental variable (n by p_Z matrix)
#' @param tau Quantile (number between 0 and 1)
#' @param O_neg, O_pos Indices for residuals whose sign is fixed to be negative
#'  and positive, respectively (vectors)
#' @param M A large number that bounds the absolute value of the residuals
#' (a positive number); defaults to 10
#' @param params Gurobi parameters, see \link{https://www.gurobi.com/documentation/9.1/refman/parameter_descriptions.html}
iqr_milp <- function(Y,
                     X,
                     D,
                     Z,
                     tau,
                     O_neg,
                     O_pos,
                     M,
                     params = list(TimeLimit = 300,
                                   FeasibilityTol = 1e-6,
                                   OutputFlag = 0)) {

  # Get dimensions of data
  n <- length(Y)
  n_D <- nrow(D)
  n_X <- nrow(X)
  n_Z <- nrow(Z)
  p_D <- ncol(D)
  p_X <- ncol(X)
  p_Z <- ncol(Z)

  # Ensure that number of observations is the same
  stopifnot(all.equal(n, n_D))
  stopifnot(all.equal(n, n_X))
  stopifnot(all.equal(n, n_Z))

  # Create vector of 1s
  ones <- rep(1, n)

  # Find indices of observations whose residuals are not fixed
  O_neg <- sort(O_neg)
  O_pos <- sort(O_pos)
  O <- c(O_neg, O_pos)        # indices of fixed residuals
  I <- setdiff(seq_len(n), O) # indices of free residuals
  n_O <- length(O)            # number of fixed residuals
  n_I <- length(I)            # number of free residuals
  if (!is.null(O)) {
    # If a residual is positive, then the associated k must be 1, which means
    # the dual variable, a, must also be 1. Accordingly, the associated l must
    # be 0.
    # If a residual is negative, then the associated l must be 1, which means
    # the dual variable, a, must also be 1. Accordingly, the associated k must
    # be 0.

    a_fix <- rep(NA_real_, n)
    a_fix[O_pos] <- 1
    a_fix[O_neg] <- 0
    a_O <- a_fix[O]     # binary vector of fixed dual variables, length = n_O

    k_fix <- rep(NA_real_, n)
    k_fix[O_pos] <- 1
    k_fix[O_neg] <- 0
    k_O <- k_fix[O]

    l_fix <- rep(NA_real_, n)
    l_fix[O_pos] <- 0
    l_fix[O_neg] <- 1
    l_O <- l_fix[O]
  }

  # Decision variables in order from left/top to right/bottom:
  # 1. beta_X
  # 2. beta_Z_plus
  # 3. beta_Z_minus ; note: abs(beta_Z) = beta_Z_plus + beta_Z_minus
  # 4. beta_D
  # 5. u
  # 6. v
  # 7. a
  # 8. k
  # 9. l

  # Objective: minimize absolute value of \beta_Z
  obj <- c(rep(0, p_X), # beta_X
           rep(1, p_Z), # beta_Z_plus
           rep(1, p_Z), # beta_Z_minus
           rep(0, p_D), # beta_D
           rep(0, n),   # u
           rep(0, n),   # v
           rep(0, n),   # a
           rep(0, n),   # k
           rep(0, n))   # l

  # Primal Feasibility Constraint (11)
  A_pf <- cbind(X,                  # beta_X
                Z,                  # beta_Z_plus
                -Z,                 # beta_Z_minus
                D,                  # beta_D
                diag(1, nrow = n),  # u
                -diag(1, nrow = n), # v
                diag(0, nrow = n),  # a
                diag(0, nrow = n),  # k
                diag(0, nrow = n))  # l
  b_pf <- Y
  sense_pf <- rep("=", n)

  # Dual Feasibility Constraint (13) and (14)
  A_df_X <- cbind(diag(0, nrow = n),  # beta_X
                  diag(0, nrow = n),  # beta_Z_plus
                  diag(0, nrow = n),  # beta_Z_minus
                  diag(0, nrow = n),  # beta_D
                  diag(0, nrow = n),  # u
                  diag(0, nrow = n),  # v
                  t(X),               # a
                  diag(0, nrow = n),  # k
                  diag(0, nrow = n))  # l
  b_df_X <- (1 - tau) * t(X) %*% ones
  sense_df_X <- rep("=", p_X)


  A_df_Z <- cbind(diag(0, nrow = n),  # beta_X
                  diag(0, nrow = n),  # beta_Z_plus
                  diag(0, nrow = n),  # beta_Z_minus
                  diag(0, nrow = n),  # beta_D
                  diag(0, nrow = n),  # u
                  diag(0, nrow = n),  # v
                  t(Z),               # a
                  diag(0, nrow = n),  # k
                  diag(0, nrow = n))  # l
  b_df_Z <- (1 - tau) * t(Z) %*% ones
  sense_df_Z <- rep("=", p_Z)

  # Complementary Slackness (16) and (17)
  A_cs_uk <- cbind(diag(0, nrow = n),       # beta_X
                   diag(0, nrow = n),       # beta_Z_plus
                   diag(0, nrow = n),       # beta_Z_minus
                   diag(0, nrow = n),       # beta_D
                   diag(1, nrow = n),       # u
                   diag(0, nrow = n),       # v
                   diag(0, nrow = n),       # a
                   -M * diag(1, nrow = n),  # k
                   diag(0, nrow = n))       # l
  b_cs_uk <- rep(0, n)
  sense_cs_uk <- rep("<=", n)

  A_cs_vl <- cbind(diag(0, nrow = n),       # beta_X
                   diag(0, nrow = n),       # beta_Z_plus
                   diag(0, nrow = n),       # beta_Z_minus
                   diag(0, nrow = n),       # beta_D
                   diag(0, nrow = n),       # u
                   diag(1, nrow = n),       # v
                   diag(0, nrow = n),       # a
                   diag(0, nrow = n),       # k
                   -M * diag(1, nrow = n))  # l
  b_cs_vl <- rep(0, n)
  sense_cs_vl <- rep("<=", n)

  A_cs_ak <- cbind(diag(0, nrow = n),   # beta_X
                   diag(0, nrow = n),   # beta_Z_plus
                   diag(0, nrow = n),   # beta_Z_minus
                   diag(0, nrow = n),   # beta_D
                   diag(0, nrow = n),   # u
                   diag(0, nrow = n),   # v
                   diag(1, nrow = n),   # a
                   -diag(1, nrow = n),  # k
                   diag(0, nrow = n))   # l
  b_cs_ak <- rep(0, n)
  sense_cs_ak <- rep(">=", n)

  A_cs_al <- cbind(diag(0, nrow = n), # beta_X
                   diag(0, nrow = n), # beta_Z_plus
                   diag(0, nrow = n), # beta_Z_minus
                   diag(0, nrow = n), # beta_D
                   diag(0, nrow = n), # u
                   diag(0, nrow = n), # v
                   diag(1, nrow = n), # a
                   diag(0, nrow = n), # k
                   diag(1, nrow = n)) # l
  b_cs_al <- rep(0, n)
  sense_cs_al <- rep("<=", n)

  # Non-negativity and Boundedness Constraints (12) and (15)
  lb <- c(rep(-Inf, p_X), # beta_X
          rep(0, p_Z),    # beta_Z_plus
          rep(0, p_Z),    # beta_Z_minus
          rep(-Inf, p_D), # beta_D
          rep(0, n),      # u
          rep(0, n),      # v
          rep(0, n),      # a
          rep(0, n),      # k
          rep(0, n))      # l
  ub <- c(rep(Inf, p_X),  # beta_X
          rep(Inf, p_Z),  # beta_Z_plus
          rep(Inf, p_Z),  # beta_Z_minus
          rep(Inf, p_D),  # beta_D
          rep(Inf, n),    # u
          rep(Inf, n),    # v
          rep(1, n),      # a
          rep(1, n),      # k
          rep(1, n))      # l

  # Integrality Constraint (see vtype) (18)
  vtype <- c(rep("C", p_X), # beta_X
             rep("C", p_Z), # beta_Z_plus
             rep("C", p_Z), # beta_Z_minus
             rep("C", p_D), # beta_D
             rep("C", n),   # u
             rep("C", n),   # v
             rep("B", n),   # a
             rep("B", n),   # k
             rep("B", n))   # l

  # Pre-processing: fix residuals of outliers

}
