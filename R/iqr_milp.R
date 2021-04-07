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
#' @param tau Quantile (number between 0 and 1)
#' @param O_neg,O_pos Indices for residuals whose sign is fixed to be negative
#'  and positive, respectively (vectors)
#' @param M A large number that bounds the absolute value of the residuals
#' (a positive number); defaults to 10
#' @param TimeLimit Maximum time (in seconds) spent on a linear program;
#'  defaults to 300, will be appended to \code{params}
#' @param projection If TRUE (default), project D on the space spanned by X and
#'  Z to construct the vector of functions of transformed instruments; else,
#'  let Z be the instruments for endogenous variables
#' @param params Gurobi parameters, see \url{https://www.gurobi.com/documentation/9.1/refman/parameter_descriptions.html}
#' @param quietly If TRUE (default), sends messages during execution (boolean)
#'
#' @return A named list of
#'  \enumerate{
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
#'    \item beta_Phi: coefficients on instruments (beta_Phi_plus - beta_Phi_minus)
#'    \item resid: residuals (u - v)
#'    \item objval: value of objective function (absolute value of beta_Phi)
#'  }
iqr_milp <- function(Y,
                     X,
                     D,
                     Z,
                     tau,
                     O_neg,
                     O_pos,
                     M,
                     TimeLimit = 300,
                     projection = TRUE,
                     params = list(FeasibilityTol = 1e-6,
                                   OutputFlag = 0),
                     quietly = TRUE) {

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

  if (projection) {
    # Obtain fitted values from projecting D on space spanned by X and Z
    XZ <- cbind(X, Z)
    proj_matrix <- solve(t(XZ) %*% XZ) %*% t(XZ) %*% D
    Phi <- XZ %*% proj_matrix
  } else {
    Phi <- Z
  }
  p_Phi <- ncol(Phi)

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

  # Objective: minimize absolute value of \beta_Phi
  obj <- c(rep(0, p_X), # beta_X
           rep(1, p_Phi), # beta_Phi_plus
           rep(1, p_Phi), # beta_Phi_minus
           rep(0, p_D), # beta_D
           rep(0, n),   # u
           rep(0, n),   # v
           rep(0, n),   # a
           rep(0, n),   # k
           rep(0, n))   # l

  # Primal Feasibility Constraint (11)
  A_pf <- cbind(X,                  # beta_X
                Phi,                  # beta_Phi_plus
                -Phi,                 # beta_Phi_minus
                D,                  # beta_D
                diag(1, nrow = n),  # u
                -diag(1, nrow = n), # v
                diag(0, nrow = n),  # a
                diag(0, nrow = n),  # k
                diag(0, nrow = n))  # l
  b_pf <- Y
  sense_pf <- rep("=", n)

  # Dual Feasibility Constraint (13) and (14)
  A_df_X <- cbind(diag(0, nrow = p_X),  # beta_X
                  diag(0, nrow = p_X),  # beta_Phi_plus
                  diag(0, nrow = p_X),  # beta_Phi_minus
                  diag(0, nrow = p_X),  # beta_D
                  diag(0, nrow = p_X),  # u
                  diag(0, nrow = p_X),  # v
                  t(X),                 # a
                  diag(0, nrow = p_X),  # k
                  diag(0, nrow = p_X))  # l
  b_df_X <- (1 - tau) * t(X) %*% ones
  sense_df_X <- rep("=", p_X)


  A_df_Phi <- cbind(diag(0, nrow = p_Phi),  # beta_X
                    diag(0, nrow = p_Phi),  # beta_Phi_plus
                    diag(0, nrow = p_Phi),  # beta_Phi_minus
                    diag(0, nrow = p_Phi),  # beta_D
                    diag(0, nrow = p_Phi),  # u
                    diag(0, nrow = p_Phi),  # v
                    t(Phi),             # a
                    diag(0, nrow = p_Phi),  # k
                    diag(0, nrow = p_Phi))  # l
  b_df_Phi <- (1 - tau) * t(Phi) %*% ones
  sense_df_Phi <- rep("=", p_Phi)

  # Complementary Slackness (16) and (17)
  A_cs_uk <- cbind(diag(0, nrow = n),       # beta_X
                   diag(0, nrow = n),       # beta_Phi_plus
                   diag(0, nrow = n),       # beta_Phi_minus
                   diag(0, nrow = n),       # beta_D
                   diag(1, nrow = n),       # u
                   diag(0, nrow = n),       # v
                   diag(0, nrow = n),       # a
                   -M * diag(1, nrow = n),  # k
                   diag(0, nrow = n))       # l
  b_cs_uk <- rep(0, n)
  sense_cs_uk <- rep("<=", n)

  A_cs_vl <- cbind(diag(0, nrow = n),       # beta_X
                   diag(0, nrow = n),       # beta_Phi_plus
                   diag(0, nrow = n),       # beta_Phi_minus
                   diag(0, nrow = n),       # beta_D
                   diag(0, nrow = n),       # u
                   diag(1, nrow = n),       # v
                   diag(0, nrow = n),       # a
                   diag(0, nrow = n),       # k
                   -M * diag(1, nrow = n))  # l
  b_cs_vl <- rep(0, n)
  sense_cs_vl <- rep("<=", n)

  A_cs_ak <- cbind(diag(0, nrow = n),   # beta_X
                   diag(0, nrow = n),   # beta_Phi_plus
                   diag(0, nrow = n),   # beta_Phi_minus
                   diag(0, nrow = n),   # beta_D
                   diag(0, nrow = n),   # u
                   diag(0, nrow = n),   # v
                   diag(1, nrow = n),   # a
                   -diag(1, nrow = n),  # k
                   diag(0, nrow = n))   # l
  b_cs_ak <- rep(0, n)
  sense_cs_ak <- rep(">=", n)

  A_cs_al <- cbind(diag(0, nrow = n), # beta_X
                   diag(0, nrow = n), # beta_Phi_plus
                   diag(0, nrow = n), # beta_Phi_minus
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

  # Integrality Constraint (see vtype) (18)
  vtype <- c(rep("C", p_X), # beta_X
             rep("C", p_Phi), # beta_Phi_plus
             rep("C", p_Phi), # beta_Phi_minus
             rep("C", p_D), # beta_D
             rep("C", n),   # u
             rep("C", n),   # v
             rep("B", n),   # a
             rep("B", n),   # k
             rep("B", n))   # l

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

    A_pp_a <- c(rep(0, p_X),  # beta_X
                rep(0, p_Phi),  # beta_Phi_plus
                rep(0, p_Phi),  # beta_Phi_minus
                rep(0, p_D),  # beta_D
                rep(0, n),    # u
                rep(0, n),    # v
                fixed,        # a
                rep(0, n),    # k
                rep(0, n))    # l
    b_a_fixed <- rep(0, n)
    b_a_fixed[O_pos] <- 1
    b_a_fixed[O_neg] <- 0
    b_pp_a <- c(rep(0, p_X),  # beta_X
                rep(0, p_Phi),  # beta_Phi_plus
                rep(0, p_Phi),  # beta_Phi_minus
                rep(0, p_D),  # beta_D
                rep(0, n),    # u
                rep(0, n),    # v
                b_a_fixed,    # a
                rep(0, n),    # k
                rep(0, n))    # l
    sense_pp_a <- rep("=", n)

    A_pp_k <- c(rep(0, p_X),  # beta_X
                rep(0, p_Phi),  # beta_Phi_plus
                rep(0, p_Phi),  # beta_Phi_minus
                rep(0, p_D),  # beta_D
                rep(0, n),    # u
                rep(0, n),    # v
                rep(0, n),    # a
                fixed,        # k
                rep(0, n))    # l
    b_k_fixed <- rep(0, n)
    b_k_fixed[O_pos] <- 1
    b_k_fixed[O_neg] <- 0
    b_pp_k <- c(rep(0, p_X),  # beta_X
                rep(0, p_Phi),  # beta_Phi_plus
                rep(0, p_Phi),  # beta_Phi_minus
                rep(0, p_D),  # beta_D
                rep(0, n),    # u
                rep(0, n),    # v
                rep(0, n),    # a
                b_k_fixed,    # k
                rep(0, n))    # l
    sense_pp_k <- rep("=", n)

    A_pp_l <- c(rep(0, p_X),  # beta_X
                rep(0, p_Phi),  # beta_Phi_plus
                rep(0, p_Phi),  # beta_Phi_minus
                rep(0, p_D),  # beta_D
                rep(0, n),    # u
                rep(0, n),    # v
                rep(0, n),    # a
                rep(0, n),    # k
                fixed)        # l
    b_l_fixed <- rep(0, n)
    b_l_fixed[O_pos] <- 0
    b_l_fixed[O_neg] <- 1
    b_pp_l <- c(rep(0, p_X),  # beta_X
                rep(0, p_Phi),  # beta_Phi_plus
                rep(0, p_Phi),  # beta_Phi_minus
                rep(0, p_D),  # beta_D
                rep(0, n),    # u
                rep(0, n),    # v
                rep(0, n),    # a
                rep(0, n),    # k
                b_l_fixed)    # l
    sense_pp_l <- rep("=", n)

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

  # Putting it all together
  iqr <- list()
  iqr$obj <- obj
  iqr$A <- rbind(A_pf,    # Primal Feasibility
                 A_df_X,  # Dual Feasibility - X
                 A_df_Phi,  # Dual Feasibility - Phi
                 A_cs_uk, # Complementary Slackness - u and k
                 A_cs_vl, # Complementary Slackness - v and l
                 A_cs_ak, # Complementary Slackness - a and k
                 A_cs_al, # Complementary Slackness - a and l
                 A_pp_a,  # Pre-processing - fixing a
                 A_pp_k,  # Pre-processing - fixing k
                 A_pp_l)  # Pre-processing - fixing l
  iqr$rhs <- rbind(b_pf,    # Primal Feasibility
                   b_df_X,  # Dual Feasibility - X
                   b_df_Phi,  # Dual Feasibility - Phi
                   b_cs_uk, # Complementary Slackness - u and k
                   b_cs_vl, # Complementary Slackness - v and l
                   b_cs_ak, # Complementary Slackness - a and k
                   b_cs_al, # Complementary Slackness - a and l
                   b_pp_a,  # Pre-processing - fixing a
                   b_pp_k,  # Pre-processing - fixing k
                   b_pp_l)  # Pre-processing - fixing l
  iqr$sense <- rbind(sense_pf,    # Primal Feasibility
                     sense_df_X,  # Dual Feasibility - X
                     sense_df_Phi,  # Dual Feasibility - Phi
                     sense_cs_uk, # Complementary Slackness - u and k
                     sense_cs_vl, # Complementary Slackness - v and l
                     sense_cs_ak, # Complementary Slackness - a and k
                     sense_cs_al, # Complementary Slackness - a and l
                     sense_pp_a,  # Pre-processing - fixing a
                     sense_pp_k,  # Pre-processing - fixing k
                     sense_pp_l)  # Pre-processing - fixing l
  iqr$lb <- lb
  iqr$ub <- ub
  iqr$vtype <- vtype
  iqr$modelsense <- "min"
  params$TimeLimit <- TimeLimit
  result <- gurobi::gurobi(iqr, params)

  # Return results
  status <- result$status
  msg <- paste("Status of IQR program:", status)
  send_note_if(msg, !quietly, message)  # Print status of program if !quietly

  out <- list() # Initialize list of results to return
  out$result <- result
  out$status <- status
  if (status %in% c("OPTIMAL", "SUBOPTIMAL")) {
    answer <- result$answer
    out$beta_X <- answer[1:p_X]
    out$beta_Phi_plus <- answer[(p_X + 1):(p_X + p_Phi)]
    out$beta_Phi_minus <- answer[(p_X + p_Phi + 1):(p_X + 2*p_Phi)]
    out$beta_D <- answer[(p_X + 2*p_Phi + 1):(p_X + 2*p_Phi + p_D)]
    out$u <- answer[(p_X + 2*p_Phi + p_D + 1):(p_X + 2*p_Phi + p_D + n)]
    out$v <- answer[(p_X + 2*p_Phi + p_D + n + 1):(p_X + 2*p_Phi + p_D + 2*n)]
    out$a <- answer[(p_X + 2*p_Phi + p_D + 2*n + 1):(p_X + 2*p_Phi + p_D + 3*n)]
    out$k <- answer[(p_X + 2*p_Phi + p_D + 3*n + 1):(p_X + 2*p_Phi + p_D + 4*n)]
    out$l <- answer[(p_X + 2*p_Phi + p_D + 4*n + 1):(p_X + 2*p_Phi + p_D + 5*n)]

    out$beta_Phi <- out$beta_Phi_plus - out$beta_Phi_minus
    out$resid <- out$u - out$v
    out$objval <- answer$objval
  }

  return(out)
}
