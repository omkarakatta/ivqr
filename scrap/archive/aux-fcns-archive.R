### compute_xi_i -------------------------

# TODO: incorporate this function into first_approach* and
# find_subsample_in_polytope
#' @param h Active basis in terms of the data provided
#' @param Y Dependent variable (vector of length n)
#' @param X Exogenous variable (including constant vector) (n by p_X matrix)
#' @param D Endogenous variable (n by p_D matrix)
#' @param Z Instrumental variable (n by p_Z matrix)
#' @param Phi Transformation of X and Z to be used in the program;
#'  defaults to the linear projection of D on X and Z (matrix with n rows)
#' @param tau Quantile (numeric)
#' @param beta_D_proposal Coefficients on the endogeneous variables (vector of
#'  length p_D); if NULL, use \code{h_to_beta} function and the \code{h}
#'  argument to determine \code{beta_D_proposal}
#' @param beta_X_proposal Coefficients on the exogeneous variables (vector of
#'  length p_D); if NULL, use \code{h_to_beta} function and the \code{h}
#'
#' @return Matrix of dimension p by (n-p) of xi_i's not including the i's in the active basis
compute_xi_i <- function(h,
                         Y, X, D, Z, Phi = linear_projection(D, X, Z),
                         tau,
                         beta_D_proposal = NULL,
                         beta_X_proposal = NULL) {
  n <- nrow(Y)
  p <- length(h)

  # get beta_X_proposal and beta_D_proposal
  if (is.null(beta_D_proposal) | is.null(beta_X_proposal)) {
    coef <- h_to_beta(h, Y = Y, X = X, D = D, Phi = Phi)
    if (is.null(beta_D_proposal)) {
      beta_D_proposal <- coef$beta_D
    }
    if (is.null(beta_X_proposal)) {
      beta_X_proposal <- coef$beta_X
    }
  }

  Y_tilde <- Y - D %*% beta_D_proposal
  design <- cbind(X, Phi)
  designh_inv <- solve(design[h, , drop = FALSE])

  s_i <- vector("list", length = n)
  for (i in seq_len(n)) {
    if (is.element(i, h)) {
      s_i[[i]] <- 0 # if index i is in active basis, set xi to be 0
    } else {
      # NOTE: beta_Phi should be 0
      const <- (tau - as.numeric(Y_tilde[i] - X[i, ] %*% beta_X_proposal < 0))
      s_i[[i]] <- const * design[i, ] %*% designh_inv
    }
  }
  s <- t(do.call(rbind, s_i))
  xi_mat <- s[, setdiff(seq_len(n), h), drop = FALSE]
  stopifnot(nrow(xi_mat) == p)
  stopifnot(ncol(xi_mat) == n - p)

  xi_mat
}

