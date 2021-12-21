### Propose h -- "Algorithm 3" -- Algorithm A -------------------------

#' @param beta_proposal Vector of proposed coefficients beta_D and beta_X
#'  (vector of length p_D + p_X)
#' @param beta_hat Vector of beta_D and beta_X coefficients from IQR point
#'  estimate (vector of length p_D + p_X)
#' @param varcov_mat Asymptotic variance-covariance matrix; see
#'  \code{wald_varcov}
#'
#' @return density of limiting distribution of estimator evaluated at
#'  \code{beta_proposal}
density_wald <- function(beta_hat, beta_proposal, varcov_mat) {
  mvnfast::dmvn(beta_proposal, mu = beta_hat, sigma = varcov_mat)
}

#' @param residuals Residuals from IQR point estimation (vector of length n)
#' @param p_design Dimension of design matrix in QR (p_X + p_Phi) (numeric)
#' @param theta Hyperparameter (numeric)
#'
#' @return Named list
#'  \enumerate{
#'    \item \code{draws}: vector of length `n` with 1 if the index is in active
#'      basis and 0 otherwise
#'    \item \code{h_star}: Indices of the active basis
#'  }
propose_active_basis <- function(residuals, p_design, theta = 1) {
  weights <- exp(-1 * theta * residuals^2) # pointwise-operation
  n <- length(residuals)
  # 1 ball of each color in n urns;
  # probability of picking ball in urn i is equal to weights[i]
  # returns a vector with 1 if selected and 0 otherwise 
  # DONE: use `rmultinom`
  # draws <- BiasedUrn::rMFNCHypergeo(nran = 1, m = rep(1, n), n = p_design,
  #                                   odds = weights, precision = 1e-7)
  # Q: is using the while loop okay?
  while_bool <- TRUE
  while (while_bool) {
    # we pick p_design balls from length(weights) bins
    draws <- as.numeric(rmultinom(n = 1, size = p_design, prob = weights))
    # NOTE: possible to pick two balls from the same bin, hence the while loop will ensure all balls are from different bins
    # break out of while loop if we pick at most 1 ball from each bin
    if (identical(sort(unique(draws)), c(0, 1))) {
      while_bool <- FALSE
    }
  }
  h_star <- which(draws == 1)
  list(
    draws = draws,
    h_star = h_star # TODO: do we even need to store this?
  )
}

#' @param active_basis_draws Vector of length n with 1 if index is in active
#'  basis and 0 otherwise; see \code{propose_active_basis}
#' @param residuals Residuals from IQR point estimation (vector of length n)
#' @param p_design Dimension of design matrix in QR (p_X + p_Phi) (numeric)
#' @param theta Hyperparameter (numeric)
#'
#' @return density of Hypergeometric distribution evaluated at
#'  \code{density_active_basis}
density_active_basis <- function(active_basis_draws, residuals, p_design, theta = 1) {
  weights <- exp(-1 * theta * residuals^2)
  # n <- length(residuals)
  # # DONE: replace with product of weights
  # BiasedUrn::dMWNCHypergeo(x = active_basis_draws, m = rep(1, n), n = p_design,
  #                          odds = weights, precision = 1e-7)
  prod(weights[active_basis_draws == 1])
}

#' @param iterations Number of iterations through the MCMC algorithm
#' @param beta_X Coefficients on the exogenous variables (vector of length
#'  p_X); obtained from IQR point estimate
#' @param beta_D Coefficients on the endogeneous variables (vector of length
#'  p_D); obtained from IQR point estimate
#' @param Y Dependent variable (vector of length n)
#' @param X Exogenous variable (including constant vector) (n by p_X matrix)
#' @param D Endogenous variable (n by p_D matrix)
#' @param Z Instrumental variable (n by p_Z matrix)
#' @param Phi Transformation of X and Z to be used in the program;
#'  defaults to the linear projection of D on X and Z (matrix with n rows)
#' @param tau Quantile (numeric)
#' @param alpha Used for Hall and Sheather bandwidth; defaults to 0,1
#' @param theta Hyperparameter (numeric)
#' @param psi Hyperparameter; Coefficient on the variance-covariance matrix; defaults to 1 (numeric)
#' @param varcov_mat variance-covariance matrix of beta_D and beta_X
#' @param discard_burnin If TRUE (default), discard the first set of draws that are equivalent to the IQR MILP estimate
#'
#' @return A named list, each with a data frame:
#'  \enumerate{
#'    \item "beta": data frame where the rows are beta_D and beta_X while the
#'      columns are each iteration in the MCMC (excluding the burn-in period if \code{discard_burnin} is TRUE)
#'    \item "h": data frame where the rows are each entry in the active basis
#'      while the columns are each iteration in the MCMC (excluding the burn-in
#'      period if \code{discard_burnin} is TRUE)
# TODO: right now, I feed in beta_D to figure out h; maybe I should feed in h to figure out beta_D and beta_X? Also, maybe I should create beta_to_h(beta_D, ...) function
mcmc_active_basis <- function(iterations,
                              beta_X, # beta_X IQR point estimates
                              beta_D, # beta_D from IQR point estimates
                              Y, X, D, Z, Phi = linear_projection(D, X, Z),
                              tau,
                              alpha = 0.1,
                              theta = 1,
                              psi = 1, # coefficient on the varcov matrix
                              varcov_mat = NULL,
                              discard_burnin = TRUE) {
  qr <- run_concentrated_qr(beta_D = beta_D, Y = Y, X = X, D = D, Phi = Phi, tau = tau)
  residuals <- qr$residuals
  dual <- qr$dual
  initial_basis <- which(dual > 0 & dual < 1)
  u_vec <- runif(n = iterations) # create uniform RV for acceptance/rejection rule
  initial_draws <- vector("double", length = nrow(Y))
  initial_draws[initial_basis] <- 1
  beta_hat <- c(beta_D, beta_X)
  draws_current <- initial_draws
  beta_current <- beta_hat
  p_design <- ncol(X) + ncol(Phi)
  if (is.null(varcov_mat)) {
    varcov_mat <- wald_varcov(
      resid = residuals,
      alpha = alpha,
      tau = tau,
      D = D,
      X = X,
      Phi = Phi,
      psi = psi
    )$varcov
  }

  result <- vector("list", iterations) # preallocate space to store coefficients
  result_h <- vector("list", iterations) # preallocate space to store active basis
  result_prob <- vector("double", iterations) # store wald density of coefficients
  h_current <- initial_basis # TODO: refactor `draws` and `h`
  for (i in seq_len(iterations)) {
    u <- u_vec[[i]]
    h_proposal <- propose_active_basis(residuals, p_design = p_design, theta = theta)
    draws_proposal <- h_proposal$draws
    beta_proposal_full <- h_to_beta(h = h_proposal$h_star, Y = Y, X = X, D = D, Phi = Phi)
    beta_proposal <- c(beta_proposal_full$beta_D, beta_proposal_full$beta_X)
    # TODO: come up with better variable names
    wald_proposal <- density_wald(beta_hat = beta_hat,
                                  beta_proposal = beta_proposal,
                                  varcov_mat)
    geom_proposal <- density_active_basis(active_basis_draws = draws_proposal,
                                          residuals = residuals,
                                          p_design = p_design,
                                          theta = theta)
    if (i == 1) {
      wald_current <- density_wald(beta_hat = beta_hat,
                                   beta_proposal = beta_current,
                                   varcov_mat)
      geom_current <- density_active_basis(active_basis_draws = draws_current,
                                           residuals = residuals,
                                           p_design = p_design,
                                           theta = theta)
    }
    a <- wald_proposal / wald_current * geom_current / geom_proposal
    if (u < a) { # accept
      beta_current <- beta_proposal
      draws_current <- draws_proposal
      h_current <- h_proposal$h_star
      wald_current <- wald_proposal
      geom_current <- geom_proposal
    }
    result[[i]] <- beta_current
    result_h[[i]] <- h_current
    result_prob[[i]] <- wald_current
  }
  # each row is a coefficient, each column is one iteration of MCMC
  result_df <- do.call(cbind, result)
  rownames(result_df) <- c(paste0("beta_D", seq_len(ncol(D))), paste0("beta_X", seq_len(ncol(X))))
  result_h_df <- do.call(cbind, result_h)
  if (discard_burnin) {
    # find where stationary distribution begins
    stationary_begin <- min(which(result_df[1, ] != beta_hat[1]))
    # remove burn-in period
    result_df <- result_df[, stationary_begin:ncol(result_df)]
    result_h_df <- result_h_df[, stationary_begin:ncol(result_h_df)]
    result_prob <- result_prob[stationary_begin:length(result_prob)]
  }
  list("beta" = result_df, "h" = result_h_df, "prob" = result_prob)
}
