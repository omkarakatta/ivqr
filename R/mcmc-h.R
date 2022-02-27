# Helpers ----------------------------------------------------------------------

#' Compute variance-covariance matrix of naive endogeneous QR
#'
#' @param psi Coefficient on variance-covariance matrix
#'
#' @return variance-covariance matrix of beta_D and beta_X coefficients when
#' running a naive endogeneous QR of Y on D and X without an intercept
naive_varcov <- function(Y, D, X, tau, psi = 1) {
  qr_naive <- quantreg::rq(Y ~ D + X - 1, tau = tau)
  psi * summary(qr_naive, covariance = TRUE)$cov
}

#' Evaluate Asymptotic Distribution of IQR Estimator
#'
#' @param beta_star Vector of proposed coefficients beta_D and beta_X
#'  (vector of length p_D + p_X)
#' @param beta_opt Vector of beta_D and beta_X coefficients from IQR point
#'  estimate (vector of length p_D + p_X)
#' @param varcov_mat Asymptotic variance-covariance matrix; see
#'  `wald_varcov` and `naive_varcov`
#'
#' @return density of limiting distribution of estimator evaluated at
#'  \code{beta_star}
density_wald <- function(beta_star, beta_opt, varcov_mat) {
  mvnfast::dmvn(beta_star, mu = beta_opt, sigma = varcov_mat)
}

#' Compute weights on indices to propose active basis
#'
#' Indices whose residuals are small are weighted more than indices whose
#' residuals are large.
#'
#' @param residuals_opt Residuals from IQR-QR problem
#' @param theta Tuning parameter
#'
#' @return weight for each index
weight_indices <- function(residuals_opt, theta) {
  exp(-1 * theta * residuals_opt^2)
}

#' Compute unnormalized proposal density of active basis
#'
#' @param weights Weights on each index, see `weight_indices`
#' @param h Indices in the active basis
#'
#' @return unnormalized proposal density of active basis
proposal_h <- function(weights, h) {
  prod(weights[h])
}

#' Propose active basis
#'
#' Indices whose observations have small residuals will be more likely to
#' belong to the active basis, while indices whose observations have extreme
#' residuals will be less likely to belong to the active basis.
#'
#' @param weights Vector of weights for each observation
#' @param theta Tuning parameter
#' @param p Dimension of design matrix in IQR-QR (p_X + p_Phi)
#' @param num_draws Number of draws from proposal distribution of active basis
#'
#' @return indices in the active basis
draw_proposal_h <- function(weights,
                            p,
                            num_draws,
                            label_bool = FALSE,
                            label_skip = 5,
                            existing_h = c()) {
  # Choose `p` balls from `n` urns, where `n := length(residuals)`.
  # We are more likely to pick balls from urns with larger entries in `weights`.
  # It's possible to pick two balls from the same bin.
  # Hence, I will keep drawing `p` balls until all balls are drawn from
  # different bins.
  tmp <- function(x) identical(as.numeric(sort(unique(x))), c(0, 1))
  draws_raw <- rmultinom(n = num_draws, size = p, prob = weights)
  valid <- apply(draws_raw, 2, tmp)
  if (sum(valid) == 0) {
    num_draws_remaining <- num_draws
  } else {
    draws_filter <- draws_raw[, valid, drop = FALSE]
    draws_h <- t(apply(draws_filter, 2, function(draw) which(draw == 1)))
    num_draws_remaining <- num_draws - nrow(draws_h)
    existing_h <- rbind(existing_h, draws_h)
  }
  if (num_draws_remaining > 0) {
    draw_proposal_h(weights, p, num_draws_remaining,
                    label_bool, label_skip,
                    existing_h)
  } else {
    return(existing_h)
  }
}

#' Evaluate target density in the MCMC Sampler of the proposal coefficients
#'
#' See `density_wald`.
#'
#' @inheritParams density_wald
#'
#' @return density of target distribution in MCMC Sampler of the proposal
#'  distribution for the coefficients
target_h <- function(...) {
  density_wald(...)
}


# Main -------------------------------------------------------------------------

#' MH Sampler of the Proposal Distribution of Coefficients/Active Bases
#'
#' The proposal distribution Q_beta is the asymptotic distribution of the IQR
#' estimator limited to the finite support of possible solutions enumerated by
#' the active basis.
#' This MH sampler returns draws from the proposal distribution by
#'
#'  1. Proposing an active basis `h_star` according to Q_h, which puts more
#'      weight on indices with smaller residuals
#'  2. Computing coefficients `beta_star` (see `h_to_beta`)
#'  3. Computing the acceptance probability
#'  4. Accept/Reject in the usual MH manner
#'
#' @param iterations Number of iterations
#' @param initial_h Initial active basis indices
#' @param initial_beta_D,initial_beta_X Initial coefficients, e.g., IQR point
#'  estimate
#' @param residuals_opt Residuals from IQR-QR problem
#' @param varcov_mat Variance-covariance matrix
#' @param theta Tuning parameter
#' @param Y,X,D,Phi Data
#' @param discard_burnin If TRUE, remove first few samples that have the same
#'  coefficients
#'
#' @return named list
#'  1. `beta`: data frame where each row is a vector of coefficients (one
#'        row per iteration)
#'  2. `h`: data frame where each row is a vector of indices in the active
#'        basis (one row per iteration)
#'  3. `record`: binary vector, 1 if that iteration's proposal was accepted
#'  4. `stationary_begin`: all iterations prior to this number were discarded
mcmc_h <- function(
  iterations,
  h_opt,
  beta_D_opt,
  beta_X_opt,
  residuals_opt,
  theta,
  varcov_mat,
  Y, X, D, Phi,
  discard_burnin,
  manual_burnin = 1,
  label_function = function(idx) {
    print(paste("MCMC-H IDX:", idx, "/", iterations))
  },
  label_skip = floor(iterations / 5),
  label_bool = TRUE
) {
  # preliminaries
  p <- length(beta_D_opt) + length(beta_X_opt)
  beta_opt <- c(beta_D_opt, beta_X_opt)
  beta_current <- beta_opt
  h_current <- h_opt
  weights <- weight_indices(residuals_opt, theta)
  log_Q_current <- log(proposal_h(weights, h_current))
  log_P_current <- log(target_h(beta_current, beta_opt, varcov_mat))

  # pre-allocate results
  result_beta <- vector("list", iterations)
  result_h <- vector("list", iterations)
  result_record <- vector("double", iterations)
  result_P <- vector("double", iterations)

  # run MCMC
  u_vec <- runif(n = iterations)
  for (mcmc_idx in seq_len(iterations)) {
    record <- 0
    log_u <- log(u_vec[[mcmc_idx]])


    if (label_bool && mcmc_idx %% label_skip == 0) label_function(mcmc_idx)

    # Step 1: Propose active basis
    # TODO: can we vectorize this as we do with u_vec?
    # `rmultinom` is not vectorized, is it?
    h_star <- as.numeric(draw_proposal_h(weights, p, num_draws = 1)[1, ])
    log_Q_star <- log(proposal_h(weights, h_star))

    # Step 2: Compute coefficients
    beta_full <- h_to_beta(h = h_star, Y = Y, X = X, D = D, Phi = Phi)
    beta_star <- c(beta_full$beta_D, beta_full$beta_X)
    log_P_star <- log(target_h(beta_star, beta_opt, varcov_mat))

    # Step 3: Compute acceptance probability
    log_acc_prob <- log_P_star - log_P_current + log_Q_current - log_Q_star

    # Step 4: Accept/Reject
    if (log_u < log_acc_prob) {
      beta_current <- beta_star
      h_current <- h_star
      log_P_current <- log_P_star
      log_Q_current <- log_Q_star
      record <- 1
    }
    result_beta[[mcmc_idx]] <- beta_current
    result_h[[mcmc_idx]] <- h_current
    result_record[[mcmc_idx]] <- record
    result_P[[mcmc_idx]] <- exp(log_P_current)
  }

  # each row is the information returned by a single iteration of the MCMC
  beta_df <- do.call(rbind, result_beta)
  colnames(beta_df) <- c(
    paste0("beta_D", seq_len(ncol(D))),
    paste0("beta_X", seq_len(ncol(X)))
  )
  h_df <- do.call(rbind, result_h)

  # find where coefficients differ from `beta_opt`, which is the initial beta
  if (discard_burnin) {
    newcoef_burnin <- min(which(beta_df[1, ] != beta_opt[1]))
  } else {
    newcoef_burnin <- 1
  }

  # remove burn-in
  stationary_begin <- max(manual_burnin, newcoef_burnin)
  beta_df <- beta_df[stationary_begin:nrow(beta_df), ]
  h_df <- h_df[stationary_begin:nrow(h_df), ]
  result_record <- result_record[stationary_begin:length(result_record)]
  result_P <- result_P[stationary_begin:length(result_P)]

  list(
    "beta" = beta_df,
    "h" = h_df,
    "record" = result_record,
    "stationary_begin" = stationary_begin,
    "target_density" = result_P
  )
}
