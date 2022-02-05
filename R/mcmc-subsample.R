# measures norm of the FOC violation in each of the `p` entries of xi_vec
# - tau < x < 1 - tau => - tau - x < 0 & x - (1 - tau) < 0
# => max{...} = 0 => we satisfy FOC
foc_violation <- function(
  h,
  subsample,
  tau,
  xi_mat,
  params # l_norm
) {
  # p by m matrix multiplied by m-vector of 1's
  subsample_indices <- which(subsample > 0)
  xi_vec <- xi_mat[, subsample_indices] %*% rep(1, length(subsample_indices))
  left <- -1 * tau - xi_vec
  right <- xi_vec - (1 - tau)
  violation <- pmax(left, right, rep(0, length(h))) # p by 1
  sum(abs(violation) ^ params$l_norm) ^ (1 / params$l_norm)
}

# transform distance
exp_dist <- function(distance, params) {
  exp(-params$gamma * distance^params$l_power)
}

# run random walk for main and auxiliary variables
rwalk_subsample <- function(
  h,
  Y, X, D, Phi,
  tau,
  initial_subsample, # rounded version of continuous center to FOC
  iterations,
  h_alt = NULL,
  distance_function = foc_violation,
  distance_params,
  transform_function = exp_dist,
  transform_params # parameters to pass to `transform_function`
) {
  # Q: proposal of subsamples/aux variables negate each other, right?

  # Subsamples -----------------------------------------------------------------

  # preliminaries
  xi_mat <- compute_foc_conditions(h,
    Y = Y,
    X = X,
    D = D,
    Phi = Phi,
    tau = tau
  )
  if (!is.null(h_alt)) {
    xi_mat_alt <- compute_foc_conditions(
      h_alt,
      Y = Y,
      X = X,
      D = D,
      Phi = Phi,
      tau = tau
    )
  }

  # Initialize MCMC
  D_current <- initial_subsample
  dist_current <- distance_function(
    h = h,
    subsample = D_current,
    tau = tau,
    xi_mat = xi_mat,
    params = distance_params
  )
  P_current <- transform_function(dist_current, transform_params)
  membership_current <- all.equal(dist_current, 0)

  if (!is.null(h_alt)) {
    dist_alt_current <- distance_function(
      h = h_alt,
      subsample = D_current,
      tau = tau,
      xi_mat = xi_mat_alt,
      params = distance_params
    )
    P_alt_current <- transform_function(dist_alt_current, transform_params)
    membership_alt_current <- all.equal(dist_alt_current, 0)
  } else {
    dist_alt_current <- NA
    P_alt_current <- NA
    membership_alt_current <- NA
  }

  # collect draws from uniform distribution
  u_vec <- runif(n = iterations)

  # pre-allocate results
  result_record <- vector("double", iterations) # accept or reject?
  result_P <- vector("double", iterations) # Q(x | beta_star)
  result_distance <- vector("double", iterations) # FOC(beta_star) violation
  result_membership <- vector("double", iterations) # x \in FOC(beta_star)
  result_P_alt <- vector("double", iterations) # Q(x | beta_hat)
  result_distance_alt <- vector("double", iterations) # FOC(beta_hat) violation
  result_membership_alt <- vector("double", iterations) # x \in FOC(beta_hat)

  for (mcmc_idx in seq_len(iterations)) {
    record <- 0
    u <- u_vec[[mcmc_idx]]

    # Get proposals
    ones <- setdiff(which(D_current == 1), h)
    zeros <- setdiff(which(D_current == 0), h)
    one_to_zero <- sample(ones, 1, replace = FALSE)
    zero_to_one <- sample(zeros, 1, replace = FALSE)
    D_star <- D_current
    D_star[one_to_zero] <- 0
    D_star[zero_to_one] <- 1

    # Compute P_star
    dist_star <- distance_function(
      h = h,
      subsample = D_star,
      tau = tau,
      xi_mat = xi_mat,
      params = distance_params
    )
    P_star <- transform_function(dist_star, transform_params)

    # Compute acceptance probabilities and accept/reject
    acc_prob <- P_star / P_current
    if (u < acc_prob) {
      D_current <- D_star
      P_current <- P_star
      dist_current <- dist_star
      record <- 1
      membership_current <- all.equal(dist_current, 0)

      if (!is.null(h_alt)) {
        dist_alt_current <- distance_function(
          h = h_alt,
          subsample = D_current,
          tau = tau,
          xi_mat = xi_mat_alt,
          params = distance_params
        )
        P_alt_current <- transform_function(dist_alt_current, transform_params)
        membership_alt_current <- all.equal(dist_alt_current, 0)
      }
    }
    result_record[[mcmc_idx]] <- record
    result_P[[mcmc_idx]] <- P_current
    result_distance[[mcmc_idx]] <- dist_current
    result_distance_alt[[mcmc_idx]] <- dist_alt_current
    result_membership[[mcmc_idx]] <- membership_current
    result_membership_alt[[mcmc_idx]] <- membership_alt_current
    result_P_alt[[mcmc_idx]] <- P_alt_current
  }

  list(
    record = result_record,
    P = result_P,
    distance = result_distance,
    membership = result_membership,
    P_alt = result_P_alt,
    distance_alt = result_distance_alt,
  )
}
