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
  h_aux,
  initial_aux, # initial subsample for aux variables
  distance_function = foc_violation,
  distance_params,
  transform_function = exp_dist,
  transform_params # parameters to pass to `transform_function`
) {
  # Q: proposal of subsamples/aux variables negate each other...right?

  # preliminaries
  xi_mat <- compute_foc_conditions(h, Y = Y, X = X, D = D, Phi = Phi, tau = tau)
  xi_mat_aux <- compute_foc_conditions(h_aux, Y = Y, X = X, D = D, Phi = Phi, tau = tau)

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

  aux_current <- initial_aux
  dist_aux_current <- distance_function(
    h = h_aux,
    subsample = aux_current,
    tau = tau,
    xi_mat = xi_mat_aux,
    params = distance_params
  )
  Q_aux_current <- transform_function(dist_aux_current, transform_params)
  membership_aux_current <- all.equal(dist_aux_current, 0)
  tmp <- distance_function(
    h = h,
    subsample = aux_current,
    tau = tau,
    xi_mat = xi_mat,
    params = distance_params
  )
  P_aux_current <- transform_function(tmp, transform_params)

  # collect draws from uniform distribution
  u_vec <- runif(n = iterations)
  u_vec_aux <- runif(n = iterations)

  # pre-allocate results
  result_record <- vector("double", iterations)
  result_P <- vector("double", iterations)
  result_distance <- vector("double", iterations)
  result_membership <- vector("double", iterations)

  result_record_aux <- vector("double", iterations)
  result_Q_aux <- vector("double", iterations)
  result_P_aux <- vector("double", iterations)
  result_distance_aux <- vector("double", iterations)
  result_membership_aux <- vector("double", iterations)

  for (mcmc_idx in seq_len(iterations)) {
    record <- 0
    record_aux <- 0
    u <- u_vec[[mcmc_idx]]
    u_aux <- u_vec_aux[[mcmc_idx]]

    # Get proposals
    ones <- setdiff(which(D_current == 1), h)
    zeros <- setdiff(which(D_current == 0), h)
    one_to_zero <- sample(ones, 1, replace = FALSE)
    zero_to_one <- sample(zeros, 1, replace = FALSE)
    D_star <- D_current
    D_star[one_to_zero] <- 0
    D_star[zero_to_one] <- 1

    ones <- setdiff(which(aux_current == 1), h_aux)
    zeros <- setdiff(which(aux_current == 0), h_aux)
    one_to_zero <- sample(ones, 1, replace = FALSE)
    zero_to_one <- sample(zeros, 1, replace = FALSE)
    aux_star <- aux_current
    aux_star[one_to_zero] <- 0
    aux_star[zero_to_one] <- 1

    # Compute P_star
    dist_star <- distance_function(
      h = h,
      subsample = D_star,
      tau = tau,
      xi_mat = xi_mat,
      params = distance_params
    )
    P_star <- transform_function(dist_star, transform_params)

    dist_aux_star <- distance_function(
      h = h_aux,
      subsample = aux_star,
      tau = tau,
      xi_mat = xi_mat_aux,
      params = distance_params
    )
    Q_aux_star <- transform_function(dist_aux_star, transform_params)

    # Compute acceptance probabilities and accept/reject
    acc_prob <- P_star / P_current
    if (u < acc_prob) {
      D_current <- D_star
      P_current <- P_star
      dist_current <- dist_star
      record <- 1
      membership_current <- all.equal(dist_current, 0)
    }
    result_record[[mcmc_idx]] <- record
    result_P[[mcmc_idx]] <- P_current
    result_distance[[mcmc_idx]] <- dist_current
    result_membership[[mcmc_idx]] <- membership_current

    acc_prob_aux <- Q_aux_star / Q_aux_current
    if (u_aux < acc_prob_aux) {
      aux_current <- aux_star
      Q_aux_current <- Q_aux_star
      dist_aux_current <- dist_aux_star
      record_aux <- 1
      membership_aux_current <- all.equal(dist_current, 0)

      # compute P(x | beta_star)
      tmp <- distance_function(
        h = h,
        subsample = aux_current,
        tau = tau,
        xi_mat = xi_mat,
        params = distance_params
      )
      P_aux_current <- transform_function(tmp, transform_params)
    }
    result_record_aux[[mcmc_idx]] <- record_aux
    result_Q_aux[[mcmc_idx]] <- Q_aux_current
    result_P_aux[[mcmc_idx]] <- P_aux_current
    result_distance_aux[[mcmc_idx]] <- dist_aux_current
    result_membership_aux[[mcmc_idx]] <- membership_aux_current
  }

  list(
    record = result_record,
    P = result_P,
    distance = result_distance,
    membership = result_membership,
    record_aux = result_record_aux,
    Q_aux = result_Q_aux,
    P_aux = result_P_aux,
    distance_aux = result_distance_aux,
    membership_aux = result_membership_aux,
  )
}
