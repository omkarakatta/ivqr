# Helpers ----------------------------------------------------------------------

# measures norm of the FOC violation in each of the `p` entries of xi_vec
# - tau < x < 1 - tau => - tau - x < 0 & x - (1 - tau) < 0
# => max{...} = 0 => we satisfy FOC
foc_violation <- function(
  h,
  subsample,
  tau,
  xi_mat,
  xi_vec = NULL,
  minus_index = NULL,
  plus_index = NULL,
  params, # l_norm
  reference_subsample = NULL
) {
  # p by m matrix multiplied by m-vector of 1's
  if (is.null(xi_vec)) {
    xi_vec <- xi_mat %*% as.matrix(subsample, ncol = 1)
  } else {
    xi_vec <- xi_vec - xi_mat[, minus_index] + xi_mat[, plus_index]
  }

  left <- -1 * tau - xi_vec
  right <- xi_vec - (1 - tau)
  violation <- pmax(left, right, rep(0, length(h))) # p by 1
  list(
    dist = sum(abs(violation) ^ params$l_norm) ^ (1 / params$l_norm),
    xi_vec = xi_vec
  )
}

dist_to_reference_subsample <- function(
  h,
  subsample,
  tau,
  xi_mat,
  xi_vec = NULL,
  minus_index = NULL,
  plus_index = NULL,
  params, # l_norm
  reference_subsample = NULL
) {
  list(
    dist = sum(subsample != reference_subsample),
    xi_vec = xi_vec
  )
}

# transform distance
# params$min_prob: smallest acceptable value for the probability (set to -Inf
# (or any negative number) to make this ineffective)
# returns a named list containing:
# 1. result: transformed distance with correction
# 2. mult_correction: multiplicative correction
# 3. dist_transform: transformed distance without correction
exp_dist <- function(distance, params) {
  d <- exp(-1 * params$gamma * distance^params$l_power)
  list(
    result = max(d, params$min_prob),
    mult_correction = 1,
    dist_transform = d
  )
}

# params$min_prob: smallest acceptable value for the probability (set to -Inf
# (or any negative number) to make this ineffective)
# params$n: number of observations in original data
# params$m: number of observations in subsample
# params$p: number of observations in active basis
# 1. result: transformed distance with correction
# 2. mult_correction: multiplicative correction
# 3. dist_transform: transformed distance without correction
exp_dist_correction <- function(distance, params) {
  n <- params$n
  m <- params$m
  p <- params$p
  sharing <- m - distance / 2
  correction <- choose(m - p, sharing - p) * choose(n - m, m - sharing)
  d <- exp(-1 * params$gamma * distance^params$l_power)
  list(
    result = max(1 / correction * d, params$min_prob),
    mult_correction = 1 / correction,
    dist_transform = d
  )
}

# Main -------------------------------------------------------------------------

# run random walk for main and auxiliary variables
#' @param label_function Function of the MCMC index; used with `label_frequency`
#' @param label_skip Call `label_function` every `label_skip` iterations
#' @param label_bool If TRUE, use `label_function` and `label_skip` to keep
#'  track of MCMC
#' @param profile_bool If TRUE, save elapsed time of each step in rwalk
#'  subsample
#' @param save_subsamples If TRUE, save subsamples in a list
rwalk_subsample <- function(
  h,
  Y, X, D, Phi,
  tau,
  initial_subsample, # rounded version of continuous center to FOC
  iterations,
  h_alt = NULL,
  distance_function = foc_violation, # foc_violation, dist_to_reference_subsample #nolint
  distance_params,
  transform_function = exp_dist, # returns named list with `result` key
  transform_params,
  label_function = function(idx) {
    print(paste("RWALK IDX:", idx, "/", iterations))
  },
  label_skip = floor(iterations / 5),
  label_bool = TRUE,
  profile_bool = FALSE,
  save_subsamples = FALSE,
  reference_subsample = NULL
) {
  # Q: proposal of subsamples/aux variables negate each other, right?

  initial_subsample <- as.integer(initial_subsample)
  if (!is.null(reference_subsample)) {
    reference_subsample <- as.integer(reference_subsample)
  }

  if (profile_bool) time <- list()

  # Subsamples -----------------------------------------------------------------

  if (profile_bool) overall_start_time <- Sys.time()

  # preliminaries
  if (profile_bool) start_time <- Sys.time()
  xi_mat <- compute_foc_conditions(
    h,
    Y = Y,
    X = X,
    D = D,
    Phi = Phi,
    tau = tau
  )
  if (profile_bool) time$xi_mat <- difftime(Sys.time(), start_time,
                                            units = "secs")
  if (profile_bool) start_time <- Sys.time()
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
  if (profile_bool) time$xi_mat_alt <- difftime(Sys.time(), start_time,
                                                units = "secs")

  # Initialize MCMC
  if (profile_bool) start_time <- Sys.time()
  D_current <- initial_subsample
  # ensure subsample includes active basis (useful if `initial_subsample` was
  # obtained with `foc_center` where `h` was NULL)
  D_current[h] <- 1
  switch_to_zero <- sum(D_current) - sum(initial_subsample)
  D_current[alt_sample(setdiff(which(D_current == 1), h), switch_to_zero)] <- 0
  dist_current_info <- distance_function(
    h = h,
    subsample = D_current,
    tau = tau,
    xi_mat = xi_mat,
    params = distance_params,
    reference_subsample = reference_subsample
  )
  dist_current <- dist_current_info$dist
  xi_vec_current <- dist_current_info$xi_vec
  # log_P_current <- log(transform_function(dist_current, transform_params)$result) #nolint
  log_P_current <- -Inf # this will ensure we move in the first step # Q: is this okay?
  current_transformation <- list(result = NA,
                                 mult_correction = NA,
                                 dist_transform = NA)
  if (identical(distance_function, foc_violation)) {
    foc_violation_dist <- dist_current
    membership_current <- isTRUE(all.equal(foc_violation_dist, 0))
  } else {
    foc_violation_dist <- foc_violation(
      h = h,
      subsample = D_current,
      tau = tau,
      xi_mat = xi_mat,
      params = distance_params
    )$dist
    membership_current <- isTRUE(all.equal(foc_violation_dist, 0))
  }

  if (profile_bool) time$initialization <- difftime(Sys.time(), start_time,
                                                    units = "secs")

  if (profile_bool) start_time <- Sys.time()
  if (!is.null(h_alt)) {
    dist_alt_current <- distance_function(
      h = h_alt,
      subsample = D_current,
      tau = tau,
      xi_mat = xi_mat_alt,
      params = distance_params,
      reference_subsample = reference_subsample
    )$dist
    P_alt_current <- transform_function(dist_alt_current, transform_params)$result
    if (identical(distance_function, foc_violation)) {
      membership_alt_current <- isTRUE(all.equal(dist_alt_current, 0))
    } else {
      membership_alt_current <- isTRUE(all.equal(foc_violation(
        h = h_alt,
        subsample = D_current,
        tau = tau,
        xi_mat = xi_mat_alt,
        params = distance_params
      )$dist, 0))
    }
  } else {
    dist_alt_current <- NA
    P_alt_current <- NA
    membership_alt_current <- NA
  }
  if (profile_bool) time$initialization_alt <- difftime(Sys.time(), start_time,
                                                        units = "secs")

  n <- length(initial_subsample)
  m <- sum(initial_subsample)
  p <- length(h)

  # collect draws from uniform distribution
  u_vec <- runif(n = iterations)

  # pre-allocate results
  result_record <- vector("double", iterations) # accept or reject?
  result_P <- vector("double", iterations) # Q(x | beta_star)
  result_distance <- vector("double", iterations) # given by distance_function
  result_distance_star <- vector("double", iterations) # proposed subsample dist
  result_membership <- vector("double", iterations) # x \in FOC(beta_star)
  result_foc_violation <- vector("double", iterations) # FOC violation
  result_P_alt <- vector("double", iterations) # Q(x | beta_hat)
  result_distance_alt <- vector("double", iterations) # FOC(beta_hat) violation
  result_membership_alt <- vector("double", iterations) # x \in FOC(beta_hat)
  result_D <- vector("list", iterations) # subsample

  result_P_star <- vector("double", iterations) # Q(x_star | beta_star)
  result_proposed_dir <- vector("character", iterations) # closer, farther, same
  result_star_transformation <- vector("list", iterations)
  result_current_transformation <- vector("list", iterations)

  ones_reference <- setdiff(which(reference_subsample == 1), h)
  zeros_reference <- setdiff(which(reference_subsample == 0), h)

  if (profile_bool) time$iterations <- vector("list", iterations)

  for (mcmc_idx in seq_len(iterations)) {
    record <- 0
    log_u <- log(u_vec[[mcmc_idx]])

    if (label_bool && mcmc_idx %% label_skip == 0) label_function(mcmc_idx)

    ones_current <- setdiff(which(D_current == 1), h)
    zeros_current <- setdiff(which(D_current == 0), h)

    if (!is.null(reference_subsample)) {
      common_ones <- intersect(ones_reference, ones_current) #nolint
      common_zeros <- intersect(zeros_reference, zeros_current) #nolint
      different_ones <- intersect(zeros_reference, ones_current) #nolint
      different_zeros <- intersect(ones_reference, zeros_current) #nolint
      stopifnot(length(common_ones) + length(common_zeros) + length(different_ones) + length(different_zeros) + length(h) == n)
      sharing <- length(intersect(which(reference_subsample == 1), which(D_current == 1))) #nolint
      stopifnot(sharing >= p) # we must always share active basis indices
      stopifnot(all(h %in% intersect(which(reference_subsample == 1), which(D_current == 1)))) # we must always share active basis indices #nolint
      type_vec <- c()
      weight_vec <- c()
      num_closer <- (m - sharing)^2
      num_farther <- (sharing - p) * (n - 2 * m + sharing)
      num_same <- (m - sharing) * (n - 2 * (m + sharing) - p)
      if (num_closer > 0) {
        type_vec <- append(type_vec, "closer")
        weight_vec <- append(weight_vec, 1)
      }
      if (num_farther > 0) {
        type_vec <- append(type_vec, "farther")
        weight_vec <- append(weight_vec, 1)
      }
      if (num_same > 0) {
        type_vec <- append(type_vec, "same")
        weight_vec <- append(weight_vec, 1)
      }
    }

    if (profile_bool) while_start_time <- Sys.time()

    while_counter <- 0
    while_bool <- TRUE
    while (while_bool) {
      while_counter <- while_counter + 1

      if (profile_bool) start_time <- Sys.time()
      # Get proposals
      if (!is.null(reference_subsample)) {
        type <- alt_sample(type_vec, 1, prob = weight_vec)
        result_proposed_dir[[mcmc_idx]] <- type
        if (type == "closer") {
          one_to_zero <- alt_sample(different_ones, 1)
          zero_to_one <- alt_sample(different_zeros, 1)
          log_Q_star <- -log(num_closer)
          log_Q_current <- -log((sharing + 1 - p) * (n - 2 * m + sharing + 1))
        } else if (type == "farther") {
          one_to_zero <- alt_sample(common_ones, 1)
          zero_to_one <- alt_sample(common_zeros, 1)
          log_Q_star <- -log(num_farther)
          log_Q_current <- -log((m - (sharing - 1))^2)
        } else if (type == "same") {
          one_to_zero <- alt_sample(ones_current, 1)
          zero_to_one <- ifelse(one_to_zero %in% common_ones, alt_sample(different_zeros, 1), alt_sample(common_zeros, 1)) #nolint
          log_Q_star <- -log(num_same)
          log_Q_current <- -log(num_same)
        }
      } else {
        type <- NA
        result_proposed_dir[[mcmc_idx]] <- NA
        one_to_zero <- alt_sample(ones_current, 1)
        zero_to_one <- alt_sample(zeros_current, 1)
        log_Q_star <- 0
        log_Q_current <- 0
      }
      D_star <- D_current
      stopifnot(D_current[one_to_zero] == 1)
      stopifnot(D_current[zero_to_one] == 0L)
      D_star[one_to_zero] <- 0L
      D_star[zero_to_one] <- 1L
      stopifnot(sum(D_current) == m)
      stopifnot(sum(D_star) == m)
      stopifnot(sum(D_current != D_star) == 2) # ensure D_star is one-step neighbor
      if (profile_bool) {
        time$iterations[[mcmc_idx]]$get_D_star <- difftime(Sys.time(),
                                                           start_time,
                                                           units = "secs")
      }

      if (profile_bool) start_time <- Sys.time()
      # Compute P_star
      dist_star_info <- distance_function(
        h = h,
        subsample = D_star,
        tau = tau,
        xi_mat = xi_mat,
        xi_vec = xi_vec_current,
        minus_index = one_to_zero,
        plus_index = zero_to_one,
        params = distance_params,
        reference_subsample = reference_subsample
      )
      dist_star <- dist_star_info$dist
      xi_vec_star <- dist_star_info$xi_vec
      sharing_star <- m - dist_star / 2
      star_transformation <- transform_function(dist_star, transform_params)
      log_P_star <- log(star_transformation$result)
      if (profile_bool) {
        time$iterations[[mcmc_idx]]$compute_dist_star <- difftime(Sys.time(),
                                                                  start_time,
                                                                  units = "secs") #nolint
      }
      # go to next iteration of while loop if log_P_star is infinite
      if (is.infinite(log_P_star)) {
        next
      }

      if (profile_bool) start_time <- Sys.time()
      # Compute acceptance probabilities and accept/reject
      log_acc_prob <- log_P_star - log_P_current + log_Q_current - log_Q_star
      accept_reject_bool <- log_u < log_acc_prob
      if (profile_bool) {
        time$iterations[[mcmc_idx]]$compute_acc_prob <- difftime(Sys.time(),
                                                                 start_time,
                                                                 units = "secs")
      }

      # exit while loop if there are no numerical issues
      if (!is.na(accept_reject_bool)) {
        while_bool <- FALSE
      }
    } # exit while loop

    if (profile_bool) {
      time$iterations[[mcmc_idx]]$while_loop <- difftime(Sys.time(),
                                                         while_start_time,
                                                         units = "secs")
      time$iterations[[mcmc_idx]]$while_counter <- while_counter
    }

    if (accept_reject_bool) {
      D_current <- D_star
      log_P_current <- log_P_star
      dist_current <- dist_star
      xi_vec_current <- xi_vec_star
      current_transformation <- star_transformation
      record <- 1
      if (identical(distance_function, foc_violation)) {
        foc_violation_dist <- dist_current
        membership_current <- isTRUE(all.equal(foc_violation_dist, 0))
      } else {
        foc_violation_dist <- foc_violation(
          h = h,
          subsample = D_current,
          tau = tau,
          xi_mat = xi_mat,
          params = distance_params
        )$dist
        membership_current <- isTRUE(all.equal(foc_violation_dist, 0))
      }

      if (!is.null(h_alt)) {
        dist_alt_current <- distance_function(
          h = h_alt,
          subsample = D_current,
          tau = tau,
          xi_mat = xi_mat_alt,
          params = distance_params,
          reference_subsample = reference_subsample
        )$dist
        P_alt_current <- transform_function(dist_alt_current,
                                            transform_params)$result
        if (identical(distance_function, foc_violation)) {
          membership_alt_current <- isTRUE(all.equal(dist_alt_current, 0))
        } else {
          membership_alt_current <- isTRUE(all.equal(foc_violation(
            h = h_alt,
            subsample = D_current,
            tau = tau,
            xi_mat = xi_mat_alt,
            params = distance_params
          )$dist, 0))
        }
      }
    }
    result_record[[mcmc_idx]] <- record
    result_P[[mcmc_idx]] <- exp(log_P_current)
    result_P_star[[mcmc_idx]] <- exp(log_P_star)
    result_distance[[mcmc_idx]] <- dist_current
    result_current_transformation[[mcmc_idx]] <- current_transformation
    result_distance_star[[mcmc_idx]] <- dist_star
    result_star_transformation[[mcmc_idx]] <- star_transformation
    result_distance_alt[[mcmc_idx]] <- dist_alt_current
    result_membership[[mcmc_idx]] <- membership_current
    result_membership_alt[[mcmc_idx]] <- membership_alt_current
    result_foc_violation[[mcmc_idx]] <- foc_violation_dist
    result_P_alt[[mcmc_idx]] <- P_alt_current
    if (save_subsamples) {
      result_D[[mcmc_idx]] <- D_current
    }
  } # exit for loop

  if (profile_bool) time$overall <- difftime(Sys.time(), overall_start_time,
                                             units = "secs")

  if (!profile_bool) time <- NA

  if (!save_subsamples) result_D <- NA

  list(
    record = result_record,
    P = result_P,
    P_star = result_P_star,
    proposed_dir = result_proposed_dir,
    distance = result_distance,
    distance_star = result_distance_star,
    star_transformation = result_star_transformation,
    current_transformation = result_current_transformation,
    membership = result_membership,
    foc_violation = result_foc_violation,
    subsamples = result_D,
    P_alt = result_P_alt,
    distance_alt = result_distance_alt,
    profile = time
  )
}
