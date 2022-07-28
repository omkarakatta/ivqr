### Propose first subsample -- "First Approach" -------------------------

### first_approach -------------------------

#' Propose subsamples
#'
#' Propose observations, one at a time, until we have enough to create a subsample.
#' Note that observations inside the active basis will be inside the final subsample.
#'
#' @param Y Dependent variable (vector of length n)
#' @param X Exogenous variable (including constant vector) (n by p_X matrix)
#' @param D Endogenous variable (n by p_D matrix)
#' @param Z Instrumental variable (n by p_Z matrix)
#' @param tau Quantile (numeric)
#' @param h Indices of active basis (vector of length p_X + p_Phi)
#' @param subsample_size Size of subsample (numeric at most n)
#' @param beta_D_proposal Coefficients on the endogeneous variables (vector of
#'  length p_D); if NULL, use \code{h_to_beta} function and the \code{h}
#'  argument to determine \code{beta_D_proposal}
#' @param beta_X_proposal Coefficients on the exogeneous variables (vector of
#'  length p_D); if NULL, use \code{h_to_beta} function and the \code{h}
#'  argument to determine \code{beta_X_proposal}
#' @param gamma,l_norm,l_power Hyperparameters
#'
#' @return Named list
#'  \enumerate{
#'    \item \code{subsample_set}: set of indices in subsample
#'    \item \code{subsample_weights}: set of weights for each observation in subsample
#'    \item \code{prob}: probability of proposing \code{subsample_set} (except
#'      for observations already in active basis)
#'    \item \code{log_prob}: log of \code{prob}
#'    \item \code{xi}: xi vector for this subsample
#'  }
# Q: Should the beta_*_proposal correspond to the same coefficients as h_to_beta(h)? If so, I don't even need the beta_*_proposal in the arguments. I can just use the h_to_beta(h) to get these coefficients...right?
first_approach <- function(Y, X, D, Z, Phi = linear_projection(D, X, Z), tau,
                           h, subsample_size,
                           beta_D_proposal = NULL, beta_X_proposal = NULL, 
                           gamma = 1, l_norm = 1, l_power = l_norm) {
  n <- length(Y)
  p_design <- length(h)
  subsample_set <- h # store indices of subsample
  # Q: is this correct? do we need to run a regression to get beta_X_proposal?
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
  s <- do.call(rbind, s_i)

  ones <- matrix(1, nrow = 1, ncol = length(subsample_set))
  sum_across_subsample_set <- ones %*% s[subsample_set, , drop = FALSE]
  subsample_weights <- vector("double", subsample_size - length(h))

  # TODO: remove `alt_subsample_weights` after debugging
  # alt_subsample_weights <- vector("double", subsample_size - length(h))

  # we have length(h) observations in subsample; we need subsample_size - length(h) more
  for (j in seq_len(subsample_size - length(h))) {
    choices <- setdiff(seq_len(nrow(Y)), subsample_set)
    n_choices <- length(choices)
    s_remaining <- s[choices, , drop = FALSE]

    # each row is one observation that we can choose from;
    # the value is the weight in the exponent
    ones <- matrix(1, nrow = n_choices, ncol = 1)
    sum_remaining <- kronecker(ones, sum_across_subsample_set) + s_remaining

    # for each row, apply e^(-gamma * (l-norm^l))
    raw_weights <- apply(sum_remaining, 1, function(x) {
      tmp <- sum(abs(x)^l_norm) ^ (1 / l_norm)
      exp(-gamma * tmp^l_power)
    })
    # FIX: what if raw_weights are 0? then total_weights = 0, so then we get 0 / 0!!!
    # We still raise an error in rmultinom.
    # But if we were to try repeating this in the main MCMC, we would still the same weights...
    # So, we need to repeat the mcmc_active_basis to get a new proposal of the active basis and coefficients?
    total_weights <- sum(unlist(raw_weights))
    weights <- raw_weights / total_weights

    # alt_weights <- apply(sum_remaining, 1, function(x) {
    #   # - tau < x < 1 - tau => - tau - x < 0 & x - (1 - tau) < 0
    #   # => max{...} = 0 => we satisfy FOC
    #   # Put differently:
    #   # x is a 1 by p vector;
    #   # So, x - (t - tau) and (- tau - x) are both 1 by p vectors.
    #   # If each entry of these 1 by p vectors are negative, we satisfy FOC conditions
    #   # Fix: `max` must be used element-wise!!! # Q: ask GP if this is right
    #   left <- - tau - x
    #   right <- x - (1 - tau)
    #   max_result <- vector("double", length(x))
    #   for (entry in seq_along(left)) {
    #     left_entry <- left[[entry]]
    #     right_entry <- right[[entry]]
    #     max_result[[entry]] <- max(left_entry, right_entry, 0)
    #   }
    #   tmp <- sum(abs(max_result)^l_norm) ^ (1 / l_norm)
    #   exp(-gamma * tmp^l_power)
    # })

    # choose 1 element in a single vector of size `length(weights)` to be 1
    winner <- tryCatch({
      list(
        status = "OKAY",
        answer = which(rmultinom(n = 1, size = 1, prob = weights) == 1)
      )
    }, error = function(e) {
      list(
        status = "ERROR",
        status_message = e,
        sum_remaining = sum_remaining,
        problem_weights = weights # return problematic weights
      )
    })
    if (winner$status == "ERROR") {
      return(winner)
    }
    winner <- winner$answer
    subsample_weights[[j]] <- weights[winner]
    # alt_subsample_weights[[j]] <- alt_weights[winner]
    new_observation <- choices[winner]
    # TODO: store new_observation in new vector; then append to subsample_set after for loop
    subsample_set <- c(subsample_set, new_observation)
    sum_across_subsample_set <- sum_across_subsample_set + s[new_observation, , drop = FALSE]
  }
  prob <- prod(subsample_weights)
  log_prob <- sum(log(subsample_weights))
  list(
    status = "OKAY",
    status_message = "OKAY",
    prob = prob, # return unnormalized probability of creating subsample
    log_prob = log_prob, # return log of unnormalized probability of creating subsample
    subsample_set = subsample_set, # return set of indices to create subsample!
    # alt_subsample_weights = alt_subsample_weights, # return alternative weights
    subsample_weights = subsample_weights, # return weights for each observation in subsample
    # s = s, # return each entry of term of xi for all observations, not just those in the subsample
    xi = sum_across_subsample_set # return xi object # TODO: double-check that this is indeed the xi object
  )
}

### first_approach_v2 -------------------------

#' Propose subsamples
#'
#' Propose `subsample_size - (p_X + p_Phi)` observations all at once to create a subsample.
#' Note that observations inside the active basis will be inside the final subsample.
#' So, the final subsample will be of size `subsample_size`.
#'
#' @param Y Dependent variable (vector of length n)
#' @param X Exogenous variable (including constant vector) (n by p_X matrix)
#' @param D Endogenous variable (n by p_D matrix)
#' @param Z Instrumental variable (n by p_Z matrix)
#' @param tau Quantile (numeric)
#' @param h Indices of active basis (vector of length p_X + p_Phi)
#' @param subsample_size Size of subsample (numeric at most n)
#' @param beta_D_proposal Coefficients on the endogeneous variables (vector of
#'  length p_D); if NULL, use \code{h_to_beta} function and the \code{h}
#'  argument to determine \code{beta_D_proposal}
#' @param beta_X_proposal Coefficients on the exogeneous variables (vector of
#'  length p_D); if NULL, use \code{h_to_beta} function and the \code{h}
#'  argument to determine \code{beta_X_proposal}
#' @param gamma,l_norm,l_power Hyperparameters
#'
#' @return Named list
#'  \enumerate{
#'    \item \code{subsample_set}: set of indices in subsample
#'    \item \code{subsample_weights}: set of weights for each observation in subsample
#'    \item \code{prob}: probability of proposing \code{subsample_set} (except
#'      for observations already in active basis)
#'    \item \code{log_prob}: log of \code{prob}
#'    \item \code{xi}: xi vector for this subsample
#'  }
first_approach_v2 <- function(Y, X, D, Z, Phi = linear_projection(D, X, Z), tau,
                              h, subsample_size,
                              beta_D_proposal = NULL, beta_X_proposal = NULL, 
                              gamma = 1, l_norm = 1, l_power = l_norm) {
  n <- length(Y)
  p_design <- length(h)
  subsample_set <- h # store indices of subsample
  # Q: is this correct? do we need to run a regression to get beta_X_proposal?
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
  s <- do.call(rbind, s_i)

  ones <- matrix(1, nrow = 1, ncol = length(subsample_set))
  sum_across_subsample_set <- ones %*% s[subsample_set, , drop = FALSE] # this should be 0
  subsample_weights <- vector("double", subsample_size - length(h))

  choices <- setdiff(seq_len(nrow(Y)), subsample_set)
  n_choices <- length(choices)
  s_remaining <- s[choices, , drop = FALSE]

  # each row is one observation that we can choose from;
  # the value is the weight in the exponent
  ones <- matrix(1, nrow = n_choices, ncol = 1)
  sum_remaining <- kronecker(ones, sum_across_subsample_set) + s_remaining
  raw_weights <- apply(sum_remaining, 1, function(x) {
    # left <- - tau - x
    # right <- x - (1 - tau)
    # max_result <- vector("double", length(x))
    # for (entry in seq_along(left)) {
    #   left_entry <- left[[entry]]
    #   right_entry <- right[[entry]]
    #   max_result[[entry]] <- max(left_entry, right_entry, 0)
    # }
    # tmp <- sum(abs(max_result)^l_norm) ^ (1 / l_norm)
    # exp(-gamma * tmp^l_power)
    tmp <- sum(abs(x)^l_norm) ^ (1 / l_norm)
    exp(-gamma * tmp^l_power)
    # 0.99 ^ (-gamma * tmp^l_power)
  })

  # FIX: what if raw_weights are 0? then total_weights = 0, so then we get 0 / 0!!!
  # We still raise an error in rmultinom.
  # But if we were to try repeating this in the main MCMC, we would still the same weights...
  # So, we need to repeat the mcmc_active_basis to get a new proposal of the active basis and coefficients?
  total_weights <- sum(unlist(raw_weights))
  weights <- raw_weights / total_weights

  # choose remaining observations
  winner <- tryCatch({
    # make sure we don't draw any observation more than once
    while_bool <- TRUE
    while_counter <- 0
    while (while_bool) {
      while_counter <- while_counter + 1
      # print(while_counter)
      draws <- as.numeric(rmultinom(
        n = 1, size = subsample_size - length(h), prob = weights
      ))
      if (identical(sort(unique(draws)), c(0, 1))) {
        while_bool <- FALSE
      }
    }
    list(
      status = "OKAY",
      answer = which(draws == 1)
    )
  }, error = function(e) {
    list(
      status = "ERROR",
      status_message = e,
      sum_remaining = sum_remaining,
      problem_weights = weights # return problematic weights
    )
  })
  if (winner$status == "ERROR") {
    return(winner)
  }

  winner <- winner$answer
  new_observations <- choices[winner]
  subsample_set <- c(subsample_set, new_observations)
  sum_across_subsample_set <- sum_across_subsample_set + matrix(1, nrow = 1, ncol = length(new_observations)) %*% s[new_observations, , drop = FALSE]

  subsample_weights <- weights[winner]
  prob <- prod(subsample_weights)
  log_prob <- sum(log(subsample_weights))

  list(
    status = "OKAY",
    status_message = "OKAY",
    prob = prob, # return unnormalized probability of creating subsample
    log_prob = log_prob, # return log of unnormalized probability of creating subsample
    subsample_set = subsample_set, # return set of indices to create subsample!
    subsample_weights = subsample_weights, # return weights for each observation in subsample
    xi = sum_across_subsample_set # return xi object # TODO: double-check that this is indeed the xi object
  )
}

### first_approach_v4 -------------------------

#' Propose subsamples
#'
#' Propose observations, one at a time, until we have enough to create a subsample.
#' Note that observations inside the active basis will be inside the final subsample.
#' Unlike the original \code{first_approach}, I am replacing the exponential
#' function with the reciprocal.
#'
#' @param Y Dependent variable (vector of length n)
#' @param X Exogenous variable (including constant vector) (n by p_X matrix)
#' @param D Endogenous variable (n by p_D matrix)
#' @param Z Instrumental variable (n by p_Z matrix)
#' @param tau Quantile (numeric)
#' @param h Indices of active basis (vector of length p_X + p_Phi)
#' @param subsample_size Size of subsample (numeric at most n)
#' @param beta_D_proposal Coefficients on the endogeneous variables (vector of
#'  length p_D); if NULL, use \code{h_to_beta} function and the \code{h}
#'  argument to determine \code{beta_D_proposal}
#' @param beta_X_proposal Coefficients on the exogeneous variables (vector of
#'  length p_D); if NULL, use \code{h_to_beta} function and the \code{h}
#'  argument to determine \code{beta_X_proposal}
#' @param gamma,l_norm,l_power Hyperparameters
#'
#' @return Named list
#'  \enumerate{
#'    \item \code{subsample_set}: set of indices in subsample
#'    \item \code{subsample_weights}: set of weights for each observation in subsample
#'    \item \code{prob}: probability of proposing \code{subsample_set} (except
#'      for observations already in active basis)
#'    \item \code{log_prob}: log of \code{prob}
#'    \item \code{xi}: xi vector for this subsample
#'  }
first_approach_v4 <- function(Y, X, D, Z, Phi = linear_projection(D, X, Z), tau,
                           h, subsample_size,
                           beta_D_proposal = NULL, beta_X_proposal = NULL, 
                           gamma = 1, l_norm = 1, l_power = l_norm) {
  n <- length(Y)
  p_design <- length(h)
  subsample_set <- h # store indices of subsample
  # Q: is this correct? do we need to run a regression to get beta_X_proposal?
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
  s <- do.call(rbind, s_i)

  ones <- matrix(1, nrow = 1, ncol = length(subsample_set))
  sum_across_subsample_set <- ones %*% s[subsample_set, , drop = FALSE]
  subsample_weights <- vector("double", subsample_size - length(h))

  # TODO: remove `alt_subsample_weights` after debugging
  # alt_subsample_weights <- vector("double", subsample_size - length(h))

  # we have length(h) observations in subsample; we need subsample_size - length(h) more
  for (j in seq_len(subsample_size - length(h))) {
    choices <- setdiff(seq_len(nrow(Y)), subsample_set)
    n_choices <- length(choices)
    s_remaining <- s[choices, , drop = FALSE]

    # each row is one observation that we can choose from;
    # the value is the weight in the exponent
    ones <- matrix(1, nrow = n_choices, ncol = 1)
    sum_remaining <- kronecker(ones, sum_across_subsample_set) + s_remaining

    # for each row, apply e^(-gamma * (l-norm^l))
    raw_weights <- apply(sum_remaining, 1, function(x) {
      tmp <- sum(abs(x)^l_norm) ^ (1 / l_norm)
      gamma / tmp^l_power
      # exp(-gamma * tmp^l_power)
    })
    # FIX: what if raw_weights are 0? then total_weights = 0, so then we get 0 / 0!!!
    # We still raise an error in rmultinom.
    # But if we were to try repeating this in the main MCMC, we would still the same weights...
    # So, we need to repeat the mcmc_active_basis to get a new proposal of the active basis and coefficients?
    total_weights <- sum(unlist(raw_weights))
    weights <- raw_weights / total_weights

    # alt_weights <- apply(sum_remaining, 1, function(x) {
    #   # - tau < x < 1 - tau => - tau - x < 0 & x - (1 - tau) < 0
    #   # => max{...} = 0 => we satisfy FOC
    #   # Put differently:
    #   # x is a 1 by p vector;
    #   # So, x - (t - tau) and (- tau - x) are both 1 by p vectors.
    #   # If each entry of these 1 by p vectors are negative, we satisfy FOC conditions
    #   # Fix: `max` must be used element-wise!!! # Q: ask GP if this is right
    #   left <- - tau - x
    #   right <- x - (1 - tau)
    #   max_result <- vector("double", length(x))
    #   for (entry in seq_along(left)) {
    #     left_entry <- left[[entry]]
    #     right_entry <- right[[entry]]
    #     max_result[[entry]] <- max(left_entry, right_entry, 0)
    #   }
    #   tmp <- sum(abs(max_result)^l_norm) ^ (1 / l_norm)
    #   exp(-gamma * tmp^l_power)
    # })

    # choose 1 element in a single vector of size `length(weights)` to be 1
    winner <- tryCatch({
      list(
        status = "OKAY",
        answer = which(rmultinom(n = 1, size = 1, prob = weights) == 1)
      )
    }, error = function(e) {
      list(
        status = "ERROR",
        status_message = e,
        sum_remaining = sum_remaining,
        problem_weights = weights # return problematic weights
      )
    })
    if (winner$status == "ERROR") {
      return(winner)
    }
    winner <- winner$answer
    subsample_weights[[j]] <- weights[winner]
    # alt_subsample_weights[[j]] <- alt_weights[winner]
    new_observation <- choices[winner]
    # TODO: store new_observation in new vector; then append to subsample_set after for loop
    subsample_set <- c(subsample_set, new_observation)
    sum_across_subsample_set <- sum_across_subsample_set + s[new_observation, , drop = FALSE]
  }
  prob <- prod(subsample_weights)
  log_prob <- sum(log(subsample_weights))
  list(
    status = "OKAY",
    status_message = "OKAY",
    prob = prob, # return unnormalized probability of creating subsample
    log_prob = log_prob, # return log of unnormalized probability of creating subsample
    subsample_set = subsample_set, # return set of indices to create subsample!
    # alt_subsample_weights = alt_subsample_weights, # return alternative weights
    subsample_weights = subsample_weights, # return weights for each observation in subsample
    # s = s, # return each entry of term of xi for all observations, not just those in the subsample
    xi = sum_across_subsample_set # return xi object # TODO: double-check that this is indeed the xi object
  )
}

### random_walk_subsample -------------------------

#' @param initial_subsample A vector of length n with m 1's and n-m 0's; 1
#'  means that the corresponding index is in the subsample
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
#'  argument to determine \code{beta_X_proposal}
#' @param iter How many iterations in the for loop?
#' @param gamma,l_norm,l_power Tuning parameters
#' @param k How many observations do we want to remove/add to the subsample in
#'  the random walk? Larger k means larger jumps; only valid if \code{k_method}
#'  is "constant"
#' @param k_method If "constant" (default), we set \code{k} to be user-specified
#' @param distance_method How do we compute the distance to the FOC polytope?
#'  If it is 1 (default), we use the norm sum of the xi's.
#'  If it is 2, we use the normed difference of the current subsample in the
#'  random walk and \code{reference_subsample}.
#'  If it is 3, we use the transport map idea! TODO: describe this!!!
#'  If it is "simple_random_walk", then we run a simple random walk without any
#'  MCMC-related business. # TODO: test this
#'  If it is 4, we use the transport map idea with violation of the FOC
#'  condition.! TODO: describe this!!!
#'  If it is 5, we use the violation of the FOC constraints.
#' @param reference_subsample If \code{distance_method} is 2, we compare the
#'  subsamples in the random walk to this reference subsample to determine the
#'  distance; only valid if \code{distance_method} is 2; default is NULL; The
#'  intention was that this reference would be a subsample inside the global
#'  FOC polytope.
#' @param transform_method If "exp", use the exponential as the target
#'  distribution
#' @param s_i Matrix of dimension n - p that contains the xi_i_opts that are
#'  mapped from the xi_i_star according to the optimal transport map; only valid for distance_method == 3
#' @param seed For replicability
#' @param save_memory Set to TRUE to save memory by avoiding storage of subsamples
# TODO: test that we always accept the subsample if distance_method = "simple_random_walk"
# TODO: create a separate function just to do the simple random walk without
# any of the extra bells and whistles of MCMC
random_walk_subsample <- function(initial_subsample,
                                  h,
                                  Y, X, D, Z, Phi = linear_projection(D, X, Z),
                                  tau,
                                  beta_D_proposal = NULL,
                                  beta_X_proposal = NULL,
                                  iter = 1000,
                                  gamma = 1, l_norm = 1, l_power = 1,
                                  k,
                                  k_method = "constant",
                                  distance_method = 1,
                                  transform_method = "exp",
                                  reference_subsample = NULL, # for distance_method = 2
                                  s_i,
                                  seed = Sys.date(),
                                  save_memory = FALSE) {
  # k must be smaller than the number of observations in subsample minus the
  # observations in the active basis
  stopifnot(sum(initial_subsample) - length(h) > k)
  # ensure observations in active basis are in the initial subsample
  # stopifnot(all(initial_subsample[h] == 1)) # if initial_subsample is rounded,
  # we don't want to check for this.
  if (distance_method == 2) {
    stopifnot(!is.null(reference_subsample))
    stopifnot(length(reference_subsample) == length(initial_subsample))
  }

  # get beta_X_proposal and beta_D_proposal
  if (distance_method == 1 | distance_method == 5) {
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
  }

  n <- length(initial_subsample)
  m <- sum(initial_subsample)
  current_subsample <- initial_subsample
  current_prob <- 1

  # define how we transform the distance into a weight/unnormalized probability
  if (transform_method == "exp") {
    transform_function <- function(x) {
      exp(-gamma * x^l_power)
    }
  } else if (transform_method == "rec") {
    transform_function <- function(x) {
      gamma / (x^l_power)
    }
  }

  # compute distance of initial_subsample
  if (distance_method == 1) {
    # TODO: is it possible to use compute_xi_i?
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
    s <- do.call(rbind, s_i)
    distance_function <- function(x) {
      sum(abs(x) ^ l_norm) ^ (1 / l_norm)
    }
    curr_s <- s[current_subsample == 1, ]
    current_distance <- distance_function(matrix(1, nrow = 1, ncol = nrow(curr_s)) %*% curr_s)
    current_distance_prob <- transform_function(current_distance)
  } else if (distance_method == 2) {
    # find norm of global subsample and current subsample
    distance_function <- function(x) {
      sum(abs(x - reference_subsample)^l_norm) ^ (1 / l_norm)
    }
    current_distance <- distance_function(current_subsample)
    current_distance_prob <- transform_function(current_distance)
  } else if (distance_method == "simple_random_walk") {
    current_distance <- NA
    current_distance_prob <- 1
  } else if (distance_method == 3) {
    distance_function <- function(x) {
      # sum the column vectors of x
      sum_cols <- x %*% rep(1, ncol(x))
      # take norm of sum of s_i's
      sum(abs(sum_cols)^l_norm) ^ (1 / l_norm)
    }
    okay <- setdiff(seq_len(n), h)
    ones_current <- which(current_subsample[okay] == 1)
    s_i_current <- s_i[, ones_current]
    # compute norm of sum of s_i_current
    current_distance <- distance_function(s_i_current)
    # transform (e.g., exponentiate) "distance"
    current_distance_prob <- transform_function(current_distance)
  } else if (distance_method == 4) {
    distance_function <- function(x) {
      # sum the column vectors of x
      sum_cols <- x %*% rep(1, ncol(x))
      left <- - tau - sum_cols
      right <- sum_cols - (1 - tau)
      max_result <- vector("double", length(sum_cols))
      for (entry in seq_along(left)) {
        left_entry <- left[[entry]]
        right_entry <- right[[entry]]
        max_result[[entry]] <- max(left_entry, right_entry, 0)
      }
      tmp <- sum(abs(max_result)^l_norm) ^ (1 / l_norm)
      exp(-gamma * tmp^l_power) # TODO: is this correct? shouldn't this be transform_function?
      warning("need to double-check distance_method == 4")
    }
    okay <- setdiff(seq_len(n), h)
    ones_current <- which(current_subsample[okay] == 1)
    s_i_current <- s_i[, ones_current]
    # compute norm of sum of s_i_current
    current_distance <- distance_function(s_i_current)
    # transform (e.g., exponentiate) "distance"
    current_distance_prob <- transform_function(current_distance)
  } else if (distance_method == 5) {
    # TODO: is it possible to use compute_xi_i?
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
    s <- do.call(rbind, s_i) # n rows
    stopifnot(ncol(s) == length(h))
    stopifnot(nrow(s) == n)

    distance_function <- function(x) {
      sum_cols <- rep(1, m) %*% x # TODO: double-check this
      left <- - tau - sum_cols
      right <- sum_cols - (1 - tau)
      max_result <- vector("double", length(h))
      for (entry in seq_along(h)) {
        left_entry <- left[[entry]]
        right_entry <- right[[entry]]
        # - tau < x < 1 - tau => - tau - x < 0 & x - (1 - tau) < 0
        # => max{...} = 0 => we satisfy FOC
        max_result[[entry]] <- max(left_entry, right_entry, 0)
      }
      tmp <- sum(abs(max_result)^l_norm) ^ (1 / l_norm)
    }

    curr_s <- s[current_subsample == 1, ]
    current_distance <- distance_function(curr_s)
    current_distance_prob <- transform_function(current_distance)
  }

  set.seed(seed)
  u_vec <- runif(iter)
  out_subsample <- vector("list", iter)
  out_distance <- vector("list", iter)
  out_distance_prob <- vector("list", iter)
  out_a_log <- vector("double", iter)
  out_record <- vector("double", iter)
  for (i in seq_len(iter)) {
    u <- u_vec[[i]]

    # choose k
    # if (k_method == "random") {
    #   # Q: how do we choose k randomly?
    #   # Q: how do we compute proposal_prob? choose(m-p, k) * choose(n - m, k)?
    # }
    if (k_method == "constant") {
      # ensure proposal_prob / current_prob = 1
      proposal_prob <- current_prob
    }

    # random walk: switch k 1's with k 0's
    ones <- setdiff(which(current_subsample == 1), h)
    # print(which(current_subsample == 1))
    # print(ones)
    # print(h)
    zeros <- setdiff(which(current_subsample == 0), h)
    switch_ones_to_zeros <- sample(x = ones, size = k)
    switch_zeros_to_ones <- sample(x = zeros, size = k)
    proposal_subsample <- current_subsample
    proposal_subsample[switch_ones_to_zeros] <- 0
    proposal_subsample[switch_zeros_to_ones] <- 1

    # compute distance of proposal subsample
    if (distance_method == 1) {
      proposal_s <- s[proposal_subsample == 1, ]
      proposal_distance <- distance_function(matrix(1, nrow = 1, ncol = nrow(proposal_s)) %*% proposal_s)
      proposal_distance_prob <- transform_function(proposal_distance)
    } else if (distance_method == 2) {
      proposal_distance <- distance_function(proposal_subsample)
      proposal_distance_prob <- transform_function(proposal_distance)
    } else if (distance_method == "simple_random_walk") {
      proposal_distance <- NA
      proposal_distance_prob <- 1
    } else if (distance_method == 3) {
      okay <- setdiff(seq_len(n), h)
      ones_proposal <- which(proposal_subsample[okay] == 1)
      s_i_proposal <- s_i[, ones_proposal]
      # compute norm of sum of s_i_proposal
      proposal_distance <- distance_function(s_i_proposal)
      # transform (e.g., exponentiate) "distance"
      proposal_distance_prob <- transform_function(proposal_distance)
    } else if (distance_method == 4) {
      okay <- setdiff(seq_len(n), h)
      ones_proposal <- which(proposal_subsample[okay] == 1)
      s_i_proposal <- s_i[, ones_proposal]
      # compute norm of sum of s_i_proposal
      proposal_distance <- distance_function(s_i_proposal)
      # transform (e.g., exponentiate) "distance"
      proposal_distance_prob <- transform_function(proposal_distance)
    } else if (distance_method == 5) {
      proposal_s <- s[proposal_subsample == 1, ]
      proposal_distance <- distance_function(proposal_s)
      proposal_distance_prob <- transform_function(proposal_distance)
    }

    # compute acceptance probability
    # print(proposal_distance_prob) # DEBUG:
    # print(current_distance_prob) # DEBUG:
    # print(proposal_prob) # DEBUG:
    # print(current_prob) # DEBUG:
    # print(i) # DEBUG:

    accept_bool <- tryCatch({
      a_log <- log(proposal_distance_prob) - log(current_distance_prob) + log(proposal_prob) - log(current_prob)
      out_a_log[[i]] <- a_log
      bool <- log(u) < a_log
      # print(bool) # DEBUG:
      stopifnot(is.logical(bool) & !is.na(bool))
      list(
        status = "OKAY",
        bool = bool
      )
    }, error = function(e) {
      list(
        status = "ERROR",
        status_message = e,
        bool = bool,
        a_log = a_log,
        s_i_current = s_i_current,
        s_i_proposal = s_i_proposal,
        proposal_distance = proposal_distance,
        current_distance = current_distance,
        proposal_distance_prob = proposal_distance_prob,
        current_distance_prob = current_distance_prob
      )
    })
    if (accept_bool$status == "ERROR") {
      return(accept_bool)
    }

    accept_bool <- accept_bool$bool
    if (accept_bool) { # accept
      current_subsample <- proposal_subsample
      current_distance <- proposal_distance
      current_distance_prob <- proposal_distance_prob
      out_record[[i]] <- 1
    } else {
      out_record[[i]] <- 0
    }
    if (!save_memory) {
      out_subsample[[i]] <- current_subsample
    }
    out_distance[[i]] <- current_distance
    out_distance_prob[[i]] <- current_distance_prob
  }

  # TODO: compute foc_membership?

  out <- list(
    status = "OKAY",
    a_log = out_a_log,
    subsample = ifelse(save_memory, "`save_memory` is `TRUE`", out_subsample),
    distance = out_distance,
    distance_prob = out_distance_prob, # use in main MCMC
    record = out_record
  )

  # save foc membership if `distance_method == 5`
  if (distance_method == 5) {
    out$foc_membership <- sapply(out_distance, function(i) {
      isTRUE(all.equal(i, 0))
    })
  }

  out
}
