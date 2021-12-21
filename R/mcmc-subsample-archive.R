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
