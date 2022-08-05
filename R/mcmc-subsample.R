#' Algorithm 3: Random Walk
#'
#' Run an MCMC that samples subsamples from the subsample polytope.
#'
#' We generate a sample of subsamples by running an MCMC where the proposal is a
#' random walk on the subsample polytope and the target achieves uniformity in
#' the distance away from the anchor.
#'
#' See ~/projects/5_IVQR_GP/notes/mcmc-psuedocode/alg3.pdf for more details.
#'
#' @param h Indices of active basis
#' @param Y,X,D,Phi Data
#' @param tau Quantile
#' @param anchor Reference subsample that is considered to be the "origin"
#' @param initial Initial subsample from which we start the random walk
#' @param iterations Number of iterations
#' @param gamma Tuning parameter on target probability
#' @param zero_prob_tol smallest possible value for the target probability
#'
#' @return
rwalk_subsample <- function(
  h,
  Y, X, D, Phi,
  tau,
  anchor,
  initial,
  iterations,
  gamma,
  zero_prob_tol # TODO: need to add this to the note
) {
  # Setup-level Information, e.g., size of problem, foc conditions, etc. -------
  anchor <- as.integer(anchor)
  initial <- as.integer(initial)
  n <- length(anchor) # number of observations in original data set
  m <- sum(anchor)    # number of observations in subsample
  p <- length(h)      # number of indices in active basis
  xi_mat <- compute_foc_conditions(
    h,
    Y = Y,
    X = X,
    D = D,
    Phi = Phi,
    tau = tau
  )

  # Sanity Checks --------------------------------------------------------------

  # Subsamples can only contain 0's and 1's.
  stopifnot(sort(unique(anchor)) == c(0L, 1L))
  stopifnot(sort(unique(initial)) == c(0L, 1L))
  # Subsamples must be n-long vectors with m 1's.
  stopifnot(length(initial) == n)
  stopifnot(sum(initial) == m)
  # Subsamples must have 1's in the active basis indices.
  stopifnot(unique(anchor[h]) == 1)
  stopifnot(unique(initial[h]) == 1)

  # Information about Anchor ---------------------------------------------------

  anchor_ones <- setdiff(which(anchor == 1L), h)
  anchor_zeros <- which(anchor == 0L)

  # Pre-allocate Results -------------------------------------------------------
  out_record <- vector("double", iterations)
  # out_violation_norm <- # TODO: dimensions?
  out_membership <- vector("logical", iterations)
  # acceptance record
  # distance-to-anchor of proposed subsamples D⋆
  # distance-to-anchor of accepted subsamples Dᵢ
  # foc membership of proposed subsamples D⋆
  # foc membership of accepted subsamples Dᵢ
  # |S(D⋆)| -- #subsamples whose distance-to-anchor is same as that of D⋆
  # |S(Dᵢ)| -- #subsamples whose distance-to-anchor is same as that of Dᵢ
  # |C₁(Dᵢ)| -- #one-step neighbors of Dᵢ that are closer to anchor than Dᵢ
  # |S₁(Dᵢ)| -- #one-step neighbors of Dᵢ that are same to anchor than Dᵢ
  # |F₁(Dᵢ)| -- #one-step neighbors of Dᵢ that are farther to anchor than Dᵢ
  # |C₁(D⋆)| -- #one-step neighbors of D⋆ that are closer to anchor than D⋆
  # |S₁(D⋆)| -- #one-step neighbors of D⋆ that are same to anchor than D⋆
  # |F₁(D⋆)| -- #one-step neighbors of D⋆ that are farther to anchor than D⋆
  # Proposed direction, denoted d (either "closer", "same", or "farther")
  # |d₁(Dᵢ)| -- #one-step neighbors of Dᵢ that were in the direction d
  # |¬d₁(D⋆)| -- #one-step neighbors of D⋆ that were in the direction ¬d

  # MCMC Initialization --------------------------------------------------------

  # We don't need to initialize `violation_info$violation_norm` because we know
  # we will accept in the first iteration.
  # Hence, we will compute the `violation_info$violation_norm` once we accept.
  current <- initial
  log_P_current <- -Inf # accept the proposal in the first iteration
  violation_info <- list()
  violation_info$xi_vec <- NULL

  # uniform draws
  log_u_vec <- log(runif(n = iterations))

  # Run MCMC -------------------------------------------------------------------
  for (mcmc_idx in seq_len(iterations)) {
    # Draw from uniform distribution.
    log_u <- log_u_vec[[mcmc_idx]]

    # Collect information about current subsample.
    current_ones <- setdiff(which(current == 1L), h)
    current_zeros <- which(current == 0L)

    # Characterize indices based on their values in `current` and `anchor`.
    # Here is a table to describe them:
    # | name              | value in current  | value in anchor |
    # |:------------------|:-----------------:|:---------------:|
    # | common_ones       |        1          |        1        |
    # | common_zeros      |        0          |        0        |
    # | different_ones    |        1          |        0        |
    # | different_zeros   |        0          |        1        |
    common_ones <- intersect(current_ones, anchor_ones)
    common_zeros <- intersect(current_zeros, anchor_zeros)
    different_ones <- intersect(current_ones, anchor_zeros)
    different_zeros <- intersect(current_zeros, anchor_ones)
    # Aside from the active basis, there should be m - p indices that are 1's.
    stopifnot(length(common_ones) + length(different_ones) == m - p)
    # There should be n - m indices that are 0's.
    stopifnot(length(common_zeros) + length(different_zeros) == n - m)

    # Define the number of indices that are 1 in both `current` and `anchor`.
    r_current <- length(common_ones) + p

    # Compute |C₁(Dᵢ)|, |S₁(Dᵢ)|, and |F₁(Dᵢ)|, i.e., number of one-step
    # neighbors of current subsample that are closer to, same distance away
    # from, and farther away from `anchor` as the current subsample.
    onestep_current <- onestep_cardinalities(n, m, p, r_current)

    # Choose a direction uniformly at random.
    d <- sample(
      names(onestep_current),
      1,
      prob = as.numeric(onestep_current > 0)
    )

    # Choose a one-step neighbor of `current` uniformly at random in the
    # direction `d`.
    # Compute number `r` of shared observations between proposal and anchor.
    if (d == "closer") {
      one_to_zero <- alt_sample(different_ones, 1)
      zero_to_one <- alt_sample(different_ones, 1)
      r_star <- r_current + 1
    } else if (d == "same") {
      one_to_zero <- alt_sample(current_ones, 1)
      zero_to_one <- ifelse(
        one_to_zero %in% common_ones,
        alt_sample(different_zeros, 1),
        alt_sample(common_zeros, 1)
      )
      r_star <- r_current
    } else if (d == "farther") {
      one_to_zero <- alt_sample(common_ones, 1)
      zero_to_one <- alt_sample(common_zeros, 1)
      r_star <- r_current - 1
    }
    # Construct proposal subsample.
    star <- current
    star[one_to_zero] <- 0L
    star[zero_to_one] <- 1L
    # Compute proposal probability of D⋆ given Dᵢ, i.e., log(Q(D⋆|Dᵢ))
    log_Q_star <- log_Q(onestep_current, d)
    # Compute proposal probability of D⋆ given Dᵢ, i.e., Q(Dᵢ|D⋆)
    onestep_star <- onestep_cardinalities(n, m, p, r_star)
    log_Q_current <- log_Q(onestep_star, opposite_direction(d))

    # Construct target probability of proposal.
    # P(D⋆|Dᵦ) = 1/S(D⋆) * exp{-γ * ‖D⋆ - Dᵦ‖}
    log_P_star <- max(
      log_P(gamma, star, anchor, n, m, p, r_star),
      zero_prob_tol
    )

    # Construct acceptance probability ratio.
    log_acc_prob <- log_P_star - log_P_current + log_Q_current - log_Q_star

    # If we accept, update `current` information:
    # `record`
    # `current`
    # `log_P_current`
    # `violation_info`
    # `membership`
    record <- 0
    if (log_u < log_acc_prob) {
      record <- 1
      current <- star
      log_P_current <- log_P_star
      violation_info <- foc_violation(
        h,
        current,
        tau,
        xi_mat,
        xi_vec = violation_info$xi_vec, # incumbent value of xi_vec
        minus_index = one_to_zero,
        plus_index = zero_to_one
      )
      membership <- isTRUE(all.equal(violation_info$violation_norm, 0))
    }
    out_record[[mcmc_idx]] <- record
    out_membership <- membership
  }
  list(
    record = out_record,
    membership = out_membership
  )
}

# Helpers ----------------------------------------------------------------------

# @param n Number of observations in original data set
# @param m Number of observations in subsample
# @param p Number of indices in active basis
# @param r Number of indices that are 1 in both subsample and anchor
onestep_cardinalities <- function(n, m, p, r) {
  # All subsamples must be 1 in the active basis.
  # Hence, number of shared 1's must be at least `p`.
  stopifnot(r >= p)
  c(
    "closer" = (m - r)^2,
    "same" = (m - r) * (n - 2 * (m + r) - p),
    "farther" = (r - p) * (n - 2 * m + r)
  )
}

# @param n Number of observations in original data set
# @param m Number of observations in subsample
# @param p Number of indices in active basis
# @param r Number of indices that are 1 in both subsample and anchor
global_cardinality <- function(n, m, p, r) {
  # All subsamples must be 1 in the active basis.
  # Hence, number of shared 1's must be at least `p`.
  stopifnot(r >= p)
  choose(m - p, r - p) * choose(n - m, m - r)
}

# Obtain "opposite" direction.
opposite_direction <- function(direction) {
  if (direction == "closer") {
    return("farther")
  } else if (direction == "same") {
    return("same")
  } else if (direction == "farther") {
    return("farther")
  } else {
    stop(paste("Unrecognized direction, opposite is not defined:", direction))
  }
}

# Obtain proposal density in logs.
log_Q <- function(onestep, direction) {
  -log(onestep[[direction]]) - log(sum(as.numeric(onestep > 0)))
}

# Obtain target density in logs.
log_P <- function(gamma, star, anchor, n, m, p, r) {
  # The (Euclidean) distance between two subsamples is simply the number of
  # indices whose values are different across the two subsamples.
  distance <- sum(star != anchor)
  -log(global_cardinality(n, m, p, r)) - gamma * distance
}

# Measure norm of the FOC violation in each of the `p` entries of `xi_vec`.
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
  l_norm = 2
) {
  if (is.null(xi_vec)) {
    # p by m matrix multiplied by m-vector of 1's
    xi_vec <- xi_mat %*% as.matrix(subsample, ncol = 1)
  } else {
    xi_vec <- xi_vec - xi_mat[, minus_index] + xi_mat[, plus_index]
  }
  left <- -1 * tau - xi_vec
  right <- xi_vec - (1 - tau)
  violation <- pmax(left, right, rep(0, length(h))) # p by 1
  list(
    violation_norm = sum(abs(violation)^l_norm) ^ (1 / l_norm), # Q: dimensions?
    xi_vec = xi_vec
  )
}
