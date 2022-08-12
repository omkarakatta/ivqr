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

  # Pre-allocate Results -------------------------------------------------------

  # At each iteration, we will need to keep track of the information for three
  # types of subsamples:
  #   1. incumbent subsample
  #   2. proposed subsample
  #   3. accepted subsample
  # If we reject the proposal, then the accepted subsample is the same as the
  # incumbent subsample.
  # If we accept, then the accepted subsample is the same as the proposed
  # subsample, and the proposed subsample becomes the new incumbent.
  #
  # Note that info for incumbent subsample will be calculated in the previous
  # iteration but saved in the current iteration. Meanwhile, info for the
  # proposed subsample will be calculated and saved in the middle of the
  # iteration. Finally, the info for the accepted subsample will be saved at the
  # end of the iteration based on whether we accept or reject.
  #
  # For each of these types of subsamples, we need to store:
  #   1. distance-to-anchor
  #   2. FOC membership
  #   3. S, i.e., #subsamples with the same distance-to-anchor
  #   4. one-step info
  #     - C₁, i.e., #one-step neighbors that are closer to anchor
  #     - S₁, i.e., #one-step neighbors that are same distance away from anchor
  #     - F₁, i.e., #one-step neighbors that are farther away from anchor
  #
  # In addition to storing information about these subsamples, we'll also store
  # the following at each iteration:
  #   1. proposed direction
  #   2. whether we accept or reject the proposal
  #   3. acceptance-probability ratio components, i.e.,
  #      log_Q_star
  #      log_Q_incumbent
  #      log_P_star
  #      log_P_incumbent
  #      log_acc_prob

  out_incumbent <- vector("list", iterations)
  out_proposal <- vector("list", iterations)
  out_accepted <- vector("list", iterations)
  out_info <- vector("list", iterations)

  # MCMC Initialization --------------------------------------------------------

  # We don't need to initialize `violation_info$violation_norm` because we know
  # we will accept in the first iteration.
  # Hence, we will compute the `violation_info$violation_norm` once we accept.
  incumbent <- initial
  log_P_incumbent <- -Inf # accept the proposal in the first iteration
  violation_info <- list()
  violation_info$xi_vec <- NULL

  # uniform draws
  log_u_vec <- log(runif(n = iterations))

  # Run MCMC -------------------------------------------------------------------
  for (mcmc_idx in seq_len(iterations)) {
    # Draw from uniform distribution.
    log_u <- log_u_vec[[mcmc_idx]]

    # Collect information about incumbent subsample.
    # After the first iteration, the incumbent information will be updated based
    # on whether we accept (i.e., proposal becomes incumbent) or reject (i.e.,
    # incumbent remains the same).
    if (mcmc_idx == 1) {
      incumbent_indices <- characterize_indices(
        subsample = incumbent,
        anchor = anchor,
        h = h
      )
      r_incumbent <- incumbent_indices$sharing
      incumbent_info <- subsample_information(
        subsample = incumbent,
        anchor = anchor,
        h = h,
        r = r_incumbent,
        tau = tau,
        xi_mat = xi_mat
      )
    }
    out_incumbent[[mcmc_idx]] <- incumbent_info

    # Choose a direction uniformly at random.
    d <- sample(
      names(incumbent_info$onestep),
      1,
      prob = as.numeric(incumbent_info$onestep > 0)
    )
    out_info[[mcmc_idx]]$d <- d

    # Choose a one-step neighbor of `incumbent` uniformly at random in the
    # direction `d`.
    # Compute number `r` of shared observations between proposal and anchor.
    if (d == "closer") {
      one_to_zero <- alt_sample(incumbent_indices$different_ones, 1)
      zero_to_one <- alt_sample(incumbent_indices$different_zeros, 1)
      r_star <- r_incumbent + 1
    } else if (d == "same") {
      incumbent_ones <- setdiff(which(incumbent == 1L), h)
      one_to_zero <- alt_sample(incumbent_ones, 1)
      zero_to_one <- ifelse(
        one_to_zero %in% incumbent_indices$common_ones,
        alt_sample(incumbent_indices$different_zeros, 1),
        alt_sample(incumbent_indices$common_zeros, 1)
      )
      r_star <- r_incumbent
    } else if (d == "farther") {
      one_to_zero <- alt_sample(incumbent_indices$common_ones, 1)
      zero_to_one <- alt_sample(incumbent_indices$common_zeros, 1)
      r_star <- r_incumbent - 1
    }
    # Construct proposal subsample.
    star <- incumbent
    star[one_to_zero] <- 0L
    star[zero_to_one] <- 1L
    stopifnot(sum(star) == m)

    proposal_info <- subsample_information(
      subsample = star,
      anchor = anchor,
      h = h,
      r = r_star,
      tau = tau,
      xi_mat = xi_mat
    )
    out_proposal[[mcmc_idx]] <- proposal_info

    # Compute proposal probability of D⋆ given Dᵢ, i.e., log(Q(D⋆|Dᵢ))
    log_Q_star <- log_Q(incumbent_info$onestep, d)
    # Compute proposal probability of D⋆ given Dᵢ, i.e., Q(Dᵢ|D⋆)
    log_Q_incumbent <- log_Q(proposal_info$onestep, opposite_direction(d))

    # Construct target probability of proposal.
    # P(D⋆|Dᵦ) = 1/S(D⋆) * exp{-γ * ‖D⋆ - Dᵦ‖}
    log_P_star <- max(
      log_P(gamma, star, anchor, proposal_info$global_cardinality),
      zero_prob_tol
    )

    # Construct acceptance probability ratio.
    log_acc_prob <- log_P_star - log_P_incumbent + log_Q_incumbent - log_Q_star

    out_info[[mcmc_idx]]$log_acc_prob <- log_acc_prob
    out_info[[mcmc_idx]]$log_Q_star <- log_Q_star
    out_info[[mcmc_idx]]$log_Q_incumbent <- log_Q_incumbent
    out_info[[mcmc_idx]]$log_P_star <- log_P_star
    out_info[[mcmc_idx]]$log_P_incumbent <- log_P_incumbent

    # If we don't accept,
    #   - let the `record` be 0
    #   - let `accepted` info be given by `incumbent_info`
    # Else, if we accept,
    #   - let the `record` be 1
    #   - let `accepted` info be given by `proposal_info`
    #   - replace `incumbent` info with `proposal` info
    out_info[[mcmc_idx]]$record <- 0
    out_accepted[[mcmc_idx]] <- incumbent_info
    if (log_u < log_acc_prob) {
      out_info[[mcmc_idx]]$record <- 1
      out_accepted[[mcmc_idx]] <- proposal_info

      # for next iteration, update the `incumbent` with `proposal` details
      incumbent <- star
      log_P_incumbent <- log_P_star
      incumbent_info <- proposal_info
      incumbent_indices <- characterize_indices(
        subsample = incumbent,
        anchor = anchor,
        h = h
      )
      r_incumbent <- r_star
    }
  }

  # Return results -------------------------------------------------------------

  list(
    incumbent = out_incumbent,
    proposal = out_proposal,
    accepted = out_accepted,
    info = out_info
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
log_P <- function(gamma, star, anchor, global_cardinality) {
  # The (Euclidean) distance between two subsamples is simply the number of
  # indices whose values are different across the two subsamples.
  distance <- sum(star != anchor)
  -log(global_cardinality) - gamma * distance
}

# Measure norm of the FOC violation in each of the `p` entries of `xi_vec`.
# - tau < x < 1 - tau => - tau - x < 0 & x - (1 - tau) < 0
# => max{...} = 0 => we satisfy FOC
#
# `xi_mat` is p by n matrix containing the information needed to evaluate
# whether a subsample satisfies the FOC conditions. The sum of the `m` columns
# gives us a p by 1 vector called `xi_vec`. Each component of `xi_vec` must be
# between -tau and 1-tau.
#
# If two subsamples are one-step neighbors of each other, then there is only a
# single observation in the subsample A that is not in subsample B (called
# `minus_index`, and there is
# a single observation in subsample B (called `plus_index`) that is not in A.
# So, given `xi_vec` for subsample A, we can compute `xi_vec` for B by
# subtracting the column of `xi_mat` corresponding to `minus_index` and adding
# the column of `xi_mat` corresponding to `plus_index`.
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
    # p by n matrix multiplied by m-vector of 1's
    xi_vec <- xi_mat %*% as.matrix(subsample, ncol = 1)
  } else {
    xi_vec <- xi_vec - xi_mat[, minus_index] + xi_mat[, plus_index]
  }
  left <- -1 * tau - xi_vec
  right <- xi_vec - (1 - tau)
  violation <- pmax(left, right, rep(0, length(h))) # p by 1
  list(
    violation_norm = sum(abs(violation)^l_norm) ^ (1 / l_norm), # length is 1
    xi_vec = xi_vec
  )
}

# Characterize indices by comparing `subsample` and `anchor`
#' @param subsample An n-long vector with m 1's (including in active basis); we
#'  want to collect information about this subsample
#' @param anchor Anchor subsample
#' @param h Active basis indices
characterize_indices <- function(
  subsample,
  anchor,
  h
) {
  # PERF: maybe anchor_* can be moved outside so we don't have to recompute
  # them? Some of these variables are "constant".
  n <- length(anchor)
  m <- sum(anchor)
  p <- length(h)
  anchor_ones <- setdiff(which(anchor == 1L), h)
  anchor_zeros <- which(anchor == 0L)
  ones <- setdiff(which(subsample == 1L), h)
  zeros <- which(subsample == 0L)

  # Characterize indices based on their values in `incumbent` and `anchor`.
  # Note that active basis indices are not included when describing *_ones.
  # Here is a table to describe them:
  # | name              | value in subsample  | value in anchor |
  # |:------------------|:-------------------:|:---------------:|
  # | common_ones       |        1            |        1        |
  # | common_zeros      |        0            |        0        |
  # | different_ones    |        1            |        0        |
  # | different_zeros   |        0            |        1        |
  common_ones <- intersect(ones, anchor_ones)
  common_zeros <- intersect(zeros, anchor_zeros)
  different_ones <- intersect(ones, anchor_zeros)
  different_zeros <- intersect(zeros, anchor_ones)
  # Aside from the active basis, there should be m - p indices that are 1's.
  stopifnot(length(common_ones) + length(different_ones) == m - p)
  # There should be n - m indices that are 0's.
  stopifnot(length(common_zeros) + length(different_zeros) == n - m)

  # Return results
  list(
    common_ones = common_ones,
    common_zeros = common_zeros,
    different_ones = different_ones,
    different_zeros = different_zeros,
    sharing = length(common_ones) + p
  )
}

# Collect information about a subsample
#' @param subsample An n-long vector with m 1's (including in active basis); we
#'  want to collect information about this subsample
#' @param anchor Anchor subsample
#' @param h Active basis indices
#' @param xi_mat See `compute_foc_conditions`
#' @param tau Quantile
#'
#' @return Named list:
#'  1. distance: distance of subsample to anchor
#'  2. onestep: See `onestep_cardinalities()`
#'  3. global_cardinality: See `global_cardinality()`
#'  4. foc_violation_norm
#'  5. membership
subsample_information <- function(
  subsample,
  anchor,
  h,
  r,
  xi_mat,
  tau
) {
  n <- length(subsample)
  m <- sum(subsample)
  p <- length(h)

  distance <- sum(subsample != anchor)

  # Compute |C₁(Dᵢ)|, |S₁(Dᵢ)|, and |F₁(Dᵢ)|, i.e., number of one-step
  # neighbors of incumbent subsample that are closer to, same distance away
  # from, and farther away from `anchor` as the incumbent subsample.
  onestep <- onestep_cardinalities(n, m, p, r)

  # Compute |S(Dᵢ)|, i.e., number of subsamples that are the same distance away
  # from anchor as Dᵢ
  global <- global_cardinality(n, m, p, r)

  # FOC violation
  foc_info <- foc_violation(h, subsample, tau, xi_mat)
  foc_violation_norm <- foc_info$violation_norm
  membership <- isTRUE(all.equal(foc_violation_norm, 0))

  # Return results as a named list
  list(
    "distance_to_anchor" = distance,
    "onestep" = onestep,
    "global_cardinality" = global,
    "foc_violation_norm" = foc_violation_norm,
    "membership" = membership
  )
}
