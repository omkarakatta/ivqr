### one-step neighbors -------------------------

# TODO: document
#' Check membership of one-step neighbors
#'
#' @param subsample An n-long vector with m ones
#' @param reference An n-long vector from which we measure the distance to the
#'  one-step neighbors
#' @param h A vector of indices that are the active basis
#' @param Y,X,D,Phi data
#' @param tau Quantile
#' @param MEMBERSHIP_FCN function for checking membership
#' @param ... arguments for STATUS
#'
#' @return A named list with two elements
#'  1. distance: vector with distances between neighbors and \code{reference}
#'  2. status: vector indicating whether neighbor is inside polytope
onestep <- function(subsample, reference,
                    h, Y, X, D, Phi, tau,
                    MEMBERSHIP_FCN = foc_membership_v3,
                    ...) {
  stopifnot(subsample[h] == 1)
  ones <- setdiff(which(subsample == 1), h)
  zeros <- which(subsample == 0)
  status_vec <- vector("double", length(ones) * length(zeros))
  distance_vec <- vector("double", length(ones) * length(zeros))
  counter <- 0
  for (one_to_zero in ones) {
    for (zero_to_one in zeros) {
      counter <- counter + 1
      neighbor <- subsample
      neighbor[one_to_zero] <- 0
      neighbor[zero_to_one] <- 1
      distance <- sum((neighbor - reference)^2)^(0.5)
      sub_ind <- which(neighbor == 1)
      membership_info <- MEMBERSHIP_FCN(
        h = which(sub_ind %in% h),
        Y_subsample = Y[sub_ind, , drop = FALSE],
        X_subsample = X[sub_ind, , drop = FALSE],
        D_subsample = D[sub_ind, , drop = FALSE],
        Phi_subsample = Phi[sub_ind, , drop = FALSE],
        tau = tau,
        ...
      )
      status_vec[[counter]] <- as.integer(membership_info$status)
      distance_vec[[counter]] <- distance
    }
  }
  list(
       distance = distance_vec,
       status = status_vec
  )
}

### exhaustive_membership -------------------------

# check all possible subsamples to see if it is inside FOC polytope
# returns list of:
# 1. status_vec: logical vector, TRUE when the subsample is inside FOC polytope
# 2. subsample_list: list of n-vectors with m 1's with entries in h
exhaustive_membership <- function(
  h, n, m,
  Y, X, D, Phi, tau,
  MEMBERSHIP_FCN = foc_membership_v3,
  ...
) {
  p <- length(h)

  subsample_template <- rep(0, n)
  subsample_template[h] <- 1
  possible_indices <- which(subsample_template == 0)
  tmp <- expand.grid(rep(list(possible_indices), m - p))
  keep_rows <- apply(tmp, 1, function(row) {
    row <- as.numeric(row)
    # columns must be in strictly monotonic order
    # example: c(1,1) => FALSE
    # example: c(2,1) => FALSE
    # example: c(2,3) => TRUE
    identical(sort(unique(row)), row)
  })
  subsample_indices_mat <- tmp[keep_rows, ]
  num_subsamples <- choose(n - p, m - p)
  # stopifnot(nrow(subsample_indices_mat) == num_subsamples)

  status_vec <- vector("double", num_subsamples)
  subsample_list <- vector("list", num_subsamples)
  for (i in seq_len(num_subsamples)) {
    new_indices <- as.numeric(subsample_indices_mat[i, ])
    new_subsample <- subsample_template
    new_subsample[new_indices] <- 1
    subsample_list[[i]] <- new_subsample
    stopifnot(sum(new_subsample) == m)
    sub_ind <- which(new_subsample == 1)
    status_vec[i] <- MEMBERSHIP_FCN(
      h = which(sub_ind %in% curr_h),
      Y_subsample = Y[sub_ind, , drop = FALSE],
      X_subsample = X[sub_ind, , drop = FALSE],
      D_subsample = D[sub_ind, , drop = FALSE],
      Phi_subsample = Phi[sub_ind, , drop = FALSE],
      tau = tau,
      ...
    )$status
  }

  list(
    status_vec = status_vec,
    subsample_list = subsample_list
  )
}

