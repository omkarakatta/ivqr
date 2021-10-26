### Title: Propose subsamples by walking along the simplex
###
### Description: One part of the MCMC sampler of the subsampling distribution of
### the IQR estimator is the proposal of subsamples. We could propose the
### subsamples uniformly at random. However, the FOC polytope containing
### subsamples is very small, so we likely won't propose any subsamples whose
### solution is the proposed beta_D. We can try to propose subsamples with
### a distribution that has greater density on the subsamples in the FOC
### polytope. However, since the FOC polytope is specific to the beta_D, the
### proposal distribution in the numerator (which is using the proposed beta_D)
### and the proposal distribution in the denominator (which is using the current
### beta_D) will not have the same normalization constant. Therefore, we need to
### compute the normalization constant, which may not be easy. Alternatively,
### we can come up with a method whose normalization constant is invariant to
### the beta_D.
### The method we code up here does something like this. We first identify a
### subsample that is inside the FOC polytope. See
### `find_subsample_in_polytope()`. This subsample will be called
### `center_subsample` (i.e., the subsample that is in the "center" of the FOC
### polytope). Next, we start walking along the simplex of subsamples. The way
### we walk is important. At our current subsample we ask for the `n_directions`
### that bring us closer to DC. Among these `n_directions`, we choose 1 at
### random. While this description sweeps a lot of details under the rug, the
### key to this idea is that the way we propose the subsamples are the same
### regardless of the beta_D. The directions we choose will be different, but
### since we choose the direction randomly from the same number of directions,
### the normalization constant will be the same in the numerator and denominator
### of the acceptance rejection probability.
###
### Author: Omkar A. Katta

### choose_entering_directions -------------------------
#' @param current_subsample A binary vector of length `n` with `m` 1's
#' @param last_entrant The observation that was last to enter the
#'  `current_subsample`.
#' @param n_directions Number of directions from which we want to choose
#' @param center_subsample A binary vector of length `n` with `m` 1's that is
#'  inside the FOC polytope of interest
# TODO: finish documenting
get_directions <- function(current_subsample,
                           last,
                           n_directions,
                           center_subsample,
                           mode = "enter") {
  # partition observations
  in_current <- which(current_subsample == 1)
  in_center <- which(center_subsample == 1)
  out_current <- which(current_subsample == 0)
  out_center <- which(center_subsample == 0)

  if (tolower(mode) == "enter") {
    in_center_out_current <- sort(setdiff(in_center, in_current))
    out_current_out_center <- sort(intersect(out_current, out_center))
    vec_list <- list(in_center_out_current, out_current_out_center)
  } else if (tolower(mode) == "exit") {
    in_current_out_center <- sort(setdiff(in_current, in_center))
    in_current_in_center <- sort(intersect(in_center, in_current))
    vec_list <- list(in_current_out_center, in_current_in_center)
  } else {
    stop("Wrong value for `mode`")
  }
  get_valid_directions(vec_list, n_directions, last)
}

# get directions from vec_list until we run out
# recursive function
# TODO: document
get_valid_directions <- function(vec_list, n_dir, last, dir = c()) {
  vec <- vec_list[[1]]
  if (length(vec) >= n_dir) {
    valid_dir <- vec[vec > last]
    remaining <- n_dir - length(valid_dir)
    if (remaining <= 0) {
      dir <- c(dir, valid_dir[seq_len(n_dir)])
    } else {
      dir <- c(dir, valid_dir, vec[seq_len(remaining)])
    }
  } else {
    dir <- c(dir, vec)
    n_dir <- n_dir - length(dir)
    vec_list <- vec_list[2:length(vec_list)]
    dir <- get_valid_directions(vec_list, n_dir, last, dir)
  }
  dir
}
