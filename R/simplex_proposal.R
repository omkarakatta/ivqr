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
# maybe there is some recursion that I can take advantage of...
get_entering_directions <- function(current_subsample,
                                    last_entrant,
                                    n_directions,
                                    center_subsample) {
  # partition observations
  in_current <- which(current_subsample == 1)
  in_center <- which(center_subsample == 1)
  out_current <- which(current_subsample == 0)
  out_center <- which(center_subsample == 0)
  in_current_out_center <- sort(setdiff(in_current, in_center))
  in_center_out_current <- sort(setdiff(in_center, in_current))
  in_current_in_center <- sort(intersect(in_center, in_current))
  out_current_out_center <- sort(intersect(out_current, out_center))

  # preallocate vector of directions
  # dir <- vector("double", n_directions)

  vec <- in_center_out_current
  n_dir <- n_directions
  dir <- c()

  if (length(vec) >= n_dir) {
    valid_dir <- vec[vec > last_entrant]
    remaining <- n_dir - length(valid_dir)
    if (remaining <= 0) {
      dir <- c(dir, valid_dir[seq_len(n_dir)])
    } else {
      dir <- c(dir, valid_dir, vec[seq_len(remaining)])
    }
  } else {
    dir <- c(dir, vec)
    vec <- out_current_out_center # update choices
    n_dir <- n_dir - length(dir) # update how many we need to choose
    if (length(vec) >= n_dir) {
      valid_dir <- vec[vec > last_entrant]
      remaining <- n_dir - length(valid_dir)
      if (remaining <= 0) {
        dir <- c(dir, valid_dir[seq_len(n_dir)])
      } else {
        dir <- c(dir, valid_dir, vec[seq_len(remaining)])
      }
    }
  }

  dir
}
