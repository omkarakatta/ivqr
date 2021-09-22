### Title: Pseudo-Marginal Monte Carlo Markov Chain Sampling of the Subsampling Distribution
###
### Description: See ~/BFI/5_IVQR_GP/fromGuillaume/pseudi_marginal_mcmc_subsampling_v2.pdf
### This file contains functions that will be helpful for sampling from the
### subsampling distribution of the IQR estimator. The subsampling distribution
### closely mimics the true distribution of the IQR estimator as seen in
### ~/BFI/ivqr/scrap/sim/bootstrap/subsample.R. Creating the subsampling
### distribution is quite expensive: for each subsample, we need to solve the
### IQR problem. Instead, let's try developing the subsampling distribution
### (i.e., generate samples from the subsampling distribution) without having
### to solve the problem.
###
### Bird's Eye View -------------------------
###
### Here is a bird's-eye-view of how we can sample from the subsampling distribution.
### 1. propose a candidate solution to the IQR problem according to a distribution `Q`
###    (see Algorithm A below)
### 2. determine the probability `P` that the proposal is the solution to an IQR problem of a subsample
###    (see Algorithm B below)
### 3. construct the acceptance probability using `P` and `Q`; be sure to correct for the intensity of the proposals
### 4. accept or reject the proposal, and repeat the process
### 5. the set of accepted proposals = samples from the subsampling distribution
###
### This idea is not yet complete. We need to answer the following questions:
### A. How do we propose the candidate solution? That is, what is `Q`?
### B. And how do we determine `P`?
###
### Let's tackle Question A: -------------------------
###
### First, recall that the solution to the IQR problem is also a solution to a QR problem.
### Since we know there is a semi-closed form solution to the QR problem in
### terms of the active basis, there is a way to express the solution to the
### IQR problem in terms of the active basis. This motivates the sampling of
### the IQR estimator in the following way:
###
### Algorithm A: proposing active basis and obtaining candidate solutions as a by-product
###   1. propose active basis according to some proposal distribution
###   2. compute solution to QR problem using semi-closed form expression
###   3. determine probability `phi`; this is the almost the limit distribution of the IQR estimator 
###   4. use proposal distribution and `phi` to determine acceptance probability;
###      remember to correct for the frequency of proposal
###   5. accept or reject the QR solution to come up with candidate solutions in Step 1. of the birds-eye-view
###
### There are two questions that remain:
###   A1. How do we propose the active basis?
###
### To answer A1, let's construct some weights for each of the n observations in our original data set.
### These weights will depend on the residuals obtained from solving the IQR problem using the full data set.
### The larger the residuals are in magnitude, the smaller the weights.
### Then, define Q to be a multivariate noncentral Fisher hypergeometric distribution.
### That is, the probability that the proposed active basis is comprised of
### h_{1}...h_{n} indices is proportional to the product of the corresponding
### weights for h_{1},...,h_{n}.
### Using these n indices, we can compute the solution to the IQR problem using
### the closed-form solution to the QR problem.
###
### As mentioned in the algorithm, `phi` is the almost limit distribution of the IQR estimator.
### As we've seen in the past, using the limit distribution yields really large
### confidence intervals. It makes for a good proposal distribution, which is
### why we're using it in Algorithm A.
###
### I've been writing that `phi` is "almost" the limit distribution of the IQR estimator.
### We want `phi` to be the limit distribution AFTER we restrict it to the
### support of the coefficients that solve the QR.
### Thus, if we let N denote the density of the Wald distribution and
### `\hat{\beta}` denote the result of step 2 of Algorithm A, we compute `phi`
### as follows:
### `N(\hat{\beta}) / \sum_{i \in \support(\beta)} N(i)`.
### The denominator is the sum of the limit distribution evaluated at each point in the support.
###
### At the end of Algorithm A, we will have what we need to proceed to Step 2 of the bird's-eye-view algorithm.
###
### Now, let's tackle Question B: -------------------------
###
### We actually don't need the probability that the proposed solution from Step
### 1 is the solution to an IQR problem of a subsample. We actually just need
### an unbiased estimate of this probability.
### Below is an algorithm that achieves this. Note that this algorithm replaces
### steps 2-4 of the bird's-eye-view algorithm, not just stepj 2.
###
### Algorithm B:
###   1. store the candidate solution from the end of Step 1 of the
###      bird's-eye-view algorithm (a.k.a. end of Step 5 of Algorithm A above)
###      and all other info e.g. density of the proposal
###   2. sample `M` subsamples from the space of possible subsamples according to a distribution `R`
###   3. for each subsample, determine the probability that the candidate
###      solution is the IQR estimator for that subsample
###   4. construct the acceptance probabilities using the probability of
###      proposing the candidate solution, probability of proposing the subsamples,
###      and the probability that the candidate solution is indeed the IQR estimator
###      for the subsample (we'll write this as an average across each of the
###      subsamples); correct for the intensity of the proposals
###   5. accept or reject the candidate solution and the subsamples
###
### Observe that the acceptance probability in Step 4 uses information from the
### previous steps of Algorithm B and the information from Algorithm A.
###
### Now the following questions remain:
### B1. How do we sample the subsamples and compute the probability of actually sampling those subsamples?
###
### One approach would be to sample the subsamples uniformly from the space of all possible subsamples.
### That is, we randomly choose `m` of the `n` observations uniformly at random
### to create a single subsample, and repeat the process times to create a
### total of `M` subsamples. The probability that each of these subsamples is
### realized is `1 / (n choose m)`.
### There is a major drawback with this strategy.
### None of these `M` subsamples are likely going to yield an IQR problem whose
### solution is the candidate solution from Step 1 of Algorithm B.
### This is because the space of subsamples that satisfy the FOCs of the IQR
### problem at the candidate solution is extremely small relative the space of
### all possible subsamples.
### So, we will constantly be rejecting our candidate solution.
### What we would like to have instead is a way to sample subsamples that are
### more likely to satisfy the FOCs of the IQR problem at the candidate
### solution; BUT we will need to correct for this sort of aggressive proposal behavior.
###
### We have two approaches to accomplish this.
### First, we can build each subsample point-by-point.
### Second, we start with an initial sample (perhaps using the first approach).
### Then, we try to move closer to the space of the subsamples that satisfy the
### FOCs of the IQR problem using the simplex algorithm. Along the path we take
### to get closer to the desired space, we can choose `N` points and uses these
### as subsamples. Some of these may lie inside the desired space.
###
### Author: Omkar A. Katta

### Preliminaries -------------------------

#' Verify membership of a data set in the FOC conditions
#'
#' @param h Indices of the active basis written in terms of the subsample data [p-dimensional vector]
#' @param X_subsample Covariates in subsample [m by p_X matrix]
#' @param D_subsample Endogeneous variables in subsample [m by p_D matrix]
#' @param Phi_subsample Transformed instruments in subsample [m by p_Phi]
#' @param y_subsample Outcome vector in the subsample [m by 1 matrix]
#' @param tau Quantile [numeric]
#' @param beta_D Coefficients on the endogeneous variable; ideally obtained from \code{h} [p_D by 1 matrix]
# TODO: create a function that figures out what beta_D from h (from GP)
foc_membership <- function(h, D_subsample, Phi_subsample, X_subsample, y_subsample, tau, beta_D) {
  # check dimensions
  m <- nrow(y_subsample)
  # TODO: consider removing this when we are done to improve speed
  stopifnot(nrow(D_subsample) == m)
  stopifnot(nrow(Phi_subsample) == m)
  stopifnot(nrow(X_subsample) == m)

  p <- ncol(X_subsample) + ncol(Phi_subsample) # we concentrate out D_subsample
  stopifnot(length(h) == p)

  # compute relevant data
  Dh <- D_subsample[h, , drop = FALSE]
  Phih <- Phi_subsample[h, , drop = FALSE]
  Xh <- X_subsample[h, , drop = FALSE]
  yh <- y_subsample[h]
  design <- cbind(Phi_subsample, X_subsample) # design matrix
  designh <- design[h, , drop = FALSE] # design matrix

  # compute b(h)
  bh <- solve(designh) %*% yh # bh is a p by 1 matrix
  stopifnot(nrow(bh) == p)
  stopifnot(ncol(bh) == 1)

  # compute resid_subsample
  # TODO: figure out beta_D from h if beta_D is not specified
  stopifnot(nrow(beta_D) == p_D)
  stopifnot(ncol(beta_D) == 1)
  y_tilde_subsample <- y_subsample - D_subsample %*% beta_D
  reg <- quantreg::rq(y_tilde_subsample ~ X_subsample + Phi_subsample)
  resid_subsample <- reg$residuals

  # create indices that are not in h
  # TODO: the indices of h need to be written in terms of the subsample.
  # e.g.:
  # suppose the 5th observation in the full data is the smallest index in the active basis.
  # further suppose the 5th observation is the first one present in the subsample.
  # then, this means that `1` should be present in `h` since the first observation in the subsample corresponds to the 5th observation in full data.
  noth <- setdiff(seq_len(m), h)

  # compute xi(h, D)
  resid_noth <- matrix(resid_subsample[noth], ncol = 1) # (m - p) by 1 matrix
  ind_mat <- diag(tau - as.numeric(resid_noth < 0)) # (m - p) by (m - p) matrix
  design_noth <- design[noth, , drop = FALSE] # (m - p) by p matrix
  summand <- ind_mat %*% design_noth
  sum_summand <- matrix(1, nrow = 1, ncol = length(noth))
  xi <- t(sum_summand %*% solve(Xh))
  stopifnot(nrow(xi) == p)
  stopifnot(ncol(xi) == 1)

  # # TODO: vectorize
  # # compute xi(h,D)
  # xi <- 0
  # for (i in seq_along(noth)) {
  #   index <- noth[[i]] # if index is in h, then residual is 0; otherwise, residual is not 0 (in general position)
  # #   ind <- y_subsample[[index]] - design[[index, , drop = FALSE]] %*% bh < 0
  #   ind <- resid_subsample[[index]] < 0
  #   summand <- (tau - ind) %*% X_subsample[i, , drop = FALSE] %*% solve(Xh)
  #   xi <- xi + summand
  # }

  # HACK: this is old code that iterated through all `m` observations; this is deprecated
  # BUG: the `ind` should be the residuals of the QR of interest.
  # # compute xi(h,D)
  # xi <- 0 # initialize sum
  # # iterate through all `m` observations of the subsample
  # # TODO: sum over indices that aren't in `h`
  # for (i in seq_along(resid_subsample)) {
  #   if (!all.equal(resid_subsample[[i]], 0)) {
  #     # compute indicator
  #     ind <- y_subsample[[i]] - X_subsample[i, , drop = FALSE] %*% bh < 0 # scalar
  #     # compute summand
  #     summand <- (tau - ind) %*% X_subsample[i, , drop = FALSE] %*% solve(Xh)
  #     xi <- xi + summand
  #   }
  # }

  # check if xi satisfies FOC inequality
  left <- matrix(-1 * tau %*% rep(1, p), ncol = 1)
  right <- matrix((1 - tau) %*% rep(1, p), ncol = 1)
  all((left <= xi) & (xi <= right)) # returns TRUE if both are true, FALSE otherwise
}
# TODO: Sanity Check -- get subsample, solve for the beta, and then pass the corresponding h and the subsample into foc_membership
# TODO: Sanity Check -- check for the false case
