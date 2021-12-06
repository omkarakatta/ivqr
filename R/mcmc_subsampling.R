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

### h_to_beta -------------------------

#' Find the coefficients given the active basis
#'
#' Given the active basis, find the coefficients that solve the IQR problem.
#'
#' @param h Active basis in terms of the data provided
#' @param Y Dependent variable (vector of length n)
#' @param X Exogenous variable (including constant vector) (n by p_X matrix)
#' @param D Endogenous variable (n by p_D matrix)
#' @param Z Instrumental variable (n by p_Z matrix)
#' @param Phi Transformation of X and Z to be used in the program;
#'  defaults to the linear projection of D on X and Z (matrix with n rows)
#'
#' @return A named list of coefficients
#' \enumerate{
#'   \item \code{beta_D}: coefficients on the endogeneous variables
#'   \item \code{beta_X}: coefficients on the exogenous variables
#'   \item \code{beta_Phi}: coefficients on the transformed instruments
#' }
#'
#' @family mcmc_subsampling
h_to_beta <- function(h, Y, X, D, Z, Phi = linear_projection(D, X, Z)) {
  # dimensions
  p_Phi <- ncol(Phi)
  p_D <- ncol(D)
  stopifnot(p_Phi == p_D)
  p_X <- ncol(X)
  p <- p_Phi + p_X

  design <- cbind(X, Phi)
  designh <- design[h, ]
  a <- solve(designh, Y[h])
  B <- solve(designh, D[h, , drop = FALSE])

  # Solve system from the lower rows, corresponding to beta_Phi = 0
  a_Phi <- a[(p_X + 1):p]
  B_Phi <- B[(p_X + 1):p, ]

  beta_D <- solve(B_Phi, a_Phi)
  beta_X_Phi <- a - B %*% beta_D
  beta_X <- beta_X_Phi[1:p_X]
  beta_Phi <- beta_X_Phi[(p_X + 1):p]

  list("beta_D" = beta_D,
       "beta_X" = beta_X,
       "beta_Phi" = beta_Phi)
}

### foc_membership -------------------------
#' Verify membership of a data set in the FOC conditions
#'
#' From the active basis, we derive the FOC conditions. If a data set satisfies
#' the FOC conditions, then we know that the coefficients obtained from this
#' active basis solve the IQR problem for this data set.
#'
#' @param h Indices of the active basis written in terms of the subsample data
#'  (p-dimensional vector)
#' @param Y_subsample Outcome vector in the subsample (m by 1 matrix)
#' @param X_subsample Covariates in subsample (m by p_X matrix)
#' @param D_subsample Endogeneous variables in subsample (m by p_D matrix)
#' @param Phi_subsample Transformed instruments in subsample (m by p_Phi)
#' @param tau Quantile (numeric)
#' @param beta_D Coefficients on the endogeneous variable; ideally obtained
#'  from \code{h} (p_D by 1 matrix)
#'
#' @return Named list
#'  \enumerate{
#'    \item \code{status}: TRUE if subsample satisfies FOC conditions; FALSE
#'      otherwise
#'    \item \code{xi}: vector that must be contained within -tau and 1-tau to
#'      satisfy FOC conditions
#'  }
#'
#' @family mcmc_subsampling
foc_membership <- function(
  h,
  Y_subsample,
  X_subsample,
  D_subsample,
  Phi_subsample,
  tau,
  beta_D = h_to_beta(h,
                     Y = Y_subsample,
                     X = X_subsample,
                     D = D_subsample,
                     Phi = Phi_subsample)$beta_D
) {
  # TODO: the indices of h need to be written in terms of the subsample.
  # e.g.:
  # suppose the 5th observation in the full data is the smallest index in the
  # active basis. further suppose the 5th observation is the first one present
  # in the subsample. then, this means that `1` should be present in `h` since
  # the first observation in the subsample corresponds to the 5th observation
  # in full data.

  # check dimensions
  m <- nrow(Y_subsample)
  # TODO: remove this when we are done to improve speed
  stopifnot(nrow(D_subsample) == m)
  stopifnot(nrow(Phi_subsample) == m)
  stopifnot(nrow(X_subsample) == m)

  p <- ncol(X_subsample) + ncol(Phi_subsample) # we concentrate out D_subsample
  stopifnot(length(h) == p)

  p_D <- ncol(D_subsample)

  # compute relevant data
  # Dh <- D_subsample[h, , drop = FALSE]
  # Phih <- Phi_subsample[h, , drop = FALSE]
  # Xh <- X_subsample[h, , drop = FALSE]
  yh <- Y_subsample[h]
  design <- cbind(X_subsample, Phi_subsample) # design matrix
  designh <- design[h, , drop = FALSE] # design matrix

  # compute b(h)
  bh <- solve(designh) %*% yh # bh is a p by 1 matrix
  stopifnot(nrow(bh) == p)
  stopifnot(ncol(bh) == 1)

  # compute resid_subsample
  stopifnot(nrow(beta_D) == p_D)
  stopifnot(ncol(beta_D) == 1)
  y_tilde_subsample <- Y_subsample - D_subsample %*% beta_D
  # assume beta_Phi is 0 => doesn't show up in residuals
  reg <- quantreg::rq(y_tilde_subsample ~ X_subsample - 1,
                      tau = tau)
  resid_subsample <- reg$residuals

  # create indices that are not in h
  noth <- setdiff(seq_len(m), h)

  # compute xi(h, D)
  resid_noth <- matrix(resid_subsample[noth], ncol = 1) # (m - p) by 1 matrix
  ind_mat <- diag(tau - as.numeric(resid_noth < 0)) # (m - p) by (m - p) matrix
  design_noth <- design[noth, , drop = FALSE] # (m - p) by p matrix
  summand <- ind_mat %*% design_noth
  sum_summand <- matrix(1, nrow = 1, ncol = length(noth)) %*% summand
  xi <- t(sum_summand %*% solve(designh))
  stopifnot(nrow(xi) == p)
  stopifnot(ncol(xi) == 1)

  # DEBUG: see each entry of sum in xi object
  s_i <- vector("list", length(noth))
  for (i in seq_len(nrow(summand))) {
    s_i[[i]] <- summand[i, ] %*% solve(designh)
  }
  s <- do.call("rbind", s_i)

  # check if xi satisfies FOC inequality
  left <- matrix(-1 * tau %*% rep(1, p), ncol = 1)
  right <- matrix((1 - tau) %*% rep(1, p), ncol = 1)
  list(
    status = all((left <= xi) & (xi <= right)), # returns TRUE if both are true, else FALSE
    s = s, # return entries of xi object
    xi = xi,
    left = left,
    right = right
  )
}
# TODO: Sanity Check -- check for the false case

#' Verify membership of a data set in the FOC conditions
#'
#' If a data set satisfies the FOC conditions with respect to some active
#' basis, then a quantile regression of the concentrated-out outcome using the
#' endogeneous coefficients given by the active basis on the covariates and the
#' transformed instruments should be 0.
#'
#' @param h Indices of the active basis written in terms of the subsample data
#'  (p-dimensional vector)
#' @param Y_subsample Outcome vector in the subsample (m by 1 matrix)
#' @param X_subsample Covariates in subsample (m by p_X matrix)
#' @param D_subsample Endogeneous variables in subsample (m by p_D matrix)
#' @param Phi_subsample Transformed instruments in subsample (m by p_Phi)
#' @param tau Quantile (numeric)
#' @param beta_D Coefficients on the endogeneous variable; ideally obtained
#'  from \code{h} (p_D by 1 matrix)
#'
#' @return Named list
#'  \enumerate{
#'    \item \code{status}: TRUE if subsample satisfies FOC conditions; FALSE
#'      otherwise
#'    \item \code{beta_Phi}: coefficients on transfomed instruments when
#'      running a QR of Y - D %*% beta_D on X and Phi, where beta_D is given by
#'      \code{h}
#'    \item \code{norm}: L1 norm of \code{beta_Phi}
#'  }
#'
#' @family mcmc_subsampling
foc_membership_v2 <- function(
  h,
  Y_subsample,
  X_subsample,
  D_subsample,
  Phi_subsample,
  tau,
  beta_D = h_to_beta(h,
                     Y = Y_subsample,
                     X = X_subsample,
                     D = D_subsample,
                     Phi = Phi_subsample)$beta_D,
  tolerance = 1e-9
) {
  # run a regression of Y_tilde ~ X and Phi using just the subsample
  Y_tilde_subsample <- Y_subsample - D_subsample %*% beta_D
  qr <- quantreg::rq(Y_tilde_subsample ~ Phi_subsample + X_subsample - 1, tau = tau)
  beta_Phi <- coef(qr)[seq_len(ncol(Phi))]
  beta_Phi_norm <- sum(abs(beta_Phi))
  # check if L1 norm of the beta_Phi's are less than some tolerance
  status <- beta_Phi_norm < tolerance
  list(
    beta_Phi = beta_Phi,
    norm = beta_Phi_norm,
    status = status
  )
}

# TODO: document
foc_membership_v3 <- function(
  h,
  Y_subsample,
  X_subsample,
  D_subsample,
  Phi_subsample,
  tau
) {

  m <- nrow(Y_subsample)
  p <- ncol(X_subsample) + ncol(Phi_subsample) # we concentrate out D_subsample
  stopifnot(length(h) == p)

  # compute relevant data
  yh <- Y_subsample[h]
  design <- cbind(X_subsample, Phi_subsample) # design matrix
  designh <- design[h, , drop = FALSE] # design matrix

  # compute b(h)
  bh <- solve(designh) %*% yh # bh is a p by 1 matrix
  stopifnot(nrow(bh) == p)
  stopifnot(ncol(bh) == 1)

  # compute resid_subsample
  beta <- h_to_beta(h, Y = Y_subsample, X = X_subsample, D = D_subsample, Phi =
                    Phi_subsample)
  beta_D <- beta$beta_D
  beta_X <- beta$beta_X
  resid_subsample <- Y_subsample - D_subsample %*% beta_D - X_subsample %*% beta_X

  # create indices that are not in h
  noth <- setdiff(seq_len(m), h)

  # compute xi(h, D)
  resid_noth <- matrix(resid_subsample[noth], ncol = 1) # (m - p) by 1 matrix
  ind_mat <- diag(tau - as.numeric(resid_noth < 0)) # (m - p) by (m - p) matrix
  design_noth <- design[noth, , drop = FALSE] # (m - p) by p matrix
  summand <- ind_mat %*% design_noth
  sum_summand <- matrix(1, nrow = 1, ncol = length(noth)) %*% summand
  xi <- t(sum_summand %*% solve(designh))
  stopifnot(nrow(xi) == p)
  stopifnot(ncol(xi) == 1)

  # DEBUG: see each entry of sum in xi object
  s_i <- vector("list", length(noth))
  for (i in seq_len(nrow(summand))) {
    s_i[[i]] <- summand[i, ] %*% solve(designh)
  }
  s <- do.call("rbind", s_i)

  # check if xi satisfies FOC inequality
  left <- matrix(-1 * tau %*% rep(1, p), ncol = 1)
  right <- matrix((1 - tau) %*% rep(1, p), ncol = 1)
  list(
    status = all((left <= xi) & (xi <= right)), # returns TRUE if both are true, else FALSE
    s = s, # return entries of xi object
    xi = xi,
    left = left,
    right = right
  )
}

#' @param beta_D Endogenous coefficients (vector of length p_D)
#' @param Y Dependent variable (vector of length n)
#' @param X Exogenous variable (including constant vector) (n by p_X matrix)
#' @param D Endogenous variable (n by p_D matrix)
#' @param Z Instrumental variable (n by p_Z matrix)
#' @param Phi Transformation of X and Z to be used in the program;
#'  defaults to the linear projection of D on X and Z (matrix with n rows)
#' @param tau Quantile (numeric)
#'
#' @return Results of QR(Y - D %*% beta_D ~ X + Phi) for a given \code{tau}
run_concentrated_qr <- function(beta_D, Y, X, D, Z, Phi = linear_projection(D, X, Z), tau) {
  quantreg::rq(Y - D %*% beta_D ~ X + Phi - 1, tau = tau)
}


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
      distance <- sum((neighbor - unrounded_center)^2)^(0.5)
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

### Propose h -- "Algorithm 3" -- Algorithm A -------------------------

# TODO: check email for any TODOs
# TODO: check reference code for any TODOs
# TODO: check how these functions work *together* to make it easy to use in terms of input/output flows

#' @param beta_proposal Vector of proposed coefficients beta_D and beta_X
#'  (vector of length p_D + p_X)
#' @param beta_hat Vector of beta_D and beta_X coefficients from IQR point
#'  estimate (vector of length p_D + p_X)
#' @param varcov_mat Asymptotic variance-covariance matrix; see
#'  \code{wald_varcov}
#'
#' @return density of limiting distribution of estimator evaluated at
#'  \code{beta_proposal}
density_wald <- function(beta_hat, beta_proposal, varcov_mat) {
  mvnfast::dmvn(beta_proposal, mu = beta_hat, sigma = varcov_mat)
}

#' @param residuals Residuals from IQR point estimation (vector of length n)
#' @param p_design Dimension of design matrix in QR (p_X + p_Phi) (numeric)
#' @param theta Hyperparameter (numeric)
#'
#' @return Named list
#'  \enumerate{
#'    \item \code{draws}: vector of length `n` with 1 if the index is in active
#'      basis and 0 otherwise
#'    \item \code{h_star}: Indices of the active basis
#'  }
propose_active_basis <- function(residuals, p_design, theta = 1) {
  weights <- exp(-1 * theta * residuals^2) # pointwise-operation
  n <- length(residuals)
  # 1 ball of each color in n urns;
  # probability of picking ball in urn i is equal to weights[i]
  # returns a vector with 1 if selected and 0 otherwise 
  # DONE: use `rmultinom`
  # draws <- BiasedUrn::rMFNCHypergeo(nran = 1, m = rep(1, n), n = p_design,
  #                                   odds = weights, precision = 1e-7)
  # Q: is using the while loop okay?
  while_bool <- TRUE
  while (while_bool) {
    # we pick p_design balls from length(weights) bins
    draws <- as.numeric(rmultinom(n = 1, size = p_design, prob = weights))
    # NOTE: possible to pick two balls from the same bin, hence the while loop will ensure all balls are from different bins
    # break out of while loop if we pick at most 1 ball from each bin
    if (identical(sort(unique(draws)), c(0, 1))) {
      while_bool <- FALSE
    }
  }
  h_star <- which(draws == 1)
  list(
    draws = draws,
    h_star = h_star # TODO: do we even need to store this?
  )
}

#' @param active_basis_draws Vector of length n with 1 if index is in active
#'  basis and 0 otherwise; see \code{propose_active_basis}
#' @param residuals Residuals from IQR point estimation (vector of length n)
#' @param p_design Dimension of design matrix in QR (p_X + p_Phi) (numeric)
#' @param theta Hyperparameter (numeric)
#'
#' @return density of Hypergeometric distribution evaluated at
#'  \code{density_active_basis}
density_active_basis <- function(active_basis_draws, residuals, p_design, theta = 1) {
  weights <- exp(-1 * theta * residuals^2)
  # n <- length(residuals)
  # # DONE: replace with product of weights
  # BiasedUrn::dMWNCHypergeo(x = active_basis_draws, m = rep(1, n), n = p_design,
  #                          odds = weights, precision = 1e-7)
  prod(weights[active_basis_draws == 1])
}

#' @param iterations Number of iterations through the MCMC algorithm
#' @param beta_X Coefficients on the exogenous variables (vector of length
#'  p_X); obtained from IQR point estimate
#' @param beta_D Coefficients on the endogeneous variables (vector of length
#'  p_D); obtained from IQR point estimate
#' @param Y Dependent variable (vector of length n)
#' @param X Exogenous variable (including constant vector) (n by p_X matrix)
#' @param D Endogenous variable (n by p_D matrix)
#' @param Z Instrumental variable (n by p_Z matrix)
#' @param Phi Transformation of X and Z to be used in the program;
#'  defaults to the linear projection of D on X and Z (matrix with n rows)
#' @param tau Quantile (numeric)
#' @param alpha Used for Hall and Sheather bandwidth; defaults to 0,1
#' @param theta Hyperparameter (numeric)
#' @param psi Hyperparameter; Coefficient on the variance-covariance matrix; defaults to 1 (numeric)
#' @param varcov_mat variance-covariance matrix of beta_D and beta_X
#' @param discard_burnin If TRUE (default), discard the first set of draws that are equivalent to the IQR MILP estimate
#'
#' @return A named list, each with a data frame:
#'  \enumerate{
#'    \item "beta": data frame where the rows are beta_D and beta_X while the
#'      columns are each iteration in the MCMC (excluding the burn-in period if \code{discard_burnin} is TRUE)
#'    \item "h": data frame where the rows are each entry in the active basis
#'      while the columns are each iteration in the MCMC (excluding the burn-in
#'      period if \code{discard_burnin} is TRUE)
# TODO: right now, I feed in beta_D to figure out h; maybe I should feed in h to figure out beta_D and beta_X? Also, maybe I should create beta_to_h(beta_D, ...) function
mcmc_active_basis <- function(iterations,
                              beta_X, # beta_X IQR point estimates
                              beta_D, # beta_D from IQR point estimates
                              Y, X, D, Z, Phi = linear_projection(D, X, Z),
                              tau,
                              alpha = 0.1,
                              theta = 1,
                              psi = 1, # coefficient on the varcov matrix
                              varcov_mat = NULL,
                              discard_burnin = TRUE) {
  qr <- run_concentrated_qr(beta_D = beta_D, Y = Y, X = X, D = D, Phi = Phi, tau = tau)
  residuals <- qr$residuals
  dual <- qr$dual
  initial_basis <- which(dual > 0 & dual < 1)
  u_vec <- runif(n = iterations) # create uniform RV for acceptance/rejection rule
  initial_draws <- vector("double", length = nrow(Y))
  initial_draws[initial_basis] <- 1
  beta_hat <- c(beta_D, beta_X)
  draws_current <- initial_draws
  beta_current <- beta_hat
  p_design <- ncol(X) + ncol(Phi)
  if (is.null(varcov_mat)) {
    varcov_mat <- wald_varcov(
      resid = residuals,
      alpha = alpha,
      tau = tau,
      D = D,
      X = X,
      Phi = Phi,
      psi = psi
    )$varcov
  }

  result <- vector("list", iterations) # preallocate space to store coefficients
  result_h <- vector("list", iterations) # preallocate space to store active basis
  result_prob <- vector("double", iterations) # store wald density of coefficients
  h_current <- initial_basis # TODO: refactor `draws` and `h`
  for (i in seq_len(iterations)) {
    u <- u_vec[[i]]
    h_proposal <- propose_active_basis(residuals, p_design = p_design, theta = theta)
    draws_proposal <- h_proposal$draws
    beta_proposal_full <- h_to_beta(h = h_proposal$h_star, Y = Y, X = X, D = D, Phi = Phi)
    beta_proposal <- c(beta_proposal_full$beta_D, beta_proposal_full$beta_X)
    # TODO: come up with better variable names
    wald_proposal <- density_wald(beta_hat = beta_hat,
                                  beta_proposal = beta_proposal,
                                  varcov_mat)
    geom_proposal <- density_active_basis(active_basis_draws = draws_proposal,
                                          residuals = residuals,
                                          p_design = p_design,
                                          theta = theta)
    if (i == 1) {
      wald_current <- density_wald(beta_hat = beta_hat,
                                   beta_proposal = beta_current,
                                   varcov_mat)
      geom_current <- density_active_basis(active_basis_draws = draws_current,
                                           residuals = residuals,
                                           p_design = p_design,
                                           theta = theta)
    }
    a <- wald_proposal / wald_current * geom_current / geom_proposal
    if (u < a) { # accept
      beta_current <- beta_proposal
      draws_current <- draws_proposal
      h_current <- h_proposal$h_star
      wald_current <- wald_proposal
      geom_current <- geom_proposal
    }
    result[[i]] <- beta_current
    result_h[[i]] <- h_current
    result_prob[[i]] <- wald_current
  }
  # each row is a coefficient, each column is one iteration of MCMC
  result_df <- do.call(cbind, result)
  rownames(result_df) <- c(paste0("beta_D", seq_len(ncol(D))), paste0("beta_X", seq_len(ncol(X))))
  result_h_df <- do.call(cbind, result_h)
  if (discard_burnin) {
    # find where stationary distribution begins
    stationary_begin <- min(which(result_df[1, ] != beta_hat[1]))
    # remove burn-in period
    result_df <- result_df[, stationary_begin:ncol(result_df)]
    result_h_df <- result_h_df[, stationary_begin:ncol(result_h_df)]
    result_prob <- result_prob[stationary_begin:length(result_prob)]
  }
  list("beta" = result_df, "h" = result_h_df, "prob" = result_prob)
}

### Propose first subsample -- "First Approach" -------------------------

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
                                  seed = Sys.date()) {
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
  if (distance_method == 1) {
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
      exp(-gamma * tmp^l_power)
    }
    okay <- setdiff(seq_len(n), h)
    ones_current <- which(current_subsample[okay] == 1)
    s_i_current <- s_i[, ones_current]
    # compute norm of sum of s_i_current
    current_distance <- distance_function(s_i_current)
    # transform (e.g., exponentiate) "distance"
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
    out_subsample[[i]] <- current_subsample
    out_distance[[i]] <- current_distance
    out_distance_prob[[i]] <- current_distance_prob
  }

  # TODO: compute foc_membership?

  list(
    status = "OKAY",
    a_log = out_a_log,
    subsample = out_subsample,
    distance = out_distance,
    distance_prob = out_distance_prob, # use in main MCMC
    record = out_record
  )
}

# find_center_subsample_polytope -------------------------

# The center of the subsample simplex, {D \in [0,1]^n s.t. sum(D_i) = m}, should
# be rep(m/n, length = n). This function creates a program that verifies this.
# Examples:
# tmp <- find_center_subsample_polytope(n = 4, m = 2)
# tmp$center # should be roughly m/n = 1/2
# Q: The log variant of the program showing up as infeasible. Why?
find_center_subsample_polytope <- function(
  n, m,
  gencontype = "power",
  a = ifelse(gencontype == "power", 0.5, exp(1)),
  params = list()
) {
  vertices_index <- combn(n, m)
  num_vertices <- ncol(vertices_index)
  stopifnot(num_vertices == choose(n, m))
  n_ones <- rep(0, n)
  vertices <- t(apply(vertices_index, 2, function(i) {
    n_ones[i] <- 1
    n_ones
  }))
  stopifnot(ncol(vertices) == n)
  stopifnot(nrow(vertices) == num_vertices)
  vertices_i <- c(t(vertices)) # collapse by rows
  stopifnot(length(vertices_i) == num_vertices * n)
  slack_i <- diag(1, length(vertices_i))
  zero_mat <- diag(0, length(vertices_i))
  dv_i <- do.call("rbind", rep(list(diag(1, nrow = n)), length = num_vertices))
  define_slack <- cbind(dv_i,      # coordinates of center
                        -slack_i,  # positive part of slack
                        slack_i,   # negative part of slack
                        zero_mat,  # f(pos slack), where f = sqrt by default
                        zero_mat)  # f(neg slack), where f = sqrt by default
  num_decision_vars <- ncol(define_slack)

  sum_A <- c(rep(1, n), rep(0, num_decision_vars - n))
  sum_rhs <- m

  model <- list()
  model$A <- rbind(define_slack, sum_A)
  model$rhs <- c(vertices_i, sum_rhs)
  model$sense <- "="
  model$modelsense <- "max"
  model$lb <- rep(0, num_decision_vars)
  # model$ub <- c(rep(1, n), rep(Inf, num_decision_vars - n))
  model$ub <- c(rep(1, num_decision_vars)) # slack in each index must be <= 1
  model$obj <- c(
    rep(0, n),
    rep(0, 2 * length(vertices_i)),
    rep(1, 2 * length(vertices_i))
  )

  xstart <- n
  ystart <- n + 2 * length(vertices_i)
  gencon <- vector("list", 2 * length(vertices_i))
  for (j in seq_len(2 * length(vertices_i))) {
    yvar <- ystart + j
    xvar <- xstart + j
    gencon[[j]]$xvar <- xvar
    gencon[[j]]$yvar <- yvar
    gencon[[j]]$a <- a
  }
  if (gencontype == "log") {
    model$genconloga <- gencon
  } else if (gencontype == "power") {
    model$genconpow <- gencon
  }

  sol <- gurobi::gurobi(model, params)

  list(
    model = model,
    sol = sol,
    center = sol$x[seq_len(n)]
  )
}

# find_subsample_in_polytope -------------------------

# goal: given an active basis, find a subsample inside the polytope
# Q: is it possible for there to be multiple solutions? In which case, maybe this doesn't provide a deterministic map from the active basis to the subsample, depending on how Gurobi breaks ties.
#' Find subsample in the FOC polytope
#'
#' It is useful to start our random walk at a subsample that is inside the FOC
#' polytope. To find this polytope, we solve a linear program that aims to
#' minimize the absolute L-1 norm of the weighted average of the xi_i objects.
#' The weights are inside the unit interval and sum to the size of our
#' subsample. Most of the weights will be integral, but some will be between 0
#' and 1. We will round the non-integral solutions so they are integra. If the
#' weight for an observaion is 1, then that observation belongs in the
#' subsample. Otherwise, the observation doesn't belong in the subsample. Due
#' to the rounding, it isn't guaranteed that the resulting subsample belongs in
#' the FOC polytope, but it shouldn't be too far away from the polytope.
#' Note that this program excludes the observations in the active basis -- we
#' must add them back in after-the-fact.
#'
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
#' @param subsample_size Size of subsample
#' @param params Named list of parameters to send to Gurobi
#' @param type Either "C" or "B" to denote if our "omega" variables are
#'  continuous or binary; If "C", then the solution is the continuous center
#'  of our FOC polytope; If "B", then the solution is the actual center
# We need something inside the polytope. This continuous method with our
# arbitrary rounding doesn't guarantee this. Our workaround:
# 1. exhaustive search over non-integral solutions
# 2. solve the program as an MILP
# TODO: solve the program as an MILP and see if the result is inside the FOC
# polytope
# TODO: code up improved linear program to find DC (ask GP for the note)
find_subsample_in_polytope <- function(
  h,
  Y, X, D, Z, Phi = linear_projection(D, X, Z),
  tau,
  beta_D_proposal = NULL,
  beta_X_proposal = NULL,
  subsample_size,
  params = list(OutputFlag = 0),
  type = "C"
) {

  n <- nrow(Y)
  p <- length(h)

  # get beta_X_proposal and beta_D_proposal

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
  s <- t(do.call(rbind, s_i))
  xi_mat <- s[, setdiff(seq_len(n), h)]
  stopifnot(nrow(xi_mat) == p)
  stopifnot(ncol(xi_mat) == n - p)

  # Decision variables in order from left/top to right/bottom:
  # 1. omega --- (n - p) by 1
  # 2. xi --- (p by 1)
  # 3. ximinus --- (p by 1)
  # 4. xiplus --- (p by 1)

  num_omega <- n - p
  num_xi <- p
  num_ximinus <- num_xi
  num_xiplus <- num_xi
  num_decision_vars <- num_omega + num_xi + num_ximinus + num_xiplus

  model <- list()
  model$modelsense <- "min"
  model$lb <- c(
    rep(0, num_omega),
    rep(-Inf, num_xi),
    rep(0, num_ximinus),
    rep(0, num_xiplus)
  )
  model$ub <- c(
    rep(1, num_omega),
    rep(Inf, num_xi),
    rep(Inf, num_ximinus),
    rep(Inf, num_xiplus)
  )
  model$obj <- c(
    rep(0, num_omega),
    rep(0, num_xi),
    rep(1, num_ximinus),
    rep(1, num_xiplus)
  )
  model$vtype <- c(
    rep(type, num_omega),
    rep("C", num_decision_vars - num_omega)
  )

  # define xiplus and ximinus
  tmp <- diag(1, nrow = num_xi)
  zero_mat <- diag(0, nrow = num_omega)
  xisign_A <- expand_matrix(
    cbind(tmp, -tmp, tmp), # TODO: maybe the final minus should be a plus?
    newrow = num_xi,
    newcol = num_decision_vars,
    row_direction = "bottom",
    col_direction = "left"
  )
  xisign_sense <- rep("=", num_xi)
  xisign_rhs <- rep(0, num_xi)

  # define xi
  tmp <- cbind(xi_mat, diag(-1, nrow = num_xi))
  xi_A <- expand_matrix(
    tmp,
    newrow = num_xi,
    newcol = num_decision_vars,
    row_direction = "bottom",
    col_direction = "right"
  )
  xi_sense <- rep("=", num_xi)
  xi_rhs <- rep(0, num_xi)

  # sum of omega
  omega_A <- matrix(c(
    rep(1, num_omega), rep(0, num_decision_vars - num_omega)
  ), nrow = 1)
  omega_sense <- "="
  omega_rhs <- subsample_size - p

  model$A <- rbind(xisign_A, xi_A, omega_A)
  model$rhs <- c(xisign_rhs, xi_rhs, omega_rhs)
  model$sense <- c(xisign_sense, xi_sense, omega_sense)

  sol <- gurobi::gurobi(model, params)
  status <- sol$status
  omega <- sol$x[seq_len(num_omega)]

  if (status != "OPTIMAL") {
    return(list(
      status = status,
      status_message = paste("Gurobi status:", status),
      model = model,
      sol = sol,
      omega = omega
    ))
  }

  # turn continuous solution into an integral one
  if (type == "C") {
    current_sum <- length(which(omega == 1)) # how many integral 1's do we have?
    remaining <- subsample_size - p - current_sum # how many need to be switched?
    to_be_rounded <- omega > 0 & omega < 1 # which can we switch?
    # NOTE: `rank` works better than `order` when there are ties
    max_indices <- rank(-omega, ties.method = "random") # 1 = largest
    # switch the largest numbers that aren't 1
    switch <- which(max_indices <= current_sum + remaining & max_indices > current_sum)
    omega_mod <- omega
    omega_mod[switch] <- 1
    omega_mod[which(omega_mod < 1)] <- 0
    omega_mod <- round(omega_mod, 0)
    stopifnot(all.equal(sum(omega_mod), subsample_size - p))
  } else if (type == "B") {
    omega_mod <- omega
  }

  list(
    model = model,
    sol = sol,
    status = status,
    omega = omega,
    omega_mod = omega_mod,
    xi_mat = xi_mat
  )
}

### find_chebyschev_center -------------------------
# TODO: document
find_chebyschev_center <- function(
  h,
  Y, X, D, Z, Phi = linear_projection(D, X, Z),
  tau,
  beta_D_proposal = NULL,
  beta_X_proposal = NULL,
  subsample_size,
  params = list(OutputFlag = 0),
  type = "C",
  l_norm = 2
) {

  n <- nrow(Y)
  p <- length(h)

  # get beta_X_proposal and beta_D_proposal

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
  s <- t(do.call(rbind, s_i))
  xi_mat <- s[, setdiff(seq_len(n), h)]
  stopifnot(nrow(xi_mat) == p)
  stopifnot(ncol(xi_mat) == n - p)

  # Decision variables in order from left/top to right/bottom:
  # 1. omega --- (n - p) by 1
  # 2. r --- 1 by 1

  num_omega <- n - p
  num_r <- 1
  num_decision_vars <- num_omega + num_r

  model <- list()
  model$modelsense <- "max"
  model$lb <- rep(0, num_decision_vars)
  model$ub <- c(
    rep(1, num_omega),
    rep(Inf, num_r)
  )
  model$obj <- c(
    rep(0, num_omega),
    rep(1, num_r)
  )
  model$vtype <- c(
    rep(type, num_omega),
    rep("C", num_r)
  )

  # sum of omega
  omega_A <- matrix(c(
    rep(1, num_omega), rep(0, num_decision_vars - num_omega)
  ), nrow = 1)
  omega_sense <- "="
  omega_rhs <- subsample_size - p

  # FOC boundary
  right_A <- vector("list", p)
  left_A <- vector("list", p)
  for (j in seq_len(p)) {
    xi_j <- xi_mat[j, ]
    xi_j_norm <- sum(abs(xi_j)^l_norm) ^ (1 / l_norm)
    right_A[[j]] <- c(xi_mat[j, ], xi_j_norm)
    left_A[[j]] <- c(-xi_mat[j, ], xi_j_norm)
  }
  right_A <- do.call(rbind, right_A)
  left_A <- do.call(rbind, left_A)
  foc_A <- rbind(right_A, left_A)
  foc_sense <- rep("<=", 2 * p)
  foc_rhs <- c(rep(1-tau, p), rep(tau, p))

  # constraints
  model$A <- rbind(omega_A, foc_A)
  model$sense <- c(omega_sense, foc_sense)
  model$rhs <- c(omega_rhs, foc_rhs)

  sol <- gurobi::gurobi(model, params)
  status <- sol$status
  omega <- sol$x[seq_len(num_omega)]

  if (status != "OPTIMAL") {
    return(list(
      status = status,
      status_message = paste("Gurobi status:", status),
      model = model,
      sol = sol,
      omega = omega,
      xi_mat = xi_mat
    ))
  }

  # turn continuous solution into an integral one
  if (type == "C") {
    current_sum <- length(which(omega == 1)) # how many integral 1's do we have?
    remaining <- subsample_size - p - current_sum # how many need to be switched?
    to_be_rounded <- omega > 0 & omega < 1 # which can we switch?
    # NOTE: `rank` works better than `order` when there are ties
    max_indices <- rank(-omega, ties.method = "random") # 1 = largest
    # switch the largest numbers that aren't 1
    switch <- which(max_indices <= current_sum + remaining & max_indices > current_sum)
    omega_mod <- omega
    omega_mod[switch] <- 1
    omega_mod[which(omega_mod < 1)] <- 0
    omega_mod <- round(omega_mod, 0)
    stopifnot(all.equal(sum(omega_mod), subsample_size - p))
  } else if (type == "B") {
    omega_mod <- omega
  }

  list(
    model = model,
    sol = sol,
    status = status,
    omega = omega,
    omega_mod = omega_mod,
    xi_mat = xi_mat
  )
}

### find_center_repellent -------------------------
# TODO: document
find_center_repellent <- function(
  h,
  Y, X, D, Z, Phi = linear_projection(D, X, Z),
  tau,
  beta_D_proposal = NULL,
  beta_X_proposal = NULL,
  subsample_size,
  params = list(OutputFlag = 0),
  type = "C",
  gencontype = "power", # "power" or "log" or "max"
  a = ifelse(gencontype == "power", 0.5, exp(1)),
  simplex_repel = TRUE, # repel away from facets of the simplex
  foc_repel = TRUE # repel away from the FOC conditions
) {

  n <- nrow(Y)
  p <- length(h)

  if (!simplex_repel & !foc_repel) {
    stop("At least one of `foc_repel` and `simplex_repel` must be TRUE. Both can't be false.")
  }
  if (simplex_repel & ((n - p) - 1) <= 0) {
    stop("facet_center denominator must be positive i.e. we need n - p - 1 > 0")
  }

  xi_mat <- compute_xi_i(
    h = h,
    Y = Y, X = X, D = D, Z = Z, Phi = Phi,
    tau = tau,
    beta_D_proposal = beta_D_proposal,
    beta_X_proposal = beta_X_proposal
  )

  # Decision variables in order from left/top to right/bottom:
  # 1. omega --- (n - p) by 1
  # 2. right_slack --- p by 1
  # 3. left_slack --- p by 1
  # 4. right_slack_transformed --- p by 1
  # 5. left_slack_transformed --- p by 1
  # 6. simplex_slack --- 2n by 1
  # 7. simplex_slack_transformed --- 2n by 1
  # NOTE: when gencontype = "max", we have one additional decision variable
  # which is the max of the slack variables; we aren't considering this variable
  # in `num_decision_vars`
  num_omega <- n - p
  num_right_slack <- p
  num_left_slack <- p
  num_right_slack_transform <- num_right_slack
  num_left_slack_transform <- num_left_slack
  num_simplex_slack <- 2 * (n - p)
  num_simplex_slack_transform <- 2 * (n - p)
  num_decision_vars <- num_omega + num_right_slack + num_left_slack +
    num_right_slack_transform + num_left_slack_transform +
    num_simplex_slack + num_simplex_slack_transform

  model <- list()
  model$modelsense <- ifelse(gencontype == "max", "min", "max")
  model$lb <- c(
    rep(0, num_omega),
    rep(0, num_right_slack),
    rep(0, num_left_slack),
    rep(-Inf, num_right_slack_transform),
    rep(-Inf, num_right_slack_transform),
    rep(0, num_simplex_slack),
    rep(-Inf, num_simplex_slack_transform)
  )
  model$ub <- c(
    rep(1, num_omega),
    rep(Inf, num_right_slack + num_left_slack),
    rep(Inf, num_right_slack_transform + num_left_slack_transform),
    rep(Inf, num_simplex_slack),
    rep(Inf, num_simplex_slack_transform)
  )
  model$obj <- c(
    rep(0, num_omega + num_right_slack + num_left_slack),
    rep(1, num_right_slack_transform + num_left_slack_transform),
    rep(0, num_simplex_slack),
    rep(1, num_simplex_slack_transform)
  )
  model$vtype <- c(
    rep(type, num_omega),
    rep("C", num_decision_vars - num_omega)
  )

  # sum of omega
  omega_A <- matrix(c(
    rep(1, num_omega), rep(0, num_decision_vars - num_omega)
  ), nrow = 1)
  omega_sense <- "="
  omega_rhs <- subsample_size - p

  # FOC slack variables
  right_A <- vector("list", p)
  left_A <- vector("list", p)
  for (j in seq_len(p)) {
    xi_j <- xi_mat[j, ]
    zeros <- rep(0, p)
    ones <- zeros
    ones[j] <- 1
    right_A[[j]] <- c(xi_j, ones, zeros, zeros, zeros)
    left_A[[j]] <- c(-xi_j, zeros, ones, zeros, zeros)
  }
  right_A <- do.call(rbind, right_A)
  left_A <- do.call(rbind, left_A)
  foc_A <- rbind(right_A, left_A)
  foc_A <- expand_matrix(
    foc_A, newrow = nrow(foc_A), newcol = num_decision_vars,
    row_direction = "bottom", col_direction = "right"
  )
  foc_sense <- rep("=", 2 * p)
  foc_rhs <- c(rep(1 - tau, p), rep(tau, p))

  # transform FOC slack variables
  xstart <- num_omega
  ystart <- xstart + num_right_slack + num_left_slack
  foc_gencon <- vector("list", num_left_slack_transform + num_right_slack_transform)
  for (j in seq_len(2 * p)) {
    yvar <- ystart + j
    xvar <- xstart + j
    foc_gencon[[j]]$xvar <- xvar
    foc_gencon[[j]]$yvar <- yvar
    foc_gencon[[j]]$a <- a
  }

  # simplex slack variables
  simplex_rhs <- vector("double", num_simplex_slack)
  simplex_lhs <- vector("list", num_simplex_slack)
  counter <- 0
  subsample_center <- rep((subsample_size - p) / (n - p), n - p)
  for (j in c(0, 1)) {
    for (i in seq_len(n - p)) {
      counter <- counter + 1
      e_ij <- rep(0, num_simplex_slack)
      e_ij[[counter]] <- 1
      facet_center <- rep(
        ifelse(j == 0, (subsample_size - p) / ((n - p) - 1), ((subsample_size - p) - 1) / ((n - p) - 1)),
        n - p
      )
      facet_center[i] <- j
      normal_vec <- subsample_center - facet_center
      normal_unit <- normal_vec / sqrt(sum(normal_vec^2))
      simplex_rhs[[counter]] <- -sum(normal_unit * facet_center)
      # if (j == 0) {
      #   normal_A <- c(normal_unit, rep(0, n))
      # } else {
      #   normal_A <- c(rep(0, n), normal_unit)
      # }
      simplex_lhs[[counter]] <- c(
        -normal_unit, # num_omega
        rep(0, num_left_slack + num_right_slack +
            num_left_slack_transform + num_right_slack_transform),
        e_ij, # num_simplex_slack
        rep(0, num_simplex_slack_transform) # num_simplex_slack_transform
      )
    }
  }
  simplex_A <- do.call(rbind, simplex_lhs)
  stopifnot(ncol(simplex_A) == num_decision_vars)
  stopifnot(nrow(simplex_A) == num_simplex_slack)
  simplex_sense <- rep("=", num_simplex_slack)

  # transform simplex slack variables
  xstart <- num_omega + num_left_slack + num_right_slack +
    num_left_slack_transform + num_right_slack_transform
  ystart <- xstart + num_simplex_slack
  simplex_gencon <- vector("list", num_simplex_slack_transform)
  for (j in seq_len(num_simplex_slack_transform)) {
    yvar <- ystart + j
    xvar <- xstart + j
    simplex_gencon[[j]]$xvar <- xvar
    simplex_gencon[[j]]$yvar <- yvar
    simplex_gencon[[j]]$a <- a
  }

  # constraints

  if (gencontype == "max") {
    # add one more decision variable for the max slack variable
    omega_A <- cbind(omega_A, rep(0, nrow(omega_A)))
    foc_A <- cbind(foc_A, rep(0, nrow(foc_A)))
    simplex_A <- cbind(simplex_A, rep(0, nrow(simplex_A)))
    model$lb <- c(model$lb, 0)
    model$ub <- c(model$ub, Inf)
    model$vtype <- c(model$vtype, "C")


    # max >= slack => -slack + max >= 0
    foc_slack_tmp <- diag(-1, num_right_slack + num_left_slack)
    max_foc_slack <- cbind(
      matrix(0, nrow = nrow(foc_slack_tmp), ncol = num_omega),
      foc_slack_tmp,
      matrix(0, nrow = nrow(foc_slack_tmp), ncol = num_right_slack_transform +
        num_left_slack_transform + num_simplex_slack +
        num_simplex_slack_transform),
      rep(1, length = nrow(foc_slack_tmp))
    )
    max_foc_slack_sense <- rep(">=", length = nrow(foc_slack_tmp))
    max_foc_slack_rhs <- rep(0, length = nrow(foc_slack_tmp))

    simplex_slack_tmp <- diag(-1, num_simplex_slack_transform)
    max_simplex_slack <- cbind(
      matrix(0, nrow = nrow(simplex_slack_tmp), ncol = num_omega +
        num_right_slack + num_left_slack + num_right_slack_transform +
        num_left_slack_transform),
      simplex_slack_tmp,
      matrix(0, nrow = nrow(simplex_slack_tmp), ncol = num_simplex_slack_transform),
      rep(1, length = nrow(simplex_slack_tmp))
    )
    max_simplex_slack_sense <- rep(">=", length = nrow(simplex_slack_tmp))
    max_simplex_slack_rhs <- rep(0, length = nrow(simplex_slack_tmp))
    model$obj <- c(
      rep(0, num_decision_vars),
      1
    )
  } else {
    max_foc_slack <- c()
    max_foc_slack_sense <- c()
    max_foc_slack_rhs <- c()
    max_simplex_slack <- c()
    max_simplex_slack_sense <- c()
    max_simplex_slack_rhs <- c()
  }

  model$A <- rbind(omega_A, foc_A, simplex_A)
  model$sense <- c(omega_sense, foc_sense, simplex_sense)
  model$rhs <- c(omega_rhs, foc_rhs, simplex_rhs)
  gencon <- append(foc_gencon, simplex_gencon)

  if (!foc_repel) { # no foc, only simplex
    model$A <- rbind(omega_A, simplex_A, max_simplex_slack)
    model$sense <- c(omega_sense, simplex_sense, max_simplex_slack_sense)
    model$rhs <- c(omega_rhs, simplex_rhs, max_simplex_slack_rhs)
    gencon <- simplex_gencon
    model$obj <- c(
      rep(0, num_omega + num_right_slack + num_left_slack),
      rep(0, num_right_slack_transform + num_left_slack_transform),
      rep(0, num_simplex_slack),
      rep(1, num_simplex_slack_transform)
    )
    if (gencontype == "max") {
      model$obj <- c(
        rep(0, num_omega + num_right_slack + num_left_slack),
        rep(0, num_right_slack_transform + num_left_slack_transform),
        rep(0, num_simplex_slack),
        rep(0, num_simplex_slack_transform),
        1
      )
    }
  }
  if (!simplex_repel) { # no simplex, only foc
    model$A <- rbind(omega_A, foc_A, max_foc_slack)
    model$sense <- c(omega_sense, foc_sense, max_foc_slack_sense)
    model$rhs <- c(omega_rhs, foc_rhs, max_foc_slack_rhs)
    gencon <- foc_gencon
    model$obj <- c(
      rep(0, num_omega + num_right_slack + num_left_slack),
      rep(1, num_right_slack_transform + num_left_slack_transform),
      rep(0, num_simplex_slack),
      rep(0, num_simplex_slack_transform)
    )
    if (gencontype == "max") {
      model$obj <- c(
        rep(0, num_omega + num_right_slack + num_left_slack),
        rep(0, num_right_slack_transform + num_left_slack_transform),
        rep(0, num_simplex_slack),
        rep(0, num_simplex_slack_transform),
        1
      )
    }
  }

  if (gencontype == "log") {
    model$genconloga <- gencon
  } else if (gencontype == "power") {
    model$genconpow <- gencon
  }

  sol <- gurobi::gurobi(model, params)
  status <- sol$status
  omega <- sol$x[seq_len(num_omega)]

  if (status != "OPTIMAL") {
    return(list(
      status = status,
      status_message = paste("Gurobi status:", status),
      model = model,
      sol = sol,
      omega = omega,
      xi_mat = xi_mat
    ))
  }

  # turn continuous solution into an integral one
  if (type == "C") {
    current_sum <- length(which(omega == 1)) # how many integral 1's do we have?
    remaining <- subsample_size - p - current_sum # how many need to be switched?
    to_be_rounded <- omega > 0 & omega < 1 # which can we switch?
    # NOTE: `rank` works better than `order` when there are ties
    max_indices <- rank(-omega, ties.method = "random") # 1 = largest
    # switch the largest numbers that aren't 1
    switch <- which(max_indices <= current_sum + remaining & max_indices > current_sum)
    omega_mod <- omega
    omega_mod[switch] <- 1
    omega_mod[which(omega_mod < 1)] <- 0
    omega_mod <- round(omega_mod, 0)
    stopifnot(all.equal(sum(omega_mod), subsample_size - p))
  } else if (type == "B") {
    omega_mod <- omega
  }

  list(
    model = model,
    sol = sol,
    status = status,
    omega = omega,
    omega_mod = omega_mod,
    xi_mat = xi_mat
  )
}

### compute_xi_i -------------------------
# TODO: incorporate this function into first_approach* and
# find_subsample_in_polytope
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
#'
#' @return Matrix of dimension p by (n-p) of xi_i's not including the i's in the active basis
compute_xi_i <- function(h,
                         Y, X, D, Z, Phi = linear_projection(D, X, Z),
                         tau,
                         beta_D_proposal = NULL,
                         beta_X_proposal = NULL) {
  n <- nrow(Y)
  p <- length(h)

  # get beta_X_proposal and beta_D_proposal
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
  s <- t(do.call(rbind, s_i))
  xi_mat <- s[, setdiff(seq_len(n), h), drop = FALSE]
  stopifnot(nrow(xi_mat) == p)
  stopifnot(ncol(xi_mat) == n - p)

  xi_mat
}

### ot -------------------------
#' Compute OT between identical uniform distributions
#'
#' Compute transport map between U(1,...,n-p) and U(1,...,n-p) where C(i, j) =
#' norm difference of pre[, i] and post[, j] where i and j are between 1 and
#' n-p.
#'
#' We use Gurobi to solve the OT problem if \code{method} is "gurobi".
#' We use the transport package to solve the OT problem if \code{method} is "transport".
ot <- function(pre, post, params = list(OutputFlag = 0), method = "gurobi") {
  n_minus_p <- ncol(pre)
  stopifnot(n_minus_p == ncol(post))

  # prelims
  num_decision_vars <- n_minus_p^2

  # compute cost matrix
  c_ij <- matrix(0, nrow = n_minus_p, ncol = n_minus_p)
  for (col in seq_len(n_minus_p)) {
    # get diagonals
    c_ij[col, col] <- sum(abs(pre[, col] - post[, col]))
    for (i in seq_len(n_minus_p - col)) {
      row <- i + col
      # create a lower-triangular matrix with the cost
      c_ij[row, col] <- sum(abs(pre[, row] - post[, col]))
    }
  }
  # fill out the cost matrix
  c_ij[upper.tri(c_ij)] <- c_ij[lower.tri(c_ij)]

  if (tolower(method) == "gurobi") {
    # create constraints
    const_pre <- vector("list", length = n_minus_p)
    const_post <- vector("list", length = n_minus_p)
    for (i in seq_len(n_minus_p)) {
      zeros_left <- matrix(0, nrow = n_minus_p, ncol = i - 1)
      ones <- matrix(1, nrow = n_minus_p, ncol = 1)
      zeros_right <- matrix(0, nrow = n_minus_p, ncol = n_minus_p - i)
      a_mat <- cbind(zeros_left, ones, zeros_right)
      const_pre[[i]] <- c(a_mat)
      const_post[[i]] <- c(t(a_mat))
    }
    a_mat_pre <- do.call(rbind, const_pre)
    a_mat_post <- do.call(rbind, const_post)
    a_mat <- rbind(a_mat_pre, a_mat_post)
    stopifnot(ncol(a_mat) == num_decision_vars)
    stopifnot(nrow(a_mat) == n_minus_p * 2)

    # create program
    model <- list()
    model$obj <- c(c_ij) # turn into vector (go down each column)
    model$A <- a_mat
    model$sense <- rep("=", length = n_minus_p * 2)
    model$rhs <- rep(1, length = n_minus_p * 2)
    model$vtype <- rep("C", num_decision_vars)

    sol <- gurobi::gurobi(model, params)
    status <- sol$status
    t_ij <- matrix(sol$x, nrow = n_minus_p, ncol = n_minus_p)
    map <- apply(t_ij, 1, function(x) which(x == 1))
  } else if (tolower(method) == "transport") {
    unif <- rep(1, length = n_minus_p)
    sol <- transport::transport(unif, unif, costm = c_ij)
    model <- NA
    status <- NA
    map <- sol$to
  }

  list(
    model = model, # gurobi-specific
    status = status, # gurobi-specific
    sol = sol,
    c_ij = c_ij,
    map = map
  )
}
