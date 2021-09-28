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
  B <- solve(designh, D[h, ])

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
#' @return TRUE if the subsample satisfies FOC conditions; FALSE otherwise
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
  reg <- quantreg::rq(y_tilde_subsample ~ X_subsample + Phi_subsample - 1,
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

  # check if xi satisfies FOC inequality
  left <- matrix(-1 * tau %*% rep(1, p), ncol = 1)
  right <- matrix((1 - tau) %*% rep(1, p), ncol = 1)
  all((left <= xi) & (xi <= right)) # returns TRUE if both are true, else FALSE
}
# TODO: Sanity Check -- check for the false case

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
  # TODO: use `rmultinom`
  draws = BiasedUrn::rMFNCHypergeo(nran = 1, m = rep(1, n), n = p_design,
                                   odds = weights, precision = 1e-7)
  h_star = which(draws == 1)
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
  n <- length(residuals)
  # TODO: replace with product of weights
  BiasedUrn::dMWNCHypergeo(x = active_basis_draws, m = rep(1, n), n = p_design,
                           odds = weights, precision = 1e-7)
}

#' @param iterations Number of iterations through the MCMC algorithm
#' @param initial_beta_hat Initial value for c(beta_D, beta_X), e.g., coefficients from IQR point estimate
#' @param initial_draws Initial value for active basis, represented as a vector
#'  of 0's and 1's (1 if index is in active basis, 0 otherwise) (vector of
#'  length n) # Q: does this have to be the active basis associated with
#'  initial_beta_hat?
#' @param residuals Residuals from IQR point estimation (vector of length n)
#' @param theta Hyperparameter (numeric)
# TODO: update documentation
mcmc_active_basis <- function(iterations,
                              beta_X, # beta_X IQR point estimates
                              beta_D, # beta_D from IQR point estimates
                              Y, X, D, Z, Phi = linear_projection(D, X, Z),
                              tau,
                              alpha = 0.1,
                              theta = 1,
                              discard_burnin = TRUE) {
  qr <- run_concentrated_qr(beta_D = beta_D, Y = Y, X = X, D = D, Phi = Phi, tau = tau)
  residuals <- qr$residuals
  dual <- qr$dual
  initial_basis <- which(dual > 0 & dual < 1)
  initial_draws <- seq_len(nrow(Y))
  initial_draws[initial_basis] <- 1
  beta_hat <- c(beta_D, beta_X)
  draws_current <- initial_draws
  beta_current <- beta_hat
  p_design <- ncol(X) + ncol(Phi)
  u_vec <- runif(n = iterations) # create uniform RV for acceptance/rejection rule
  varcov_mat <- wald_varcov(
    resid = residuals,
    alpha = alpha,
    tau = tau,
    D = D,
    X = X,
    Phi = Phi
  )$varcov

  result <- vector("list", iterations) # preallocate space to store coefficients
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
    wald_current <- density_wald(beta_hat = beta_hat,
                                 beta_proposal = beta_current,
                                 varcov_mat) # TODO: presumably, we would have already computed this...right?
    geom_proposal <- density_active_basis(active_basis_draws = draws_proposal,
                                          residuals = residuals,
                                          p_design = p_design,
                                          theta = theta)
    geom_current <- density_active_basis(active_basis_draws = draws_current,
                                         residuals = residuals,
                                         p_design = p_design,
                                         theta = theta)
    a <- wald_proposal / wald_current * geom_current / geom_proposal
    if (u < a) { # accept
      beta_current <- beta_proposal
      draws_current <- draws_proposal
    }
    result[[i]] <- beta_current # accept => we save beta_proposal, otherwise we save beta_current
  }
  # each row is a coefficient, each column is one iteration of MCMC
  result_df <- do.call(cbind, result)
  if (discard_burnin) {
    # find where stationary distribution begins
    stationary_begin <- min(which(result_df[1, ] != beta_hat[1]))
    # remove burn-in period
    result_df <- result_df[, stationary_begin:ncol(result_df)]
  }
  result_df
}
