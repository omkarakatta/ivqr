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
  reg <- quantreg::rq(y_tilde_subsample ~ X_subsample + Phi_subsample - 1,
                      tau = tau)
  resid_subsample <- reg$residuals
  # TODO: I need to compute the residuals after *assuming* beta_Phi = 0

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
    xi = xi
  )
}
# TODO: Sanity Check -- check for the false case

# TODO: document
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
  # run a regression of Y_tilde ~ X and Phi
  # using just the subsample
  # check if L1 norm of the beta_Phi's are less than some tolerance
  # with default tolerance being 1e-9

  Y_tilde_subsample <- Y_subsample - D_subsample %*% beta_D
  qr <- quantreg::rq(Y_tilde_subsample ~ Phi_subsample + X_subsample - 1, tau = tau)
  beta_Phi <- coef(qr)[seq_len(ncol(Phi))]
  beta_Phi_norm <- sum(abs(beta_Phi))
  status <- beta_Phi_norm < tolerance
  list(
    beta_Phi = beta_Phi,
    norm = beta_Phi_norm,
    status = status
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
                              discard_burnin = TRUE) {
  qr <- run_concentrated_qr(beta_D = beta_D, Y = Y, X = X, D = D, Phi = Phi, tau = tau)
  residuals <- qr$residuals
  dual <- qr$dual
  initial_basis <- which(dual > 0 & dual < 1)
  initial_draws <- vector("double", length = nrow(Y))
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
#'  arugment to determine \code{beta_X_proposal}
#' @param gamma,l Hyperparameters
#'
#' @return Named list
#'  \enumerate{
#'    \item \code{subsample_set}: set of indices in subsample
#'    \item \code{prob}: probability of proposing \code{subsample_set} (except
#'      for observations already in active basis)
#'    \item \code{log_prob}: log of \code{prob}
#'  }
# Q: Should the beta_*_proposal correspond to the same coefficients as h_to_beta(h)? If so, I don't even need the beta_*_proposal in the arguments. I can just use the h_to_beta(h) to get these coefficients...right?
first_approach <- function(Y, X, D, Z, Phi = linear_projection(D, X, Z), tau,
                           h, subsample_size,
                           beta_D_proposal = NULL, beta_X_proposal = NULL, 
                           gamma = 1, l = 2) {
  n <- length(Y)
  p_design <- length(h)
  subsample_set <- h # store indices of subsample
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
  alt_subsample_weights <- vector("double", subsample_size - length(h))
  # we have length(h) observations in subsample; we need subsample_size - length(h) more
  for (j in seq_len(subsample_size - length(h))) {
    choices <- setdiff(seq_len(nrow(Y)), subsample_set)
    n_choices <- length(choices)
    s_remaining <- s[choices, , drop = FALSE]

    # each row is one observation that we can choose from; the value is the weight in the exponent
    ones <- matrix(1, nrow = n_choices, ncol = 1)
    sum_remaining <- kronecker(ones, sum_across_subsample_set) + s_remaining

    # for each row, apply e^(-gamma * (l-norm^l))
    weights <- apply(sum_remaining, 1, function(x){
      exp(-gamma * sum(abs(x)^l))
    })

    alt_weights <- apply(sum_remaining, 1, function(x){
      # - tau < x < 1 - tau => - tau - x < 0 & x - (1 - tau) < 0
      # => max{...} = 0 => we satisfy FOC
      # Put differently:
      # x is a 1 by p vector;
      # So, x - (t - tau) and (- tau - x) are both 1 by p vectors.
      # If each entry of these 1 by p vectors are negative, we satisfy FOC conditions
      # Fix: `max` must be used element-wise!!! # Q: ask GP if this is right
      left <- - tau - x
      right <- x - (1 - tau)
      max_result <- vector("double", length(x))
      for (entry in seq_along(left)) {
        left_entry <- left[[entry]]
        right_entry <- right[[entry]]
        max_result[[entry]] <- max(left_entry, right_entry, 0)
      }
      exp(-gamma * sum(abs(max_result))^l)
    })

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
        problem_weights = weights # return problematic weights
      )
    })
    if (winner$status == "ERROR") {
      return(winner)
    }
    winner <- winner$answer
    subsample_weights[[j]] <- weights[winner]
    alt_subsample_weights[[j]] <- alt_weights[winner]
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
    alt_subsample_weights = alt_subsample_weights, # return alternative weights
    subsample_weights = subsample_weights, # return weights for each observation in subsample
    xi = sum_across_subsample_set, # return xi object # TODO: double-check that this is indeed the xi object
    s = s # return each entry of term of xi for all observations, not just those in the subsample
  )
}
