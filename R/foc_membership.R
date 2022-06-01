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
                     Phi = Phi_subsample)$beta_D,
  beta_X = h_to_beta(h,
                     Y = Y_subsample,
                     X = X_subsample,
                     D = D_subsample,
                     Phi = Phi_subsample)$beta_X
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
  resid_subsample <- y_tilde_subsample - X_subsample %*% beta_X

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

### foc_membership_v2 -------------------------

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

### foc_membership_v3 -------------------------

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
