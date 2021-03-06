% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcmc-h.R
\name{mcmc_h}
\alias{mcmc_h}
\title{MH Sampler of the Proposal Distribution of Coefficients/Active Bases}
\usage{
mcmc_h(
  iterations,
  h_opt,
  beta_D_opt,
  beta_X_opt,
  residuals_opt,
  theta,
  varcov_mat,
  Y,
  X,
  D,
  Phi,
  discard_burnin,
  manual_burnin = 1,
  unique_beta_quota = rep(0, ncol(D)),
  largest_min = rep(Inf, ncol(D)),
  smallest_max = rep(-Inf, ncol(D)),
  max_iterations = Inf,
  label_function = function(idx, total = iterations) {
     print(paste("MCMC-H IDX:",
    idx, "/", total))
 },
  label_skip = floor(iterations/5),
  label_bool = TRUE,
  always_accept = FALSE
)
}
\arguments{
\item{iterations}{Number of iterations}

\item{residuals_opt}{Residuals from IQR-QR problem}

\item{theta}{Tuning parameter}

\item{varcov_mat}{Variance-covariance matrix}

\item{Y, X, D, Phi}{Data}

\item{discard_burnin}{If TRUE, remove first few samples that have the same
coefficients}

\item{unique_beta_quota}{Vector of length p_D; minimum number of unique beta
(defaults to 0)}

\item{largest_min}{Vector of length p_D; sample from MCMC must contain
beta_D below corresponding value of \code{largest_min}; defaults to Inf}

\item{smallest_max}{Vector of length p_D; sample from MCMC must contain
beta_D above corresponding value of \code{smallest_max}; defaults to -Inf}

\item{max_iterations}{Total number of iterations before stopping the
program; defaults to Inf, but set to a lower number when using \code{largest_min}
or \code{smallest_max}}

\item{always_accept}{If FALSE (default), sample from \code{proposal_h}. Else,
sample from \code{target_h}.}

\item{initial_h}{Initial active basis indices}

\item{initial_beta_D, initial_beta_X}{Initial coefficients, e.g., IQR point
estimate}
}
\value{
named list
\enumerate{
\item \code{beta}: data frame where each row is a vector of coefficients (one
row per iteration)
\item \code{h}: data frame where each row is a vector of indices in the active
basis (one row per iteration)
\item \code{record}: binary vector, 1 if that iteration's proposal was accepted
\item \code{stationary_begin}: all iterations prior to this number were discarded
}
}
\description{
The proposal distribution Q_beta is the asymptotic distribution of the IQR
estimator limited to the finite support of possible solutions enumerated by
the active basis.
This MH sampler returns draws from the proposal distribution by
}
\details{
\enumerate{
\item Proposing an active basis \code{h_star} according to Q_h, which puts more
weight on indices with smaller residuals
\item Computing coefficients \code{beta_star} (see \code{h_to_beta})
\item Computing the acceptance probability
\item Accept/Reject in the usual MH manner
}
}
