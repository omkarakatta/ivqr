% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/warmstart.R
\name{compute_warmstart}
\alias{compute_warmstart}
\title{Compute warm-starting solutions}
\usage{
compute_warmstart(
  Y,
  X,
  D,
  Z,
  Phi = linear_projection(D, X, Z),
  tau,
  method = NULL
)
}
\arguments{
\item{Y}{Dependent variable (vector of length n)}

\item{X}{Exogenous variable (including constant vector) (n by p_X matrix)}

\item{D}{Endogenous variable (n by p_D matrix)}

\item{Z}{Instrumental variable (n by p_Z matrix)}

\item{Phi}{Transformation of X and Z to be used in the program;
defaults to the linear projection of D on X and Z (n by p_D matrix)}

\item{tau}{Quantile (numeric between 0 and 1)}

\item{method}{Determines how we compute the warm-starting solution;
defaults to NULL, see 'Details'}
}
\value{
Either NULL or a vector of values to warm-start the program; the
output is meant to be fed into the \code{start} parameter of \code{iqr_milp}
}
\description{
The result of this function should be fed into the \code{start} argument of
\code{iqr_milp}.
}
\details{
If \code{method} is NULL, we won't warm-start the program at all, i.e., we
return NULL.

If \code{method} contains the phrase "naive3", we construct the
warm-starting solution with two quantile regressions.
First, we run a preliminary quantile regression of Y on D and X to obtain
\code{beta_D}. Then, we concentrate out D from Y according to \code{beta_D}.
After regressing the concentrated-out Y on Phi and X, we obtain
\code{beta_Phi} and \code{beta_X}. Using this information, we can construct
our decision vector.

If \code{method} contains the phrase "full", then we not only specify the
coefficients in the decision vector, we also specify the dual variables,
residuals, and the binary variables.
}
