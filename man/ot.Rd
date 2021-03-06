% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aux-fcns.R
\name{ot}
\alias{ot}
\title{Compute OT between identical uniform distributions}
\usage{
ot(pre, post, params = list(OutputFlag = 0), method = "gurobi")
}
\description{
Compute transport map between U(1,...,n-p) and U(1,...,n-p) where C(i, j) =
norm difference of pre\link{, i} and post\link{, j} where i and j are between 1 and
n-p.
}
\details{
We use Gurobi to solve the OT problem if \code{method} is "gurobi".
We use the transport package to solve the OT problem if \code{method} is "transport".
}
