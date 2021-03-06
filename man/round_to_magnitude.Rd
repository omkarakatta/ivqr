% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aux-fcns.R
\name{round_to_magnitude}
\alias{round_to_magnitude}
\title{Round up the nearest order of magnitude}
\usage{
round_to_magnitude(x)
}
\arguments{
\item{x}{Number to be rounded (numeric)}
}
\value{
Nearest order of magnitude larger than \code{x}
}
\description{
For example, if \code{x} is 12, this function returns 100.
If \code{x} is 0.12, this function returns 0.1.
}
