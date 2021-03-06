% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aux-fcns.R
\name{parse_mult_sols}
\alias{parse_mult_sols}
\title{Parse multiple solution (.sol) files}
\usage{
parse_mult_sols(sol_dir, expr = "sol$", value = "list")
}
\arguments{
\item{sol_dir}{Path to directory with .sol files}

\item{expr}{Expression to limit files in \code{sol_dir}; defaults to
searching for "sol" extension}

\item{value}{How should the information be returned? defaults to "list",
accepts "data.frame" if the number of decision veariables are the same
across solution files}
}
\description{
Obtain solution from multiple .sol files.
}
