% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mat-fcns.R
\name{block_diagonal}
\alias{block_diagonal}
\title{Create a block diagonal matrix}
\usage{
block_diagonal(mat_list, fill = 0)
}
\arguments{
\item{mat_list}{The list of matrices in the order they appear from top-right to
bottom-left in the final matrix}

\item{fill}{The value of entries in the off-block entries; defaults to 0}
}
\value{
A block-diagonal matrix
}
\description{
Combine a list of matrices into a block diagonal matrix
}
\details{
The matrices provided in \code{mat_list} will be arranged as block matrices
with the first listed matrix at the top-right and the last matrix at the
bottom left.
The entries that are not in each block are filled with \code{fill}, which
defaults to 0.
}
\examples{
\dontrun{
block_diagonal(diag(1, 2), diag(1, 3)) #~ this is identical to `diag(1, 5)`
}
}
\seealso{
\code{\link{expand_matrix}}

Other matrix operations: 
\code{\link{expand_matrix}()}
}
\concept{matrix operations}
