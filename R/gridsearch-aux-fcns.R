### Meta -------------------------
###
### Title: Auxiliary Functions for Grid Search Procedure
###
### Description: Grid search can be decomposed into a series of steps:
### 1. Create a grid
### 2. Evaluate objective on a grid
### 3. Find coordinate(s) with minimum objective
### 4. Optionally create a new grid centered on the coordinate(s) from the
###   previous step and repeat steps 2-4.
###
### Author: Omkar A. Katta
###

### get_initial_beta_D -------------------------
#' Obtain naive coefficient estimates for \code{beta_D}
#'
#' Naively regress \code{Y} on \code{D}, \code{Z}, and \code{X}, and return
#' the coefficients on \code{D}.
#'
#' @param Y Outcome matrix (more precisely, a column vector)
#' @param D Matrix of endogenous variables
#' @param Z Matrix of instruments
#' @param X Matrix of covariates (incl. column of 1's for the intercept)
#' @param ... Arguments to be passed to \code{quantreg::rq()}
#'
#' @return A vector of coefficients on the endogenous variable
get_initial_beta_D <- function(Y, D, Z, X, ...) {
  qr <- quantreg::rq(Y ~ D + Z + X - 1, ...)
  stats::coef(qr)[seq_len(ncol(D))]
}

### center_out -------------------------
#' Create a univariate grid of values from the center out
#'
#' Given some center, extend in the positive and negative directions to create
#' a grid of values
#'
#' @param center Output of \code{get_initial_beta_D} or a number to center
#'  the grid (numeric)
#' @param increment Granularity of the grid (numeric)
#' @param length Number of grid values to the right and left of the grid
#'  (numeric)
#'
#' @return A vector of grid values
center_out <- function(center, increment, length) {
  pos <- seq(from = center, by = increment, length.out = length + 1)
  neg <- seq(to = center, by = increment, length.out = length + 1)
  grid <- unique(c(neg, pos))
  grid
}
