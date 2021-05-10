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
#' @param Y Dependent variable (vector of length n)
#' @param X Exogenous variable (including constant vector) (n by p_X matrix)
#' @param D Endogenous variable (n by p_D matrix)
#' @param Z Instrumental variable (n by p_Z matrix)
#' @param tau Quantile of interest (numeric between 0 and 1)
#' @param ... Arguments to be passed to \code{quantreg::rq()}
#'
#' @return A vector of coefficients on the endogenous variable
get_initial_beta_D <- function(Y, X, D, Z, tau, ...) {
  msg <- "`tau` is meant to be a single numeric."
  send_note_if(msg, length(tau) > 1, warning)
  qr <- quantreg::rq(Y ~ D + Z + X - 1, tau = tau, ...)
  stats::coef(qr)[seq_len(ncol(D))]
}

### center_out_uni -------------------------
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
center_out_uni <- function(center, increment, length) {
  pos <- seq(from = center, by = increment, length.out = length + 1)
  neg <- seq(to = center, by = increment, length.out = length + 1)
  grid <- unique(c(neg, pos))
  grid
}

### center_out_grid -------------------------
#' Create multidimensional grid of values
#'
#' Apply \code{center_out_uni} element-wise to each argument of this function
#' and take all possible combinations of coordinate grid values to create a
#' data frame of grid coordinates
#'
#' @param center Output of \code{get_initial_beta_D}; values that act as
#'  the center of each axis of the grid (numeric vector of length p_D)
#' @param increment Granularity of the grid in each axis
#'  (numeric vector of length p_D)
#' @param length Number of grid values to the right and left of each axis
#'  (numeric vector of length p_D)
#'
#' @return A data frame with each row acting as the coordinates of the grid.
#'  If \code{center} is a named vector, the columns of the data frame will
#'  share the names of the vector. If \code{center} is not a named vector,
#'  the columns of the data frame are Var1, Var2, etc.
center_out_grid <- function(center, increment, length) {
  stopifnot(length(center) == length(increment))
  stopifnot(length(center) == length(length))
  grid_list <- vector("list", length(center))
  for (i in seq_along(center)) {
    vec <- center_out_uni(center[[i]], increment[[i]], length[[i]])
    grid_list[[i]] <- vec
  }
  grid <- expand.grid(grid_list)
  if (!is.null(names(center))) {
    colnames(grid) <- names(center)
  }
  grid
}

### get_iqr_objective -------------------------
#' Get value of IQR objective given coefficients on endogeneous variable
#'
#' Obtain the sum of the absolute value of the coefficients on the instruments
#' in a quantile regression of \code{Y} after concentrating out \code{D}
#' according to \code{beta_D} on \code{Z} and \code{X}.
#'
#' @param beta_D Vector of coefficients on the endogenous variable
#'  (numeric vector)
#' @param Y Dependent variable (vector of length n)
#' @param X Exogenous variable (including constant vector) (n by p_X matrix)
#' @param D Endogenous variable (n by p_D matrix)
#' @param Z Instrumental variable (n by p_Z matrix)
#' @param tau Quantile of interest (numeric between 0 and 1)
#' @param ... Arguments to be passed to \code{quantreg::rq()}
#'
#' @return A named list of two entries:
#'  \enumerate{
#'    \item \code{beta_Z}: named vector of coefficients on the instruments
#'    \item \code{tau}: quantile of interest
#'    \item \code{obj}: sum of absolute value of \code{beta_Z}, i.e., value of
#'      IQR objective
#'  }
get_iqr_objective <- function(beta_D, Y, X, D, Z, tau, ...) {
  msg <- "`tau` is meant to be a single numeric."
  send_note_if(msg, length(tau) > 1, stop, call. = FALSE)
  qr <- quantreg::rq(Y - D %*% beta_D ~ Z + X - 1, tau = tau, ...)

  # `as.data.frame(coef(qr))` returns a data frame with each row corresponding
  # to a coefficient and each column corresponding to a quantile.
  beta_Z <- as.data.frame(coef(qr))[seq_len(ncol(Z)), ]
  names(beta_Z) <- colnames(Z)
  list(beta_Z = beta_Z, tau = tau, obj = sum(abs(beta_Z)))
}

### get_iqr_objective_grid -------------------------
#' Compute IQR objective given grid of coefficients on endogeneous variables
#'
#' For each set of \code{beta_D} suggested by \code{grid}, compute the
#' sum of the absolute values of \code{beta_Z}
#'
#'
#'
get_iqr_objective_grid <- function(grid,
                                   Y,
                                   X,
                                   D,
                                   Z,
                                   tau,
                                   update = round(nrow(grid) / 20),
                                   ...) {
  msg <- "`update` is meant to be a single number less than `nrow(grid)`."
  send_note_if(msg, update > nrow(grid), stop, call. = FALSE)
  beta_Z_coef <- vector("list", length = nrow(grid))
  objective <- vector("double", length = nrow(grid))
  for (i in seq_len(nrow(grid))) {
    if (i %% update == 0) {
      msg <- paste("Grid:", i, "out of", nrow(grid))
      print(msg)
    }
    beta_D_vec <- as.numeric(grid[i, ])
    result <- get_iqr_objective(beta_D_vec, Y, X, D, Z, tau, ...)
    beta_Z_coef[[i]] <- result$beta_Z
    objective[[i]] <- result$obj
  }
  beta_Z_coef <- do.call(rbind, beta_Z_coef)
  cbind(grid, beta_Z_coef, objective)
}
