### Meta -------------------------
###
### Title: Grid Search Procedure
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
  # TODO: Use colnames(D) for name of coef vector if colnames(D) is not NULL
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
  beta_Z <- as.data.frame(stats::coef(qr))[seq_len(ncol(Z)), ]
  if (!is.null(colnames(Z))) {
    names(beta_Z) <- colnames(Z)
  }
  list(beta_Z = beta_Z, tau = tau, obj = sum(abs(beta_Z)))
}

### gridsearch -------------------------
#' Compute IQR objective given grid of coefficients on endogeneous variables
#'
#' For each set of \code{beta_D} suggested by \code{grid}, compute the
#' sum of the absolute values of \code{beta_Z}
#'
#' This code is not run in parallel.
#'
#' @param grid Data frame with p_D columns where each row is a set of
#'  \code{beta_D} coefficients; Output of \code{center_out_grid}
#' @param Y Dependent variable (vector of length n)
#' @param X Exogenous variable (including constant vector) (n by p_X matrix)
#' @param D Endogenous variable (n by p_D matrix)
#' @param Z Instrumental variable (n by p_Z matrix)
#' @param tau Quantile of interest (numeric between 0 and 1)
#' @param update Send progress after every \code{update} values of grid search;
#'  defaults ot \code{round(nrow(grid)) / 20}, i.e., 20 updates (numeric)
#' @param ... Arguments to be passed to \code{quantreg::rq()}
#'
#' @return A data frame of dimension \code{nrow(grid)} by p_D + p_Z + 1 where
#'  each row corresponds to one set of coordinate values on the grid, the
#'  corresponding values for the instrument coefficient, and the resulting
#'  IQR objective
gridsearch <- function(grid,
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
  result_with_min_obj <- cbind(iteration = seq_len(nrow(grid)),
                               grid, beta_Z_coef, objective) %>%
    dplyr::mutate(min_obj = objective == min(objective))

  # print argmin of grid search
  print(result_with_min_obj[result_with_min_obj$min_obj, ])

  # return results of grid evaluations
  return(invisible(result_with_min_obj))
}

### gridsearch_parallel -------------------------
#' Compute IQR objective given grid of coefficients on endogeneous variables
#'
#' For each set of \code{beta_D} suggested by \code{grid}, compute the
#' sum of the absolute values of \code{beta_Z}
#'
#' This code is run in parallel.
#' To store results as they are found, specify a directory with \code{log_dir}
#' to store intermediate CSV files. Then, use \code{concatenate_csvs} to
#' merge CSVs into one file, i.e.,
#' \code{cols = c("iteration", colnames(grid), colnames(Z), "objective")}
#' \code{concatenate_csv(dir = log_dir, cols = cols, header = FALSE)}
#'
#' Note: if the function terminates without any problems, the intermediary
#' CSV files will be deleted, and the results are saved as a CSV in
#' \code{log_dir} as "gridsearch_results.csv".
#'
#' @import foreach
#'
#' @param grid Data frame with p_D columns where each row is a set of
#'  \code{beta_D} coefficients; Output of \code{center_out_grid}
#' @param Y Dependent variable (vector of length n)
#' @param X Exogenous variable (including constant vector) (n by p_X matrix)
#' @param D Endogenous variable (n by p_D matrix)
#' @param Z Instrumental variable (n by p_Z matrix)
#' @param tau Quantile of interest (numeric between 0 and 1)
#' @param log_dir Path of log directory to store parallelized results;
#'  if NULL (default), log is not saved (string)
#' @param cores Number of cores to be used in parallelization process
#' @param ... Arguments to be passed to \code{quantreg::rq()}
#'
#' @return A data frame of dimension \code{nrow(grid)} by p_D + p_Z + 3 where
#'  each row corresponds to one set of coordinate values on the grid, the
#'  corresponding values for the instrument coefficient, and the resulting
#'  IQR objective as well as the row number of the \code{grid} and whether
#'  the objective is the smallest within the grid search
gridsearch_parallel <- function(grid,
                                Y,
                                X,
                                D,
                                Z,
                                tau,
                                log_dir = NULL,
                                cores = parallel::detectCores()[1] - 2,
                                ...) {
  create_log <- FALSE
  if (!is.null(log_dir)) {
    if (dir.exists(log_dir) | file.exists(log_dir)) {
      stop(paste(log_dir, "already exists. Choose a different `log_dir`."))
    } else { # don't overwrite file
      create_log <- TRUE
      dir.create(log_dir)
    }
  }

  # set up cluster
  cl <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl))

  # name of results
  cols <- c("iteration", colnames(grid), colnames(Z), "objective")

  # find IQR objective for each grid coordinate
  i <- NULL # avoid undefined global variable note in R CMD check results
  result <- foreach (i = seq_len(nrow(grid)),
           .combine = rbind,
           .export = c("get_iqr_objective")) %dopar% {
    beta_D_vec <- as.numeric(grid[i, ])
    names(beta_D_vec) <- colnames(grid)
    tmp <- get_iqr_objective(beta_D_vec, Y, X, D, Z, tau, ...)
    result <- c(i, beta_D_vec, tmp$beta_Z, objective = tmp$obj)
    names(result) <- cols

    # store results in log
    if (create_log) {
      file_path <- paste0(log_dir, "/", "iteration", i, ".csv")
      file.create(file_path)
      cat(result, sep = ",", file = file_path)
      cat("\n", sep = ",", file = file_path, append = TRUE) # add newline at EOF
    }
    result
  }

  # Find smallest IQR objective
  result_with_min_obj <- result %>%
    as.data.frame() %>%
    dplyr::mutate(min_obj = objective == min(objective))

  # save result in CSV
  if (create_log) {
    unlink(log_dir, recursive = T)
    dir.create(log_dir)
    utils::write.csv(result_with_min_obj,
                     paste0(log_dir, "/", "gridsearch_results.csv"))
  }

  # print argmin of grid search
  print(result_with_min_obj[result_with_min_obj$min_obj, ])

  # return results of grid evaluations
  return(invisible(result_with_min_obj))
}
