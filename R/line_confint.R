### Meta -------------------------
###
### Title: Construct a univariate confidence interval via a line search
###
### Description: For different values of the null \beta_j (where \beta_j could
### refer to \beta_X or \beta_D), determine whether the associated test
### statistic rejects or fails to reject the null. The values that fail to
### reject belong in the confidence interval. We search for the smallest and
### largest value of the confidence interval.
###
### Author: Omkar A. Katta
###

### line_confint -------------------------
#' Construct a univariate confidence interval via a line search
#'
#' Search for the smallest/largest value of a particular coefficient where the
#' null hypothesis of the coefficient being equal to that value is not
#' rejected.
#'
#' The suggested value of \code{beta_null} is given by the output of
#' \code{iqr_milp} or \code{preprocess_iqr_milp}.  For example, if \code{index}
#' is 2, \code{endogeneous} is TRUE, and the \code{iqr_milp} estimates are
#' stored in \code{iqr}, then \code{iqr$beta_D[2]} should be the value of
#' \code{beta_null}.
#'
#'
#' @param index Index of the coefficient of interest
#'  (numeric between 1 and p_D if \code{endogeneous} is TRUE;
#'  numeric between 1 and p_X if \code{endogeneous} is FALSE)
#' @param endogeneous If TRUE (default), \code{index} refers to the
#'  coefficients on the endogeneous variables; if FALSE, \code{index} refers to
#'  the coefficients on the exogeneous variables (boolean)
#' @param beta_null Initial value of line search; suggested value is given by
#'  \code{iqr_milp}; if NULL (default), naive quantile regression point
#'  estimate will be used (numeric)
#' @param stopping_tolerance Smallest step size before stopping the line
#'  search; if NULL (default), takes on the value of \code{width_ratio} times
#'  the width of the naive confidence interval (numeric)
#' @param width_ratio Determines the relationship between the size of naive
#'  confidence interval and default value of \code{stopping_tolerance};
#' @param step_size Initial step size of line search; if NULL(default), use
#'  half the width of the naive confidence interval (numeric)
#' @param step_rate Rate at which we decrease the step size when we cross the
#'  border of the confidence interval; defaults to 0.5 (numeric, less than 1)
#'  defaults to 1 (numeric)
#' @param cores Number of cores to be used in parallelization process; defaults
#'  to 2 (no more than 2 needed)
#' @param log_dir Path of log directory to store parallelized results;
#'  if NULL (default), log is not saved (string)
#' @param log_name Name of CSV file (including extension) to store results;
#'  defaults to "gridsearch_results.csv" (string)
#' @param remove_intermediate_csvs If TRUE, intermediate csvs created during
#'  while loop will be removed after line search is finished; defaults to FALSE
#' @param return_setup If TRUE, return the step size and naive quantile
#' regression summary without carrying out the line search; defaults to FALSE
#' @param p_val_tol,small_change_count_tol If the change in p-value is less
#'  than \code{p_val_tol} for \code{small_change_count_tol} consecutive
#'  iterations, we stop the while loop; defaults to 1e-6 and 10
#' @inheritParams test_stat
#'
#' @import foreach
#'
#' @export
line_confint <- function(index,
                         endogeneous = TRUE,
                         beta_null = NULL,
                         stopping_tolerance = NULL,
                         width_ratio = 1,
                         step_size = NULL,
                         step_rate = 0.5,
                         cores = 2,
                         log_dir = NULL,
                         log_name = "line_search.csv",
                         remove_intermediate_csvs = FALSE,
                         return_setup = FALSE,
                         p_val_tol = 1e-6,
                         small_change_count_tol = 10,
                         alpha = 0.1,
                         Y,
                         X,
                         D,
                         Z,
                         Phi = linear_projection(D, X, Z),
                         tau,
                         B = NULL,
                         orthogonalize_statistic = FALSE,
                         homoskedasticity = FALSE,
                         kernel = "Powell",
                         show_progress = TRUE,
                         FUN = preprocess_iqr_milp,
                         ...) {

  # Start clock
  clock_start <- Sys.time()
  message(paste("Clock started:", clock_start))

  out <- list() # Initialize list of results to return

  # Get dimensions of data
  n <- length(Y)
  n_D <- nrow(D)
  n_X <- nrow(X)
  n_Z <- nrow(Z)
  n_Phi <- nrow(Phi)
  p_D <- ncol(D)
  p_X <- ncol(X)

  # Ensure that number of observations is the same
  stopifnot(all.equal(n, n_D))
  stopifnot(all.equal(n, n_X))
  stopifnot(all.equal(n, n_Z))
  stopifnot(all.equal(n, n_Phi))

  # Ensure that there are some RHS variables
  stopifnot(p_D + p_X > 0)

  # Send warning if more than 2 cores are being used
  msg <- "No more than 2 cores are needed for line search."
  send_note_if(msg, cores > 2, warning)

  # Check that index is an integer that is between 1 and p_D or p_X (inclusive)
  stopifnot(index > 0) # can't be negative
  stopifnot(length(index) == 1) # we project on axis of only one coefficient
  stopifnot(round(index) == index) # must be integer
  stopifnot(index <= ifelse(endogeneous, p_D, p_X))
  out$index <- index
  out$endogeneous <- endogeneous

  # set up log
  create_log <- FALSE
  date_time <- format(Sys.time(), "%y%m%d_%H%M%S")
  out$date_time <- date_time
  if (!is.null(log_dir)) {
    # create path of log file
    # note that date and time will be prepended: yymmdd_hhmmss
    stub <- paste0("index", index, "_", ifelse(endogeneous, "endo", "exo"),
                   "_tau", tau,
                   "_kernel", ifelse(homoskedasticity, "homoskedasticity", kernel))
    log_path <- paste0(log_dir, "/", date_time, "_", stub, "_", log_name)
    if (file.exists(log_path)) {
      msg <- paste(log_path,
                   "already exists. Choose a different `log_name` or `log_dir`")
      stop(msg)
    } else {
      create_log <- TRUE
      # create directory if nonexistent
      dir.create(log_dir, showWarnings = FALSE)
    }
    message(paste("Log created"))
  }

  # Determine stopping tolerance via naive QR
  if (is.null(stopping_tolerance) | is.null(beta_null) | is.null(step_size)) {
    if (endogeneous) {
      qr <- quantreg::rq(Y ~ D + X - 1, tau = tau)
    } else {
      qr <- quantreg::rq(Y ~ X + D - 1, tau = tau)
    }
    qr_summary <- summary(qr)
    out$qr_summary <- qr_summary
    if (ncol(qr_summary$coef) == 4) {
      center <- c(qr_summary$coef[index, "Value"])
      se <- c(qr_summary$coef[index, "Std. Error"])
      crit_val <- qnorm(1 - alpha / 2)
      bounds <- c(center - crit_val * se, center + crit_val * se)
    } else if (ncol(qr_summary$coef) == 3) {
      bounds <- qr_summary$coef[index, c("lower bd", "upper bd"), drop = FALSE]
    }
    width <- width_ratio * (bounds[2] - bounds[1])
    if (is.null(stopping_tolerance)) {
      stopping_tolerance <- round_to_magnitude(width * 0.01)
    }
    if (is.null(beta_null)) {
      beta_null <- qr_summary$coef[index, 1]
    }
    if (is.null(step_size)) {
      step_size <- width / 2
    }
  }
  # return(out) # DEBUG
  stopifnot(is.numeric(stopping_tolerance))
  stopifnot(is.numeric(beta_null))
  stopifnot(is.numeric(step_size))
  stopifnot(step_size > 0) # TODO: If this happens, print the bounds; If this happens, should I stop or set the step size with another heuristic?
  out$stopping_tolerance <- stopping_tolerance
  out$beta_null <- beta_null
  out$initial_step_size <- step_size

  send_note_if(paste0("left bound: ", bounds[1]), show_progress, message)
  send_note_if(paste0("right bound: ", bounds[2]), show_progress, message)
  send_note_if(paste0("width: ", width), show_progress, message)
  send_note_if(paste0("stopping tolerance: ", stopping_tolerance),
               show_progress, message)
  send_note_if(paste0("step size: ", step_size), show_progress, message)

  if (return_setup) {
    return(out)
  }

  # Construct null hypothesis
  beta_D_null <- rep(NA, p_D)
  beta_X_null <- rep(NA, p_X)
  if (endogeneous) {
    beta_D_null[index] <- beta_null
  } else {
    beta_X_null[index] <- beta_null
  }

  # Is initial value in confidence interval?
  initial_test_stat <- test_stat(
    beta_D_null = beta_D_null,
    beta_X_null = beta_X_null,
    alpha = alpha,
    Y = Y,
    X = X,
    D = D,
    Z = Z,
    Phi = Phi,
    tau = tau,
    B = B,
    orthogonalize_statistic = orthogonalize_statistic,
    homoskedasticity = homoskedasticity,
    kernel = kernel,
    show_progress = show_progress,
    # FUN = FUN, # TODO: allow user to choose between preprocess_iqr_milp and iqr_milp
    ...
  )
  if (!is.null(initial_test_stat$ended_early)) {
    return(out)
  }
  out$initial_test_stat <- initial_test_stat
  initial_rejected <- initial_test_stat$p_val < alpha
  if (initial_rejected) {
    warning("Initial value of beta_null was rejected.")
    return(out)
    # TODO: do a line search to find a value in the confidence interval
  }

  msg <- "Computed initial test statistic"
  send_note_if(msg, show_progress, message)

  out$homoskedasticity <- homoskedasticity
  kernel <- ifelse(homoskedasticity, "homoskedasticity", kernel)
  out$kernel <- kernel
  out$alpha <- alpha
  out$cores <- cores

  # Line Search
  # set up cluster
  cl <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl))
  confint <- foreach (direction = c(-1, 1),
                      .combine = c,
                      .export = c("test_stat",
                                  "preprocess_iqr_milp",
                                  "iqr_milp")) %dopar% {
    type <- ifelse(direction == 1, "max", "min")
    step <- step_size
    current_beta <- beta_null
    beta_D_null <- rep(NA, p_D)
    beta_X_null <- rep(NA, p_X)
    current_ts_reject <- FALSE # if not FALSE, then we wouldn't do line search
    current_p_val <- initial_test_stat$p_val
    counter <- 0
    small_change_in_p_val <- 0
    while (step > stopping_tolerance &
           small_change_in_p_val < small_change_count_tol) {
      # TODO: stop searching if the change in p-value is very small
      clock_start_while <- Sys.time()
      counter <- counter + 1
      old_p_val <- current_p_val # save previous p-value
      old_ts_reject <- current_ts_reject # save previous status of test
      old_beta <- current_beta # save previous beta
      current_beta <- current_beta + step * direction # update beta
      # construct null
      if (endogeneous) {
        beta_D_null[index] <- current_beta
      } else {
        beta_X_null[index] <- current_beta
      }
      # find test stat
      ts <- test_stat(
        beta_D_null = beta_D_null,
        beta_X_null = beta_X_null,
        alpha = alpha,
        Y = Y,
        X = X,
        D = D,
        Z = Z,
        Phi = Phi,
        tau = tau,
        B = B,
        orthogonalize_statistic = orthogonalize_statistic,
        homoskedasticity = homoskedasticity,
        kernel = kernel,
        show_progress = show_progress,
        # FUN = preprocess_iqr_milp, # TODO: let user decide what this is
        ...
      )
      # determine test status
      current_p_val <- ts$p_val
      current_ts_reject <- ts$p_val < alpha
      # check if change in p-value is small
      if (abs(current_p_val - old_p_val) < p_val_tol) {
        small_change_in_p_val <- small_change_in_p_val + 1
      } else {
        small_change_in_p_val <- 0
      }
      # update direction and step size
      if (old_ts_reject != current_ts_reject) {
        direction <- -1 * direction
        step <- step * step_rate
      }

      clock_end_while <- Sys.time()
      elapsed_while <- difftime(clock_end_while,
                                clock_start_while,
                                units = "mins")

      # save results
      results <- c(type,                # min or max
                   counter,             # iteration through while loop
                   index,               # coefficient of interest
                   homoskedasticity,    # homoskedasticity
                   kernel,              # kernel (only relevant if homoskedasticity is FALSE)
                   endogeneous,         # coefficient of interest
                   current_beta,        # current value of null
                   current_p_val,       # current p-value
                   ts$test_stat,        # current test-stat
                   current_ts_reject,   # p-val < alpha => reject null?
                   elapsed_while,       # time for current iteration
                   step,                # new step
                   direction,           # new direction
                   stopping_tolerance)  # stopping tolerance
      if (create_log) {
        file_path <- paste0(log_dir, "/",
                            date_time, "_", stub, "_",
                            type, "_counter_", counter, ".csv")
        file.create(file_path)
        cat(results, sep = ",", file = file_path)
        cat("\n", sep = ",",
            file = file_path,
            append = TRUE) # add newline at EOF
      }
    }
    out$counter <- counter # TODO: this is not being saved; save min_counter and max_counter

    if (small_change_in_p_val < small_change_count_tol) {
      # Linearly interpolate between previous two values of beta
      # One of these values of beta will be inside the confidence interval, while
      # the other value is outside the confidence interval.
      pair_p_val <- c(old_p_val, current_p_val)
      ordered <- order(pair_p_val) # ordered[1] = index of smaller p-value
      pi <- (alpha - pair_p_val[ordered[1]]) /
        (pair_p_val[ordered[2]] - pair_p_val[ordered[1]])
      pair_beta <- c(old_beta, current_beta)
      beta_border <- (1 - pi) * pair_beta[ordered[1]] + pi * pair_beta[ordered[2]]
      beta_border # this is what we return at end of foreach loop
    } else {
      paste("p-val flattens for", type)
    }
  }
  flattens_min <- grepl("p flattens for min", confint)
  flattens_max <- grepl("p flattens for max", confint)
  flattens <- flattens_min + flattens_max
  if (sum(flattens) == 2) {
    # if both min and max directions have plateauing p-values,
    # then return the following value of confint
    confint <- c("p-val flattens for min", "p-val flattens for max")
    warning("p-val flattens for min and max")
  } else if (sum(flattens_min) == 1) {
    # if only the min direction has plateauing p-value,
    # then return the message first then the upper bound
    confint <- c(confint[flattens_min], confint[!flattens_min])
    warning("p-val flattens for min")
  } else if (sum(flattens_max) == 1) {
    # if only the min direction has plateauing p-value,
    # then return the lower bound first then the message
    confint <- c(confint[!flattens_max], confint[flattens_max])
    warning("p-val flattens for max")
  } else if (sum(flattens) == 0) {
    # if none of the directions has plateauing p-value,
    # then return the lower bound first then the upper bound
    confint <- sort(confint)
  } else {
    warning("unexpected output of foreach loop")
    out$confint <- confint
    return(out)
  }
  out$lower <- confint[1] # lower bound
  out$upper <- confint[2] # upper bound

  # Stop the clock
  clock_end <- Sys.time()
  elapsed_time <- difftime(clock_end, clock_start, units = "mins")
  out$time_elapsed <- elapsed_time
  message(paste("Clock stopped:", clock_end))

  # save results
  results <- data.frame(index = index,
                        endogeneous = endogeneous,
                        homoskedasticity = homoskedasticity,
                        kernel = kernel,
                        alpha = alpha,
                        lower = confint[1],
                        upper = confint[2],
                        cores = cores,
                        minutes = elapsed_time)

  if (create_log) {
    # save results in CSV
    utils::write.csv(results, log_path)
    # TODO: add option to concatenate these files instead of removing them
    if (remove_intermediate_csvs) {
      # remove files save during while loop
      # TODO: maybe let the pattern use date_time to not disturb other files
      file.remove(list.files(log_dir, pattern = "counter", full.names = T))
    }
  }

  return(out)
}
