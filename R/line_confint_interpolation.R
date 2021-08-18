### Meta -------------------------
###
### Title: Construct a univariate confidence interval via a line search
###
### Description: For different values of the null \beta_j (where \beta_j could
### refer to \beta_X or \beta_D), determine whether the associated test
### statistic rejects or fails to reject the null. The values that fail to
### reject belong in the confidence interval. We search for the smallest and
### largest value of the confidence interval.
### This code is based on ./R/line_confint.R, but this also includes an interpolation step:
### We compute the short-iqr for extreme values of beta_D, one to the left and
### another to the right of the initial milp estimate. Then, we (linearly)
### interpolate the values of the test-stat using these three values (left,
### right, and initial aka milp_estimate). Then, using the interpolation, we
### come up with our best guess for the value of the beta_D that produces the
### test-stat closest to the critical value. We then intelligently choose the
### step size and perform the line search as usual.
###
### Author: Omkar A. Katta
###

### line_confint_interpolation -------------------------
#' Construct a univariate confidence interval via a line search and interpolation
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
#' The time limit for computing the test-statistic at each step of the line
#' search code is twice as long as the time it takes to compute the initial
#' test statistic. If the time limit is too restrictive thereby ending the
#' computation of the test-statistic too early, we skip the step and move a bit
#' forward in the line search code and temporarily increase the time limit for
#' this new step.
#'
#' The interpolation procedure is as follows:
#' We compute the test-statistic for two extreme values of beta_D, one to the
#' left and another to the right of the initial milp estimate. Then, we
#' (linearly) interpolate the values of the test-stat using these three values
#' (left, right, and initial aka milp estimate). Then using the interpolation,
#' we come up with our best guess for the value of the beta_D that produces the
#' test-stat closest to the critical value. We then intelligently choose the
#' step size and perform the line search as usual.
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
#' @param initial_TimeLimit Time limit (in sec) for each iteration of
#'  preprocessing in the computation of initial test-statistic; defaults to
#'  heuristic (scalar)
#' @param initial_globalTimeLimit Time limit (in sec) for computation of
#'  initial test-statistic; defaults to Inf (scalar)
#' @param save_log If TRUE, save log files from Gurobi programs in
#' \code{log_dir}; defaults to FALSE (boolean)
#' @inheritParams test_stat
#'
#' @import foreach
#' @importFrom stats rnorm
#'
#' @export
line_confint_interpolation <- function(index,
                                       endogeneous = TRUE,
                                       beta_null = NULL,
                                       stopping_tolerance = NULL,
                                       width_ratio = 1,
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
                                       residuals = NULL,
                                       show_progress = TRUE,
                                       # FUN = preprocess_iqr_milp,
                                       initial_TimeLimit = NULL,
                                       initial_globalTimeLimit = Inf,
                                       save_log = FALSE,
                                       ...) {

  ### Prep -------------------------

  # Start clock
  clock_start <- Sys.time()
  message(paste("Clock started:", clock_start))

  out <- list() # Initialize list of results to return
  out$homoskedasticity <- homoskedasticity
  kernel <- ifelse(homoskedasticity, "homoskedasticity", kernel)
  out$kernel <- kernel
  out$alpha <- alpha

  crit_val <- stats::qnorm(1 - alpha / 2)
  out$crit_val <- crit_val

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

  ### Determine stopping tolerance via naive QR -------------------------
  if (is.null(stopping_tolerance) | is.null(beta_null)) {
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
      crit_val <- stats::qnorm(1 - alpha / 2)
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
  }
  # return(out) # DEBUG
  stopifnot(is.numeric(stopping_tolerance))
  stopifnot(is.numeric(beta_null))
  out$stopping_tolerance <- stopping_tolerance
  out$beta_null <- beta_null

  send_note_if(paste0("left bound: ", bounds[1]), show_progress, message)
  send_note_if(paste0("right bound: ", bounds[2]), show_progress, message)
  send_note_if(paste0("width: ", width), show_progress, message)
  send_note_if(paste0("stopping tolerance: ", stopping_tolerance),
               show_progress, message)

  if (return_setup) {
    return(out)
  }

  ### Initial test-stat -------------------------

  # Construct null hypothesis
  beta_D_null <- rep(NA, p_D)
  beta_X_null <- rep(NA, p_X)
  if (endogeneous) {
    beta_D_null[index] <- beta_null
  } else {
    beta_X_null[index] <- beta_null
  }

  # TODO: save initial test-stat info in CSV with suffix counter = 0
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
    residuals = residuals,
    TimeLimit = initial_TimeLimit,
    globalTimeLimit = initial_globalTimeLimit,
    LogFileName = ifelse(save_log,
                         paste0(log_dir, "/", date_time, "_initialTestStat"),
                         ""),
    # FUN = FUN, # TODO: allow user to choose between preprocess_iqr_milp and iqr_milp
    ...
  )
  if (!is.null(initial_test_stat$ended_early)) {
    return(out) # diagnostic messages printed by test_stat()
  }
  out$initial_test_stat <- initial_test_stat
  initial_rejected <- initial_test_stat$p_val < alpha
  initial_time <- as.numeric(initial_test_stat$time_elapsed * 60) # time in seconds
  if (initial_rejected) {
    warning("Initial value of beta_null was rejected.")
    return(out)
    # TODO: do a line search to find a value in the confidence interval
  }

  msg <- paste("Computed initial test statistic:", Sys.time())
  send_note_if(msg, show_progress, message)

  ### figure out left and right boundaries for interpolation -------------------------
  wald <- wald_univariate(
    center = beta_null,
    endogeneous = endogeneous,
    index = index,
    resid = residuals,
    alpha = alpha,
    tau = tau,
    D = D,
    X = X,
    Z = Z,
    Phi = Phi
  )
  left_beta <- wald$lower
  right_beta <- wald$upper
  out$wald <- wald

  ### Interpolation Exercise -------------------------

  # TODO: parallelize computation of left_test_stat and right_test_stat?

  # evaluate left boundary
  beta_D_null <- rep(NA, p_D)
  beta_X_null <- rep(NA, p_X)
  if (endogeneous) {
    beta_D_null[index] <- left_beta
  } else {
    beta_X_null[index] <- left_beta
  }
  left_test_stat <- test_stat(
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
    residuals = residuals,
    TimeLimit = initial_TimeLimit,
    globalTimeLimit = initial_globalTimeLimit,
    LogFileName = ifelse(save_log,
                         paste0(log_dir, "/", date_time, "_initialTestStat"),
                         ""),
    # FUN = FUN, # TODO: allow user to choose between preprocess_iqr_milp and iqr_milp
    ...
  )
  out$left <- left_test_stat

  # evaluate right boundary
  beta_D_null <- rep(NA, p_D)
  beta_X_null <- rep(NA, p_X)
  if (endogeneous) {
    beta_D_null[index] <- right_beta
  } else {
    beta_X_null[index] <- right_beta
  }
  right_test_stat <- test_stat(
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
    residuals = residuals,
    TimeLimit = initial_TimeLimit,
    globalTimeLimit = initial_globalTimeLimit,
    LogFileName = ifelse(save_log,
                         paste0(log_dir, "/", date_time, "_initialTestStat"),
                         ""),
    # FUN = FUN, # TODO: allow user to choose between preprocess_iqr_milp and iqr_milp
    ...
  )
  out$right <- right_test_stat

  # interpolation of left_test_stat, right_test_stat, and initial_test_stat
  interpolation_data <- data.frame(
    x = c(beta_null, left_beta, right_beta),
    y = c(initial_test_stat$test_stat, left_test_stat$test_stat, right_test_stat$test_stat)
  )
  reg <- lm(y ~ x, data = interpolation_data)
  reg_coeffs <- coef(reg)
  reg_shiftdown <- reg_coeffs
  reg_shiftdown['(Intercept)'] <- reg_coeffs['(Intercept)'] - crit_val
  beta_1 <- polyroot(coef(reg_shiftdown))
  reg_shiftup <- reg_coeffs
  reg_shiftup['(Intercept)'] <- reg_coeffs['(Intercept)'] + crit_val
  beta_2 <- polyroot(coef(reg_shiftup))
  # determine min_beta_candidate and max_beta_candidate
  min_beta_candidate <- min(beta_1, beta_2)
  max_beta_candidate <- max(beta_1, beta_2)
  # determine whether *_beta_candidate is associated with positive or negative crit_val
  # i.e., is the test-stat incr. or decr. with respect to beta_D
  # this is useful when computing the "error" from our interpolation (see min_delta_y and max_delta_y)
  min_crit_val <- ifelse(min_beta_candidate == beta_1, -crit_val, crit_val)
  max_crit_val <- ifelse(max_beta_candidate == beta_1, -crit_val, crit_val)
  stopifnot(min_crit_val == -max_crit_val)

  ### Find step size and direction for *_beta_candidate -------------------------

  # evaluate min_beta_candidate
  beta_D_null <- rep(NA, p_D)
  beta_X_null <- rep(NA, p_X)
  if (endogeneous) {
    beta_D_null[index] <- min_beta_candidate
  } else {
    beta_X_null[index] <- min_beta_candidate
  }
  min_test_stat <- test_stat(
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
    residuals = residuals,
    TimeLimit = initial_TimeLimit,
    globalTimeLimit = initial_globalTimeLimit,
    LogFileName = ifelse(save_log,
                         paste0(log_dir, "/", date_time, "_minBetaCandidate"),
                         ""),
    # FUN = FUN, # TODO: allow user to choose between preprocess_iqr_milp and iqr_milp
    ...
  )
  out$min_test_stat <- min_test_stat
  # figure out min_step_size and min_direction
  # step size is determined as follows:
  #   1. if candidate is inside confidence interval, find a nearby point outside confidence interval
  #      if candidate is outside confidence interval, find a nearby point inside confidence interval
  #   2. linearly interpolate between these two points, and find the slope
  #   3. measure the error from the crit_val and actual test stat of candidate: actual test - crit_val
  min_ts_reject <- min_test_stat$p_val < alpha
  if (min_ts_reject) {
    min_direction <- 1 # if I reject, it means I'm to the left of the lower bound
    min_interpolation_data <- data.frame(
      x = c(min_beta_candidate, beta_null), # interpolate with a point inside the confidence interval
      y = c(min_test_stat$test_stat, initial_test_stat$test_stat)
    )
  } else {
    min_direction <- -1 # if I don't reject, it means I'm to the right of the lower bound
    min_interpolation_data <- data.frame(
      x = c(min_beta_candidate, left_beta), # interpolate with the closest point outside the confidence interval
      y = c(min_test_stat$test_stat, left_test_stat$test_stat)
    )
  }
  min_reg <- lm(y ~ x, data = min_interpolation_data)
  min_slope <- coef(min_reg)['x']
  min_delta_y <- min_test_stat$test_stat - min_crit_val
  min_step_size <- abs(min_delta_y / min_slope) # TODO: what does the sign of this quantity mean?

  # evaluate max_beta_candidate
  beta_D_null <- rep(NA, p_D)
  beta_X_null <- rep(NA, p_X)
  if (endogeneous) {
    beta_D_null[index] <- max_beta_candidate
  } else {
    beta_X_null[index] <- max_beta_candidate
  }
  max_test_stat <- test_stat(
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
    residuals = residuals,
    TimeLimit = initial_TimeLimit,
    globalTimeLimit = initial_globalTimeLimit,
    LogFileName = ifelse(save_log,
                         paste0(log_dir, "/", date_time, "_maxBetaCandidate"),
                         ""),
    # FUN = FUN, # TODO: allow user to choose between preprocess_iqr_milp and iqr_milp
    ...
  )
  out$max_test_stat <- max_test_stat
  # figure out max_step_size and max_direction
  # step size is determined as follows:
  #   1. if candidate is inside confidence interval, find a nearby point outside confidence interval
  #      if candidate is outside confidence interval, find a nearby point inside confidence interval
  #   2. linearly interpolate between these two points, and find the slope
  #   3. measure the error from the crit_val and actual test stat of candidate: actual test - crit_val
  max_ts_reject <- max_test_stat$p_val < alpha
  if (max_ts_reject) {
    max_direction <- -1 # if I reject, it means I'm to the right of the upper bound
    max_interpolation_data <- data.frame(
      x = c(max_beta_candidate, beta_null), # interpolate with a point inside the confidence interval
      y = c(max_test_stat$test_stat, initial_test_stat$test_stat)
    )
  } else {
    max_direction <- 1 # if I don't reject, it means I'm to the left of the upper bound
    max_interpolation_data <- data.frame(
      x = c(max_beta_candidate, right_beta), # interpolate with the closest point outside the confidence interval
      y = c(max_test_stat$test_stat, right_test_stat$test_stat)
    )
  }
  max_reg <- lm(y ~ x, data = max_interpolation_data)
  max_slope <- coef(max_reg)['x']
  max_delta_y <- max_test_stat$test_stat - max_crit_val
  max_step_size <- abs(max_delta_y / max_slope) # TODO: what does the sign of this quantity mean?


  ### Line Search -------------------------

  msg <- paste("Starting Line Search", Sys.time())
  send_note_if(msg, show_progress, message)

  # set up cluster
  out$cores <- cores
  cl <- parallel::makeCluster(cores, outfile = paste0(log_dir, "/", date_time, "_outfile.txt")) # TODO: let user decide where and when to save the outfile
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl))
  confint <- foreach (type = c("min", "max"),
                      .combine = c,
                      .export = c("test_stat",
                                  "send_note_if",
                                  "preprocess_iqr_milp")) %dopar% {
    # type <- ifelse(direction == 1, "max", "min")
    step <- ifelse(type == "min", min_step_size, max_step_size) # initial step size
    direction <- ifelse(type == "min", min_direction, max_direction) # initial direction
    current_beta <- ifelse(type == "min", min_beta_candidate, max_beta_candidate)
    beta_D_null <- rep(NA, p_D)
    beta_X_null <- rep(NA, p_X)
    current_ts_reject <- NA
    current_p_val <- NA
    counter <- 0
    small_change_in_p_val <- 0
    # set the time limit to be 2 * time limit of initial test-stat computation
    time_limit <- initial_time * 2
    init <- TRUE # first iteration of while loop shouldn't move from current_beta
    ts <- list(); ts$ended_early <- NULL # initialize ts$ended_early to be NULL
    # TODO: stop the while loop if the width of the confidence interval is very large (think about how to present the results)
    while (step > stopping_tolerance &
           small_change_in_p_val < small_change_count_tol) {
      clock_start_while <- Sys.time() # start the clock for current step
      counter <- counter + 1 # update counter
      send_note_if(paste("Type & Counter:", type, counter), show_progress, message) # DEBUG: consider removing this
      old_beta <- current_beta # save previous beta
      if (!is.null(ts$ended_early)) {
        message(paste(type, counter, "PREVIOUSLY ENDED EARLY")) # DEBUG: remove later
        # If the previous test-stat computation ended early, then we'll take a
        # tiny step in the same direction with a small perturbation and then
        # continue the line search.
        # This only applies *after* we try one iteration of the while loop
        # (else, ts$ended_early is undefined)
        perturb <- stats::rnorm(1, mean = 0, sd = step / 5)
        current_beta <- current_beta + step * direction * 0.5 + perturb
      } else {
        message(paste(type, counter, "PREVIOUSLY DID NOT END EARLY")) # DEBUG: remove later
        # If the previous test-stat computation successfully ran, then we'll
        # take a step in the specified direction and then continue the line
        # search.
        old_p_val <- current_p_val # save previous p-value if previous step wasn't skipped
        old_ts_reject <- current_ts_reject # save previous status of test if previous step wasn't skipped
        if (!init) {
          current_beta <- current_beta + step * direction # update beta
        } else {
          # don't change current_beta and flip init to FALSE
          init <- FALSE
        }
      }
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
        residuals = residuals,
        TimeLimit = time_limit,
        globalTimeLimit = time_limit,
        LogFileName = ifelse(save_log,
                             paste0(log_dir, "/", date_time, "_", type, "_step", counter),
                             ""),
        # FUN = preprocess_iqr_milp, # TODO: let user decide what this is
        ...
      )
      # If the test-stat takes too long, we end it early and skip the step.
      # We then move forward in the line search by taking a smaller step in the
      # same direction with a slight perturbation (see above where we define
      # the next step).
      # Otherwise, we record the results, check if the p-value flattens, and
      # change the direction if need be.
      # message(paste(type, counter, "did not end early:", is.null(ts$ended_early))) # DEBUG: remove later
      if (!is.null(ts$ended_early)) {
        message(paste(type, counter, "ENDED EARLY")) # DEBUG: remove later
        time_limit <- time_limit * 2 # double the time limit temporarily
        current_p_val <- "skipped"
        current_ts_reject <- "skipped"
        current_ts <- "skipped"
      } else {
        message(paste(type, counter, "DID NOT END EARLY")) # DEBUG: remove later
        time_limit <- initial_time * 2 # reset time limit
        # determine test status
        current_p_val <- ts$p_val
        current_ts_reject <- ts$p_val < alpha
        current_ts <- ts$test_stat
        # check if change in p-value is small
        if (abs(current_p_val - old_p_val) < p_val_tol) {
          small_change_in_p_val <- small_change_in_p_val + 1
        } else {
          small_change_in_p_val <- 0
        }
        message(paste(type, counter, "Update p-val count")) # DEBUG: remove later
        # update direction and step size
        if (current_ts_reject) {
          # if we reject, then we want to move to the right if we're searching
          # for the lower bound or move to the left if we're searching for the
          # upper bound.
          direction <- ifelse(type == "min", 1, -1)
        } else {
          # if we fail to reject, then we want to move the left if we're
          # searching for the lower bound or move to the right if we're
          # searching for the upper bound.
          direction <- ifelse(type == "min", -1, 1)
        }
        step <- step * step_rate # TODO: determine step rate!!!
        ### OLD CODE, TODO: remove
        ### if (old_ts_reject != current_ts_reject) {
        ###   direction <- -1 * direction
        ###   step <- step * step_rate
        ### }
        message(paste(type, counter, "Update direction (if necessary)")) # DEBUG: remove later
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
                   current_ts,          # current test-stat
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
        message(paste(type, counter, "SAVE RESULTS")) # DEBUG: remove later
      }
    } # end of while loop

    message(paste(type, "leave while loop:", TRUE)) # DEBUG: remove later
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
  } # end of for loop

  message(paste("leave foreach loop:", TRUE)) # DEBUG: remove later
  flattens_min <- grepl("p-val flattens for min", confint)
  flattens_max <- grepl("p-val flattens for max", confint)
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
    if (remove_intermediate_csvs) {
      # remove files saved during while loop
      file.remove(list.files(log_dir, pattern = date_time, full.names = T))
    } else {
      concatenate_csvs(
        log_dir,
        cols = c("type",
                 "counter",
                 "index",
                 "homoskedasticity",
                 "kernel",
                 "endogeneous",
                 "current_beta",
                 "current_p_val",
                 "test_stat",
                 "current_ts_reject",
                 "elapsed_while",
                 "new_step",
                 "new_direction",
                 "stopping_tolerance"),
        pattern = paste0(date_time, "_", stub, "_", "min", "_counter_"),
        merged_name = paste0(date_time, "_", stub, "_", "min", ".csv"),
        remove_after_merge = TRUE,
        header = FALSE
      )
      concatenate_csvs(
        log_dir,
        cols = c("type",
                 "counter",
                 "index",
                 "homoskedasticity",
                 "kernel",
                 "endogeneous",
                 "current_beta",
                 "current_p_val",
                 "test_stat",
                 "current_ts_reject",
                 "elapsed_while",
                 "new_step",
                 "new_direction",
                 "stopping_tolerance"),
        pattern = paste0(date_time, "_", stub, "_", "max", "_counter_"),
        merged_name = paste0(date_time, "_", stub, "_", "max", ".csv"),
        remove_after_merge = TRUE,
        header = FALSE
      )
    }
  }

  return(out)
}
