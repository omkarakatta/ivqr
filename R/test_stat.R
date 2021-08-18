### Meta -------------------------
###
### Title: Create the test statistic and p-value under strong identification
###
### Description: Under the null hypothesis that beta_D,J and beta_X,J is a
### specified vector, we can create a test statistic under the assumption of
### homoskedasticity or heteroskedasticity. If |J| + |K| = 1, then we have a
### test statistic that is a linear form. If |J| + |K| > 1, then we have a test
### statistic that is a quadratic form. The test statistic, in a nutshell,
### evaluates the dual feasibility constraint constructed with objects that are
### the result of a short IQR problem where we concentrate out D_J and X_K from
### Y with the null values of beta_D,J and beta_X,J.
###
### Author: Omkar A. Katta
###


### test_stat -------------------------
#' Compute the test statistic and p-value under strong identification
#'
#' Under user-specified null hypotheses (multivariate or univariate) on the
#' endogeneous and exogeneous variables' coefficients, compute the test
#' statistic and p-value associated with a specified \code{alpha}-level under
#' homoskedasticity or heteroskedasticity.
#'
#' For example, if \code{beta_D_null = c(NA, 0, NA)} and \code{beta_X_null =
#' c(1, NA)}, then \eqn{\beta_{D,2} = 0} and \eqn{\beta_{X,1} = 1} is the null
#' hypothesis.
#' Note that in this example, p_D is 3 and p_X is 2, which corresponds to the
#' lengths of the coefficient vectors.
#'
#' Based on \code{beta_D_null}, we find the vector J, which contains the
#' indices of the coefficients specified under the null.
#' We can then define D_J and Phi_J as the columns of D and Phi whose indices
#' are in J.  Similarly, we can define D_J_minus and Phi_J_minus as the columns
#' of D and Phi whose indices are not in J.  If J is the empty set, then D_J
#' and Phi_J should be thought of as being empty matrices without any
#' dimensions.  If J specified all indices from 1 to p_D, then D_J_minus and
#' Phi_J_minus should be thought of as being empty.  These conventions apply to
#' the vector K, which contains the indices of the coefficients on the
#' endogeneous variable specified uner the null.
#'
#' If the test is univariate (i.e., |J| + |K| = 1), then the test statistic
#' has a linear form. Otherwise, the test statistic has a quadratic form.
#'
#' @param beta_D_null Vector of coefficients on the endogeneous variable under
#'  the null hypothesis; if a coefficient is not specified under the null, let
#'  the corresponding entry be NA (vector of length p_D with |J| non-zero
#'  entries)
#' @param beta_X_null Vector of coefficients on the covariates under
#'  the null hypothesis; if a coefficient is not specified under the null, let
#'  the corresponding entry be NA (vector of length p_X with |K| non-zero
#'  entries)
#' @param alpha Alpha level of the test; defaults to 0.1; only used when
#'  \code{homoskedasticity} is FALSE (numeric between 0 and 1)
#' @param Y Dependent variable (vector of length n)
#' @param X Exogenous variable (including constant vector) (n by p_X matrix)
#' @param D Endogenous variable (n by p_D matrix)
#' @param Z Instrumental variable (n by p_Z matrix)
#' @param Phi Transformation of X and Z to be used in the program;
#'  defaults to the linear projection of D on X and Z (matrix with n rows)
#' @param tau Quantile (number between 0 and 1)
#' @param B Matrix that enters numerator of test statistic
#'  (n by |J| + |K| matrix); If NULL (default), this matrix is (Phi_J, X_K),
#'  where Phi_J and X_K are the columns of Phi and X with indices J and K;
#' @param orthogonalize_statistic If TRUE, \eqn{\tilde{B}} will be used in
#'  numerator of test statistic; defaults to FALSE; for advanced users only
#' @param homoskedasticity If TRUE, assume density of error at 0 is constant;
#'  defaults to FALSE (boolean)
#' @param kernel Only active if \code{homoskedasticity} is FALSE; either
#'  "Powell" (default) to use the Powell estimator or
#'  "Gaussian" to use a Gaussian kernel; only used when
#'  \code{homoskedasticity} is FALSE
#' @param a_hat Vector (n by 1) of dual variables; if NULL (default), use
#'  dual-variables from short-iqr regression
#' @param residuals Residuals from IQR MILP program; if NULL (default), use
#'  residuals from short-iqr regression
#' @param show_progress If TRUE (default), sends progress messages during
#'  execution (boolean); also passed to \code{preprocess_iqr_milp}
#' @param print_results If TRUE (default), print the test-statistic, p-value,
#'  and alpha level (boolean)
#' @param ... Arguments passed to \code{preprocess_iqr_milp}
test_stat <- function(beta_D_null,
                      beta_X_null,
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
                      a_hat = NULL,
                      residuals = NULL,
                      kernel = "Powell",
                      show_progress = TRUE,
                      print_results = TRUE,
                      # FUN = preprocess_iqr_milp,
                      ...) {

  # Start clock
  clock_start <- Sys.time()
  msg <- paste("Clock started:", clock_start)
  send_note_if(msg, show_progress, message)

  out <- list() # Initialize list of results to return

  # Check arguments
  # if (!identical(FUN, preprocess_iqr_milp) & !identical(FUN, iqr_milp)) {
  #   stop("`FUN` should either be `preprocess_iqr_milp` or `iqr_milp`.")
  # }
  kernel <- tolower(kernel)
  if (kernel != "powell" & kernel != "gaussian" & !homoskedasticity) {
    stop(paste0("`kernel` should either be 'Powell' or 'Gaussian', not ",
                kernel))
  }

  # Get dimensions of data
  n <- length(Y)
  p_D <- ncol(D)
  p_X <- ncol(X)

  # Check dimensions of vectors under null hypothesis
  stopifnot(length(beta_D_null) == p_D)
  stopifnot(length(beta_X_null) == p_X)
  stopifnot(alpha < 1 & alpha > 0)

  out$beta_D_null <- beta_D_null
  out$beta_X_null <- beta_X_null
  out$alpha <- alpha

  send_note_if("Obtained data", show_progress, message)

  # Get indices as a boolean vector: TRUE => specified under null, FALSE => o/w
  J <- !is.na(beta_D_null)
  K <- !is.na(beta_X_null)
  cardinality_J <- sum(J)
  cardinality_K <- sum(K)
  beta_D_null_nontrivial <- beta_D_null[J]
  beta_X_null_nontrivial <- beta_X_null[K]

  if (cardinality_J + cardinality_K == 0) {
    stop("Invalid null: beta_D_null and beta_X_null cannot be all NA values.")
  }

  # Construct matrices
  D_J <- D[, J, drop = FALSE]
  D_J_minus <- D[, !J, drop = FALSE]
  X_K <- X[, K, drop = FALSE]
  X_K_minus <- X[, !K, drop = FALSE]
  # Z_J <- Z[, J, drop = FALSE] # not used
  Z_J_minus <- Z[, !J, drop = FALSE]
  Phi_J <- Phi[, J, drop = FALSE]
  Phi_J_minus <- Phi[, !J, drop = FALSE]

  send_note_if("Constructed basic matrices", show_progress, message)

  # Concentrate out D_J and X_K
  tmp <- Y
  if (cardinality_J > 0) {
    tmp <- tmp - D_J %*% beta_D_null_nontrivial
  }
  if (cardinality_K > 0) {
    tmp <- tmp - X_K %*% beta_X_null_nontrivial
  }
  Y_tilde <- tmp

  send_note_if("Concentrated out Y", show_progress, message)

  # Obtain \hat{a} and residuals via short-iqr regression
  if (is.null(a_hat) | is.null(residuals)) {
    if (ncol(D_J_minus) == 0) {
      # If there are no endogeneous variables, return quantile regression results:
      msg <- paste("p_D is 0 -- running QR instead of IQR MILP...")
      send_note_if(msg, TRUE, warning)
      qr <- quantreg::rq(Y_tilde ~ X_K_minus - 1, tau = tau)
      short_iqr <- qr
      FUN <- "qr"
    } else {
      short_iqr <- preprocess_iqr_milp(
        Y = Y_tilde,
        X = X_K_minus,
        D = D_J_minus,
        Z = Z_J_minus, # not really important since we specify Phi
        Phi = Phi_J_minus,
        tau = tau,
        show_progress = show_progress,
        ...
      )
      FUN <- "preprocess_iqr_milp"
    }

    send_note_if("Computed short-IQR solution", show_progress, message)

    if (FUN == "preprocess_iqr_milp") {
      short_iqr_result <- short_iqr$final_fit
    } else if (FUN == "qr") {
      short_iqr_result <- short_iqr
    }
    out$short_iqr <- short_iqr_result
  }

  # residuals for estimating variance matrix aka Psi
  # short-iqr residuals are used to construct a_hat (and consequently b_hat)
  # `residuals` argument is used to construct Psi matrix (aka variance of residuals)
  # If the `residuals` argument is not provided, we use the short-iqr residuals.
  if (is.null(residuals)) {
    msg <- "`residuals` not provided for variance estimation; using short-iqr residuals instead"
    warning(msg)
    resid <- short_iqr_result$resid
    out$resid <- resid
  } else {
    resid <- residuals
    out$resid <- resid
  }

  # dual variables for computing test-stat
  if (is.null(a_hat)) {
    # use dual variables from short-iqr regression
    # To get a_hat under full-vector inference (i.e., cardinality_J == p_D),
    # we run a quantile regression, not the usual IQR MILP procedure.
    # Hence, we don't need to print "status" or "objval" if we are in
    # the full-vector setting because quantile regression doesn't return such
    # results.
    if (cardinality_J != p_D) {
      if (short_iqr_result$status != "OPTIMAL") {
        warning(paste("Short IQR Status:", short_iqr_result$status))
      }
      if (short_iqr_result$status == "TIME_LIMIT" || short_iqr_result$objval != 0) {
        message(paste("Short IQR Objective Value:", short_iqr_result$objval))
        out$ended_early <- TRUE
        return(out)
      }
      a_hat <- short_iqr_result$a
    } else {
      # In full-vector inference, we use the dual variables of quantile
      # regression to define a_hat.
      a_hat <- short_iqr_result$dual
    }
    out$a_hat <- a_hat
  } else {
    # user-specified a_hat vector
    out$a_hat <- a_hat
  }

  # Obtain \hat{b}
  b_hat <- a_hat - (1 - tau)

  # Create B
  if (is.null(B)) {
    B <- cbind(Phi_J, X_K)
  }
  stopifnot(nrow(B) == n)
  stopifnot(ncol(B) == cardinality_J + cardinality_K)
  out$B <- B

  send_note_if("Created `B`", show_progress, message)

  # Create \tilde{B} depending on homoskedasticity or heteroskedasticity
  B_minus <- cbind(Phi_J_minus, X_K_minus)
  stopifnot(nrow(B_minus) == n)
  stopifnot(ncol(B_minus) == p_D - cardinality_J + p_X - cardinality_K)
  C_minus <- cbind(D_J_minus, X_K_minus)
  stopifnot(nrow(C_minus) == n)
  stopifnot(ncol(C_minus) == p_D - cardinality_J + p_X - cardinality_K)

  if (homoskedasticity) {
    Psi <- diag(1, nrow = n)
  } else {
    # Hall and Sheather (1988) bandwidth
    tmp_a <- n ^ (1 / 3)
    tmp_b <- stats::qnorm(1 - 0.5 * alpha) ^ (2 / 3)
    tmp_c <- 1.5 * (stats::dnorm(stats::qnorm(tau)) ^ 2)
    tmp_d <- 2 * (stats::qnorm(tau) ^ 2) + 1
    hs <- tmp_a * tmp_b * ((tmp_c / tmp_d) ^ (1 / 3))
    out$hs <- hs
    if (kernel == "powell") {
      bw <- hs
      # Note that the 1 / (2 * n * bw) is negated in the formula for B_tilde
      Psi <- diag(as.numeric(abs(resid) < bw), nrow = n, ncol = n)
    } else if (kernel == "gaussian") {
      bw <- hs
      # Note that the 1 / (n * bw) is negated in the formula for B_tilde
      Psi <- diag(stats::dnorm(resid / bw), nrow = n, ncol = n)
    } else {
      stop(
       "Let `homoskedasticity` be TRUE or choose an appropriate `kernel`."
      )
    }
  }
  G_minus <- Psi %*% C_minus
  B_tilde <- B - B_minus %*% solve(t(G_minus) %*% B_minus) %*% t(G_minus) %*% B

  send_note_if("Created `B_tilde`", show_progress, message)

  out$Psi <- Psi
  out$B_tilde <- B_tilde
  out$homoskedasticity <- homoskedasticity
  out$kernel <- ifelse(homoskedasticity, "homoskedasticity", kernel)

  # Construct S_n with orthogonalize_statistic
  if (orthogonalize_statistic) {
    S_n <- n ^ (-1 / 2) * t(B_tilde) %*% b_hat
  } else {
    S_n <- n ^ (-1 / 2) * t(B) %*% b_hat
  }
  stopifnot(nrow(S_n) == cardinality_J + cardinality_K)
  stopifnot(ncol(S_n) == 1)

  out$S_n <- S_n

  send_note_if("Created `S_n`", show_progress, message)

  # Construct L_n or Q_n depending on |J| + |K| == or != 1
  # Compute p-value
  if (cardinality_J + cardinality_K == 1) {
    denom <- sqrt(tau * (1 - tau) * t(B_tilde) %*% B_tilde / n)
    out$denom <- denom
    test_stat <- S_n / denom
    p_val <- 2 * (1 - stats::pnorm(abs(test_stat)))
  } else {
    M_n <- (1 / n) * t(B_tilde) %*% B_tilde
    out$M_n <- M_n
    test_stat <- t(S_n) %*% solve(M_n) %*% S_n / (tau * (1 - tau))
    p_val <- 1 - stats::pchisq(test_stat, df = cardinality_J + cardinality_K)
  }

  send_note_if("Computed test statistic and p-value", show_progress, message)

  if (print_results) {
    print(paste("Alpha:", alpha))
    print(paste("Test Statistic:", test_stat))
    print(paste("p-value:", p_val))
  }

  stopifnot(nrow(test_stat) == 1)
  stopifnot(ncol(test_stat) == 1)
  out$test_stat <- as.numeric(test_stat)
  out$p_val <- p_val

  # Stop the clock
  clock_end <- Sys.time()
  elapsed_time <- difftime(clock_end, clock_start, units = "mins")
  out$time_elapsed <- elapsed_time

  if (print_results) {
    print(paste("Time Elapsed (mins):", elapsed_time))
  }

  return(out)
}
