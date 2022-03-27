# We have several notions of the "center of FOC polytope".
# `foc_center_method` is a wrapper for these notions, most of which are in
# `R/foc-center-archive.R`.
# `foc_center` is the prevailing notion that we use.

# foc_center_method ------------------------------------------------------------

foc_center_method <- function(FUN, ...) {
  FUN(...)
}

# find_rounded_center ----------------------------------------------------------

#' @param omega_mod See `omega_mod` from output of `foc_center`
find_rounded_center <- function(omega_mod, n, h) {
  rounded_center_indices <- setdiff(seq_len(n), h)[which(omega_mod == 1)]
  rounded_center <- vector("integer", n)
  rounded_center[rounded_center_indices] <- 1
  rounded_center[h] <- 1
  rounded_center
}

# find_continuous_center -------------------------------------------------------

find_continuous_center <- function(omega, n, h) {
  cts_center_indices <- setdiff(seq_len(n), h)[which(omega > 0)]
  cts_center <- vector("integer", n)
  cts_center[cts_center_indices] <- omega[which(omega > 0)]
  cts_center[h] <- 1
  cts_center
}

# foc_center -------------------------------------------------------------------

# If beta_D is NULL, we will compute it from `h`
# `h` needs to be in terms of the data arguments
foc_center <- function(
  h,
  beta_D = NULL,
  beta_X = NULL,
  subsample_size,
  Y, X, D, Phi,
  tau,
  params = list(OutputFlag = 0)
) {
  n <- nrow(Y)
  p <- length(h)

  if (n - p - 1 <= 0) {
    stop("`facet_center` denominator must be positive; we need n - p - 1 > 0.")
  }

  tmp <- compute_foc_conditions(h, beta_D, beta_X,
                                Y = Y, X = X, D = D, Phi = Phi, tau = tau)
  xi_mat <- tmp[, setdiff(seq_len(n), h), drop = FALSE]

  # Decision variables in order from left/top to right/bottom
  # 1. omega --- (n - p) by 1
  # 2. right_slack --- p by 1
  # 3. left_slack --- p by 1
  # 4. simplex_slack --- 2(n - p) by 1
  # 5. max_slack --- 1 by 1
  num <- list() # keep track of dimensions
  num$omega <- n - p
  num$right_slack <- p
  num$left_slack <- p
  num$simplex_slack <- 2 * (n - p)
  num$max_slack <- 1
  num$decision_vars <- sum(unlist(num))

  # Set Up Model
  model <- list()
  model$modelsense <- "min"
  model$lb <- rep(0, num$decision_vars)
  model$ub <- c(
    rep(1, num$omega),
    rep(Inf, num$decision_vars - num$omega)
  )
  model$obj <- c(
    rep(0, num$decision_vars - num$max_slack),
    rep(1, num$max_slack)
  )
  model$vtype <- rep("C", num$decision_vars)

  # Constraint 1: sum of omega is m-p (*_omega)
  A_omega <- matrix(c(
    rep(1, num$omega), rep(0, num$decision_vars - num$omega)
  ), nrow = 1)
  sense_omega <- "="
  rhs_omega <- subsample_size - p

  # Constraint 2: FOC slack variables (2p constraints) (*_foc)
  A_right <- vector("list", p)
  A_left <- vector("list", p)
  ones <- diag(1, p)
  zeros <- diag(0, p)
  A_right <- do.call(cbind, list(xi_mat, ones, zeros))
  A_left <- do.call(cbind, list(-1 * xi_mat, zeros, ones))
  A_foc <- rbind(A_right, A_left)
  A_foc <- expand_matrix(
    A_foc, newrow = 2 * p, newcol = num$decision_vars,
    row_direction = "bottom", col_direction = "right"
  )
  sense_foc <- rep("=", 2 * p)
  rhs_foc <- c(rep(1 - tau, p), rep(tau, p))

  # Constraint 3: simplex slack variables (2*(n-p) constraints) (*_simplex)
  subsample_center <- rep((subsample_size - p) / (n - p), n - p)
  rhs_simplex <- vector("double", num$simplex_slack)
  lhs_simplex <- vector("list", num$simplex_slack)
  counter <- 0
  zeros <- rep(0, num$simplex_slack)
  for (j in c(0, 1)) {
    for (i in seq_len(n - p)) {
      counter <- counter + 1
      facet_center <- rep((subsample_size - p - j) / (n - p - 1), n - p)
      facet_center[i] <- j
      normal_vec <- subsample_center - facet_center
      normal_unit <- normal_vec / sqrt(sum(normal_vec^2))
      rhs_simplex[[counter]] <- -1 * sum(normal_unit * facet_center)

      e_ij <- zeros
      e_ij[[counter]] <- 1
      lhs_simplex[[counter]] <- c(
        -1 * normal_unit,
        rep(0, num$left_slack + num$right_slack),
        e_ij,
        rep(0, num$max_slack)
      )
    }
  }
  A_simplex <- do.call(rbind, lhs_simplex)
  sense_simplex <- rep("=", num$simplex_slack)

  # Constraint 4: max of all slack variables (*_max)
  model$genconmax <- vector("list", num$max_slack)
  model$genconmax[[1]]$resvar <- num$decision_vars
  model$genconmax[[1]]$vars <- num$omega +
    seq_len(num$right_slack + num$left_slack + num$simplex_slack)

  # Combine constraints
  model$A <- rbind(A_omega, A_foc, A_simplex)
  model$sense <- c(sense_omega, sense_foc, sense_simplex)
  model$rhs <- c(rhs_omega, rhs_foc, rhs_simplex)

  # Solve program
  sol <- gurobi::gurobi(model, params)
  status <- sol$status
  omega <- sol$x[seq_len(num$omega)]

  if (status != "OPTIMAL") {
    return(list(
      status = status,
      status_message = paste("Gurobi status:", status),
      model = model,
      sol = sol,
      omega = omega,
      xi_mat = xi_mat
    ))
  }

  # Convert continuous solution to integral solution
  current_sum <- length(which(omega == 1)) # how many integral 1's do we have
  max_indices <- rank(-omega, ties.method = "random") # first is largest

  # switch the largest numbers that aren't 1
  # remaining non-integral solutions become 0
  switch <- which(max_indices <= subsample_size - p & max_indices > current_sum)
  omega_mod <- omega
  omega_mod[switch] <- 1
  omega_mod[which(omega_mod < 1)] <- 0
  omega_mod <- round(omega_mod, 0) # ensure they are integral

  # find rounded and continuous center
  rounded_center <- find_rounded_center(omega_mod, n, h)
  cts_center <- find_continuous_center(omega, n, h)

  list(
    model = model,
    sol = sol,
    status = status,
    omega = omega,
    omega_mod = omega_mod,
    xi_mat = xi_mat,
    rounded_center = rounded_center,
    cts_center = cts_center
  )
}
