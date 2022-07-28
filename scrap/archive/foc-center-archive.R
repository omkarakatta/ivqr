# see `~/BFI/5_IVQR_GP/fromGuillaume/mcmc_subsampling/center_polytope/`

# find_center_subsample_polytope -------------------------

# The center of the subsample simplex, {D \in [0,1]^n s.t. sum(D_i) = m}, should
# be rep(m/n, length = n). This function creates a program that verifies this.
# Examples:
# tmp <- find_center_subsample_polytope(n = 4, m = 2)
# tmp$center # should be roughly m/n = 1/2
# Q: The log variant of the program showing up as infeasible. Why?
find_center_subsample_polytope <- function(
  n, m,
  gencontype = "power",
  a = ifelse(gencontype == "power", 0.5, exp(1)),
  params = list()
) {
  vertices_index <- combn(n, m)
  num_vertices <- ncol(vertices_index)
  stopifnot(num_vertices == choose(n, m))
  n_ones <- rep(0, n)
  vertices <- t(apply(vertices_index, 2, function(i) {
    n_ones[i] <- 1
    n_ones
  }))
  stopifnot(ncol(vertices) == n)
  stopifnot(nrow(vertices) == num_vertices)
  vertices_i <- c(t(vertices)) # collapse by rows
  stopifnot(length(vertices_i) == num_vertices * n)
  slack_i <- diag(1, length(vertices_i))
  zero_mat <- diag(0, length(vertices_i))
  dv_i <- do.call("rbind", rep(list(diag(1, nrow = n)), length = num_vertices))
  define_slack <- cbind(dv_i,      # coordinates of center
                        -slack_i,  # positive part of slack
                        slack_i,   # negative part of slack
                        zero_mat,  # f(pos slack), where f = sqrt by default
                        zero_mat)  # f(neg slack), where f = sqrt by default
  num_decision_vars <- ncol(define_slack)

  sum_A <- c(rep(1, n), rep(0, num_decision_vars - n))
  sum_rhs <- m

  model <- list()
  model$A <- rbind(define_slack, sum_A)
  model$rhs <- c(vertices_i, sum_rhs)
  model$sense <- "="
  model$modelsense <- "max"
  model$lb <- rep(0, num_decision_vars)
  # model$ub <- c(rep(1, n), rep(Inf, num_decision_vars - n))
  model$ub <- c(rep(1, num_decision_vars)) # slack in each index must be <= 1
  model$obj <- c(
    rep(0, n),
    rep(0, 2 * length(vertices_i)),
    rep(1, 2 * length(vertices_i))
  )

  xstart <- n
  ystart <- n + 2 * length(vertices_i)
  gencon <- vector("list", 2 * length(vertices_i))
  for (j in seq_len(2 * length(vertices_i))) {
    yvar <- ystart + j
    xvar <- xstart + j
    gencon[[j]]$xvar <- xvar
    gencon[[j]]$yvar <- yvar
    gencon[[j]]$a <- a
  }
  if (gencontype == "log") {
    model$genconloga <- gencon
  } else if (gencontype == "power") {
    model$genconpow <- gencon
  }

  sol <- gurobi::gurobi(model, params)

  list(
    model = model,
    sol = sol,
    center = sol$x[seq_len(n)]
  )
}


# find_subsample_in_polytope -------------------------

# goal: given an active basis, find a subsample inside the polytope
# Q: is it possible for there to be multiple solutions? In which case, maybe this doesn't provide a deterministic map from the active basis to the subsample, depending on how Gurobi breaks ties.
#' Find subsample in the FOC polytope
#'
#' It is useful to start our random walk at a subsample that is inside the FOC
#' polytope. To find this polytope, we solve a linear program that aims to
#' minimize the absolute L-1 norm of the weighted average of the xi_i objects.
#' The weights are inside the unit interval and sum to the size of our
#' subsample. Most of the weights will be integral, but some will be between 0
#' and 1. We will round the non-integral solutions so they are integra. If the
#' weight for an observaion is 1, then that observation belongs in the
#' subsample. Otherwise, the observation doesn't belong in the subsample. Due
#' to the rounding, it isn't guaranteed that the resulting subsample belongs in
#' the FOC polytope, but it shouldn't be too far away from the polytope.
#' Note that this program excludes the observations in the active basis -- we
#' must add them back in after-the-fact.
#'
#' @param h Active basis in terms of the data provided
#' @param Y Dependent variable (vector of length n)
#' @param X Exogenous variable (including constant vector) (n by p_X matrix)
#' @param D Endogenous variable (n by p_D matrix)
#' @param Z Instrumental variable (n by p_Z matrix)
#' @param Phi Transformation of X and Z to be used in the program;
#'  defaults to the linear projection of D on X and Z (matrix with n rows)
#' @param tau Quantile (numeric)
#' @param beta_D_proposal Coefficients on the endogeneous variables (vector of
#'  length p_D); if NULL, use \code{h_to_beta} function and the \code{h}
#'  argument to determine \code{beta_D_proposal}
#' @param beta_X_proposal Coefficients on the exogeneous variables (vector of
#'  length p_D); if NULL, use \code{h_to_beta} function and the \code{h}
#'  argument to determine \code{beta_X_proposal}
#' @param subsample_size Size of subsample
#' @param params Named list of parameters to send to Gurobi
#' @param type Either "C" or "B" to denote if our "omega" variables are
#'  continuous or binary; If "C", then the solution is the continuous center
#'  of our FOC polytope; If "B", then the solution is the actual center
# We need something inside the polytope. This continuous method with our
# arbitrary rounding doesn't guarantee this. Our workaround:
# 1. exhaustive search over non-integral solutions
# 2. solve the program as an MILP
# TODO: solve the program as an MILP and see if the result is inside the FOC
# polytope
# TODO: code up improved linear program to find DC (ask GP for the note)
find_subsample_in_polytope <- function(
  h,
  Y, X, D, Z, Phi = linear_projection(D, X, Z),
  tau,
  beta_D_proposal = NULL,
  beta_X_proposal = NULL,
  subsample_size,
  params = list(OutputFlag = 0),
  type = "C"
) {

  n <- nrow(Y)
  p <- length(h)

  # get beta_X_proposal and beta_D_proposal

  if (is.null(beta_D_proposal) | is.null(beta_X_proposal)) {
    coef <- h_to_beta(h, Y = Y, X = X, D = D, Phi = Phi)
    if (is.null(beta_D_proposal)) {
      beta_D_proposal <- coef$beta_D
    }
    if (is.null(beta_X_proposal)) {
      beta_X_proposal <- coef$beta_X
    }
  }
  Y_tilde <- Y - D %*% beta_D_proposal
  design <- cbind(X, Phi)
  designh_inv <- solve(design[h, , drop = FALSE])

  s_i <- vector("list", length = n)
  for (i in seq_len(n)) {
    if (is.element(i, h)) {
      s_i[[i]] <- 0 # if index i is in active basis, set xi to be 0
    } else {
      # NOTE: beta_Phi should be 0
      const <- (tau - as.numeric(Y_tilde[i] - X[i, ] %*% beta_X_proposal < 0))
      s_i[[i]] <- const * design[i, ] %*% designh_inv
    }
  }
  s <- t(do.call(rbind, s_i))
  xi_mat <- s[, setdiff(seq_len(n), h)]
  stopifnot(nrow(xi_mat) == p)
  stopifnot(ncol(xi_mat) == n - p)

  # Decision variables in order from left/top to right/bottom:
  # 1. omega --- (n - p) by 1
  # 2. xi --- (p by 1)
  # 3. ximinus --- (p by 1)
  # 4. xiplus --- (p by 1)

  num_omega <- n - p
  num_xi <- p
  num_ximinus <- num_xi
  num_xiplus <- num_xi
  num_decision_vars <- num_omega + num_xi + num_ximinus + num_xiplus

  model <- list()
  model$modelsense <- "min"
  model$lb <- c(
    rep(0, num_omega),
    rep(-Inf, num_xi),
    rep(0, num_ximinus),
    rep(0, num_xiplus)
  )
  model$ub <- c(
    rep(1, num_omega),
    rep(Inf, num_xi),
    rep(Inf, num_ximinus),
    rep(Inf, num_xiplus)
  )
  model$obj <- c(
    rep(0, num_omega),
    rep(0, num_xi),
    rep(1, num_ximinus),
    rep(1, num_xiplus)
  )
  model$vtype <- c(
    rep(type, num_omega),
    rep("C", num_decision_vars - num_omega)
  )

  # define xiplus and ximinus
  tmp <- diag(1, nrow = num_xi)
  zero_mat <- diag(0, nrow = num_omega)
  xisign_A <- expand_matrix(
    cbind(tmp, -tmp, tmp), # TODO: maybe the final minus should be a plus?
    newrow = num_xi,
    newcol = num_decision_vars,
    row_direction = "bottom",
    col_direction = "left"
  )
  xisign_sense <- rep("=", num_xi)
  xisign_rhs <- rep(0, num_xi)

  # define xi
  tmp <- cbind(xi_mat, diag(-1, nrow = num_xi))
  xi_A <- expand_matrix(
    tmp,
    newrow = num_xi,
    newcol = num_decision_vars,
    row_direction = "bottom",
    col_direction = "right"
  )
  xi_sense <- rep("=", num_xi)
  xi_rhs <- rep(0, num_xi)

  # sum of omega
  omega_A <- matrix(c(
    rep(1, num_omega), rep(0, num_decision_vars - num_omega)
  ), nrow = 1)
  omega_sense <- "="
  omega_rhs <- subsample_size - p

  model$A <- rbind(xisign_A, xi_A, omega_A)
  model$rhs <- c(xisign_rhs, xi_rhs, omega_rhs)
  model$sense <- c(xisign_sense, xi_sense, omega_sense)

  sol <- gurobi::gurobi(model, params)
  status <- sol$status
  omega <- sol$x[seq_len(num_omega)]

  if (status != "OPTIMAL") {
    return(list(
      status = status,
      status_message = paste("Gurobi status:", status),
      model = model,
      sol = sol,
      omega = omega
    ))
  }

  # turn continuous solution into an integral one
  if (type == "C") {
    current_sum <- length(which(omega == 1)) # how many integral 1's do we have?
    remaining <- subsample_size - p - current_sum # how many need to be switched?
    to_be_rounded <- omega > 0 & omega < 1 # which can we switch?
    # NOTE: `rank` works better than `order` when there are ties
    max_indices <- rank(-omega, ties.method = "random") # 1 = largest
    # switch the largest numbers that aren't 1
    switch <- which(max_indices <= current_sum + remaining & max_indices > current_sum)
    omega_mod <- omega
    omega_mod[switch] <- 1
    omega_mod[which(omega_mod < 1)] <- 0
    omega_mod <- round(omega_mod, 0)
    stopifnot(all.equal(sum(omega_mod), subsample_size - p))
  } else if (type == "B") {
    omega_mod <- omega
  }

  list(
    model = model,
    sol = sol,
    status = status,
    omega = omega,
    omega_mod = omega_mod,
    xi_mat = xi_mat
  )
}

### find_chebyschev_center -------------------------

# TODO: document
find_chebyschev_center <- function(
  h,
  Y, X, D, Z, Phi = linear_projection(D, X, Z),
  tau,
  beta_D_proposal = NULL,
  beta_X_proposal = NULL,
  subsample_size,
  params = list(OutputFlag = 0),
  type = "C",
  l_norm = 2
) {

  n <- nrow(Y)
  p <- length(h)

  # get beta_X_proposal and beta_D_proposal

  if (is.null(beta_D_proposal) | is.null(beta_X_proposal)) {
    coef <- h_to_beta(h, Y = Y, X = X, D = D, Phi = Phi)
    if (is.null(beta_D_proposal)) {
      beta_D_proposal <- coef$beta_D
    }
    if (is.null(beta_X_proposal)) {
      beta_X_proposal <- coef$beta_X
    }
  }
  Y_tilde <- Y - D %*% beta_D_proposal
  design <- cbind(X, Phi)
  designh_inv <- solve(design[h, , drop = FALSE])

  s_i <- vector("list", length = n)
  for (i in seq_len(n)) {
    if (is.element(i, h)) {
      s_i[[i]] <- 0 # if index i is in active basis, set xi to be 0
    } else {
      # NOTE: beta_Phi should be 0
      const <- (tau - as.numeric(Y_tilde[i] - X[i, ] %*% beta_X_proposal < 0))
      s_i[[i]] <- const * design[i, ] %*% designh_inv
    }
  }
  s <- t(do.call(rbind, s_i))
  xi_mat <- s[, setdiff(seq_len(n), h)]
  stopifnot(nrow(xi_mat) == p)
  stopifnot(ncol(xi_mat) == n - p)

  # Decision variables in order from left/top to right/bottom:
  # 1. omega --- (n - p) by 1
  # 2. r --- 1 by 1

  num_omega <- n - p
  num_r <- 1
  num_decision_vars <- num_omega + num_r

  model <- list()
  model$modelsense <- "max"
  model$lb <- rep(0, num_decision_vars)
  model$ub <- c(
    rep(1, num_omega),
    rep(Inf, num_r)
  )
  model$obj <- c(
    rep(0, num_omega),
    rep(1, num_r)
  )
  model$vtype <- c(
    rep(type, num_omega),
    rep("C", num_r)
  )

  # sum of omega
  omega_A <- matrix(c(
    rep(1, num_omega), rep(0, num_decision_vars - num_omega)
  ), nrow = 1)
  omega_sense <- "="
  omega_rhs <- subsample_size - p

  # FOC boundary
  right_A <- vector("list", p)
  left_A <- vector("list", p)
  for (j in seq_len(p)) {
    xi_j <- xi_mat[j, ]
    xi_j_norm <- sum(abs(xi_j)^l_norm) ^ (1 / l_norm)
    right_A[[j]] <- c(xi_mat[j, ], xi_j_norm)
    left_A[[j]] <- c(-xi_mat[j, ], xi_j_norm)
  }
  right_A <- do.call(rbind, right_A)
  left_A <- do.call(rbind, left_A)
  foc_A <- rbind(right_A, left_A)
  foc_sense <- rep("<=", 2 * p)
  foc_rhs <- c(rep(1-tau, p), rep(tau, p))

  # constraints
  model$A <- rbind(omega_A, foc_A)
  model$sense <- c(omega_sense, foc_sense)
  model$rhs <- c(omega_rhs, foc_rhs)

  sol <- gurobi::gurobi(model, params)
  status <- sol$status
  omega <- sol$x[seq_len(num_omega)]

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

  # turn continuous solution into an integral one
  if (type == "C") {
    current_sum <- length(which(omega == 1)) # how many integral 1's do we have?
    remaining <- subsample_size - p - current_sum # how many need to be switched?
    to_be_rounded <- omega > 0 & omega < 1 # which can we switch?
    # NOTE: `rank` works better than `order` when there are ties
    max_indices <- rank(-omega, ties.method = "random") # 1 = largest
    # switch the largest numbers that aren't 1
    switch <- which(max_indices <= current_sum + remaining & max_indices > current_sum)
    omega_mod <- omega
    omega_mod[switch] <- 1
    omega_mod[which(omega_mod < 1)] <- 0
    omega_mod <- round(omega_mod, 0)
    stopifnot(all.equal(sum(omega_mod), subsample_size - p))
  } else if (type == "B") {
    omega_mod <- omega
  }

  list(
    model = model,
    sol = sol,
    status = status,
    omega = omega,
    omega_mod = omega_mod,
    xi_mat = xi_mat
  )
}

### find_center_repellent -------------------------

# TODO: document
find_center_repellent <- function(
  h,
  Y, X, D, Z, Phi = linear_projection(D, X, Z),
  tau,
  beta_D_proposal = NULL,
  beta_X_proposal = NULL,
  subsample_size,
  params = list(OutputFlag = 0),
  type = "C",
  gencontype = "power", # "power" or "log" or "max"
  a = ifelse(gencontype == "power", 0.5, exp(1)),
  simplex_repel = TRUE, # repel away from facets of the simplex
  foc_repel = TRUE # repel away from the FOC conditions
) {

  n <- nrow(Y)
  p <- length(h)

  if (!simplex_repel & !foc_repel) {
    stop("At least one of `foc_repel` and `simplex_repel` must be TRUE. Both can't be false.")
  }
  if (simplex_repel & ((n - p) - 1) <= 0) {
    stop("facet_center denominator must be positive i.e. we need n - p - 1 > 0")
  }

  xi_mat <- compute_xi_i(
    h = h,
    Y = Y, X = X, D = D, Z = Z, Phi = Phi,
    tau = tau,
    beta_D_proposal = beta_D_proposal,
    beta_X_proposal = beta_X_proposal
  )

  # Decision variables in order from left/top to right/bottom:
  # 1. omega --- (n - p) by 1
  # 2. right_slack --- p by 1
  # 3. left_slack --- p by 1
  # 4. right_slack_transformed --- p by 1
  # 5. left_slack_transformed --- p by 1
  # 6. simplex_slack --- 2n by 1
  # 7. simplex_slack_transformed --- 2n by 1
  # NOTE: when gencontype = "max", we have one additional decision variable
  # which is the max of the slack variables; we aren't considering this variable
  # in `num_decision_vars`
  num_omega <- n - p
  num_right_slack <- p
  num_left_slack <- p
  num_right_slack_transform <- num_right_slack
  num_left_slack_transform <- num_left_slack
  num_simplex_slack <- 2 * (n - p)
  num_simplex_slack_transform <- 2 * (n - p)
  num_decision_vars <- num_omega + num_right_slack + num_left_slack +
    num_right_slack_transform + num_left_slack_transform +
    num_simplex_slack + num_simplex_slack_transform

  model <- list()
  model$modelsense <- ifelse(gencontype == "max", "min", "max")
  model$lb <- c(
    rep(0, num_omega),
    rep(0, num_right_slack),
    rep(0, num_left_slack),
    rep(-Inf, num_right_slack_transform),
    rep(-Inf, num_right_slack_transform),
    rep(0, num_simplex_slack),
    rep(-Inf, num_simplex_slack_transform)
  )
  model$ub <- c(
    rep(1, num_omega),
    rep(Inf, num_right_slack + num_left_slack),
    rep(Inf, num_right_slack_transform + num_left_slack_transform),
    rep(Inf, num_simplex_slack),
    rep(Inf, num_simplex_slack_transform)
  )
  model$obj <- c(
    rep(0, num_omega + num_right_slack + num_left_slack),
    rep(1, num_right_slack_transform + num_left_slack_transform),
    rep(0, num_simplex_slack),
    rep(1, num_simplex_slack_transform)
  )
  model$vtype <- c(
    rep(type, num_omega),
    rep("C", num_decision_vars - num_omega)
  )

  # sum of omega
  omega_A <- matrix(c(
    rep(1, num_omega), rep(0, num_decision_vars - num_omega)
  ), nrow = 1)
  omega_sense <- "="
  omega_rhs <- subsample_size - p

  # FOC slack variables
  right_A <- vector("list", p)
  left_A <- vector("list", p)
  for (j in seq_len(p)) {
    xi_j <- xi_mat[j, ]
    zeros <- rep(0, p)
    ones <- zeros
    ones[j] <- 1
    right_A[[j]] <- c(xi_j, ones, zeros, zeros, zeros)
    left_A[[j]] <- c(-xi_j, zeros, ones, zeros, zeros)
  }
  right_A <- do.call(rbind, right_A)
  left_A <- do.call(rbind, left_A)
  foc_A <- rbind(right_A, left_A)
  foc_A <- expand_matrix(
    foc_A, newrow = nrow(foc_A), newcol = num_decision_vars,
    row_direction = "bottom", col_direction = "right"
  )
  foc_sense <- rep("=", 2 * p)
  foc_rhs <- c(rep(1 - tau, p), rep(tau, p))

  # transform FOC slack variables
  xstart <- num_omega
  ystart <- xstart + num_right_slack + num_left_slack
  foc_gencon <- vector("list", num_left_slack_transform + num_right_slack_transform)
  for (j in seq_len(2 * p)) {
    yvar <- ystart + j
    xvar <- xstart + j
    foc_gencon[[j]]$xvar <- xvar
    foc_gencon[[j]]$yvar <- yvar
    foc_gencon[[j]]$a <- a
  }

  # simplex slack variables
  simplex_rhs <- vector("double", num_simplex_slack)
  simplex_lhs <- vector("list", num_simplex_slack)
  counter <- 0
  subsample_center <- rep((subsample_size - p) / (n - p), n - p)
  for (j in c(0, 1)) {
    for (i in seq_len(n - p)) {
      counter <- counter + 1
      e_ij <- rep(0, num_simplex_slack)
      e_ij[[counter]] <- 1
      facet_center <- rep(
        ifelse(j == 0, (subsample_size - p) / ((n - p) - 1), ((subsample_size - p) - 1) / ((n - p) - 1)),
        n - p
      )
      facet_center[i] <- j
      normal_vec <- subsample_center - facet_center
      normal_unit <- normal_vec / sqrt(sum(normal_vec^2))
      simplex_rhs[[counter]] <- -sum(normal_unit * facet_center)
      # if (j == 0) {
      #   normal_A <- c(normal_unit, rep(0, n))
      # } else {
      #   normal_A <- c(rep(0, n), normal_unit)
      # }
      simplex_lhs[[counter]] <- c(
        -normal_unit, # num_omega
        rep(0, num_left_slack + num_right_slack +
            num_left_slack_transform + num_right_slack_transform),
        e_ij, # num_simplex_slack
        rep(0, num_simplex_slack_transform) # num_simplex_slack_transform
      )
    }
  }
  simplex_A <- do.call(rbind, simplex_lhs)
  stopifnot(ncol(simplex_A) == num_decision_vars)
  stopifnot(nrow(simplex_A) == num_simplex_slack)
  simplex_sense <- rep("=", num_simplex_slack)

  # transform simplex slack variables
  xstart <- num_omega + num_left_slack + num_right_slack +
    num_left_slack_transform + num_right_slack_transform
  ystart <- xstart + num_simplex_slack
  simplex_gencon <- vector("list", num_simplex_slack_transform)
  for (j in seq_len(num_simplex_slack_transform)) {
    yvar <- ystart + j
    xvar <- xstart + j
    simplex_gencon[[j]]$xvar <- xvar
    simplex_gencon[[j]]$yvar <- yvar
    simplex_gencon[[j]]$a <- a
  }

  # constraints

  if (gencontype == "max") {
    # add one more decision variable for the max slack variable
    omega_A <- cbind(omega_A, rep(0, nrow(omega_A)))
    foc_A <- cbind(foc_A, rep(0, nrow(foc_A)))
    simplex_A <- cbind(simplex_A, rep(0, nrow(simplex_A)))
    model$lb <- c(model$lb, 0)
    model$ub <- c(model$ub, Inf)
    model$vtype <- c(model$vtype, "C")


    # max >= slack => -slack + max >= 0
    foc_slack_tmp <- diag(-1, num_right_slack + num_left_slack)
    max_foc_slack <- cbind(
      matrix(0, nrow = nrow(foc_slack_tmp), ncol = num_omega),
      foc_slack_tmp,
      matrix(0, nrow = nrow(foc_slack_tmp), ncol = num_right_slack_transform +
        num_left_slack_transform + num_simplex_slack +
        num_simplex_slack_transform),
      rep(1, length = nrow(foc_slack_tmp))
    )
    max_foc_slack_sense <- rep(">=", length = nrow(foc_slack_tmp))
    max_foc_slack_rhs <- rep(0, length = nrow(foc_slack_tmp))

    simplex_slack_tmp <- diag(-1, num_simplex_slack_transform)
    max_simplex_slack <- cbind(
      matrix(0, nrow = nrow(simplex_slack_tmp), ncol = num_omega +
        num_right_slack + num_left_slack + num_right_slack_transform +
        num_left_slack_transform),
      simplex_slack_tmp,
      matrix(0, nrow = nrow(simplex_slack_tmp), ncol = num_simplex_slack_transform),
      rep(1, length = nrow(simplex_slack_tmp))
    )
    max_simplex_slack_sense <- rep(">=", length = nrow(simplex_slack_tmp))
    max_simplex_slack_rhs <- rep(0, length = nrow(simplex_slack_tmp))
    model$obj <- c(
      rep(0, num_decision_vars),
      1
    )
  } else {
    max_foc_slack <- c()
    max_foc_slack_sense <- c()
    max_foc_slack_rhs <- c()
    max_simplex_slack <- c()
    max_simplex_slack_sense <- c()
    max_simplex_slack_rhs <- c()
  }

  model$A <- rbind(omega_A, foc_A, simplex_A)
  model$sense <- c(omega_sense, foc_sense, simplex_sense)
  model$rhs <- c(omega_rhs, foc_rhs, simplex_rhs)
  gencon <- append(foc_gencon, simplex_gencon)

  if (!foc_repel) { # no foc, only simplex
    model$A <- rbind(omega_A, simplex_A, max_simplex_slack)
    model$sense <- c(omega_sense, simplex_sense, max_simplex_slack_sense)
    model$rhs <- c(omega_rhs, simplex_rhs, max_simplex_slack_rhs)
    gencon <- simplex_gencon
    model$obj <- c(
      rep(0, num_omega + num_right_slack + num_left_slack),
      rep(0, num_right_slack_transform + num_left_slack_transform),
      rep(0, num_simplex_slack),
      rep(1, num_simplex_slack_transform)
    )
    if (gencontype == "max") {
      model$obj <- c(
        rep(0, num_omega + num_right_slack + num_left_slack),
        rep(0, num_right_slack_transform + num_left_slack_transform),
        rep(0, num_simplex_slack),
        rep(0, num_simplex_slack_transform),
        1
      )
    }
  }
  if (!simplex_repel) { # no simplex, only foc
    model$A <- rbind(omega_A, foc_A, max_foc_slack)
    model$sense <- c(omega_sense, foc_sense, max_foc_slack_sense)
    model$rhs <- c(omega_rhs, foc_rhs, max_foc_slack_rhs)
    gencon <- foc_gencon
    model$obj <- c(
      rep(0, num_omega + num_right_slack + num_left_slack),
      rep(1, num_right_slack_transform + num_left_slack_transform),
      rep(0, num_simplex_slack),
      rep(0, num_simplex_slack_transform)
    )
    if (gencontype == "max") {
      model$obj <- c(
        rep(0, num_omega + num_right_slack + num_left_slack),
        rep(0, num_right_slack_transform + num_left_slack_transform),
        rep(0, num_simplex_slack),
        rep(0, num_simplex_slack_transform),
        1
      )
    }
  }

  if (gencontype == "log") {
    model$genconloga <- gencon
  } else if (gencontype == "power") {
    model$genconpow <- gencon
  }

  sol <- gurobi::gurobi(model, params)
  status <- sol$status
  omega <- sol$x[seq_len(num_omega)]

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

  # turn continuous solution into an integral one
  if (type == "C") {
    current_sum <- length(which(omega == 1)) # how many integral 1's do we have?
    remaining <- subsample_size - p - current_sum # how many need to be switched?
    to_be_rounded <- omega > 0 & omega < 1 # which can we switch?
    # NOTE: `rank` works better than `order` when there are ties
    max_indices <- rank(-omega, ties.method = "random") # 1 = largest
    # switch the largest numbers that aren't 1
    switch <- which(max_indices <= current_sum + remaining & max_indices > current_sum)
    omega_mod <- omega
    omega_mod[switch] <- 1
    omega_mod[which(omega_mod < 1)] <- 0
    omega_mod <- round(omega_mod, 0)
    stopifnot(all.equal(sum(omega_mod), subsample_size - p))
  } else if (type == "B") {
    omega_mod <- omega
  }

  rounded_center_indices <- setdiff(seq_len(n), h)[which(omega_mod == 1)]
  rounded_center <- vector("integer", n)
  rounded_center[rounded_center_indices] <- 1
  rounded_center[h] <- 1

  cts_center_indices <- setdiff(seq_len(n), h)[which(omega > 0)]
  cts_center <- vector("integer", n)
  cts_center[cts_center_indices] <- omega[which(omega > 0)]
  cts_center[h] <- 1

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
