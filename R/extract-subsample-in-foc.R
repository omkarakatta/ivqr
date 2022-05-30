extract_subsample_in_foc <- function(
  h,
  beta_D,
  beta_X,
  subsample_size,
  Y, X, D, Phi,
  tau,
  gurobi_params = list(OutputFlag = 0)
) {
  n <- nrow(Y)
  p <- length(h)
  tmp <- compute_foc_conditions(
    h, beta_D, beta_X,
    Y = Y, X = X, D = D, Phi = Phi, tau = tau
  )
  xi_mat <- tmp[, setdiff(seq_len(n), h), drop = FALSE]

  model <- list()
  model$obj <- xi_mat[1, ]
  model$A <- rbind(
    xi_mat,
    -xi_mat,
    rep(1, n - p)
  )
  model$rhs <- c(
    rep(-tau, p),
    rep(tau - 1, p),
    subsample_size - p
  )
  model$sense <- c(
    rep(">=", p),
    rep(">=", p),
    "="
  )
  model$vtype <- rep("B", n - p)
  model$modelsense <- "min"

  result <- gurobi::gurobi(model, gurobi_params)
  if (result$status != "OPTIMAL") {
    warning("extract_subsample_in_foc() -> not optimal program")
    return(result)
  }

  subsample <- vector("double", n)
  subsample[setdiff(seq_len(n), h)] <- result$x
  subsample[h] <- 1
  membership <- foc_membership(
    h = which(which(subsample == 1) %in% h),
    Y_subsample = Y[subsample == 1, , drop = FALSE],
    X_subsample = X[subsample == 1, , drop = FALSE],
    D_subsample = D[subsample == 1, , drop = FALSE],
    Phi_subsample = Phi[subsample == 1, , drop = FALSE],
    tau = tau,
    beta_D = beta_D
  )
  if (!membership$status) {
    warning("subsample not in polytope")
  }
  list(
    model = model,
    result = result,
    subsample = subsample,
    membership = membership
  )
}
