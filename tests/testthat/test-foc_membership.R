# FOC Membership
library(dplyr)
library(quantreg)

test_that("Autor", {

  # filter to black females
  full_data <- as.data.frame(cbind(Y_educ, X_educ, D_educ, Z_educ,
                             EducFE_educ, DInt_educ, ZInt_educ)) %>%
    filter(female > 0) %>%
    filter(wh_race < 0) %>%
    filter(oth_race < 0)

  # create data
  y_mat <- as.matrix(full_data[, y_name])
  colnames(y_mat) <- y_name
  covs_remove <- grepl("Idistyr|yq|wageadjh1_8|qtremph1_8|race|female",
                       covs_name)
  covs_keep <- covs_name[!covs_remove]
  x_mat <- cbind(ones = rep(1, nrow(full_data)),
                 as.matrix(full_data[, educfe_name]),
                 as.matrix(full_data[, covs_keep]))
  z_mat <- as.matrix(full_data[, c(z_name, zint_name)])
  d_mat <- as.matrix(full_data[, c(d_name, dint_name)])
  phi_mat <- linear_projection(d_mat, x_mat, z_mat)
  new_data <- cbind(y_mat, x_mat, z_mat, d_mat, phi_mat)

  # create subsample
  n <- nrow(new_data)
  m <- 100
  set.seed(5)
  sampled_rows <- sample(x = n, size = m, replace = FALSE)
  y_subsample <- y_mat[sampled_rows, , drop = FALSE]
  x_subsample <- x_mat[sampled_rows, , drop = FALSE]
  z_subsample <- z_mat[sampled_rows, , drop = FALSE]
  d_subsample <- d_mat[sampled_rows, , drop = FALSE]
  phi_subsample <- phi_mat[sampled_rows, , drop = FALSE]
  tau <- 0.25

  # solve subsample IQR
  result <- preprocess_iqr_milp(
    Y = y_subsample,
    X = x_subsample,
    Z = z_subsample,
    D = d_subsample,
    Phi = phi_subsample,
    tau = tau,
    show_progress = FALSE,
    show_iterations = FALSE,
    quietly = FALSE
  )

  # compute dual *exactly*
  y_tilde <- y_subsample - d_subsample %*% result$final_fit$beta_D
  reg <- quantreg::rq(y_tilde ~ x_subsample + phi_subsample - 1, tau = tau)
  dual <- reg$dual
  # compute active basis
  h <- which(dual > 0 & dual < 1)

  # does the subsample solve the FOC given the active basis h?
  result <- foc_membership(
    h = h,
    Y_subsample = y_subsample,
    X_subsample = x_subsample,
    D_subsample = d_subsample,
    Phi_subsample = phi_subsample,
    tau = tau
  )
  expect_true(result)
})

test_that("Chen-Lee", {

  # create data
  n <- 10000
  sim <- chen_lee(n = n)
  Y <- sim$Y
  D <- sim$D
  X <- sim$X
  Z <- sim$Z
  Phi <- linear_projection(D, X, Z)

  # create subsample
  m <- 100
  set.seed(2)
  sampled_rows <- sort(sample(x = n, size = m, replace = FALSE))
  Y_subsample <- Y[sampled_rows, , drop = FALSE]
  X_subsample <- X[sampled_rows, , drop = FALSE]
  Z_subsample <- Z[sampled_rows, , drop = FALSE]
  D_subsample <- D[sampled_rows, , drop = FALSE]
  Phi_subsample <- Phi[sampled_rows, , drop = FALSE]
  tau <- 0.25

  # solve IQR in subsample
  result <- preprocess_iqr_milp(
    Y = Y_subsample,
    X = X_subsample,
    Z = Z_subsample,
    D = D_subsample,
    Phi = Phi_subsample,
    tau = tau,
    show_progress = FALSE,
    show_iterations = FALSE,
    quietly = FALSE
  )

  # compute dual exactly
  Y_tilde <- Y_subsample - D_subsample %*% result$final_fit$beta_D
  reg <- quantreg::rq(Y_tilde ~ X_subsample + Phi_subsample - 1, tau = tau)
  dual <- reg$dual
  # compute active basis
  h <- which(dual > 0 & dual < 1)

  # does the subsample solve the FOC given the active basis h?
  result <- foc_membership(
    h = h,
    Y_subsample = Y_subsample,
    X_subsample = X_subsample,
    D_subsample = D_subsample,
    Phi_subsample = Phi_subsample,
    tau = tau
  )

  expect_true(result)
})
