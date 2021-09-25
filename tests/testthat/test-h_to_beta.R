get_beta_h <- function(n, tau) {
  # create simulation
  sim <- chen_lee(n = n)
  Y <- sim$Y
  D <- sim$D
  X <- sim$X
  Z <- sim$Z
  # compute answer
  result <- iqr_milp(
    Y = Y,
    D = D,
    X = X,
    Z = Z,
    tau = tau,
    M = 10
  )
  # collect answer
  beta_D <- result$beta_D
  beta_Phi <- result$beta_Phi
  beta_X <- result$beta_X
  # use h_to_beta
  Y_tilde <- Y - D %*% beta_D
  qr <- quantreg::rq(Y_tilde ~ X + linear_projection(D, X, Z) - 1, tau = tau)
  h <- which(qr$dual > 0 & qr$dual < 1)
  beta_h <- h_to_beta(h, Y = Y, X = X, D = D, Z = Z)
  list(
    beta_h = beta_h,
    beta_D = beta_D,
    beta_Phi = beta_Phi,
    beta_X = beta_X
  )
}

test_that("Chen Lee, tau = 0.5", {
  result <- get_beta_h(n = 100, tau = 0.5)
  expect_equal(
    result$beta_h$beta_D,
    result$beta_D
  )
  expect_equal(
    result$beta_h$beta_Phi,
    result$beta_Phi
  )
  expect_equal(
    result$beta_h$beta_X,
    result$beta_X
  )
})

test_that("Chen Lee, tau = 0.6", {
  result <- get_beta_h(n = 100, tau = 0.6)
  expect_equal(
    result$beta_h$beta_D,
    result$beta_D
  )
  expect_equal(
    result$beta_h$beta_Phi,
    result$beta_Phi
  )
  expect_equal(
    result$beta_h$beta_X,
    result$beta_X
  )
})

test_that("Chen Lee, tau = 0.25", {
  result <- get_beta_h(n = 100, tau = 0.25)
  expect_equal(
    result$beta_h$beta_D,
    result$beta_D
  )
  expect_equal(
    result$beta_h$beta_Phi,
    result$beta_Phi
  )
  expect_equal(
    result$beta_h$beta_X,
    result$beta_X
  )
})
