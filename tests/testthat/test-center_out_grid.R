test_that("create multidimensional grid", {
  expect_equal(
    center_out_grid(c(1, 2), c(0.1, 1), c(1, 2)),
    expand.grid(Var1 = c(0.9, 1.0, 1.1), Var2 = c(0, 1, 2, 3, 4)),
    ignore_attr = TRUE
  )
})

test_that("named grid", {
  expect_equal(
    center_out_grid(c(a = 1, b = 2), c(0.1, 1), c(1, 2)),
    expand.grid(a = c(0.9, 1.0, 1.1), b = c(0, 1, 2, 3, 4)),
    ignore_attr = TRUE
  )
})

test_that("same length between `center`, `increment`, and `length`", {
  expect_error(
    center_out_grid(c(a = 1, b = 2), c(0.1), c(1, 2)),
  )
  expect_error(
    center_out_grid(c(a = 1, b = 2), c(0.1, 0.2), c(1)),
  )
  expect_error(
    center_out_grid(c(a = 1), c(0.1, 1), c(1, 2)),
  )
  expect_error(
    center_out_grid(c(a = 1), list(0.1, 1), c(1, 2)),
  )
})

test_that("asymmetric grid", {
  expect_equal(
    center_out_grid(c(a = 1, b = 2), list(c(1, 2), 1), list(1, 2)),
    expand.grid(a = c(0, 1, 3), b = c(0, 1, 2, 3, 4)),
    ignore_attr = TRUE
  )
})

test_that("error: long increment/length values for one or more coefficients", {
  expect_error(
    center_out_grid(c(a = 1, b = 2), list(c(1, 2, 3), 1), list(1, 2)),
    "index 1"
  )
  expect_error(
    center_out_grid(c(a = 1, b = 2), list(c(1, 2), 1), list(2, c(1, 2, 3))),
    "index 2"
  )
})
