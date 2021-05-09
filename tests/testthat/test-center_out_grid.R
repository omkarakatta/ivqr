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

test_that("same length", {
  expect_error(
    center_out_grid(c(a = 1, b = 2), c(0.1), c(1, 2)),
  )
  expect_error(
    center_out_grid(c(a = 1, b = 2), c(0.1, 0.2), c(1)),
  )
  expect_error(
    center_out_grid(c(a = 1), c(0.1, 1), c(1, 2)),
  )
})
