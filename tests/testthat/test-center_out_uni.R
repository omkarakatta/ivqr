test_that("univariate grid from center out", {
  expect_equal(
    center_out_uni(5, 1, 2),
    c(3, 4, 5, 6, 7)
  )
  expect_equal(
    center_out_uni(5, 0.1, 2),
    c(4.8, 4.9, 5.0, 5.1, 5.2)
  )
  expect_equal(
    center_out_uni(0.4, 3, 1),
    c(-2.6, 0.4, 3.4)
  )
})

test_that("asymmetric grid from center out", {
  expect_equal(
    center_out_uni(5, c(1, 2), c(2, 1)),
    c(3, 4, 5, 7)
  )
})

test_that("error: greater than 2 elements in `increment` or `length`", {
  expect_error(
    center_out_uni(5, c(1, 2, 3), 1),
    "`increment`"
  )
  expect_error(
    center_out_uni(5, 1, c(1, 2, 3)),
    "`length`"
  )
})
