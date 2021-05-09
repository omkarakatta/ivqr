test_that("Univariate grid from center out", {
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
