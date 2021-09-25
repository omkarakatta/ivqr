test_that("Ensure that p_D >= 3", {
  expect_warning(
    chen_lee(p_D = 1),
    "p_D is less than 3"
  )
  expect_warning(
    chen_lee(p_D = 2),
    "p_D is less than 3"
  )
  expect_error(
    chen_lee(p_D = 3),
    NA
  )
  expect_error(
    chen_lee(p_D = 4),
    NA
  )
})

test_that("Ensure that n is a positive integer", {
  expect_error(
    chen_lee(n = -10),
    "n must be a positive integer."
  )
  expect_error(
    chen_lee(n = 1.2),
    "n must be a positive integer."
  )
  expect_error(
    chen_lee(n = 4),
    NA
  )
})
