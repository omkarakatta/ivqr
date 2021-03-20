msg1 <- "message 1"
msg2 <- "message 2"

list0 <- list()
list1 <- list(msg1)
list2 <- list(msg1, msg2)
list2_newline <- list(msg1, "\nmessage 2")

test_that("conditionally collect note with an empty list", {
  expect_equal(
    collect_note_if(list0, msg1, condition = TRUE, newline = FALSE),
    list1
  )
  expect_equal(
    collect_note_if(list1, msg2, condition = FALSE),
    list1
  )
  expect_equal(
    collect_note_if(list1, msg2, condition = TRUE, newline = FALSE),
    list2
  )
  expect_equal(
    collect_note_if(list1, msg2, condition = TRUE, newline = TRUE),
    list2_newline
  )
})
