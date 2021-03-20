test_that("send note to console conditionally", {
  expect_message(send_note_if("I'm a message", 2 + 2 == 4, message))
  expect_message(send_note_if("I'm a hidden message", 2 + 2 != 4, message), NA)
  expect_warning(send_note_if("I'm a warning", 2 + 2 == 4, warning))
  expect_warning(send_note_if("I'm a hidden warning", 2 + 2 != 4, warning), NA)
})
