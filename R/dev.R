### Meta -------------------------
###
### Title: Developper Functions
###
### Description: Auxiliary functions that aim to make package development easier
###
### Author: Omkar A. Katta
###

### send_note_if -------------------------

#' Conditionally send note to console
#'
#' Send a message or warning if \code{condition} is TRUE
#'
#' The intended values of \code{fcn} are \code{\link[base]{message}} and
#' \code{\link[base]{warning}}.
#' If \code{condition} is TRUE, then send the message to the console according
#' to \code{fcn}.
#' Otherwise, do not send anything to the console.
#'
#' @param note A string that will be sent to the console
#' @param condition If TRUE, message will be sent to console.
#' @param fcn A function, e.g., \code{\link[base]{message}}
#'  or \code{\link[base]{warning}}
#' @param ... Additional arguments for \code{fcn}
#'
#' @return A note in the console if \code{quietly} is FALSE
#'
#' @examples
#' \dontrun{
#' send_note_if("This is a warning.", TRUE, warning)
#' }
#'
#' @family dev
send_note_if <- function(note, condition, fcn, ...) {
  if (condition) fcn(note, ...)
}

### collect_note_if ---------------------------

#' Conditionally add note to list
#'
#' Add an entry to a list if \code{condition} is TRUE
#'
#' Before this function, error messages would be sent one at a time.
#' Resultantly, users would spend time in a cycle of: run code, debug code,
#' and repeat. However, it would be more efficient if the user received all
#' possible error messages, rather than just the first error message, and then
#' start debugging. This would ease the debugging process and perhaps later
#' errors would provide insight into fixing the current error.
#'
#' This function was created with the intention of collecting error messages
#' in a list, and then printing all the error messages at the end.
#' The hope was that users can debug all at one go, rather than debug
#' incrementally.
#'
#' The goal is to pass this list to \code{stop}, \code{warning}, etc.
#' To print each message on its own line, set \code{newline} to TRUE.
#'
#' \code{msg_list} must be a list; \code{new_msg} must be a string; and
#' \code{condition} and \code{newline} must be booleans.
#'
#' @param msg_list A list of messages
#' @param new_msg A new message to add to \code{msg_list}
#' @param condition If TRUE, then add \code{new_msg} to \code{msg_list}; else,
#'  \code{msg_list} remains unchanged
#' @param newline If TRUE, then add a new line to \code{new_msg} before
#'  appending to \code{msg_list}
#'
#' @return If \code{add} is TRUE, then the result will be all the messages in
#'  \code{msg_list}. Else, the function returns \code{msg_list}.
#'
#' @family dev
collect_note_if <- function(msg_list, new_msg, condition, newline = TRUE) {
  if (condition) {
    if (newline) new_msg <- paste0("\n", new_msg)
    msg_list[[length(msg_list) + 1]] <- new_msg
  }
  msg_list
}
