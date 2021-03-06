% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dev.R
\name{collect_note_if}
\alias{collect_note_if}
\title{Conditionally add note to list}
\usage{
collect_note_if(msg_list, new_msg, condition, newline = TRUE)
}
\arguments{
\item{msg_list}{A list of messages}

\item{new_msg}{A new message to add to \code{msg_list}}

\item{condition}{If TRUE, then add \code{new_msg} to \code{msg_list}; else,
\code{msg_list} remains unchanged}

\item{newline}{If TRUE, then add a new line to \code{new_msg} before
appending to \code{msg_list}}
}
\value{
If \code{add} is TRUE, then the result will be all the messages in
\code{msg_list}. Else, the function returns \code{msg_list}.
}
\description{
Add an entry to a list if \code{condition} is TRUE
}
\details{
Before this function, error messages would be sent one at a time.
Resultantly, users would spend time in a cycle of: run code, debug code,
and repeat. However, it would be more efficient if the user received all
possible error messages, rather than just the first error message, and then
start debugging. This would ease the debugging process and perhaps later
errors would provide insight into fixing the current error.

This function was created with the intention of collecting error messages
in a list, and then printing all the error messages at the end.
The hope was that users can debug all at one go, rather than debug
incrementally.

The goal is to pass this list to \code{stop}, \code{warning}, etc.
To print each message on its own line, set \code{newline} to TRUE.

\code{msg_list} must be a list; \code{new_msg} must be a string; and
\code{condition} and \code{newline} must be booleans.
}
\seealso{
Other dev: 
\code{\link{send_note_if}()}
}
\concept{dev}
