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
### decompose_A ---------------------------
#' Decompose constraint matrix
#'
#' Print what each row of the A matrix represents
#'
#' Using the dimensions of the inputs, this function prints the constraints
#' associated with each row.
#'
#' @param Y Dependent variable (vector of length n)
#' @param X Exogenous variable (including constant vector) (n by p_X matrix)
#' @param D Endogenous variable (n by p_D matrix)
#' @param Z Instrumental variable (n by p_Z matrix)
#' @param O_neg,O_pos Indices for residuals whose sign is fixed to be negative
#'  and positive, respectively (vectors)
#' @param projection If TRUE (default), project D on the space spanned by X and
#'  Z to construct the vector of functions of transformed instruments; else,
#'  let Z be the instruments for endogenous variables (boolean)
#'
#' @return A vector whose length is the same as the number of constraints;
#'  \enumerate{
#'    \item pf: primal feasibility
#'    \item df_X: dual feasibility, X
#'    \item df_Phi: dual feasibility, Phi
#'    \item cs_uk: complementary slackness, u and k
#'    \item cs_vl: complementary slackness, v and l
#'    \item cs_ak: complementary slackness, a and k
#'    \item cs_al: complementary slackness, a and l
#'    \item pp_a: preprocessing, a
#'    \item pp_k: preprocessing, k
#'    \item pp_l: preprocessing, l
#'  }
decompose_A <- function(Y,
                        X,
                        D,
                        Z,
                        O_neg = NULL,
                        O_pos = NULL,
                        projection = TRUE) {
  # Get dimensions of data
  n <- length(Y)
  n_D <- nrow(D)
  n_X <- nrow(X)
  n_Z <- nrow(Z)
  p_D <- ncol(D)
  p_X <- ncol(X)
  p_Z <- ncol(Z)

  if (projection) {
    # Obtain fitted values from projecting D on space spanned by X and Z
    XZ <- cbind(X, Z)
    proj_matrix <- solve(t(XZ) %*% XZ) %*% t(XZ) %*% D
    Phi <- XZ %*% proj_matrix
  } else {
    Phi <- Z
  }
  p_Phi <- ncol(Phi)

  # Primal Feasibility
  pf <- rep("pf", n)

  # Dual Feasibility
  df_X <- rep("df_X", p_X)
  df_Phi <- rep("df_Phi", p_Phi)

  # Complementary Slackness
  cs_uk <- rep("cs_uk", n)
  cs_vl <- rep("cs_vl", n)
  cs_ak <- rep("cs_ak", n)
  cs_al <- rep("cs_al", n)

  # Preprocessing
  O_neg <- sort(O_neg)
  O_pos <- sort(O_pos)
  O <- c(O_neg, O_pos)        # indices of fixed residuals
  if (!is.null(O)) {
    pp_a <- rep("pp_a", n)
    pp_k <- rep("pp_k", n)
    pp_l <- rep("pp_l", n)
  } else {
    pp_a <- rep("pp_a", 0)
    pp_k <- rep("pp_k", 0)
    pp_l <- rep("pp_l", 0)
  }

  c(pf, df_X, df_Phi, cs_uk, cs_vl, cs_ak, cs_al, pp_a, pp_k, pp_l)

}

### write_prm ---------------------------

#' Save Gurobi Parameters
#'
#' Create a .prm file with parameters from Gurobi
#'
#' While \code{gurobi::gurobi_write()} can save Gurobi models, it cannot save
#' the parameters. The function \code{write_prm} overcomes this limitation
#' by creating a \code{.prm} file with the parameters.
#'
#' @param params A list of Gurobi parameters
#' @param path The directory to store the resulting \code{.prm} file
#' @param filename The name of the resulting \code{.prm} file
#'
#' @return A \code{.prm} file named \code{filename} in the directory \code{path}
#'  that stores \code{params}
#'
#' @examples
#' \dontrun{
#' params <- list(NonConvex = 0)
#' write_prm(params, ".", "example")
#' file.remove("./example.prm")
#' }
#'
#' @family debugging
#' @seealso \code{\link[gurobi]{gurobi_write}}
write_prm <- function(params, path, filename) {
  names <- names(params)
  values <- sapply(params, function(x) x[1])
  lines <- paste(names, values)
  prm <- paste(lines, collapse = "\n")
  target <- paste0(path, "/", filename, ".prm")
  cat(prm, file = target)
}
