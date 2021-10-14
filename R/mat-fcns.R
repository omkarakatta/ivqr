### expand_matrix ---------------------------

#' Expand a Matrix
#'
#' Expand an m-by-n matrix to be a larger p-by-q matrix
#'
#' This function is useful to write a simple matrix whose entries can themselves
#' be grouped as a matrix.
#'
#' @param mat An m-by-n matrix to embed into a larger p-by-q matrix
#' @param newrow The number of rows of the larger matrix (p)
#' @param newcol The number of columns of the larger matrix (c)
#' @param row_direction Adds rows to "top" or "bottom" of \code{mat}
#' @param col_direction Adds columns to "left" or "right" of \code{mat}
#' @param fill The value of entries in the new rows and columns; defaults to 0
#' @param quietly If TRUE, suppresses messages; defaults to TRUE
#'
#' @return a larger matrix that embeds \code{mat}
#'
#' @examples
#' \dontrun{
#' expand_matrix(matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE),
#'               newrow = 3, newcol = 3,
#'               row_direction = "bottom",
#'               col_direction = "right")
#' }
#'
#' @family matrix operations
#' @seealso \code{\link{block_diagonal}}
expand_matrix <- function(mat, newrow, newcol,
                          row_direction, col_direction,
                          fill = 0, quietly = TRUE) {
  current_row <- nrow(mat)
  current_col <- ncol(mat)
  stopifnot(current_row <= newrow)
  stopifnot(current_col <= newcol)

  addrow <- newrow - current_row
  addcol <- newcol - current_col

  rowvalid <- ifelse(addrow == 0, FALSE, TRUE)
  colvalid <- ifelse(addcol == 0, FALSE, TRUE)

  if (rowvalid & !quietly) message("Expanding by rows")
  if (colvalid & !quietly) message("Expanding by cols")

  if (rowvalid) {
    temp <- matrix(fill, nrow = addrow, ncol = ncol(mat))
    if (!(row_direction %in% c("bottom", "top"))) {
      stop("`row_direction` can only be 'top' or 'bottom'.")
    }
    if (row_direction == "bottom") {
      mat <- rbind(mat, temp)
    } else if (row_direction == "top") {
      mat <- rbind(temp, mat)
    }
  }
  if (colvalid) {
    temp <- matrix(fill, nrow = nrow(mat), ncol = addcol)
    if (!(col_direction %in% c("left", "right"))) {
      stop("`col_direction` can only be 'right' or 'left'.")
    }
    if (col_direction == "right") {
      mat <- cbind(mat, temp)
    } else if (col_direction == "left") {
      mat <- cbind(temp, mat)
    }
  }

  mat
}

### block_diagonal ---------------------------

#' Create a block diagonal matrix
#'
#' Combine a list of matrices into a block diagonal matrix
#'
#' The matrices provided in \code{mat_list} will be arranged as block matrices
#' with the first listed matrix at the top-right and the last matrix at the
#' bottom left.
#' The entries that are not in each block are filled with \code{fill}, which
#' defaults to 0.
#'
#' @param mat_list The list of matrices in the order they appear from top-right to
#'  bottom-left in the final matrix
#' @param fill The value of entries in the off-block entries; defaults to 0
#'
#' @return A block-diagonal matrix
#'
#' @examples
#' \dontrun{
#' block_diagonal(diag(1, 2), diag(1, 3)) #~ this is identical to `diag(1, 5)`
#' }
#' @family matrix operations
#' @seealso \code{\link{expand_matrix}}
block_diagonal <- function(mat_list, fill = 0) {
  cols <- sapply(mat_list, ncol)
  left_cols <- cumsum(cols)
  total_cols <- sum(cols)
  mat_expanded <- lapply(
    seq_along(mat_list),
    function(mat_index) {
      mat <- mat_list[[mat_index]]
      mat <- expand_matrix(mat,
                           newrow = nrow(mat),
                           newcol = left_cols[[mat_index]],
                           row_direction = "top",
                           col_direction = "left",
                           fill = fill)
      mat <- expand_matrix(mat,
                           newrow = nrow(mat),
                           newcol = total_cols,
                           row_direction = "bottom",
                           col_direction = "right",
                           fill = fill)
    }
  )
  do.call(rbind, mat_expanded)
}
