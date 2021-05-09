#' @title synthex1 - Synthetic Example 1
#' @description 500 observations and 3 endogeneous variables
#' @format named list of 5 matrices:
#' \describe{
#'   \item{\code{Y}}{outcome variable, 500 by 1 matrix}
#'   \item{\code{D}}{endogeneous variables, 500 by 3 matrix}
#'   \item{\code{Z}}{instruments, 500 by 3 matrix}
#'   \item{\code{X}}{intercept, 500 by 1 matrix of 1's}
#'   \item{\code{errors}}{model error and shocks defining D, 500 by 4 matrix}
#' }
"synthex1"
