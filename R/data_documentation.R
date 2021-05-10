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

#' @title Education-relevant Data
#' @description 24684 observations with known education status
#' \describe{
#'    \item{\code{Y_educ}}{outcome variable, column matrix}
#'    \item{\code{X_educ}}{covariate matrix, 50-column matrix}
#'    \item{\code{D_educ}}{endogeneous matrix, column matrix (emp)}
#'    \item{\code{Z_educ}}{instrument matrix, column matrix (emp_pct)}
#'    \item{\code{EducFE_educ}}{education fixed effects matrix, 2-column matrix (hsgrad, hsplus)}
#'    \item{\code{DInt_educ}}{endogeneous x education matrix, 2-column matrix}
#'    \item{\code{ZInt_educ}}{instrument x education matrix, 2-column matrix}
#' }
#' @name education_autor
"Y_educ"

#' @rdname education_autor
"X_educ"

#' @rdname education_autor
"D_educ"

#' @rdname education_autor
"Z_educ"

#' @rdname education_autor
"EducFE_educ"

#' @rdname education_autor
"DInt_educ"

#' @rdname education_autor
"ZInt_educ"

