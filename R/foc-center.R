# We have several notions of the "center of FOC polytope".
# `foc_center_method` is a wrapper for these notions, most of which are in
# `R/foc-center-archive.R`.
# `foc_center` is the prevailing notion that we use.

# foc_center_method ------------------------------------------------------------

foc_center_method <- function(FUN, ...) {
  FUN(...)
}
