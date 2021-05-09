### Meta -------------------------
###
### Title: Simulate data for testing {ivqr}
###
### Description: `synthex1` is a dataset with 500 observations and 3 endogeneous
### variables. It is formatted as a list. See `help(synthex1)` for more details.
###
### Author: Cheryl Liu (modified by Omkar A. Katta)
### Original: See ChenyueLiu/ivqr/tests/testthat/test20190828.R
###

### Preliminaries -------------------------

library(MASS)
set.seed(2)
n <- 500

### Create instruments -------------------------

Z1 <- rnorm(n)
Z2 <- rnorm(n)
Z3 <- rnorm(n)
Z <- cbind(Z1, Z2, Z3)

### Create shocks and errors -------------------------

Sigma <- matrix(c(1, 0.4, 0.6, -0.2,
                  0.4, 1, 0, 0,
                  0.6, 0, 1, 0,
                  -0.2, 0, 0, 1),
                nrow = 4,
                ncol = 4,
                byrow = TRUE)
EpsV1V2V3 <- mvrnorm(n, mu = c(0, 0, 0, 0), Sigma = 0.25 * Sigma)

### Create endogeneous variables -------------------------

D1D2D3 <- matrix(NA, n, 3)
for (i in seq_len(n)) {
  D1D2D3[i, 1] <- pnorm(Z1[i] + EpsV1V2V3[i, 2])
  D1D2D3[i, 2] <- 2 * pnorm(Z2[i] + EpsV1V2V3[i, 3])
  D1D2D3[i, 3] <- 1.5 * pnorm(Z3[i] + EpsV1V2V3[i, 4])
  omitted_variable <- pnorm(Z1[i] + EpsV1V2V3[i, 2])
}

### Create intercept -------------------------

X <- matrix(1, n, 1)

D <- D1D2D3

Y <- X + D[, 1] + D[, 2] + D[, 3] + EpsV1V2V3[, 1] + omitted_variable
Y <- matrix(Y, n, 1)

### Compile and use data -------------------------

synthex1 <- list("Y" = Y, "D" = D, "Z" = Z, "X" = X, "errors" = EpsV1V2V3)

usethis::use_data(synthex1, overwrite = TRUE)
