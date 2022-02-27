#' Simulated data for GOFAR
#'
#' Simulated data with low-rank and sparse coefficient matrix.
#'
#' @usage data(simulate_gofar)
#'
#' @format A list of variables for the analysis  using GOFAR(S) and GOFAR(P):
#' \describe{
#'   \item{Y}{Generated response matrix}
#'   \item{X}{Generated predictor matrix}
#'   \item{U}{specified value of U}
#'   \item{V}{specified value of V}
#'   \item{D}{specified value of D}
#'   \item{n}{sample size}
#'   \item{Xsigma}{covariance matrix used to generate predictors in X}
#'   \item{C0}{intercept value in the coefficient matrix}
#'   \item{familygroup}{index set of the type of multivariate outcomes: "1" for Gaussian, "2" for Bernoulli, "3" for Poisson outcomes}
#' Mishra, Aditya, Dipak K. Dey, Yong Chen, and Kun Chen. Generalized co-sparse factor regression. Computational Statistics & Data Analysis 157 (2021): 107127}
"simulate_gofar"
#' @example
#' data(simulate_gofar)

