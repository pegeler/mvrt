#' @useDynLib mvrt
#' @importFrom Rcpp sourceCpp
NULL

#' Convert correlation matrix to covariance matrix
#'
#' @param R Correlation matrix
#' @param var Vector of variances
#' @export
convert_R2S <- function(R, var) {

  stopifnot(identical(
    dim(V_sqrt <- sqrt(diag(var))),
    dim(R)
  ))

  V_sqrt %*% R %*% V_sqrt
}

#' Populate a correlation matrix
#'
#' @param r Pearson's correlation to populate off-diagonals
#' @param dim Matrix dimensions
#' @export
make_cor_mat <- function (r, dim = 2) {
  out <- matrix(r, dim, dim)
  diag(out) <- 1
  out
}

#' Create a multivariate t-distributed random matrix
#'
#' @param n Number of observations
#' @param mu Vector of sample means
#' @param S Covariance matrix
#' @param df Degrees of freedom
#' @param max.norm Maximum matrix norm
#' @param type.norm Type of matrix norm passed to \code{\link[base]{norm}}
#' @importFrom stats cor cov2cor rt
#' @export
mvrtR <- function(n, mu, S, df = n - 1, max.norm = Inf, type.norm = "m") {

  g.t <- t(chol(S))

  accept_data <- FALSE

  R_ref <- cov2cor(S)

  while (!accept_data) {

    random_matrix <- matrix(rt(n*length(mu), df),nrow = length(mu))
    deviation <- g.t %*% random_matrix

    out <- t(mu + deviation)

    accept_data <- norm(R_ref - cor(out), type = type.norm) <= max.norm

  }

  out

}

