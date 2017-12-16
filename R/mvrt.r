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
#' These functions will create a multivariate t-distributed random matrix of a
#' specified mean and covariance structure.
#'
#' The \code{mvrt2} fuction differs from \code{mvrt} in that it will ensure that the
#' resulting random data has a correlation matrix within a specified criterion.
#' No element of the resulting correlation matrix will deviate from the
#' corresponding element of the correlation matrix generated from
#' specifications by more than the absolute value of \code{max_norm}. This
#' checking comes at the expense of computational speed.
#'
#' @param n Number of observations.
#' @param mu Vector of sample means.
#' @param S Covariance matrix.
#' @param df Degrees of freedom.
#' @param max_norm Maximum matrix norm in difference between target correlation
#'   matrix and the correlation matrix resulting from randomly generated data.
#'   The greatest possible norm is 2. Use care with extemely small values, as
#'   several iterations may be required to obtain such a matrix.
#' @param max_iterations Maximum number of iterations for the loop to search
#'   for a matrix meeting \code{max_norm} criterion.
#'
#' @examples
#' S <- convert_R2S(make_cor_mat(.9),2:3)
#' x <- mvrt( 30, 4:5, S, 29)
#' y <- mvrt2(30, 4:5, S, 29, .01)
#'
#' cor(x)
#' cor(y)
#' @name mvrt
NULL

#' Create a multivariate t-distributed random matrix
#'
#' @param n Number of observations
#' @param mu Vector of sample means
#' @param S Covariance matrix
#' @param df Degrees of freedom
#' @param max_norm Maximum matrix norm
#' @param max_iterations Maximum number of iterations for the loop to search
#'   for a matrix meeting \code{max_norm} criterion.
#' @param type_norm Type of matrix norm passed to \code{\link[base]{norm}}
#' @importFrom stats cor cov2cor rt
#' @export
mvrtR <- function(
  n,
  mu,
  S,
  df = n - 1,
  max_norm = 2,
  max_iterations = 1000,
  type_norm = "m")
{

  .Deprecated("mvrt2")

  g.t <- t(chol(S))

  R_ref <- cov2cor(S)

  for (iterations in seq_len(max_iterations)) {

    random_matrix <- matrix(rt(n*length(mu), df),nrow = length(mu))
    deviation <- g.t %*% random_matrix

    if (norm(R_ref - cor(t(mu + deviation)), type = type_norm) <= max_norm)
      return(t(mu + deviation))

  }

  stop(
    "Correlation structure with max norm of ", max_norm,
    " was not obtained in ", max_iterations, " iterations"
    )

}
