#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//' @rdname mvrt
//' @export
// [[Rcpp::export]]
arma::mat mvrt(int n, arma::vec mu, arma::mat S, int df=1)
{

  // Cholesky decomp and transpose covariance matrix
  arma::mat g = chol(S).t();

  // Generate the random data
  auto x_vec = rt(mu.size() * n, df);
  arma::mat x = arma::mat(&x_vec[0], mu.size(), n, false, true);

  // Give the random data covariance structure and add mean offset
  x = g * x;
  x.each_col() += mu;

  return x.t();
}

