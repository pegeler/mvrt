#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//' @rdname mvrt
//' @export
// [[Rcpp::export]]
arma::mat mvrt(int n, arma::vec mu, arma::mat S, int df=1)
{

  arma::mat x(mu.size(), n);
  arma::mat g = chol(S).t();

  for (int i=0; i < n; i++)
  {
    arma::vec rand = as<arma::vec>(rt(mu.size(),df));
    x.col(i) = mu + g * rand;
  }

  return x.t();
}
