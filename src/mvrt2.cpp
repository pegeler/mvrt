#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//' @rdname mvrt
//' @export
// [[Rcpp::export]]
arma::mat mvrt2(
    int n,
    arma::vec mu,
    arma::mat S,
    int df=1,
    double max_norm=2,
    int max_iterations=1000
  )
{

  arma::mat x(mu.size(), n);
  arma::mat g = chol(S).t();
  arma::mat V_sqrt_inv = diagmat(1 / sqrt(S.diag()));
  arma::mat target_cor = V_sqrt_inv * S * V_sqrt_inv;

  for (int iteration=0; iteration < max_iterations; iteration++)
  {

    // Generate the random data
    for (int i=0; i < n; i++)
    {
      arma::vec rand = as<arma::vec>(rt(mu.size(),df));
      x.col(i) = mu + g * rand;
    }

    // Compare to target
    arma::mat diff_matrix = abs(cor(x.t()) - target_cor);

    if (diff_matrix.max() <= max_norm) return x.t();

  }

  stop(
    "Did not generate matrix with max norm of %f in %i iterations",
    max_norm,
    max_iterations
  );

}
