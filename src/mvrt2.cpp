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

  // Cholesky decomp and transpose covariance matrix
  arma::mat g = chol(S).t();

  // Get correlation matrix of user-input S matrix
  arma::mat V_sqrt_inv = diagmat(1 / sqrt(S.diag()));
  arma::mat target_cor = V_sqrt_inv * S * V_sqrt_inv;

  for (int iteration=0; iteration < max_iterations; iteration++)
  {

    // Generate the random data
    arma::vec x_vec = as< arma::vec >( rt(mu.size() * n, df) );
    arma::mat x = arma::mat( (const double*)x_vec.begin(), mu.size(), n );

    // Give the random data covariance structure and add mean offset
    x = g * x;
    x.each_col() += mu;

    // Compare to target and retrun if meets specification
    arma::mat diff_matrix = abs(cor(x.t()) - target_cor);

    if (diff_matrix.max() <= max_norm) return x.t();

  }

  stop(
    "Did not generate matrix with max norm of %f in %i iterations",
    max_norm,
    max_iterations
  );

}
