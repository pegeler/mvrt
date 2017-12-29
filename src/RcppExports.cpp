// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// mvrt2
arma::mat mvrt2(int n, arma::vec mu, arma::mat S, int df, double max_norm, int max_iterations);
RcppExport SEXP _mvrt_mvrt2(SEXP nSEXP, SEXP muSEXP, SEXP SSEXP, SEXP dfSEXP, SEXP max_normSEXP, SEXP max_iterationsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type S(SSEXP);
    Rcpp::traits::input_parameter< int >::type df(dfSEXP);
    Rcpp::traits::input_parameter< double >::type max_norm(max_normSEXP);
    Rcpp::traits::input_parameter< int >::type max_iterations(max_iterationsSEXP);
    rcpp_result_gen = Rcpp::wrap(mvrt2(n, mu, S, df, max_norm, max_iterations));
    return rcpp_result_gen;
END_RCPP
}
// mvrt
arma::mat mvrt(int n, arma::vec mu, arma::mat S, int df);
RcppExport SEXP _mvrt_mvrt(SEXP nSEXP, SEXP muSEXP, SEXP SSEXP, SEXP dfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type S(SSEXP);
    Rcpp::traits::input_parameter< int >::type df(dfSEXP);
    rcpp_result_gen = Rcpp::wrap(mvrt(n, mu, S, df));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_mvrt_mvrt2", (DL_FUNC) &_mvrt_mvrt2, 6},
    {"_mvrt_mvrt", (DL_FUNC) &_mvrt_mvrt, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_mvrt(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
