// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// row_replicate
arma::mat row_replicate(arma::vec weight, int rn);
RcppExport SEXP _scstruc_row_replicate(SEXP weightSEXP, SEXP rnSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type weight(weightSEXP);
    Rcpp::traits::input_parameter< int >::type rn(rnSEXP);
    rcpp_result_gen = Rcpp::wrap(row_replicate(weight, rn));
    return rcpp_result_gen;
END_RCPP
}
// col_replicate
arma::mat col_replicate(arma::vec weight, int rn);
RcppExport SEXP _scstruc_col_replicate(SEXP weightSEXP, SEXP rnSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type weight(weightSEXP);
    Rcpp::traits::input_parameter< int >::type rn(rnSEXP);
    rcpp_result_gen = Rcpp::wrap(col_replicate(weight, rn));
    return rcpp_result_gen;
END_RCPP
}
// plasso_fit
arma::vec plasso_fit(arma::mat X, arma::vec y, int maxIter, double lambda, double gamma, double eps);
RcppExport SEXP _scstruc_plasso_fit(SEXP XSEXP, SEXP ySEXP, SEXP maxIterSEXP, SEXP lambdaSEXP, SEXP gammaSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type maxIter(maxIterSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(plasso_fit(X, y, maxIter, lambda, gamma, eps));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_scstruc_row_replicate", (DL_FUNC) &_scstruc_row_replicate, 2},
    {"_scstruc_col_replicate", (DL_FUNC) &_scstruc_col_replicate, 2},
    {"_scstruc_plasso_fit", (DL_FUNC) &_scstruc_plasso_fit, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_scstruc(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
