// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// multMat
NumericVector multMat(NumericMatrix m1, NumericMatrix m2);
RcppExport SEXP _GPEDM_multMat(SEXP m1SEXP, SEXP m2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type m1(m1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type m2(m2SEXP);
    rcpp_result_gen = Rcpp::wrap(multMat(m1, m2));
    return rcpp_result_gen;
END_RCPP
}
// matrix2vector
NumericVector matrix2vector(NumericMatrix m);
RcppExport SEXP _GPEDM_matrix2vector(SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type m(mSEXP);
    rcpp_result_gen = Rcpp::wrap(matrix2vector(m));
    return rcpp_result_gen;
END_RCPP
}
// innerProd
double innerProd(NumericVector x, NumericVector y);
RcppExport SEXP _GPEDM_innerProd(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(innerProd(x, y));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_GPEDM_multMat", (DL_FUNC) &_GPEDM_multMat, 2},
    {"_GPEDM_matrix2vector", (DL_FUNC) &_GPEDM_matrix2vector, 1},
    {"_GPEDM_innerProd", (DL_FUNC) &_GPEDM_innerProd, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_GPEDM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
