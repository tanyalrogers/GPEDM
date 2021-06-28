
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector multMat(NumericMatrix m1, NumericMatrix m2) {
  NumericVector m3 = m1 * m2;
  m3.attr("dim") = Dimension(m1.nrow(), m1.ncol());
  return m3;
}
// [[Rcpp::export]]
NumericVector matrix2vector(NumericMatrix m){
  m.attr("dim") = R_NilValue;
  return(m);
}
// [[Rcpp::export]]
double innerProd(NumericVector x, NumericVector y) {
  return std::inner_product(x.begin(), x.end(), y.begin(), 0.0);
}
