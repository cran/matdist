#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::mat RcppCholesky(arma::mat& X){
  arma::mat L = chol(X, "lower");
  return(L);
}
