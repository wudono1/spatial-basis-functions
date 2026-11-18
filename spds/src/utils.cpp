#include <iostream>
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
arma::mat expCorr(const arma::mat& distmat,
                  const double& psi) {

  int nrow = distmat.n_rows;
  int ncol = distmat.n_cols;
  double dd;
  arma::mat covmat(nrow, ncol);

  for (int i = 0; i < nrow; ++i) {
    for (int j = 0; j < ncol; ++j) {
      dd = distmat(i,j);
      covmat(i,j) = exp(-dd / psi);
    }
  }

  return covmat;

}

double logit(double theta, double a, double b){
  return log((theta-a)/(b-theta));
}

double logitInv(double z, double a, double b){
  return b-(b-a)/(1+exp(z));
}


