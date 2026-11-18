#ifndef UTILS_H
#define UTILS_H
#include <RcppArmadillo.h>

double logit(double theta, double a, double b);

double logitInv(double z, double a, double b);

arma::mat expCorr(const arma::mat& distmat, const double& psi);

#endif
