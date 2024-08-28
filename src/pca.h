#ifndef PACA_PACA_H
#define PACA_PACA_H

#include <pacabase.h>

Rcpp::List cpp_prcomp(const Eigen::MatrixXd& X, bool center = true, bool scale = false, 
Rcpp::Nullable<int> rank = R_NilValue, Rcpp::Nullable<double> tol  = R_NilValue, int verbosity  = 1);

#endif