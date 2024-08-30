#ifndef PACA_RPACA_H
#define PACA_RPACA_H

#include <pacabase.h>

Rcpp::List cpp_rPACA(const Eigen::MatrixXd& X, const Eigen::MatrixXd& Y,
                     int k, int niter, int batch, int rank,
                     bool normalize = false, int verbosity = 1);

#endif