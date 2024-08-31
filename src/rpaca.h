#ifndef PACA_RPACA_H
#define PACA_RPACA_H

#include <pacabase.h>

Eigen::VectorXi sampleWithoutReplacement(int n, int k, std::mt19937& gen);

Eigen::MatrixXd scaleCPP_batchmasked(const Eigen::MatrixXd& x, const Eigen::VectorXi& batch_idx);

Rcpp::List cpp_rPACA(const Eigen::MatrixXd& X, const Eigen::MatrixXd& Y,
                     int k, int niter, int batch, int rank,
                     bool normalize = false, int verbosity = 1);

#endif