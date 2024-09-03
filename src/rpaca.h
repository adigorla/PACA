#ifndef PACA_RPACA_H
#define PACA_RPACA_H

#include <pacabase.h>

Eigen::VectorXi sampleWithoutReplacement(int n, int k, std::mt19937& gen);

Eigen::MatrixXd scaleCPP_batchmasked(const Eigen::MatrixXd& x, const Eigen::VectorXi& batch_idx);

void residMatrix(Eigen::MatrixXd& X, const Eigen::VectorXd& PAC);

void singlearpaca_pvt(const Eigen::MatrixXd& X, const Eigen::MatrixXd& Y, 
                     Eigen::MatrixXd& iterPAC, Eigen::VectorXd& iterEIG, Eigen::VectorXd& klist,
                     double threshold, const int niter, const int batch,
                     const int rank, bool normalize);

Rcpp::List cpp_rPACA(const Eigen::MatrixXd& X, const Eigen::MatrixXd& Y,
                     int k, int niter, int batch, int rank,
                     bool normalize = false, int verbosity = 1);


Rcpp::List cpp_autorPACA(const Eigen::MatrixXd& X, const Eigen::MatrixXd& Y,
                     int niter, int batch, int rank,
                     double threshold = 10.0,
                     bool normalize = false, int verbosity = 1);

#endif