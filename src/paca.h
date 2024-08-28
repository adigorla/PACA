#ifndef PACA_PACA_H
#define PACA_PACA_H

#include <pacabase.h>

void selectK_pvt(const Eigen::MatrixXd& X, const Eigen::MatrixXd& Y, int& selK, bool normalize,
                 Eigen::VectorXd& eigs, Eigen::MatrixXd& A, Eigen::MatrixXd& B,
                 Eigen::MatrixXd& U, Eigen::MatrixXd& V, double thrsh);


Rcpp::List cpp_selectK(const Eigen::MatrixXd& X, const Eigen::MatrixXd& Y,
                       bool normalize = true, double threshold  = 10.0,
                       int verbosity = 1);


Rcpp::List cpp_PACA(const Eigen::MatrixXd& Xo, const Eigen::MatrixXd& Yo, int k,
                    bool normalize = true, bool retCCA = false,
                    int verbosity = 1);


Rcpp::List cpp_autoPACA(const Eigen::MatrixXd& Xo, const Eigen::MatrixXd& Yo,
                        bool normalize = true, bool retCCA = false,
                        double threshold = 10.0, int verbosity = 1);


#endif