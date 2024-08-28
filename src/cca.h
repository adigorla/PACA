#ifndef PACA_CCA_H
#define PACA_CCA_H

#include <pacabase.h>

void solveNormalEQ(const Eigen::MatrixXd& G, Eigen::LLT<Eigen::MatrixXd, Eigen::Upper>& LLtX, Eigen::LLT<Eigen::MatrixXd, Eigen::Upper>& LLtY,
                   Eigen::MatrixXd& retA, Eigen::MatrixXd& retB, Eigen::VectorXd& retRho, bool flip = false);

void doCCA_pvt(Eigen::MatrixXd& X, Eigen::MatrixXd& Y, bool normalize,
               Eigen::VectorXd& eigs, Eigen::MatrixXd& A, Eigen::MatrixXd& B,
               Eigen::MatrixXd& U, Eigen::MatrixXd& V);

Rcpp::List cpp_CCA(const Eigen::MatrixXd& Xo, const Eigen::MatrixXd& Yo,
                   bool normalize = true, int verbosity = 1);

#endif