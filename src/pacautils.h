#ifndef PACA_UTILS_H
#define PACA_UTILS_H

#include "pacabase.h"

Eigen::MatrixXd getUV(Eigen::MatrixXd& x, Eigen::MatrixXd& eigspc);

Eigen::MatrixXd UV1_calc(const Eigen::MatrixXd& UV, int k, int dUV1);

Eigen::MatrixXd correctedMat_calc(const Eigen::MatrixXd& UV, const Eigen::MatrixXd& XY,  bool flip) ; 

Eigen::VectorXd topPC_loading(Eigen::MatrixXd& A);

void adjCanCor(Eigen::MatrixXd& x, Eigen::MatrixXd& a, Eigen::MatrixXd& y, Eigen::MatrixXd& b,
               Eigen::MatrixXd& retU, Eigen::MatrixXd& retV, Eigen::VectorXd& retRho);

#endif