#include "pacabase.h"
#include "utils.h"
#include "pacautils.h"

// [[Rcpp::plugins(cpp14)]]
// [[Rcpp::depends(RcppEigen)]]
using namespace std;
using namespace Rcpp;

////////////////////////////////// Core Computation Helpers //////////////////////////////////

Eigen::MatrixXd getUV(Eigen::MatrixXd& x, Eigen::MatrixXd& eigspc) {
  Eigen::MatrixXd x_centered = centerCPP(x);
  return x_centered * eigspc;
};

Eigen::MatrixXd UV1_calc(Eigen::MatrixXd& UV, int k, int dUV1){
  Eigen::MatrixXd UV_k = UV.leftCols(k);  // Select first k columns of U
  // Calculate the column sums of squares
  Eigen::VectorXd col_sums = UV_k.colwise().squaredNorm();

  // Calculate the square root of the column sums
  Eigen::VectorXd sqrt_col_sums = col_sums.array().sqrt();

  // Replicate sqrt_col_sums across the rows
  Eigen::MatrixXd divisor = sqrt_col_sums.transpose().replicate(dUV1, 1);

  // Perform the element-wise division
  Eigen::MatrixXd result = UV_k.array() / divisor.array();

  return result;
};

Eigen::MatrixXd correctedMat_calc(Eigen::MatrixXd& UV, const Eigen::MatrixXd& XY,  bool flip = false){

  Eigen::VectorXd col_means = XY.colwise().mean();  // Calculate column means

  // Replicate the column means across the rows to form means_matrix
  Eigen::MatrixXd means_matrix = col_means.transpose().replicate(XY.rows(), 1);

  // Center X by subtracting the means_matrix
  Eigen::MatrixXd XY_centered = XY - means_matrix;

  // Perform the projection onto U_1
  Eigen::MatrixXd projection = UV * (UV.transpose() * XY_centered);

  // Subtract the projection from X_centered
  Eigen::MatrixXd XY_til = XY_centered - projection;

  // Add means_matrix back to get X_tilde
  XY_til += means_matrix;

  if (flip){
    return XY_til.transpose();
  } else{
    return XY_til;
  }
};

Eigen::VectorXd topPC_loading(Eigen::MatrixXd& A){
  A = centerCPP(A); // alt normalizeCPP
  Eigen::BDCSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);  // Perform SVD

  Eigen::VectorXd loadings = svd.matrixV().col(0);  // The first column corresponds to the loadings of the top PC // * sqrt(svd.singularValues()(0))
  Logger::LogDEBUG("tPCl stats { Vdim: ", svd.matrixV().rows(), " x ", svd.matrixV().cols(), " S(0): ", svd.singularValues()(0),
                   "Udim: ", svd.matrixU().rows(), " x ", svd.matrixU().cols()," }");
  // TODO: Seems like the sqrt(S0) is the var of the projection, is var of proj on the loading the sqrt of the singular value????=
  return loadings;
};

void adjCanCor(Eigen::MatrixXd& x, Eigen::MatrixXd& a, Eigen::MatrixXd& y, Eigen::MatrixXd& b,
               Eigen::MatrixXd& retU, Eigen::MatrixXd& retV, Eigen::VectorXd& retRho){
  {
    Timer timer("U Calc");
    retU = x*a;
  }

  {
    Timer timer("V Calc");
    retV = y*b;
  }

  {
    Timer timer("CanCor");
    retRho = (normalizeCPP(retU, false).transpose() * normalizeCPP(retV, false)).diagonal().array() / double(retU.rows() - 1);
  }
};

