#include "pacabase.h"
#include "utils.h"
#include "pacautils.h"
#include "pca.h"

// [[Rcpp::plugins(cpp14)]]
// [[Rcpp::depends(RcppEigen)]]
using namespace std;
using namespace Rcpp;

////////////////////////////////// Main //////////////////////////////////

//' @export
//' @noRd
// [[Rcpp::export(cpp_prcomp)]]
Rcpp::List cpp_prcomp(const Eigen::MatrixXd& X, bool center,
                      bool scale, Rcpp::Nullable<int> rank ,
                       Rcpp::Nullable<double> tol,
                       int verbosity) {
  Logger::SetVerbosity(verbosity);
  Timer timer("RcppEigen PCA");

  int n = X.rows();
  int p = X.cols();
  Logger::LogLOG("Performing PCA with Eigen in C++ ...");
  // Create a copy of X to avoid mutating the input
  Eigen::MatrixXd X_copy = X;

  // Center the data
  Eigen::VectorXd center_vec = Eigen::VectorXd::Zero(p);
  if (center) {
    center_vec = X_copy.colwise().mean();
    X_copy.rowwise() -= center_vec.transpose();
  }

  // Scale the data
  Eigen::VectorXd scale_vec = Eigen::VectorXd::Ones(p);
  if (scale) {
    scale_vec = (X_copy.array().square().colwise().sum() / (n - 1)).sqrt();
    for (int i = 0; i < p; ++i) {
      if (scale_vec(i) != 0) {
        X_copy.col(i) /= scale_vec(i);
      }
    }
  }

  // Determine rank
  int rank_to_use = std::min(n, p);
  if (rank.isNotNull()) {
    rank_to_use = std::min(rank_to_use, Rcpp::as<int>(rank));
  }

  // Perform SVD
  Eigen::BDCSVD<Eigen::MatrixXd> svd(X_copy, Eigen::ComputeThinU | Eigen::ComputeThinV);

  // Apply tolerance if specified
  if (tol.isNotNull()) {
    double tol_value = Rcpp::as<double>(tol);
    rank_to_use = std::min(rank_to_use,
                           (int)(svd.singularValues().array() > (tol_value * svd.singularValues()(0))).count());
  }

  // Calculate standard deviations
  Eigen::VectorXd sdev = svd.singularValues().head(rank_to_use) / std::sqrt(n - 1);
  Logger::LogLOG("PCA: DONE");
  return Rcpp::List::create(
    Rcpp::Named("sdev") = sdev,
    Rcpp::Named("rotation") = svd.matrixV().leftCols(rank_to_use),
    Rcpp::Named("x") = (X_copy * svd.matrixV().leftCols(rank_to_use)).eval(),
    Rcpp::Named("center") = center_vec,
    Rcpp::Named("scale") = scale_vec
  );
}