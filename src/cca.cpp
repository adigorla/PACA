#include "pacabase.h"
#include "utils.h"
#include "pacautils.h"
#include "cca.h"

// [[Rcpp::plugins(cpp14)]]
// [[Rcpp::depends(RcppEigen)]]
using namespace std;
using namespace Rcpp;

////////////////////////////////// Helpers //////////////////////////////////

void solveNormalEQ(const Eigen::MatrixXd& G, Eigen::LLT<Eigen::MatrixXd, Eigen::Upper>& LLtX, Eigen::LLT<Eigen::MatrixXd, Eigen::Upper>& LLtY,
                   Eigen::MatrixXd& retA, Eigen::MatrixXd& retB, Eigen::VectorXd& retRho, bool flip) {

  // 1. Use BDCSVD for better accuracy and efficiency
  Eigen::BDCSVD<Eigen::MatrixXd> svd_solver(G, Eigen::ComputeThinU | Eigen::ComputeThinV);

  // 2. Check SVD convergence
  if (!(svd_solver.matrixU().cols() > 0)) {
    Logger::LogWARN("SVD staus: FAILED");
    return;
  } else{
    Logger::LogDEBUG("SVD status : ", ((svd_solver.matrixU().cols() > 0) ? "ok" : "FAILED"));
  }

  // 3. Compute canonical correlations and vectors
  Eigen::MatrixXd U = svd_solver.matrixU();
  Eigen::MatrixXd V = svd_solver.matrixV();

  if (flip) {
    std::swap(U, V);
  }

  // 4. Solve for canonical vectors
  retA = LLtX.matrixU().solve(U);
  retB = LLtY.matrixU().solve(V);

  // // 4. Solve for canonical vectors using QR decomposition for better numerical stability
  // Eigen::HouseholderQR<Eigen::MatrixXd> qr_X(LLtX.matrixU());
  // Eigen::HouseholderQR<Eigen::MatrixXd> qr_Y(LLtY.matrixU());
  // retA = qr_X.solve(U);
  // retB = qr_Y.solve(V);

  // 5. Set canonical correlations
  retRho = svd_solver.singularValues();

  // // Log top 10 eigenvalues (or fewer if there aren't 10 available)
  // int num_top_eigenvalues = std::min(10, static_cast<int>(retRho.size()));
  // for (int i = 0; i < num_top_eigenvalues; ++i) {
  //   Logger::LogDEBUG("SVD eig : ", i, "- ", retRho(i));
  // }

};

void doCCA_pvt(Eigen::MatrixXd& X, Eigen::MatrixXd& Y, bool normalize,
               Eigen::VectorXd& eigs, Eigen::MatrixXd& A, Eigen::MatrixXd& B,
               Eigen::MatrixXd& U, Eigen::MatrixXd& V) {
  Logger::LogINFO("Running CCA ...");
  int p = X.cols();
  int n = X.rows();
  int q = Y.cols();
  Logger::LogLOG("X Shape : ", n, "x", p);
  Logger::LogLOG("Y Shape : ", n, "x", q);

  if (normalize) {
    Logger::LogLOG("Normalizing ...");
    {
      Timer timer("(In-Place) Normalization");
      X = normalizeCPP(X);
      Y = normalizeCPP(Y);
    }
    Logger::LogLOG("Done Normalizing");
  }

  // calculate the covariance matrix ; 0: unbiased estimator (1/(n-1)), 1: 2nd moment estimator (1/n)
  // Assuming normalized
  Logger::LogINFO("Estimating Covariances ...");

  Eigen::MatrixXd Cxx;
  {
    Timer timer("X Covariance");
    Cxx = covCPP(X);
  }
  Logger::LogLOG("Done with X Covariance");

  Eigen::MatrixXd Cyy;
  {
    Timer timer("Y Covariance");
    Cyy = covCPP(Y);
  }
  Logger::LogLOG("Done with Y Covariance");

  int cxr = Cxx.rows();
  int cyr = Cyy.rows();
  Logger::LogLOG("(symmetric) Cxx / Cyy Shape  : ", cxr , " / ", cyr);

  Eigen::MatrixXd Cxy;
  {
    Timer timer("Cross Covariance");
    // Cxy = (X.adjoint() * Y) / double(n - 1.0);
    Cxy = (centerCPP(X).adjoint() * centerCPP(Y)) / double(n - 1);
  }
  Logger::LogLOG("Done with XY Cross-Covariance");

  Logger::LogINFO("Solving CCA Objective ...");

  Logger::LogLOG("Performing Cholesky Decomposition ...");

  Eigen::LLT<Eigen::MatrixXd, Eigen::Upper> lltOfCX;
  Eigen::LLT<Eigen::MatrixXd, Eigen::Upper> lltOfCY;
  {
    Timer timer("Cholesky Decomposition");

    lltOfCY.compute(Cyy); // compute the Cholesky decomposition of Cyy
    if (lltOfCY.info() == Eigen::NumericalIssue){ // check for numerical issues (aka not PD)
      Logger::LogWARN("Y Covariance likely NOT PD, adding constant to it make PD");
      Cyy += Eigen::MatrixXd::Identity(cyr, cyr) * 1e-6;
      lltOfCY.compute(Cyy);
    }

    lltOfCX.compute(Cxx); // compute the Cholesky decomposition of Cxx
    if (lltOfCX.info() == Eigen::NumericalIssue){ // check for numerical issues (aka not PD)
      Logger::LogWARN("X Covariance likely NOT PD, adding constant to it make PD");
      Cxx += Eigen::MatrixXd::Identity(cxr, cxr) * 1e-6;
      lltOfCX.compute(Cxx);
    }

  }
  Logger::LogLOG("Done with Cholesky decomp");

  // // clear memory
  Cxx.resize(0, 0);
  Cyy.resize(0, 0);

  Logger::LogLOG("Performing Gamma Calculation ...");
  Eigen::MatrixXd Gamma;
  {
    // The resulting matrix should equivalent to Cxx_fi.transpose() * Cxy * Cyy_fi
    Timer timer("Gamma Calculation");

    // Solve Cxx * Xxxxy = Cxy for Xxxxy
    Eigen::MatrixXd Xxxxy = lltOfCX.matrixU().transpose().solve(Cxy);

    // Solve Cyy * Gamma = t(Xxxxy) for Gamma; t(Gamma) = Cxx^-1*Cxy*Cyy^-1*Cyx
    Gamma = lltOfCY.matrixU().transpose().solve(Xxxxy.transpose()).transpose();

    Logger::LogDEBUG("Gamma:" , Gamma.rows(), "x", Gamma.cols());
  }
  Logger::LogLOG("Done with Gamma Calculation");

  // clear memory
  Cxy.resize(0, 0);

  Logger::LogLOG("Performing SVD ...");
  {
    Timer timer("SVD");
    if (p >= q) {
      solveNormalEQ(Gamma, lltOfCX, lltOfCY,
                    A, B, eigs, false);
    } else {
      solveNormalEQ(Gamma, lltOfCX, lltOfCY,
                    A, B, eigs, true);
    }
  }
  Logger::LogLOG("Done with SVD");
  // // clear memory
  Gamma.resize(0, 0);

  adjCanCor(X, A, Y, B,
            U, V, eigs);
};


////////////////////////////////// Main //////////////////////////////////

//' @export
//' @noRd
// [[Rcpp::export(cpp_CCA)]]
Rcpp::List cpp_CCA(const Eigen::MatrixXd& Xo, const Eigen::MatrixXd& Yo,
                   bool normalize, int verbosity) {
  Logger::SetVerbosity(verbosity);
  // Create local copies of X and Y to avoid modifying the original matrices
  Eigen::MatrixXd X = Xo;
  Eigen::MatrixXd Y = Yo;

  Eigen::MatrixXd A, B, U, V;
  Eigen::VectorXd eigs;
  doCCA_pvt(X, Y, normalize, eigs, A, B, U, V);

  Logger::LogINFO("CCA: DONE");

  return Rcpp::List::create(Rcpp::Named("corr")=eigs,
                           Rcpp::Named("A")=A,
                           Rcpp::Named("B")=B,
                           Rcpp::Named("U")=U,
                           Rcpp::Named("V")=V);
};