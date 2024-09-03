#include "pacabase.h"
#include "utils.h"
#include "pacautils.h"
#include "cca.h"
#include "paca.h"

// [[Rcpp::plugins(cpp14)]]
// [[Rcpp::depends(RcppEigen)]]

using namespace std;
using namespace Rcpp;

////////////////////////////////// Helper //////////////////////////////////

void selectK_pvt(const Eigen::MatrixXd& X, const Eigen::MatrixXd& Y, int& selK, bool normalize,
                 Eigen::VectorXd& eigs, Eigen::MatrixXd& A, Eigen::MatrixXd& B,
                 Eigen::MatrixXd& U, Eigen::MatrixXd& V, double thrsh){
  // Initialize parameters
  int p = X.cols();
  // int n = X.rows();
  int q = Y.cols();
  int d = std::min(p, q);
  int min_k = 1;
  int max_k = d;
  int Ulim = max_k;
  int n_perm = std::floor(std::log2(d)/0.01);
  int counter = 0;
  int curr_k = 0;

  // Initialize variables
  Eigen::MatrixXd X_in = X;
  Eigen::MatrixXd Y_in = Y;

  if (normalize) {
    Logger::LogLOG("Normalizing ...");
    {
      Timer timer("(In-Place) Normalization");
      X_in = normalizeCPP(X_in);
      Y_in = normalizeCPP(Y_in);
    }
    Logger::LogLOG("Done Normalizing");
  }

  doCCA_pvt(X_in, Y_in, false, eigs, A, B, U, V);

  int dU1 = U.rows();
  int dV1 = V.rows();

  Logger::LogINFO("Starting K Selection binary search ...");
  // Binary Search
  Eigen::MatrixXd U1, V1;
  Eigen::MatrixXd tXtil, tYtil;
  Eigen::VectorXd Xtil_load;
  Eigen::VectorXd projected_x , projected_y, var_Xperm, var_Yperm;
  Eigen::MatrixXd Zperm, projected_Xperm, projected_Yperm;
  double var_x, var_y;

  // Synchronize R's RNG state
  GetRNGstate();

  // Use the seed from R
  // Generate a seed using R's RNG
  unsigned int seed = static_cast<unsigned int>(unif_rand() * UINT_MAX);
  Logger::LogDEBUG("setting seed :", seed);
  std::mt19937 g(seed);
  // std::random_device rd;
  // std::mt19937 g(rd());

  while (min_k < max_k){
    curr_k = std::floor((min_k + max_k)/2.0);
    counter += 1;
    if (curr_k < 2){
      break;
    }
    Logger::LogLOG("Iter:", counter, "-- K:", curr_k, ", Searching range [", min_k, ",", max_k, "]");

    U1 = UV1_calc(U, curr_k, dU1);
    V1 = UV1_calc(V, curr_k, dV1);
    // clear memory

    // Calculate new X and Y
    tXtil = correctedMat_calc(U1, X_in, true);
    tYtil = correctedMat_calc(V1, Y_in, true);
    // clear memory
    U1.resize(0, 0);
    V1.resize(0, 0);
    // M1: Porbs good till here
    // top PC loading
    Xtil_load = topPC_loading(tXtil); // compare this to other approachs

    projected_x = tXtil * Xtil_load;
    var_x = colwiseVars(projected_x)(0,0);
    projected_y = tYtil * Xtil_load;
    var_y = colwiseVars(projected_y)(0,0);
    projected_x.resize(0);
    projected_y.resize(0);

    // Fast Permutation testing
    Zperm = Xtil_load.replicate(1, n_perm);
    Xtil_load.resize(0);
    // Shuffle each column
    for (int i = 0; i < n_perm; ++i) {
      std::shuffle(Zperm.col(i).data(), Zperm.col(i).data() + Zperm.col(i).size(), g);
    }

    // Compute variances
    projected_Xperm = tXtil * Zperm;
    var_Xperm = colwiseVars(projected_Xperm);
    projected_Yperm = tYtil * Zperm;
    var_Yperm = colwiseVars(projected_Yperm);
    Zperm.resize(0, 0);
    projected_Xperm.resize(0, 0);
    projected_Yperm.resize(0, 0);

    // Count how many variances are greater than var_x2 and var_y2
    int ct_x = (var_Xperm.array() > var_x).count();
    int ct_y = (var_Yperm.array() > var_y).count();
    // clear memory
    var_Xperm.resize(0);
    var_Yperm.resize(0);

    Logger::LogDEBUG("\tPermutation results: var_x:", var_x, ", var_y:", var_y, ", ct_x:", ct_x, ", ct_y:", ct_y);

    // check conditions
    if ((ct_x < 1) && (ct_y >= 1)){
      max_k = curr_k;
    }
    else if((ct_x < 1) && (ct_y < 1)){
      if ( var_x > (thrsh * var_y)){
        max_k = curr_k;
      } else{
        min_k = curr_k + 1;
      }
    }
    else if((ct_x >= 1) && (ct_y >= 1)){
      max_k = curr_k - 1;
    }
    else {
      break;
    }
  }
  // prevents k0 from being unreasonably high
  if (curr_k >= Ulim - 1){
    curr_k = Ulim - 3;
  }

  // Synchronize the RNG state back to R
  PutRNGstate();

  selK = curr_k+1;
  Logger::LogLOG("Selected K:", selK);
  Logger::LogINFO("Done with K Selection");
};

////////////////////////////////// Main PACA functions //////////////////////////////////

//' @export
//' @noRd
// [[Rcpp::export(cpp_selectK)]]
Rcpp::List cpp_selectK(const Eigen::MatrixXd& X, const Eigen::MatrixXd& Y,
                       bool normalize, double threshold,
                       int verbosity){
 Logger::SetVerbosity(verbosity);

 Logger::LogINFO("selectK: Starting ...");
 Eigen::MatrixXd A, B, U, V;
 Eigen::VectorXd eigs;
 int selK;
 selectK_pvt(X, Y, selK, normalize, eigs, A, B, U, V, threshold);

 int dU1 = U.rows();
 Eigen::MatrixXd U0 = UV1_calc(U, selK, dU1);

 Logger::LogINFO("selectK: DONE");

 return Rcpp::List::create(Rcpp::Named("k")=selK,
                           Rcpp::Named("U0")=U0);
};


//' @export
//' @noRd
// [[Rcpp::export(cpp_PACA)]]
Rcpp::List cpp_PACA(const Eigen::MatrixXd& Xo, const Eigen::MatrixXd& Yo, int k,
                    bool normalize, bool retCCA,
                    int verbosity){
  Logger::SetVerbosity(verbosity);

  // Create local copies of X and Y to avoid modifying the original matrices
  Eigen::MatrixXd X = Xo;
  Eigen::MatrixXd Y = Yo;

  Logger::LogINFO("PACA: Starting ...");
  // perform CCA
  Eigen::MatrixXd A, B, U, V;
  Eigen::VectorXd eigs;

  //// DEBUG
  // Eigen::MatrixXd Cxx, Cyy, Cxy, Gamma;
  // doCCA_pvt(X, Y, normalize, eigs, A, B, U, V, Cxx, Cyy, Cxy, Gamma);
  ////

  doCCA_pvt(X, Y, normalize, eigs, A, B, U, V);

  Logger::LogLOG("Done with CCA");

  Logger::LogINFO("Residualizing Shared Signal ...");
  // Calculate shared components
  int dU1 = U.rows();
  Eigen::MatrixXd U0 = UV1_calc(U, k, dU1);

  // Calcute X matrix with case specific variation only
  Eigen::MatrixXd Xtil = correctedMat_calc(U0, Xo, false);

  Logger::LogINFO("PACA: DONE");

  if (retCCA){
   return Rcpp::List::create(Rcpp::Named("Xtil")=Xtil,
                             Rcpp::Named("U0")=U0,
                             Rcpp::Named("corr")=eigs,
                             Rcpp::Named("A")=A,
                             Rcpp::Named("B")=B,
                             Rcpp::Named("U")=U,
                             Rcpp::Named("V")=V);
    // Rcpp::Named("Gamma")=Gamma,
    // Rcpp::Named("Cxx")=Cxx,
    // Rcpp::Named("Cyy")=Cyy,
    // Rcpp::Named("Cxy")=Cxy)
  }
  return Rcpp::List::create(Rcpp::Named("Xtil")=Xtil,
                           Rcpp::Named("U0")=U0);
};


//' @export
//' @noRd
// [[Rcpp::export(cpp_autoPACA)]]
Rcpp::List cpp_autoPACA(const Eigen::MatrixXd& Xo, const Eigen::MatrixXd& Yo,
                        bool normalize, bool retCCA ,
                        double threshold, int verbosity){
  Logger::SetVerbosity(verbosity);

  Logger::LogINFO("autoPACA: Starting ...");
  // perform CCA + select K
  Eigen::MatrixXd A, B, U, V;
  Eigen::VectorXd eigs;

  // Estimate number of shared components
  int selK;
  selectK_pvt(Xo, Yo, selK, normalize, eigs, A, B, U, V, threshold);
  Logger::LogLOG("K_est : ", selK);


  Logger::LogINFO("Residualizing Shared Signal ...");
  // Calculate shared components
  int dU1 = U.rows();
  Eigen::MatrixXd U0 = UV1_calc(U, selK, dU1);

  // Calcute X matrix with case specific variation only
  // Eigen::MatrixXd X = Xo;
  // if (normalize){
  //   X = normalizeCPP(X);
  // }
  Eigen::MatrixXd Xtil = correctedMat_calc(U0, Xo, false);


  Logger::LogINFO("autoPACA: DONE");

  if (retCCA){
   return Rcpp::List::create(Rcpp::Named("Xtil")=Xtil,
                             Rcpp::Named("U0")=U0,
                             Rcpp::Named("k")=selK,
                             Rcpp::Named("corr")=eigs,
                             Rcpp::Named("A")=A,
                             Rcpp::Named("B")=B,
                             Rcpp::Named("U")=U,
                             Rcpp::Named("V")=V);
  }
  return Rcpp::List::create(Rcpp::Named("Xtil")=Xtil,
                           Rcpp::Named("k")=selK,
                           Rcpp::Named("U0")=U0);

};



