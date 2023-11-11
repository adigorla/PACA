// includes to make Eigen use BLAS
#define EIGEN_USE_BLAS

#include <Rcpp.h>
#include <RcppEigen.h>
#include <vector>
#include <chrono>
#include <random>
#include <algorithm>
#include <mutex>

// [[Rcpp::plugins(cpp14)]]
// [[Rcpp::depends(RcppEigen)]]
using namespace std;
using namespace Rcpp;

// TODO : autorPACA, rPACA
// TODO : how to make eigen use LAPACK?

////////////////////////////////// Loggers //////////////////////////////////

// Logger class
class Logger {
public:
  static void SetVerbosity(int level) {
    std::lock_guard<std::mutex> lock(mtx);
    verbosity = level;
  }

  template<typename... Args>
  static void LogWARN(Args... args) {
    Log(0, "WARN  : ", args...);
  }

  template<typename... Args>
  static void LogINFO(Args... args) {
    Log(1, "INFO  : ", args...);
  }

  template<typename... Args>
  static void LogLOG(Args... args) {
    Log(2, "LOG   : ", args...);
  }

  template<typename... Args>
  static void LogDEBUG(Args... args) {
    Log(3, "DEBUG : ", args...);
  }

private:
  // recurssive overloaded function to concatenate args since we cant use C++17 fold expressions
  template<typename T>
  static void concatenate_args(std::ostringstream& oss, const T& arg) {
    oss << arg << " ";
  }
  template<typename T, typename... Args>
  static void concatenate_args(std::ostringstream& oss, const T& arg, const Args&... args) {
    oss << arg << " ";
    concatenate_args(oss, args...);
  }

  template<typename... Args>
  static void Log(int level, const std::string &prefix, const Args&... args) {
    std::lock_guard<std::mutex> lock(mtx);
    if (level > verbosity) {
      return;
    }
    const auto in_time_t = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());

    std::ostringstream oss;
    oss << "[" << std::put_time(std::localtime(&in_time_t), "%Y-%m-%d %X") << "] " << prefix;

    // ((oss << args << " "), ...);
    concatenate_args(oss, args...);

    switch (level) {
    case 0:
      Rcpp::Rcout << "\033[1;33m";  // yellow
      break;
    case 1:
      // No color for INFO
      break;
    case 2:
      Rcpp::Rcout << "\033[1m\033[36m";  // Bold Cyan
      break;
    case 3:
      Rcpp::Rcout << "\033[1m\033[35m";  // Bold Magenta
      break;
    default:
      return;
    }

    Rcpp::Rcout << oss.str() << "\033[0m\n";  // reset color
  }

  static std::mutex mtx;
  static int verbosity;
};

// Timer for debugging
class Timer {
public:
  Timer(string fname){
    m_startTimept = std::chrono::high_resolution_clock::now();
    m_fname = fname;
  }
  ~Timer(){
    Stop();
  }
private:
  std::chrono::time_point<std::chrono::high_resolution_clock> m_startTimept;
  std::string m_fname;

  void Stop(){
    auto m_endTimept = std::chrono::high_resolution_clock::now();
    auto start = std::chrono::time_point_cast<std::chrono::microseconds>(m_startTimept).time_since_epoch().count();
    auto end = std::chrono::time_point_cast<std::chrono::microseconds>(m_endTimept).time_since_epoch().count();
    auto duration = end - start;
    std::ostringstream LOCAL_oss;
    LOCAL_oss << "\t{"<< m_fname <<"} took: " << duration*0.001 << "ms / " << duration << "us";

    Logger::LogDEBUG(LOCAL_oss.str());

  }
};

////////////////////////////////// Helpers //////////////////////////////////

// [[Rcpp::export]]
Eigen::MatrixXd normalizeCPP(Eigen::MatrixXd& x, bool inplace = true) {
  // col means
  Eigen::VectorXd xmeans = x.colwise().mean();
  // col sds
  Eigen::VectorXd xstd = ((x.rowwise() - xmeans.transpose()).array().square().colwise().sum() / (x.rows() - 1)).sqrt();
  // inplace
  if(!inplace){
    Eigen::MatrixXd xnew  = x;
    xnew.rowwise() -= xmeans.transpose();
    xnew.array().rowwise() /= xstd.transpose().array();
    return xnew;
  }
  x.rowwise() -= xmeans.transpose();
  x.array().rowwise() /= xstd.transpose().array();
  return x;
};

Eigen::MatrixXd centerCPP(Eigen::MatrixXd& x) {
  Eigen::VectorXd xcmeans = x.colwise().mean();
  Eigen::MatrixXd xhold(x.rows(), x.cols());
  for(int i = 0; i < x.rows(); i++) {
    xhold.row(i) = xcmeans.transpose();
  }
  return x - xhold;
};

Eigen::VectorXd colwiseVars(const Eigen::MatrixXd& matrix) {
  // Number of rows (samples)
  int rows = matrix.rows();

  // Check if there is enough data to compute variance
  if (rows < 2) {
    throw std::invalid_argument("Not enough data to compute variance");
  }

  // Calculate the means of each column
  Eigen::VectorXd mean = matrix.colwise().mean();

  // Subtract the mean from each column and square the result
  Eigen::MatrixXd centered = matrix.rowwise() - mean.transpose();
  centered = centered.array().square();

  // Compute the variance of each column
  Eigen::VectorXd variances = (centered.colwise().sum()) / (rows - 1);

  return variances;
};

Eigen::MatrixXd getUV(Eigen::MatrixXd& x, Eigen::MatrixXd& eigspc) {
  Eigen::MatrixXd x_centered = centerCPP(x);
  return x_centered * eigspc;
};

void solveNormalEQ(Eigen::MatrixXd& G, Eigen::LLT<Eigen::MatrixXd, Eigen::Upper>& LLtX, Eigen::LLT<Eigen::MatrixXd, Eigen::Upper>& LLtY,
                   Eigen::MatrixXd& retA, Eigen::MatrixXd& retB, Eigen::VectorXd& retRho, bool flip = false){
  if(flip){
    G = G.transpose();
  }

  Eigen::BDCSVD<Eigen::MatrixXd> svd_solver(G, Eigen::ComputeThinU | Eigen::ComputeThinV);
  if (!(svd_solver.matrixU().cols() > 0)){
    Logger::LogWARN("SVD status : FAILED");
  }
  Logger::LogDEBUG("SVD status : ", ((svd_solver.matrixU().cols() > 0) ? "ok" : "FAILED"));

  if(!flip){
    retA = LLtX.matrixU().solve(svd_solver.matrixU());
    retB = LLtY.matrixU().solve(svd_solver.matrixV());
    retRho = svd_solver.singularValues();
  } else{
    retA = LLtX.matrixU().solve(svd_solver.matrixV());
    retB = LLtY.matrixU().solve(svd_solver.matrixU());
    retRho = svd_solver.singularValues();
  }
};

Eigen::MatrixXd UV1_calc(Eigen::MatrixXd& UV, int k, int dUV1){
  Eigen::MatrixXd UV_k = UV.leftCols(k);  // Select first k columns of U
  Eigen::VectorXd col_sums = UV_k.array().square().colwise().sum().sqrt();  // Calculate sqrt of sum of squares of each column

  Eigen::MatrixXd ones = Eigen::MatrixXd::Ones(dUV1, 1);  // Create a dU1 x 1 matrix of ones
  Eigen::MatrixXd kronecker_product = ones * col_sums.transpose();  // This does the same thing as the kronecker product in this specific case

  Eigen::MatrixXd result = UV_k.array() / kronecker_product.array();  // Divide U_k by kronecker_product

  return result;
};

Eigen::MatrixXd correctedMat_calc(Eigen::MatrixXd& UV, Eigen::MatrixXd& XY,  bool flip = false){
  Eigen::RowVectorXd col_means = XY.colwise().mean();  // Calculate column means

  Eigen::MatrixXd centered = XY.rowwise() - col_means;  // Subtract mean from each row

  Eigen::MatrixXd XY_til = centered - UV * (UV.transpose() * centered);
  XY_til.rowwise() += col_means;  // Add back the means

  if (flip){
    return XY_til.transpose();
  } else{
    return XY_til;
  }
};

Eigen::VectorXd topPC_loading(Eigen::MatrixXd& A){
  A = centerCPP(A);
  Eigen::BDCSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);  // Perform SVD

  Eigen::VectorXd loadings = svd.matrixV().col(0);  // The first column corresponds to the loadings of the top PC
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
    Cxx = (X.adjoint() * X) / double(n - 1.0);
    // Cxx = (Cxx + Cxx.adjoint()) / 2.0;
  }
  Logger::LogLOG("Done with X Covariance");

  Eigen::MatrixXd Cyy;
  {
    Timer timer("Y Covariance");
    Cyy = (Y.adjoint() * Y) / double(n - 1.0);
    // Cyy = (Cyy + Cyy.transpose()) / 2.0;
  }
  Logger::LogLOG("Done with Y Covariance");
  int cxr = Cxx.rows();
  int cyr = Cyy.rows();
  Logger::LogLOG("(symmetric) Cxx / Cyy Shape  : ", cxr , " / ", cyr);
  Eigen::MatrixXd Cxy;
  {
    Timer timer("Cross Covariance");
    Cxy = (X.adjoint() * Y) / double(n - 1.0);
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
  // clear memory
  Cxx.resize(0, 0);
  Cyy.resize(0, 0);

  Logger::LogLOG("Performing Gamma Calculation ...");
  Eigen::MatrixXd Gamma;
  {
    Timer timer("Gamma Calculation");
    // The resulting matrix should equivalent to Cxx_fi.transpose() * Cxy * Cyy_fi
    Eigen::MatrixXd Gamma_temp = lltOfCX.matrixU().transpose().solve(Cxy);
    Gamma = lltOfCY.matrixU().transpose().solve(Gamma_temp.transpose()).transpose();
    Logger::LogDEBUG("Gamma:" , Gamma.rows(), "x", Gamma.cols());
  }
  // Test1 = Gamma;
  Logger::LogLOG("Done with Gamma Calculation");

  Logger::LogINFO("Performing SVD ...");
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

void selectK_pvt(Eigen::MatrixXd& X, Eigen::MatrixXd& Y, int& selK, bool normalize,
                 Eigen::VectorXd& eigs, Eigen::MatrixXd& A, Eigen::MatrixXd& B,
                 Eigen::MatrixXd& U, Eigen::MatrixXd& V, double thrsh){
  // Initialize parameters
  int p = X.cols();
  // int n = X.rows();
  int q = Y.cols();
  int d = std::min(p, q);
  int min_k = 1;
  int max_k = p;
  int Ulim = max_k;
  int n_perm = std::floor(std::log2(d)/0.01);
  int counter = 0;
  int curr_k = 0;

  // Initialize variables
  Eigen::MatrixXd& X_in = X;
  Eigen::MatrixXd& Y_in = Y;

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
  std::random_device rd;
  std::mt19937 g(rd());
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
  selK = curr_k+1;
  Logger::LogLOG("Selected K:", selK);
  Logger::LogINFO("Done with K Selection");
};

// TODO : add pvt function for the the sampling and iteration loop


////////////////////////////////// Mains //////////////////////////////////

int Logger::verbosity = 1;  // Default verbosity
std::mutex Logger::mtx;


//' @export
 // [[Rcpp::export]]
 Rcpp::List cpp_CCA(Eigen::MatrixXd& X, Eigen::MatrixXd& Y, bool normalize = true, int verbosity = 1) {
   Logger::SetVerbosity(verbosity);

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


//' @export
 // [[Rcpp::export]]
 Rcpp::List cpp_selectK(Eigen::MatrixXd& X, Eigen::MatrixXd& Y, bool normalize = true, double threshold  = 10.0, int verbosity = 1){
   Logger::SetVerbosity(verbosity);

   Logger::LogINFO("selectK: Starting ...");
   Eigen::MatrixXd A, B, U, V;
   Eigen::VectorXd eigs;
   int selK, dU1;
   selectK_pvt(X, Y, selK, normalize, eigs, A, B, U, V, threshold);
   dU1 = U.rows();
   Eigen::MatrixXd U0 = UV1_calc(U, selK, dU1);

   Logger::LogINFO("selectK: DONE");

   return Rcpp::List::create(Rcpp::Named("K")=selK,
                             Rcpp::Named("U0")=U0);
 };


//' @export
 // [[Rcpp::export]]
 Rcpp::List cpp_PACA(Eigen::MatrixXd& X, Eigen::MatrixXd& Y, int k, bool normalize = true, bool retCCA = false, int verbosity = 1){
   Logger::SetVerbosity(verbosity);

   Logger::LogINFO("PACA: Starting ...");
   // perform CCA
   Eigen::MatrixXd A, B, U, V;
   Eigen::VectorXd eigs;
   doCCA_pvt(X, Y, normalize, eigs, A, B, U, V);
   Logger::LogLOG("Done with CCA");

   Logger::LogINFO("Residualizing Shared Signal ...");
   int dU1 = U.rows();
   Eigen::MatrixXd U0 = UV1_calc(U, k, dU1);
   Eigen::MatrixXd Xtil = correctedMat_calc(U0, X, false);

   Logger::LogINFO("PACA: DONE");

   if (retCCA){
     return Rcpp::List::create(Rcpp::Named("Xtil")=Xtil,
                               Rcpp::Named("U0")=U0,
                               Rcpp::Named("corr")=eigs,
                               Rcpp::Named("A")=A,
                               Rcpp::Named("B")=B,
                               Rcpp::Named("U")=U,
                               Rcpp::Named("V")=V);
   }
   return Rcpp::List::create(Rcpp::Named("Xtil")=Xtil,
                             Rcpp::Named("U0")=U0);
 };


//' @export
 // [[Rcpp::export]]
 Rcpp::List cpp_autoPACA(Eigen::MatrixXd& X, Eigen::MatrixXd& Y, bool normalize = true, bool retCCA = false, double threshold = 10.0, int verbosity = 1){
   Logger::SetVerbosity(verbosity);

   Logger::LogINFO("autoPACA: Starting ...");
   // perform CCA + select K
   Eigen::MatrixXd A, B, U, V;
   Eigen::VectorXd eigs;
   int selK, dU1;
   selectK_pvt(X, Y, selK, normalize, eigs, A, B, U, V, threshold);
   Logger::LogLOG("K_est : ", selK);


   Logger::LogINFO("Residualizing Shared Signal ...");
   dU1 = U.rows();
   Eigen::MatrixXd U0 = UV1_calc(U, selK, dU1);
   Eigen::MatrixXd Xtil = correctedMat_calc(U0, X, false);


   Logger::LogINFO("autoPACA: DONE");

   if (retCCA){
     return Rcpp::List::create(Rcpp::Named("Xtil")=Xtil,
                               Rcpp::Named("U0")=U0,
                               // Rcpp::Named("K")=selK,
                               Rcpp::Named("corr")=eigs,
                               Rcpp::Named("A")=A,
                               Rcpp::Named("B")=B,
                               Rcpp::Named("U")=U,
                               Rcpp::Named("V")=V);
   }
   return Rcpp::List::create(Rcpp::Named("Xtil")=Xtil,
                             // Rcpp::Named("K")=selK,
                             Rcpp::Named("U0")=U0);

 };


// TODO : autorPACA, rPACA



