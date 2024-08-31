#include "pacabase.h"
#include "utils.h"

// [[Rcpp::plugins(cpp14)]]
// [[Rcpp::depends(RcppEigen)]]
using namespace std;
using namespace Rcpp;

////////////////////////////////// Loggers //////////////////////////////////

// Logger static member initialization
std::mutex Logger::mtx;
int Logger::verbosity = 1;

void Logger::SetVerbosity(int level) {
    std::lock_guard<std::mutex> lock(mtx);
    verbosity = level;
}

// Timer implementation
Timer::Timer(const std::string& fname) : m_fname(fname) {
    m_startTimept = std::chrono::high_resolution_clock::now();
}

Timer::~Timer() {
    Stop();
}

void Timer::Stop() {
    auto m_endTimept = std::chrono::high_resolution_clock::now();
    auto start = std::chrono::time_point_cast<std::chrono::microseconds>(m_startTimept).time_since_epoch().count();
    auto end = std::chrono::time_point_cast<std::chrono::microseconds>(m_endTimept).time_since_epoch().count();
    auto duration = end - start;
    std::ostringstream LOCAL_oss;
    LOCAL_oss << "\t{" << m_fname << "} took: " << duration * 0.001 << "ms / " << duration << "us";

    Logger::LogDEBUG(LOCAL_oss.str());
}

////////////////////////////////// General Helpers //////////////////////////////////

// [[Rcpp::export(normalizeCPP)]]
Eigen::MatrixXd normalizeCPP(Eigen::MatrixXd& x, bool inplace) {
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
  };
  x.rowwise() -= xmeans.transpose();
  x.array().rowwise() /= xstd.transpose().array();
  return x;
};

Eigen::MatrixXd scaleCPP(const Eigen::MatrixXd& x) {
    // Calculate the mean of each column
    Eigen::VectorXd means = x.colwise().mean();
    
    // Calculate the variance of each column
    Eigen::VectorXd variances = (x.rowwise() - means.transpose()).array().square().colwise().mean();
    
    // Calculate the standard deviation of each column
    Eigen::VectorXd stdDevs = variances.array().sqrt();
    
    // Create a scaled matrix
    Eigen::MatrixXd scaled = x.rowwise() - means.transpose();
    
    // Avoid division by zero
    for (int i = 0; i < stdDevs.size(); ++i) {
        if (stdDevs(i) > 1e-9) {  // Threshold to avoid division by very small numbers
            scaled.col(i).array() /= stdDevs(i);
        } else {
            // If standard deviation is very close to zero, set the column to zero
            scaled.col(i).setZero();
        }
    }
    
    return scaled;
}

Eigen::MatrixXd normalizeColEigen(const Eigen::MatrixXd& x) {
  Eigen::MatrixXd xcp = x;
  for (int i = 0; i < xcp.cols(); i++){
    xcp.col(i).normalize();
  }
  return xcp;
};

// Eigen::MatrixXd centerCPP(Eigen::MatrixXd& x) {
//   Eigen::VectorXd mean = x.colwise().mean();
//   return x.rowwise() - mean.transpose();
// }

Eigen::MatrixXd centerCPP(const Eigen::MatrixXd& x) {
  Eigen::VectorXd mean = x.colwise().mean();
  return x.rowwise() - mean.transpose();
};

Eigen::MatrixXd covCPP(const Eigen::MatrixXd& x) {
  Eigen::MatrixXd centered = centerCPP(x);
  return (centered.adjoint() * centered) / double(x.rows() - 1);
};

Eigen::MatrixXd covSymCPP(const Eigen::MatrixXd& x) {
  Eigen::MatrixXd centered = centerCPP(x);
  Eigen::MatrixXd cv_tmp = (centered.adjoint() * centered) / double(x.rows() - 1);
  Logger::LogDEBUG("SymCov");
  return (cv_tmp + cv_tmp.transpose()) / 2.0;
};

Eigen::VectorXd colwiseVars(const Eigen::MatrixXd& matrix) {
  // Number of rows (samples)
  int rows = matrix.rows();

  // Check if there is enough data to compute variance
  if (rows < 2) {
    throw std::invalid_argument("Not enough data to compute variance");
  };

  // Calculate the means of each column
  Eigen::VectorXd mean = matrix.colwise().mean();

  // Subtract the mean from each column and square the result
  Eigen::MatrixXd centered = matrix.rowwise() - mean.transpose();
  centered = centered.array().square();

  // Compute the variance of each column
  Eigen::VectorXd variances = (centered.colwise().sum()) / (rows - 1);

  return variances;
};
