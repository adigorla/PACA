#ifndef PACA_BASE_H
#define PACA_BASE_H

// includes to make Eigen use BLAS
#define EIGEN_SUPERLU_SUPPORT
#define EIGEN_USE_BLAS
// #define EIGEN_USE_LAPACKE

#include <Rcpp.h>
#include <RcppEigen.h>
#include <vector>
#include <chrono>
#include <random>
#include <algorithm>
#include <mutex>
#include <Eigen/Dense>
#include <complex>

#endif // PACA_BASE_H