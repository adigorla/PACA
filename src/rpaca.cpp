#include "pacabase.h"
#include "utils.h"
#include "pacautils.h"
#include "cca.h"
#include "pca.h"
#include "paca.h"
#include "rpaca.h"

// [[Rcpp::plugins(cpp14)]]
// [[Rcpp::depends(RcppEigen)]]

using namespace std;
using namespace Rcpp;


////////////////////////////////// Helper //////////////////////////////////
// Helper function to sample indices without replacement
// Helper function to sample indices without replacement
Eigen::VectorXi sampleWithoutReplacement(int n, int k, std::mt19937& gen) {
    std::vector<int, Eigen::aligned_allocator<int>> indices(n);
    std::iota(indices.begin(), indices.end(), 0);
    std::shuffle(indices.begin(), indices.end(), gen);
    Eigen::VectorXi result(k);
    std::copy(indices.begin(), indices.begin() + k, result.data());
    return result;
}

// TODO: write a more efficient getter/setter since we have alot of non-contiguous memory access

////////////////////////////////// Main Randomized PACA functions //////////////////////////////////

//' @export
//' @noRd
// [[Rcpp::export(cpp_rPACA)]]
Rcpp::List cpp_rPACA(const Eigen::MatrixXd& X, const Eigen::MatrixXd& Y,
                     int k, int niter, int batch, int rank,
                     bool normalize, int verbosity) {
    Logger::SetVerbosity(verbosity);
    Timer timer("rPACA");
    Logger::LogINFO("Randomized PACA with K=", k, ", subsample size=", batch, ": Starting ...");

    int m = X.rows();
    int p = X.cols();
    int q = Y.cols();

    Eigen::MatrixXd P = Eigen::MatrixXd::Zero(p, niter * rank);

    // Synchronize R's RNG state
    GetRNGstate();
    // Initialize the random number generator
    unsigned int seed = static_cast<unsigned int>(unif_rand() * UINT_MAX);
    Logger::LogDEBUG("setting seed :", seed);
    std::mt19937 gen(seed);

    // Create a vector of all row indices
    Eigen::VectorXi all_rows = Eigen::VectorXi::LinSpaced(m, 0, m - 1);

    for (int r = 0; r < niter; ++r) {
        Logger::LogINFO("Iteration ", r + 1, " of ", niter);
        Eigen::MatrixXd P_curr = Eigen::MatrixXd::Zero(p, rank);
        {
            Timer timer("Iteration " + std::to_string(r + 1));
            // Sample batch indices for X and Y
            Eigen::VectorXi inX = sampleWithoutReplacement(p, batch, gen);
            Eigen::VectorXi inY = sampleWithoutReplacement(q, batch, gen);

            // Extract batches
            Eigen::MatrixXd xBatch(m, batch);
            Eigen::MatrixXd yBatch(m, batch);
            for (int i = 0; i < batch; ++i) {
                xBatch.col(i) = X.col(inX(i));
                yBatch.col(i) = Y.col(inY(i));
            }

            // Run PACA on the batch
            Eigen::MatrixXd Xtil, U0;
            {
                Logger::LogLOG("Subsample PACA: Starting ...");
                Eigen::MatrixXd A, B, U, V;
                Eigen::VectorXd eigs;
                Eigen::MatrixXd xBatch_copy = xBatch;
                doCCA_pvt(xBatch, yBatch, normalize, eigs, A, B, U, V);
                Logger::LogLOG("Done with CCA");

                Logger::LogLOG("Residualizing Shared Signal ...");

                // Calculate shared components
                int dU1 = U.rows();
                U0 = UV1_calc(U, k, dU1);

                // Calcute X matrix with case specific variation only
                Xtil = correctedMat_calc(U0, xBatch, false);

                // Clearing Temporary Matrices with Safer Operations
                A.resize(0, 0);
                B.resize(0, 0);
                U.resize(0, 0);
                V.resize(0, 0);
                eigs.resize(0);
                xBatch_copy.resize(0, 0);
            }

            // Perform PCA on Xtil
            Eigen::MatrixXd rotation, pacapcBatch_x;
            {
                Timer timer("Subsample PCA");
                Eigen::BDCSVD<Eigen::MatrixXd> svd(Xtil.transpose(), Eigen::ComputeThinU | Eigen::ComputeThinV);
                rotation = svd.matrixV().leftCols(rank);
                pacapcBatch_x = xBatch.transpose() * svd.matrixV().leftCols(rank);
               
            }

            // Project remaining data with zeroed-out selected columns
            Eigen::MatrixXd xRemain = X;  // Copy original matrix
            // Create an index vector to select columns
            for (int i = 0; i < batch; ++i) {
                xRemain.col(inX(i)).setZero();
            }

            Eigen::MatrixXd X_tilde_remain = correctedMat_calc(U0, xRemain, false);
            P_curr = X_tilde_remain.transpose() * rotation;

            // Normalize projections
            pacapcBatch_x.colwise().normalize();
            P_curr.colwise().normalize();

            // Combine projections
            // Logger::LogDEBUG("P_curr shape: ", P_curr.rows(), "x", P_curr.cols());
            // Logger::LogDEBUG("pacapcBatch_x shape: ", pacapcBatch_x.rows(), "x", pacapcBatch_x.cols());
            for (int i = 0; i < batch; ++i) {
                P_curr.row(inX(i)) = pacapcBatch_x.col(i).transpose();
            }

            // Clear temporary matrices
            xBatch.resize(0, 0);
            yBatch.resize(0, 0);
            Xtil.resize(0, 0);
            U0.resize(0, 0);
            rotation.resize(0, 0);
            pacapcBatch_x.resize(0, 0);
            xRemain.resize(0, 0);
            X_tilde_remain.resize(0, 0);
            inX.resize(0);
            inY.resize(0);

        }

        P.block(0, r * rank, p, rank) = P_curr;

        // Clear temporary matrices
        // P_curr.setZero();
        
        Logger::LogLOG("Completed Random Sampling Iter: ", r + 1, " / ", niter);
    }
    // Synchronize the RNG state back to R
    PutRNGstate();
    // Compute final eigenvectors and eigenvalues
    Eigen::MatrixXd C = P * P.transpose();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(C);

    Logger::LogINFO("Randomized PACA: DONE");
    return Rcpp::List::create(
        Rcpp::Named("x") = eigensolver.eigenvectors().rightCols(rank),
        Rcpp::Named("eigs") = eigensolver.eigenvalues().tail(rank).reverse()
    );
}