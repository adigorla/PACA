#include "pacabase.h"
#include "utils.h"
#include "pacautils.h"
#include "cca.h"
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

Eigen::MatrixXd scaleCPP_batchmasked(const Eigen::MatrixXd& x, const Eigen::VectorXi& batch_idx) {
    Eigen::MatrixXd scaled = Eigen::MatrixXd::Zero(x.rows(), x.cols());
    
    // Create a mask for rows not in the batch
    Eigen::VectorXd mask = Eigen::VectorXd::Ones(x.rows());
    for (int i = 0; i < batch_idx.size(); ++i) {
        mask(batch_idx(i)) = 0;
    }
    // Count out-of-batch elements
    double out_batch_count = mask.sum();
    
    for (int col = 0; col < x.cols(); ++col) {
        // Extract the out-of-batch elements for this column
        Eigen::VectorXd col_out_batch = x.col(col).array() * mask.array();
        
        // Compute mean of out-of-batch elements
        double mean = col_out_batch.sum() / out_batch_count;
        
        // Compute variance of out-of-batch elements (corrected)
        double variance = ((col_out_batch.array() - mean).square() * mask.array()).sum() / (out_batch_count - 1);
        
        // Compute standard deviation
        double std_dev = std::sqrt(variance);
        
        // Scale the entire column
        if (std_dev > 1e-9) {  // Avoid division by very small numbers
            // Assuming the masked values are 0
            scaled.col(col) = ((x.col(col).array() - mean) * mask.array()) / std_dev;
            
            // scaled.col(col) = ((x.col(col).array() - mean) * mask.array()) / std_dev + x.col(col).array() * (1 - mask.array());
            // scaled.col(col) = (x.col(col).array() - mean) / std_dev;
        }
    }
    
    return scaled;
}

void residMatrix(Eigen::MatrixXd& X, const Eigen::VectorXd& PAC) {
    int n = X.rows();
    int p = X.cols();
    
    // Create V matrix
    Eigen::MatrixXd V(p, 2);
    V.col(0) = Eigen::VectorXd::Ones(p);
    V.col(1) = PAC;

    // // Calculate gamma
    // Eigen::MatrixXd gamma = V.transpose() * V;
    // Logger::LogDEBUG("gamma shape: ", gamma.rows(), "x", gamma.cols());
    // std::cout << gamma << std::endl;
    // double cond = gamma.jacobiSvd().singularValues()(0) / gamma.jacobiSvd().singularValues()(gamma.cols() - 1);
    // Logger::LogDEBUG("Condition number of gamma: ", cond);
    
    // Calculate Xtil
    Eigen::MatrixXd Xbeta = Eigen::MatrixXd::Zero(n, p);
    
    // Solve for beta and calculate Xtil for each feature
    Eigen::VectorXd beta_j(2);
    Eigen::VectorXd X_row(p);
    for (int j = 0; j < n; ++j) {
        X_row = X.row(j);
        beta_j = (V.transpose() * V).ldlt().solve(V.transpose() * X_row);
        Xbeta.row(j) = (V * beta_j).transpose();
    }
    // DEBUG
    if (!Xbeta.allFinite()) {
            throw std::runtime_error("Non-finite values detected in Xbeta");
    };

    // Update X: to residual of X
    X = X - Xbeta;
    
    // Center and scale X
    X = scaleCPP(X);
    Logger::LogDEBUG("Done w/ resid");
}


void singlearpaca_pvt(const Eigen::MatrixXd& X, const Eigen::MatrixXd& Y, 
                     Eigen::MatrixXd& iterPAC, Eigen::VectorXd& iterEIG, Eigen::VectorXd& klist,
                     double threshold, const int niter, const int batch,
                     const int rank, bool normalize) {
    // int m = X.rows();
    int p = X.cols();
    int q = Y.cols();
    

    Eigen::MatrixXd P = Eigen::MatrixXd::Zero(p, niter * rank);

    // Initialize the random number generator
    unsigned int seed = static_cast<unsigned int>(unif_rand() * UINT_MAX);
    Logger::LogDEBUG("setting seed :", seed);
    std::mt19937 gen(seed);

    for (int r = 0; r < niter; ++r) {
        Logger::LogINFO("Running randomized PACA Iteration ", r + 1, " of ", niter);
        Eigen::MatrixXd P_curr = Eigen::MatrixXd::Zero(p, rank);
        int selK;
        {
            Timer timer("Iteration " + std::to_string(r + 1));
            // Sample batch indices for X and Y
            Eigen::VectorXi Xin_idx = sampleWithoutReplacement(p, batch, gen);
            Eigen::VectorXi Yin_idx = sampleWithoutReplacement(q, batch, gen);

            // Extract batches
            Eigen::MatrixXd X_in = X(Eigen::all, Xin_idx);
            Eigen::MatrixXd Y_in = Y(Eigen::all, Yin_idx);

            // Run PACA on the batch
            Eigen::MatrixXd Xtil_in, U0;
            {
                Logger::LogLOG("Subsample PACA: Starting ...");
                Eigen::MatrixXd A, B, U, V;
                Eigen::VectorXd eigs;
                Eigen::MatrixXd X_in_copy = X_in;

                // Run PACA
                selectK_pvt(X_in, Y_in, selK, normalize, eigs, A, B, U, V, threshold);
                klist(r) = selK;

                Logger::LogINFO("Estimating within iter PACA signals ...");
                Logger::LogLOG("Residualizing Shared Signal ...");

                // Calculate shared components
                int dU1 = U.rows();
                U0 = UV1_calc(U, selK, dU1);

                // Calcute X matrix with case specific variation only
                Xtil_in = correctedMat_calc(U0, X_in, false);

                // Clearing Temporary Matrices with Safer Operations
                A.resize(0, 0);
                B.resize(0, 0);
                U.resize(0, 0);
                V.resize(0, 0);
                eigs.resize(0);
                X_in_copy.resize(0, 0);
            }

            // Perform PCA on Xtil_in
            Eigen::MatrixXd rotation, Xpac_in;
            {
                Timer timer("Subsample PCA");
                Eigen::MatrixXd XtilT_norm = scaleCPP(std::move(Xtil_in.transpose().eval()));


                Eigen::BDCSVD<Eigen::MatrixXd> svd(XtilT_norm, Eigen::ComputeThinU | Eigen::ComputeThinV);
                rotation = svd.matrixV().leftCols(rank);

                // Calculate in-batch PACA PCs
                Xpac_in = svd.matrixU().leftCols(rank);
               
            }

            // Project remaining data with zeroed-out selected columns
            Eigen::MatrixXd X_out = X;  // Copy original matrix
            // Zero out the selected columns using Eigen's indexed access
            X_out(Eigen::all, Xin_idx) = Eigen::MatrixXd::Zero(X.rows(), batch);


            // Extract out of batch paca components
            Eigen::MatrixXd Xtil_out = correctedMat_calc(U0, X_out, false);
            Eigen::MatrixXd Xpac_out = Xtil_out.transpose() * rotation;
            

            // Normalize projections
            Xpac_in = scaleCPP(Xpac_in);
            P_curr = scaleCPP_batchmasked(Xpac_out, Xin_idx);


            // Clear temporary matrices
            X_in.resize(0, 0);
            Y_in.resize(0, 0);
            X_out.resize(0, 0);
            Xtil_in.resize(0, 0);
            Xtil_out.resize(0, 0);
            U0.resize(0, 0);
            rotation.resize(0, 0);
            Xpac_out.resize(0, 0);


            // Combine projections
            P_curr(Xin_idx, Eigen::all) = Xpac_in;


            // Clear temporary matrices
            Xpac_in.resize(0, 0);
            Xin_idx.resize(0);
            Yin_idx.resize(0);

        }

        // Update P
        P.block(0, r * rank, p, rank) = P_curr;

        // Clear temporary matrices
        P_curr.resize(0, 0);
        
        Logger::LogLOG("Completed Random Sampling Iter: ", r + 1, " / ", niter);
    }

    // Calculate covariance matrix
    Eigen::MatrixXd C = covCPP(P.transpose()) ;


    // Compute final eigenvectors and eigenvalues
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(C);

    iterPAC = eigensolver.eigenvectors().rightCols(rank);
    iterEIG = eigensolver.eigenvalues().tail(rank).reverse();
    Logger::LogINFO("Randomized PACA: DONE");
};


////////////////////////////////// Main Randomized PACA functions //////////////////////////////////

//' @export
//' @noRd
// [[Rcpp::export(cpp_rPACA)]]
Rcpp::List cpp_rPACA(const Eigen::MatrixXd& X, const Eigen::MatrixXd& Y,
                     int k, int niter, int batch, int rank,
                     bool normalize, int verbosity) {
    Logger::SetVerbosity(verbosity);
    Timer timer("rPACA");
    Logger::LogINFO("Randomized PACA (with K=", k, ", subsample size=", batch, ") : Starting ...");

    // int m = X.rows();
    int p = X.cols();
    int q = Y.cols();

    Eigen::MatrixXd P = Eigen::MatrixXd::Zero(p, niter * rank);

    // Synchronize R's RNG state
    GetRNGstate();
    // Initialize the random number generator
    unsigned int seed = static_cast<unsigned int>(unif_rand() * UINT_MAX);
    Logger::LogDEBUG("setting seed :", seed);
    std::mt19937 gen(seed);

    for (int r = 0; r < niter; ++r) {
        Logger::LogINFO("Iteration ", r + 1, " of ", niter);
        Eigen::MatrixXd P_curr = Eigen::MatrixXd::Zero(p, rank);
        {
            Timer timer("Iteration " + std::to_string(r + 1));
            // Sample batch indices for X and Y
            Eigen::VectorXi Xin_idx = sampleWithoutReplacement(p, batch, gen);
            Eigen::VectorXi Yin_idx = sampleWithoutReplacement(q, batch, gen);

            // Extract batches
            // Eigen::MatrixXd X_in(m, batch);
            // Eigen::MatrixXd Y_in(m, batch);
            // for (int i = 0; i < batch; ++i) {
            //     X_in.col(i) = X.col(Xin_idx(i));
            //     Y_in.col(i) = Y.col(Yin_idx(i));
            // }
            Eigen::MatrixXd X_in = X(Eigen::all, Xin_idx);
            Eigen::MatrixXd Y_in = Y(Eigen::all, Yin_idx);

            // Run PACA on the batch
            Eigen::MatrixXd Xtil_in, U0;
            {
                Logger::LogLOG("Subsample PACA: Starting ...");
                Eigen::MatrixXd A, B, U, V;
                Eigen::VectorXd eigs;
                Eigen::MatrixXd X_in_copy = X_in;
                doCCA_pvt(X_in, Y_in, normalize, eigs, A, B, U, V);
                Logger::LogLOG("Done with CCA");

                Logger::LogLOG("Residualizing Shared Signal ...");

                // Calculate shared components
                int dU1 = U.rows();
                U0 = UV1_calc(U, k, dU1);

                // Calcute X matrix with case specific variation only
                Xtil_in = correctedMat_calc(U0, X_in, false);

                // Clearing Temporary Matrices with Safer Operations
                A.resize(0, 0);
                B.resize(0, 0);
                U.resize(0, 0);
                V.resize(0, 0);
                eigs.resize(0);
                X_in_copy.resize(0, 0);
            }

            // Perform PCA on Xtil_in
            Eigen::MatrixXd rotation, Xpac_in;
            {
                Timer timer("Subsample PCA");
                Eigen::MatrixXd XtilT_norm = scaleCPP(std::move(Xtil_in.transpose().eval()));

                // //// DEBUG
                // std::cout << "Xtil_in shape: " << Xtil_in.rows() << "x" << Xtil_in.cols() << "\n";
                // std::cout << "Xtil_Tnorm shape: " << XtilT_norm.rows() << "x" << XtilT_norm.cols() << "\n";
                // // Calculate and print means
                // Eigen::VectorXd meansxt = XtilT_norm.colwise().mean();
                // std::cout << "Xtil_Tnorm means:\n" << meansxt.head(10).transpose() << "\n\n";
                // // Calculate and print variances
                // Eigen::VectorXd variancesxt = (XtilT_norm.rowwise() - meansxt.transpose()).array().square().colwise().mean();
                // std::cout << "Xtil_Tnorm variances:\n" << variancesxt.head(10).transpose() << "\n\n";
                // ////


                Eigen::BDCSVD<Eigen::MatrixXd> svd(XtilT_norm, Eigen::ComputeThinU | Eigen::ComputeThinV);
                // std::cout << "SVD U shape: " << svd.matrixU().rows() << "x" << svd.matrixU().cols() << "\n";
                // std::cout << "SVD V shape: " << svd.matrixV().rows() << "x" << svd.matrixV().cols() << "\n";

                rotation = svd.matrixV().leftCols(rank);
                // std::cout << "rotation shape: " << rotation.rows() << "x" << rotation.cols() << "\n";
                // std::cout << "X_in BEFORE SVD (top 10 rows, left 5 cols):\n" << X_in.topLeftCorner(std::min(10, (int)X_in.rows()), std::min(5, (int)X_in.cols())) << "\n\n";

                // //// DEBUG
                // // Create a vector of indices
                // std::vector<int> indices(Xin_idx.size());
                // std::iota(indices.begin(), indices.end(), 0);

                // // Sort the indices based on the values in Xin_idx
                // std::sort(indices.begin(), indices.end(),
                //     [&Xin_idx](int i1, int i2) { return Xin_idx(i1) < Xin_idx(i2); });

                // // Take the first 5 (or less if Xin_idx has fewer than 5 elements)
                // int num_to_print = std::min(5, (int)Xin_idx.size());

                // // Create a matrix to store the selected rows for XtilT_norm (only first 10 columns)
                // int num_cols_to_print = std::min(5, (int)XtilT_norm.cols());
                // Eigen::MatrixXd selected_rows_XtilT_norm(num_to_print, num_cols_to_print);

                // // Fill the matrix with the selected rows for XtilT_norm (only first 10 columns)
                // for (int i = 0; i < num_to_print; ++i) {
                //     selected_rows_XtilT_norm.row(i) = XtilT_norm.row(Xin_idx(indices[i])).head(num_cols_to_print);
                // }
                
                // // Print the selected rows for XtilT_norm
                // std::cout << "XtilT_norm BEFORE transformation (5 rows with smallest Xin_idx values, first 5 cols):\n" 
                //         << selected_rows_XtilT_norm << "\n\n";
                // ////

                // Calculate in-batch PACA PCs
                // Eigen::MatrixXd S = svd.singularValues().head(rank).asDiagonal();
                // Xpac_in = (XtilT_norm * rotation).eval();
                Xpac_in = svd.matrixU().leftCols(rank);

                // //// DEBUG
                // std::cout << "Xpac_in shape: " << Xpac_in.rows() << "x" << Xpac_in.cols() << "\n";

                // // Calculate and print means
                // Eigen::VectorXd means = Xpac_in.colwise().mean();
                // std::cout << "Xpac_in means:\n" << means.transpose() << "\n\n";

                // // Calculate and print variances
                // Eigen::VectorXd variances = (Xpac_in.rowwise() - means.transpose()).array().square().colwise().mean();
                // std::cout << "Xpac_in variances:\n" << variances.transpose() << "\n\n";
                
                // // Create a matrix to store the selected rows for Xpac_in
                // Eigen::MatrixXd selected_rows_Xpac_in(num_to_print, rank);

                // // Fill the matrix with the selected rows for Xpac_in
                // for (int i = 0; i < num_to_print; ++i) {
                //     selected_rows_Xpac_in.row(i) = Xpac_in.row(indices[i]);
                // }

                // // Print the selected rows for Xpac_in
                // std::cout << "Xpac_in creation (5 rows with smallest Xin_idx values, all cols):\n" 
                //         << selected_rows_Xpac_in << "\n\n";
                // std::cout << "Xpac_in creation (top 10 rows, left 5 cols):\n" << Xpac_in.topLeftCorner(std::min(10, (int)Xpac_in.rows()), std::min(5, (int)Xpac_in.cols())) << "\n\n";
                // ////
               
            }

            // Project remaining data with zeroed-out selected columns
            Eigen::MatrixXd X_out = X;  // Copy original matrix
            // for (int i = 0; i < batch; ++i) {
            //     X_out.col(Xin_idx(i)).setZero();
            // }
            // Zero out the selected columns using Eigen's indexed access
            X_out(Eigen::all, Xin_idx) = Eigen::MatrixXd::Zero(X.rows(), batch);


            Eigen::MatrixXd Xtil_out = correctedMat_calc(U0, X_out, false);
            Eigen::MatrixXd Xpac_out = Xtil_out.transpose() * rotation;
            // std::cout << "Xpac_out shape: " << Xpac_out.rows() << "x" << Xpac_out.cols() << "\n";
            // std::cout << "Xpac_out BEFORE scale & Insert (top 10 rows, left 5 cols):\n" << Xpac_out.topLeftCorner(std::min(10, (int)Xpac_out.rows()), std::min(5, (int)Xpac_out.cols())) << "\n\n";

            // Normalize projections
            Xpac_in = scaleCPP(Xpac_in);
            P_curr = scaleCPP_batchmasked(Xpac_out, Xin_idx);
            // std::cout << "Xpac_out AFTER scale & BEFORE Insert (top 10 rows, left 5 cols):\n" << P_curr.topLeftCorner(std::min(10, (int)P_curr.rows()), std::min(5, (int)P_curr.cols())) << "\n\n";

            // Clear temporary matrices
            X_in.resize(0, 0);
            Y_in.resize(0, 0);
            X_out.resize(0, 0);
            Xtil_in.resize(0, 0);
            Xtil_out.resize(0, 0);
            U0.resize(0, 0);
            rotation.resize(0, 0);
            Xpac_out.resize(0, 0);
            
            // //// DEBUG
            // // Create a vector of indices
            // std::vector<int> indices(Xin_idx.size());
            // std::iota(indices.begin(), indices.end(), 0);

            // // Sort the indices based on the values in Xin_idx
            // std::sort(indices.begin(), indices.end(),
            //     [&Xin_idx](int i1, int i2) { return Xin_idx(i1) < Xin_idx(i2); });

            // // Take the first 5 (or less if Xin_idx has fewer than 5 elements)
            // int num_to_print = std::min(5, (int)Xin_idx.size());

            // // Create a matrix to store the selected rows
            // Eigen::MatrixXd selected_rows(num_to_print, Xpac_in.cols());

            // // Fill the matrix with the selected rows
            // for (int i = 0; i < num_to_print; ++i) {
            //     // selected_rows.row(i) = Xpac_in.row(indices[i]).transpose();
            //     selected_rows.row(i) = Xpac_in.row(indices[i]);
            // }
            // // Print the selected rows
            // std::cout << "Xpac_in BEFORE Insert (5 rows with smallest Xin_idx values, all cols):\n" 
            //         << selected_rows << "\n\n";
            // Logger::LogDEBUG("P_curr shape: ", P_curr.rows(), "x", P_curr.cols());
            // ////

            // Combine projections
            // for (int i = 0; i < batch; ++i) {
            //     P_curr.row(Xin_idx(i)) = Xpac_in.row(i); // Assigning rows directly
            // }
            P_curr(Xin_idx, Eigen::all) = Xpac_in;
            // std::cout << "P_curr AFTER Insert (top 10 rows, left 5 cols):\n" << P_curr.topLeftCorner(std::min(10, (int)P_curr.rows()), std::min(5, (int)P_curr.cols())) << "\n\n";
           
            // Clear temporary matrices
            Xpac_in.resize(0, 0);
            Xin_idx.resize(0);
            Yin_idx.resize(0);

        }
        // //// DEBUG
        // // Calculate and print means
        // Eigen::VectorXd meansPcrr = P_curr.colwise().mean();
        // std::cout << "P_curr means:\n" << meansPcrr.transpose() << "\n\n";

        // // Calculate and print variances
        // Eigen::VectorXd variancesPcrr = (P_curr.rowwise() - meansPcrr.transpose()).array().square().colwise().mean();
        // std::cout << "P_curr variances:\n" << variancesPcrr.transpose() << "\n\n";
        
        // // Debug print for P_curr
        // std::cout << "P_curr (top 10 rows, left 5 cols):\n" << P_curr.topLeftCorner(std::min(10, (int)P_curr.rows()), std::min(5, (int)P_curr.cols())) << "\n\n";
        // ////

        // Update P
        P.block(0, r * rank, p, rank) = P_curr;

        // Clear temporary matrices
        P_curr.resize(0, 0);
        
        Logger::LogLOG("Completed Random Sampling Iter: ", r + 1, " / ", niter);
    }
    // Debug print for P
    // std::cout << "P (top 10 rows, left 5 cols):\n" << P.topLeftCorner(std::min(10, (int)P.rows()), std::min(5, (int)P.cols())) << "\n\n";

    // Synchronize the RNG state back to R
    PutRNGstate();
    // Calculate covariance matrix
    Eigen::MatrixXd C = covCPP(P.transpose()) ;
    // //// DEBUG
    // std::cout << "C shape: " << C.rows() << "x" << C.cols() << "\n";
    // std::cout << "First 10 diagonal elements of C: " << C.diagonal().head(10).transpose() << "\n\n";
    // std::cout << "C (top 10 rows, left 5 cols):\n" << C.topLeftCorner(std::min(10, (int)C.rows()), std::min(5, (int)C.cols())) << "\n\n";
    // ////

    // Compute final eigenvectors and eigenvalues
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(C);

    Logger::LogINFO("Randomized PACA: DONE");
    return Rcpp::List::create(
        Rcpp::Named("x") = eigensolver.eigenvectors().rightCols(rank),
        Rcpp::Named("eigs") = eigensolver.eigenvalues().tail(rank).reverse()
    );
};

//' @export
//' @noRd
// [[Rcpp::export(cpp_autorPACA)]]
Rcpp::List cpp_autorPACA(const Eigen::MatrixXd& X, const Eigen::MatrixXd& Y,
                     int niter, int batch, int rank,
                     double threshold,
                     bool normalize, int verbosity) {

    Logger::SetVerbosity(verbosity);
    Timer timer("rPACA");
    Logger::LogINFO("Auto Randomized PACA (with subsample size=", batch, ") : Starting ...");

    // Create local copies of X and Y to avoid modifying the original matrices
    Eigen::MatrixXd Xcurr = X;
    Eigen::MatrixXd Ycurr = Y;

    // Create the output matrices
    Eigen::MatrixXd PACfinal = Eigen::MatrixXd::Zero(X.cols(), rank);
    Eigen::VectorXd EIGfinal = Eigen::VectorXd::Zero(rank);
    Eigen::VectorXd Klist = Eigen::VectorXd::Zero(niter);

    // Synchronize R's RNG state
    GetRNGstate();

    // Perform randomized PACA
    singlearpaca_pvt(Xcurr, Ycurr, PACfinal, EIGfinal, Klist,
                    threshold, niter, batch, rank, normalize);

    // PACfinal.col(pr) = iterPAC.col(0);
    // EIGfinal(pr) = iterEIG(0);
    // // Upadate X by residualizing out the current PAC
    // residMatrix(Xcurr, iterPAC.col(0));
    // // residMatrix(Ycurr, iterPAC.col(0));

    // Synchronize the RNG state back to R
    PutRNGstate();

    Logger::LogINFO("Auto Randomized PACA: DONE");
    return Rcpp::List::create(
        Rcpp::Named("x") = PACfinal,
        Rcpp::Named("eigs") = EIGfinal,
        Rcpp::Named("k.iter") = Klist
    );

};