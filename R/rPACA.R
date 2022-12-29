#' @title Randomized Phenotype Aware Components Analysis (rPACA)
#'
#' @description
#' This is the randomized version of the basic \code{\link{PACA}} algorithm.
#' rPACA allows us to apply PACA in regimes where p << n, i.e., in cases where the number of
#' samples is greater than the number of features.
#'
#' @param X \eqn{n_1} by \eqn{m} matrix; \cr
#'          Case (foreground) input data matrix. \cr
#'          It is recommended to normailize the feature scales as appropriate for the data modality.
#'          E.g. quantile normalization (or other comparable approaches) for RNAseq data.
#' @param Y \eqn{n_0} by \eqn{m} matrix; \cr
#'          Control (foreground) input data matrix. \cr
#'          It is recommended to normailize the feature scales as appropriate for the data modality.
#'          E.g. quantile normalization (or other comparable approaches) for RNAseq data.
#' @param k Positive integer, \eqn{k > 1}; \cr
#'          Number of, \eqn{k}, dimensions of shared variation to be removed from case data \code{X}.
#' @param nIter Positive integer; \cr
#'              Number random sampling iterations.
#' @param stepSize Positive integer; \cr
#'                 Number of cases/controls to be sampled for the training set at each iteration. \cr
#'                 Note that \eqn{k < stepSize-1} and \eqn{stepSize < min{m, n1, n0}}.
#' @param pcRank Positive integer, optional (default \eqn{2}); \cr
#'               Number of dominant principle components to be computed for the corrected case data.
#'
#' @return \code{autoPACA} returns a list containing the following components:
#' \describe{
#'    \item{x}{        \eqn{n_1} by \eqn{pcRank} matrix; \cr
#'                     the projections / scores of the randomized \emph{PACA} corrected case data.
#'    }
#'
#'    \item{eigVal}{   list of size \eqn{pcRank}; \cr
#'                     the eigen values of the returned PC components.
#'    }
#'}
#' @export
#'
#' @importFrom stats prcomp
rPACA <- function(X, Y, k , nIter = 20, stepSize = 600, pcRank = 4){
  X_pacacomb <- matrix(0, ncol=pcRank, nrow=dim(X)[2])
  X_ids <- colnames(X)
  Y_ids <- colnames(Y)

  P_comb <- matrix(0, ncol=1, nrow=dim(X)[2])
  for (r in 1:nIter){
    inX <- sample(X_ids, size=stepSize, replace=F)
    inY <- sample(Y_ids, size=stepSize, replace=F)
    xBatch <- X[,inX]
    yBatch <- Y[,inY]
    xRemain <- X[,setdiff(X_ids, inX)]

    pacaBatch <- CCA(xBatch, yBatch)

    U_1 <- pacaBatch$U[,1:k] / t(kronecker(matrix(1,1,dim(pacaBatch$U)[1]),sqrt(colSums(pacaBatch$U[,1:k]^2))))
    means_matrix <- t(kronecker(matrix(1,1,dim(xBatch)[1]),colMeans(xBatch)))
    X_centered <- xBatch - means_matrix
    X_tilde <- X_centered - (U_1 %*% (t(U_1) %*% X_centered)) + means_matrix
    rm(means_matrix, X_centered)
    pacapcBatch <- prcomp(t(X_tilde), rank. = pcRank, center = TRUE, scale. = T)

    means_matrix <- t(kronecker(matrix(1,1,dim(xRemain)[1]),colMeans(xRemain)))
    X_centered <- xRemain - means_matrix
    X_tilde_remain <- X_centered - (U_1 %*% (t(U_1) %*% X_centered)) + means_matrix
    rm(means_matrix, X_centered, U_1)
    projMat_remain <- t(X_tilde_remain) %*% pacapcBatch$rotation

    pacapcBatch <- scale(pacapcBatch$x , center = TRUE, scale = T)
    projMat_remain <- scale(projMat_remain , center = TRUE, scale = T)
    comb_proj <- rbind(pacapcBatch, projMat_remain)
    P_comb <- cbind(P_comb, comb_proj[X_ids,])
    X_pacacomb <- comb_proj[X_ids,] + X_pacacomb
    rm(comb_proj, X_tilde_remain, X_tilde)

    cat("\n\tCompleted Random Sampling Iter:", r,"/",nIter,"\n")
  }
  P_comb <- P_comb[,-1]
  Cmat <- P_comb %*% t(P_comb)
  eigPC <- eigen(Cmat)

  return(list(x = eigPC$vectors[,1:pcRank],
              eigVal = eigPC$values[,1:pcRank]))

}
