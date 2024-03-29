#' @title Randomized Phenotype Aware Components Analysis (rPACA)
#'
#' @description
#' This is the randomized version of the basic \code{\link{PACA}} algorithm.
#' rPACA allows us to apply PACA in regimes where p << n, i.e., in cases where the number of
#' samples is greater than the number of features.
#'
#' @param X \eqn{m} by \eqn{n_1} matrix; \cr
#'          Case (foreground) input data matrix. \cr
#'          It is recommended to normailize the feature scales as appropriate for the data modality.
#'          E.g. quantile normalization (or other comparable approaches) for RNAseq data.
#' @param Y \eqn{m} by \eqn{n_0} matrix; \cr
#'          Control (foreground) input data matrix. \cr
#'          It is recommended to normailize the feature scales as appropriate for the data modality.
#'          E.g. quantile normalization (or other comparable approaches) for RNAseq data.
#' @param k Positive integer, \eqn{k > 1}; \cr
#'          Number of, \eqn{k}, dimensions of shared variation to be removed from case data \code{X}.
#' @param niter Positive integer; \cr
#'              Number random sampling iterations.
#' @param batch Positive integer; \cr
#'              Number of cases/controls to be sampled for the training set at each iteration. \cr
#'              Note that \eqn{k < batch-1} and \eqn{batch < min{m, n1, n0}}.
#' @param rank Positive integer, optional (default \eqn{2}); \cr
#'             Number of dominant principle components to be computed for the corrected case data.
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
#'
#' @rdname rPACA

### TODO: SEEMS TO BE BROKEN
#       look for the code used in CONVERGE and sim analysis

paca_r <- function(X, Y, k,
                   niter = 20,
                   batch = 600,
                   rank = 5){
  X_pacacomb <- matrix(0, ncol=rank, nrow=dim(X)[2])
  X_ids <- colnames(X)
  if (is.null(X_ids)){
    X_ids <- seq(dim(X)[2])
    colnames(X) <- X_ids
  }
  Y_ids <- colnames(Y)
  if (is.null(Y_ids)){
    Y_ids <- seq(dim(Y)[2])
    colnames(Y) <- Y_ids
  }

  P_comb <- matrix(0, ncol=1, nrow=dim(X)[2])
  for (r in 1:niter){
    inX <- sample(X_ids, size=batch, replace=F)
    inY <- sample(Y_ids, size=batch, replace=F)
    xBatch <- X[,inX]
    yBatch <- Y[,inY]
    xRemain <- X[,setdiff(X_ids, inX)]

    #### NEW
    pacaBatch <- cpp_PACA(xBatch, yBatch, k, TRUE, FALSE, 0)

    pacapcBatch <- prcomp(t(pacaBatch$Xtil), rank. = rank, center = TRUE, scale. = TRUE)

    means_matrix <- t(kronecker(matrix(1,1,dim(xRemain)[1]),colMeans(xRemain)))
    X_centered <- xRemain - means_matrix
    X_tilde_remain <- X_centered - (pacaBatch$U0 %*% (t(pacaBatch$U0) %*% X_centered)) + means_matrix
    rm(means_matrix, X_centered)
    projMat_remain <- t(X_tilde_remain) %*% pacapcBatch$rotation

    pacapcBatch <- scale(pacapcBatch$x , center = TRUE, scale = T)
    projMat_remain <- scale(projMat_remain , center = TRUE, scale = T)
    comb_proj <- rbind(pacapcBatch, projMat_remain)
    P_comb <- cbind(P_comb, comb_proj[X_ids,])
    X_pacacomb <- comb_proj[X_ids,] + X_pacacomb
    rm(comb_proj, X_tilde_remain)


    ####

    # pacaBatch <- cpp_CCA(xBatch, yBatch, TRUE, 0)
    #
    # U_1 <- pacaBatch$U[,1:k] / t(kronecker(matrix(1,1,dim(pacaBatch$U)[1]),sqrt(colSums(pacaBatch$U[,1:k]^2))))
    # means_matrix <- t(kronecker(matrix(1,1,dim(xBatch)[1]),colMeans(xBatch)))
    # X_centered <- xBatch - means_matrix
    # X_tilde <- X_centered - (U_1 %*% (t(U_1) %*% X_centered)) + means_matrix
    # rm(means_matrix, X_centered)
    # pacapcBatch <- prcomp(t(X_tilde), rank. = rank, center = TRUE, scale. = T)
    #
    # means_matrix <- t(kronecker(matrix(1,1,dim(xRemain)[1]),colMeans(xRemain)))
    # X_centered <- xRemain - means_matrix
    # X_tilde_remain <- X_centered - (U_1 %*% (t(U_1) %*% X_centered)) + means_matrix
    # rm(means_matrix, X_centered, U_1)
    # projMat_remain <- t(X_tilde_remain) %*% pacapcBatch$rotation
    #
    # pacapcBatch <- scale(pacapcBatch$x , center = TRUE, scale = T)
    # projMat_remain <- scale(projMat_remain , center = TRUE, scale = T)
    # comb_proj <- rbind(pacapcBatch, projMat_remain)
    # P_comb <- cbind(P_comb, comb_proj[X_ids,])
    # X_pacacomb <- comb_proj[X_ids,] + X_pacacomb
    # rm(comb_proj, X_tilde_remain, X_tilde)

    cat("\n\tCompleted Random Sampling Iter:", r,"/",niter,"\n")
  }
  P_comb <- P_comb[,-1]
  Cmat <- P_comb %*% t(P_comb)
  eigPC <- eigen(Cmat)

  return(list(x = eigPC$vectors[,1:rank],
              eigVal = eigPC$values[1:rank]))

}
