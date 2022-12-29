#' @title Phenotype Aware Components Analysis with Automatic \eqn{k} Selection (autoPACA)
#'
#' @description
#' This is an extension to the basic \code{\link{PACA}} algorithm in which the
#' number of shared dimensions \code{k} to be removed is automatically selected.
#'
#' @param X \eqn{n_1} by \eqn{m} matrix, where \eqn{m > n_1}; \cr
#'          Case (foreground) input data matrix. \cr
#'          It is recommended to normailize the feature scales as appropriate for the data modality.
#'          E.g. quantile normalization (or other comparable approaches) for RNAseq data.
#' @param Y \eqn{n_0} by \eqn{m} matrix, where \eqn{m > n_0}; \cr
#'          Control (foreground) input data matrix. \cr
#'          It is recommended to normailize the feature scales as appropriate for the data modality.
#'          E.g. quantile normalization (or other comparable approaches) for RNAseq data.
#' @param pcRank Positive integer, optional (default \eqn{2}); \cr
#'               Number of dominant principle components to be computed for the corrected case data.
#' @param thrsh Positive real value, optional (default \eqn{10}); \cr
#'              Threshold value for the maximum ratio of variance in \emph{PACA} corrected \code{X} PCs and the variance it explain in Y
#'              which indicates the presence of residual shared variation in X.
#'
#' @return \code{autoPACA} returns a list containing the following components:
#' \describe{
#'    \item{k}{        Integer; \cr
#'                     the estimated number of shared axis of variation to be removed from the case data \code{X}.
#'    }
#'    \item{x}{        \eqn{n_1} by \eqn{pcRank} matrix; \cr
#'                     the projections / scores of the \emph{PACA} corrected case data (\code{xtil}).
#'    }
#'    \item{rotation}{  \eqn{m} by \eqn{pcRank} matrix; \cr
#'                      the rotation (eigenvectors)  of the \emph{PACA} corrected case data (\code{xtil}).
#'    }
#'
#'    \item{xtil}{     \eqn{m} by \eqn{n_1} matrix; \cr
#'                     the \emph{PACA} corrected case data, i.e., the data with the case-specific variation only. \cr
#'                     NOTE: when \eqn{residOnly = TRUE}, only the \code{xtil} matrix is returned.
#'    }
#'}
#'
#' @export
#'
#' @importFrom rsvd rpca
#' @importFrom stats qchisq
#' @importFrom stats var
#' @importFrom gtools permute
autoPACA <- function(X, Y, pcRank = 2, thrsh = 10){

  # Scale input for CCA
  stdDat <- transformCCAinput(X, Y, .center = TRUE, .scale = TRUE)

  # init parms
  n <- min(dim(stdDat$x)[2], dim(stdDat$y)[2])
  min <- 1
  max <- dim(stdDat$x)[2]
  Ulim <- max
  n_perm <- floor(log2(n)/0.01)
  nullCase <- 0
  counter <- 0

  X_in <- X_out <- stdDat$x
  Y_in <- Y_out <- stdDat$y

  res <- CCA(X_in,Y_in)

  # K selection algo Main
  while (min < max){
    k <- floor((min+max)/2)
    counter <- counter + 1
    if (k < 2){
      nullCase <- 1
      break
    }

    cat("\n\tIter: ", counter, " ; K: ", k, " ; MIN: ", min," ; MAX: ", max) #, " .... ")
    U_1 <- res$U[,1:k] / t(kronecker(matrix(1,1,dim(res$U)[1]),sqrt(colSums(res$U[,1:k]^2))))
    V_1 <- res$V[,1:k] / t(kronecker(matrix(1,1,dim(res$V)[1]),sqrt(colSums(res$V[,1:k]^2))))

    # X'_1
    means_matrix <- t(kronecker(matrix(1,1,dim(X_in)[1]),colMeans(X_in)))
    X_centered <- X_in - means_matrix
    Xtil_1 <- X_centered - (U_1 %*% (t(U_1) %*% X_centered)) + means_matrix
    rm(means_matrix, X_centered)

    # X'_2
    means_matrix <- t(kronecker(matrix(1,1,dim(X_out)[1]),colMeans(X_out)))
    X_centered <- X_out - means_matrix
    Xtil_2 <- X_centered - (U_1 %*% (t(U_1) %*% X_centered)) + means_matrix
    rm(means_matrix, X_centered)

    # Y'_2
    means_matrix <- t(kronecker(matrix(1,1,dim(Y_out)[1]),colMeans(Y_out)))
    Y_centered <- Y_out - means_matrix
    Ytil_2 <- Y_centered - (V_1 %*% (t(V_1) %*% Y_centered)) + means_matrix
    rm(means_matrix, Y_centered)

    pacapc_in <- rpca(t(Xtil_1), k = 2, center = TRUE, scale = FALSE, q = 2)
    rm(Xtil_1)
    Z1 <- pacapc_in$rotation[,1]
    rm(pacapc_in)

    var_x2 <- var(t(Xtil_2)%*%Z1)
    #df <- dim(Xtil_2)[2]
    #vx2_ci <- (df-1)*c(var_x2)/qchisq(c(.975,.025), df-1)
    var_y2 <- var(t(Ytil_2)%*%Z1)
    #df <- dim(Ytil_2)[2]
    #vy2_ci <- (df-1)*c(var_y2)/qchisq(c(.975,.025), df-1)

    #eX <- sqrt(max(eigen(t(Xtil_2)%*%Xtil_2)$values))
    #eY <- sqrt(max(eigen(t(Ytil_2)%*%Ytil_2)$values))

    ct_x <- ct_y <- 0

    # Run permuations
    for (i in 1:n_perm){
      Z <- permute(Z1)
      var_x <- var(t(Xtil_2)%*%Z)
      var_y <- var(t(Ytil_2)%*%Z)
      if(var_x > var_x2){
        ct_x <- ct_x + 1
      }
      if(var_y > var_y2){
        ct_y <- ct_y + 1
      }
    }
    rm(Xtil_2, Ytil_2)
    # not printing debugging stats for prodVer
    #cat("Xct: ", ct_x, " ; Yct: ", ct_y, " ; var(X2): ", var_x2, "; CI99 = " , vx2_ci, " ; var(Y2): ", var_y2, "; CI99 = " , vy2_ci, "; eigX2 = " , eX, "; eigY2 = " , eY,"\n")

    # check conditonts
    if ((ct_x < 1) && (ct_y >= 1)){
      max <- k
    }else if((ct_x < 1) && (ct_y < 1)){
      if ( var_x2 > (thrsh*var_y2)){
        max <- k
      } else{
        min <- k+1
      }

    }else if((ct_x >= 1) && (ct_y >= 1)){
      max <- k-1
    }else {
      nullCase <- 1
      break
    }
  }
  rm(res)

  # prevents k0 from being unresonabelly high
  if (k >= Ulim-1){
    k <- Ulim - 3
  }
  k0 <- k+1

  res <- PACA(stdDat$x, stdDat$y, k0, pcRank = pcRank, residOnly = FALSE)

  return(list(
    k = k0,
    x = res$x,
    rotation = res$rotation,
    xtil = res$xtil))

}

