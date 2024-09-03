###########################################################################################################

#' @title Randomized Phenotype Aware Components Analysis (rPACA)
#'
#' @description
#' This is the randomized version of the basic \code{\link{PACA}} algorithm.
#' rPACA allows us to apply PACA in regimes where p << n, i.e., in cases where the number of
#' samples is greater than the number of features.
#'
#' @param X \eqn{m} by \eqn{n_1} matrix; \cr
#'          Case (foreground) input data matrix. \cr
#'          It is recommended to normalize the feature scales as appropriate for the data modality.
#'          E.g. quantile normalization (or other comparable approaches) for RNAseq data.
#' @param Y \eqn{m} by \eqn{n_0} matrix; \cr
#'          Control (foreground) input data matrix. \cr
#'          It is recommended to normalize the feature scales as appropriate for the data modality.
#'          E.g. quantile normalization (or other comparable approaches) for RNAseq data.
#' @param k positive integer, optional (default: \eqn{NULL}); \cr
#'          Number of, \eqn{k}, dimensions of shared variation to be removed from case data \code{X}. \cr
#'          When \eqn{k = NULL} (default), K is automatically infered, i.e, we run auto rPACA by default.
#' @param niter Positive integer; \cr
#'              Number random sampling iterations.
#' @param batch Positive integer; \cr
#'              Number of cases/controls to be sampled for the training set at each iteration. \cr
#'              Note that \eqn{k < batch-1} and \eqn{batch < min{m, n1, n0}}.
#' @param rank Positive integer, optional (default \eqn{5}); \cr
#'             Number of dominant principle components to be computed for the corrected case data.
#' @param thrsh Positive real value, optional (default \eqn{10.0}); \cr
#'              Threshold value for the maximum ratio of variance in \emph{PACA} corrected \code{X} PCs and the variance it explain in Y
#'              which indicates the presence of residual shared variation in X.\cr
#'              Only used when \code{k = NULL}.
#' @param info Integer, optional (default: 1); \cr
#'          Verbosity level for the log generated. \cr
#'          0: Errors and warnings only \cr
#'          1: Basic informational messages \cr
#'          2: More detailed informational messages \cr
#'          3: Debug mode, all informational log is dumped
#'
#' @return \code{rpaca} returns a list containing the following components:
#' \describe{
#'    \item{x}{        \eqn{n_1} by \eqn{rank} matrix; \cr
#'                     the projections / scores of the randomized \emph{PACA} corrected case data.
#'    }
#'    \item{eigs}{   vector of size \eqn{rank}; \cr
#'                     the eigen values of the returned PACA PC components.
#'    }
#'   \item{k.iter}{  vector of size \eqn{niter}; \cr
#'                      the number of shared dimensions removed at each iteration.}
#'    }
#' @export
#'
#' @rdname rPACA
###########################################################################################################

rpaca <- function(X, Y, k=NULL, niter = 10, batch = 100, rank = 5, thrsh = 10.0, scale = FALSE, info = 1) {
  # Input shape check
  if (dim(X)[1] != dim(Y)[1]) {
    stop(sprintf(
      "RowSize X (%d) is NOT equal to RowSize Y (%d).\nInput matrices should have shape: features-by-samples (MxN) where M size should match.",
      dim(X)[1], dim(Y)[1]
    ))
  }

  if (batch >= min(c(dim(X)[1], dim(Y)[1], dim(X)[2]))) {
    stop(sprintf(
      "Batch size needs to be less than min(m, n1, n0) = %d.",
      min(c(dim(X)[1], dim(Y)[1], dim(X)[2]))
    ))
  }

  if (!is.null(k) && k >= batch - 1) {
    stop(sprintf(
      "When k is specified, k needs to be less than batch - 1 = %d.",
      batch - 1
    ))
  }

  names_list <- getNames(X, Y)

  # Call C++ function
  if (is.null(k)){
    tmp <- cpp_autorPACA(X, Y, niter, batch, rank, thrsh, scale, info)
  } else{
    tmp <- cpp_rPACA(X, Y, k, niter, batch, rank, scale, info)
    tmp$k.iter <- rep(k, niter)
  }


  # Add row and column names
  rownames(tmp$x) <- names_list$X$col
  colnames(tmp$x) <- paste0("PAC", 1:rank)
  names(tmp$eigs) <- paste0("PAC", 1:rank)
  names(tmp$k.iter) <- paste0("iter", 1:niter)


  return(tmp)
}

###########################################################################################################
