###########################################################################################################

#' @title  Phenotype Aware Components Analysis (PACA)
#'
#' @description
#' Phenotype Aware Components Analysis (PACA) is a
#' contrastive learning approach leveraging canonical correlation analysis to robustly capture weak sources of
#' subphenotypic variation. Given case-control data of any modality, PACA highlights the dominant variation in a
#' subspace that is not affected by background variation as a putative representation of phenotypic heterogeneity. We do so by
#' removing the top \code{k} components of shared variation from the cases (or foreground) \code{X}.
#' In the context of complex disease, PACA learns a gradient of variation unique to cases \code{X} in
#' a given dataset, while leveraging control samples \code{Y} for accounting for variation and imbalances of biological
#' and technical confounders between cases and controls.
#'
#' @param X \eqn{m} by \eqn{n_1} matrix, where \eqn{m > n_1}; \cr
#'          Case (foreground) input data matrix. \cr
#'
#'          Note: this input data needs to be scaled along the samples axis before being provided as input.
#'          This preprocessing can be done using the \code{\link{transformCCAinput}} function.
#' @param Y \eqn{m} by \eqn{n_0} matrix, where \eqn{m > n_0}; \cr
#'          Control (foreground) input data matrix. \cr
#'          Note: this input data needs to be scaled along the samples axis before being provided as input.
#'          This preprocessing can be done using the \code{\link{transformCCAinput}} function.
#' @param k positive integer, optional (default: \eqn{NULL}); \cr
#'          Number of, \eqn{k}, dimensions of shared variation to be removed from case data \code{X}. \cr
#'          When \eqn{k = NULL} (default), K is automatically infered, i.e, we run autoPACA by default.
#' @param scale bool, optional (default: \eqn{FALSE}); normalize (center+scale) each matrix column-wise
#' @param rank Positive integer, optional (default \eqn{2}); \cr
#'               Number of dominant principle components to be computed for the corrected case data.
#' @param thrsh Positive real value, optional (default \eqn{10}); \cr
#'              Threshold value for the maximum ratio of variance in \emph{PACA} corrected \code{X} PCs and the variance it explain in Y
#'              which indicates the presence of residual shared variation in X.
#'
#' @param ccweights bool, optional (default \eqn{FALSE}); \cr
#'                  If \eqn{TRUE}, return the \emph{PACA} corrected case data (\code{xtil}) ONLY.
#' @param info Integer, optional (default: 0); \cr
#'          Verbosity level for the log generated. \cr
#'          0: Errors and warnings only \cr
#'          1: Basic informational messages \cr
#'          2: More detailed informational messages \cr
#'          3: Debug mode, all informational log is dumped
#'
#'@return By default, \code{paca} returns a list containing the following components:
#' \describe{
#'    \item{Xtil}{     \eqn{m} by \eqn{n_1} matrix; \cr
#'                     the \emph{PACA} corrected case data, i.e., the data with the case-specific variation only.
#'    }
#'    \item{U0}{       \eqn{m} by \eqn{k} matrix; \cr
#'                     the \emph{PACA} shared components that are removed from \eqn{X}.
#'    }
#'    \item{x}{        \eqn{n_1} by \eqn{rank} matrix; \cr
#'                     the projections / scores of the \emph{PACA} corrected case data (\code{Xtil}).
#'    }
#'    \item{rotation}{  \eqn{m} by \eqn{rank} matrix; \cr
#'                      the rotation (eigenvectors)  of the \emph{PACA} corrected case data (\code{Xtil}).
#'    }
#'    \item{k}{         the number of shared components removed, int
#'    }
#'}
#'@return When \eqn{ccweights = TRUE}, \code{paca} returns a list containing the CCA direction and variates along withe the \emph{PACA} principle components:
#' \describe{
#'    \item{Xtil}{     \eqn{m} by \eqn{n_1} matrix; \cr
#'                     the \emph{PACA} corrected case data, i.e., the data with the case-specific variation only.
#'    }
#'    \item{U0}{       \eqn{m} by \eqn{k} matrix; \cr
#'                     the \emph{PACA} shared components that are removed from \eqn{X}.
#'    }
#'    \item{x}{        \eqn{n_1} by \eqn{rank} matrix; \cr
#'                     the projections / scores of the \emph{PACA} corrected case data (\code{Xtil}).
#'    }
#'    \item{rotation}{\eqn{m} by \eqn{rank} matrix; \cr
#'                      the rotation (eigenvectors)  of the \emph{PACA} corrected case data (\code{Xtil}).
#'    }
#'    \item{k}{        the number of shared components removed, int
#'    }
#'    \item{A}{       the loadings for \eqn{X}
#'    }
#'    \item{B}{       the loadings for \eqn{Y}
#'    }
#'    \item{U}{       canonical variables of \eqn{X}, calculated by column centering \eqn{X} and projecting it on \eqn{A}
#'    }
#'    \item{V}{       canonical variables of \eqn{Y}, calculated by column centering \eqn{Y} and projecting it on \eqn{B}
#'    }
#'}
#'
#' @export
#'
#' @rdname PACA
paca <- function(X, Y,
                 k = NULL,
                 scale = FALSE,
                 rank = 5,
                 thrsh = 10.0,
                 ccweights = FALSE,
                 info = 1){
  # Input shape check
  if (dim(X)[1] != dim(Y)[1]) {
    stop(sprintf(
      "RowSize X (%d) is NOT equal to RowSize Y (%d).\nInput matrices should have shape: features-by-samples (MxN) where M size should match.",
      dim(X)[1], dim(Y)[1]
    ))
  }
  if(dim(X)[1] < min(dim(X)[2], dim(Y)[2])){
    stop(sprintf(
      "Feature size is smaller than min sample size %d.\nPACA expects more samples than features, try running rPACA instead."
      , min(dim(X)[2], dim(Y)[2])
    ))
  }

  names_list <- getNames(X, Y)

  if (is.null(k)){ # run autoPACA
    tmp <- cpp_autoPACA(X, Y, scale, ccweights, thrsh, info)
    if(ccweights){
      cc_names <- paste0('CC', seq(length(tmp$corr)))
      names(tmp$corr) <- cc_names
      colnames(tmp[['A']]) <- cc_names
      colnames(tmp[['B']]) <- cc_names
      colnames(tmp[['U']]) <- cc_names
      colnames(tmp[['V']]) <- cc_names

      row.names(tmp[['A']]) <- names_list$X[['col']]
      row.names(tmp[['B']]) <- names_list$Y[['col']]
      row.names(tmp[['U']]) <- names_list$X[['row']]
      row.names(tmp[['V']]) <- names_list$Y[['row']]
    }
  } else{ # run PACA with fixed K
    print("running w/ fixed k")
    tmp <- cpp_PACA(X, Y, k, scale, ccweights, info)
    if(ccweights){
      cc_names <- paste0('CC', seq(length(tmp$corr)))
      names(tmp$corr) <- cc_names
      colnames(tmp[['A']]) <- cc_names
      colnames(tmp[['B']]) <- cc_names
      colnames(tmp[['U']]) <- cc_names
      colnames(tmp[['V']]) <- cc_names

      row.names(tmp[['A']]) <- names_list$X[['col']]
      row.names(tmp[['B']]) <- names_list$Y[['col']]
      row.names(tmp[['U']]) <- names_list$X[['row']]
      row.names(tmp[['V']]) <- names_list$Y[['row']]
      row.names(tmp[['Cxx']]) <- names_list$X[['col']]
      row.names(tmp[['Cyy']]) <- names_list$Y[['col']]
      row.names(tmp[['Cxy']]) <- names_list$X[['col']]
    }

  }

  colnames(tmp[['Xtil']]) <- names_list$X[['col']]
  row.names(tmp[['Xtil']]) <- names_list$X[['row']]
  colnames(tmp[['U0']]) <- paste0('CC', seq(dim(tmp[['U0']])[2]))
  row.names(tmp[['U0']]) <- names_list$X[['row']]

  # do PCA decomp of case specific signal
  # pca.res <- prcomp(t(tmp[['Xtil']]), rank.= rank)
  pca.res <- eigenprcomp(t(tmp[['Xtil']]), rank = rank, info = info)
  tmp$x <- pca.res$x
  tmp$rotation <- pca.res$rotation

  return(tmp)
}

###########################################################################################################

###########################################################################################################

#' @title Calculate PACA PC1 Variance
#'
#' @name paca_varPC1
#'
#' @noRd
#'
#' @usage NULL
#'
#' @param X \eqn{m} by \eqn{n_1} matrix; \cr
#'          Case (foreground) input data matrix. \cr
#'          Note: this input data needs to be scaled along the samples axis before being provided as input.
#'          This preprocessing can be done using the \code{\link{transformCCAinput}} function.
#' @param Y \eqn{m} by \eqn{n_0} matrix; \cr
#'          Control (foreground) input data matrix. \cr
#'          Note: this input data needs to be scaled along the samples axis before being provided as input.
#'          This preprocessing can be done using the \code{\link{transformCCAinput}} function.
#' @param k Positive integer, \eqn{k > 1}; \cr
#'          Number of, \eqn{k}, dimensions of shared variation to be removed from case data \code{X}.
#' @param info Integer, optional (default: 0); \cr
#'          Verbosity level for the log generated. \cr
#'          0: Errors and warnings only \cr
#'          1: Basic informational messages \cr
#'          2: More detailed informational messages \cr
#'          3: Debug mode, all informational log is dumped
#'
#' @return non-negative real value; the variance of the top \emph{PACA} PC of the input data (\code{X}).
#'
#' @importFrom stats var
#'
#' @keywords internal
#' @noRd
paca_varPC1 <- function(X, Y, k, info = 0){

  res <- cpp_PACA(X, Y, k, TRUE, FALSE, info)
  pca.res <- eigenprcomp(t(res[['Xtil']]), rank = rank, info = info)$x
  # pca.res <- rpca(t(res[['Xtil']]), k = 1, center = TRUE, scale = FALSE, q = 2)$x

  return(var(pca.res[,1]))
}

###########################################################################################################

###########################################################################################################

#' @title PACA Null Case Evaluation (paca_null)
#'
#' @name paca_null
#'
#' @description
#' This method applies a simple case/control label permutation approach to quantify the
#' statistical significance of the presence of subphenotypic variation in the cases, given a fixed \eqn{k}.
#' The procedure should be able to reject the null (no subphenotypic variation) when there is sufficiently
#' strong variation unique to the cases.
#'
#' @usage paca_null(X, Y, k, nperm = 100, info = 0)
#'
#' @param X \eqn{n_1} by \eqn{m} matrix; \cr
#'          Case (foreground) input data matrix. \cr
#'          It is recommended to normalize the feature scales as appropriate for the data modality.
#'          E.g. quantile normalization (or other comparable approaches) for RNAseq data.
#' @param Y \eqn{n_0} by \eqn{m} matrix; \cr
#'          Control (background) input data matrix. \cr
#'          It is recommended to normalize the feature scales as appropriate for the data modality.
#'          E.g. quantile normalization (or other comparable approaches) for RNAseq data.
#' @param k Positive integer, \eqn{k > 1}; \cr
#'          Number of, \eqn{k}, dimensions of shared variation to be removed from case data \code{X}.
#' @param nperm Positive integer, optional (default \eqn{100}); \cr
#'          Number of random permutations to build the empirical null distribution.
#' @param info Integer, optional (default: 0); \cr
#'          Verbosity level for the log generated. \cr
#'          0: Errors and warnings only \cr
#'          1: Basic informational messages \cr
#'          2: More detailed informational messages \cr
#'          3: Debug mode, all informational log is dumped
#'
#' @return \code{paca_null} returns a list containing the following components:
#' \describe{
#'    \item{pval}{     non-negative real value; \cr
#'                     the significance of rejecting the null hypothesis that there is no subphenotypic structure in the case data \code{X}.
#'    }
#'    \item{empVar}{   non-negative real value; \cr
#'                     the variance of the top \emph{PACA} PC of the case data (\code{xtil}).
#'    }
#'    \item{nullVars}{ list of size \eqn{nperm}; \cr
#'                     the variances of the top \emph{PACA} PC of each of the permuted null data.
#'    }
#' }
#'
#' @export
#'
#' @rdname PACAnull
paca_null <- function(X, Y, k, nperm = 100, info = 0){

  # perp data for perm
  xy <- cbind(t(X), t(Y))
  cat("\nStarting permutations...\n") # \nComb matrix dim : ", dim(xy))
  ids <- c(1:dim(xy)[2])
  colnames(xy) <- ids
  sz <- ceiling(length(ids)/2)

  # Scale input for CCA
  stdDat <- transformCCAinput(X, Y, center = TRUE, scale = TRUE)
  rm(X, Y)

  # get point stat for selected k
  empVar <- paca_varPC1(stdDat$x, stdDat$y, k, info = info)
  rm(stdDat)

  # get dist of permuted null
  nullVars <- c()
  for (i in 1:nperm){
    inCase <- sample(ids, size=sz, replace=F)
    Xs <- xy[,inCase]
    Ys <- xy[,setdiff(ids, inCase)]
    Xs <- scale(Xs, center = TRUE, scale = T)
    Ys <- scale(Ys, center = TRUE, scale = T)

    nv <- paca_varPC1(Xs, Ys, k)

    nullVars <- c(nullVars, nv)
    if (info > 1){
      if(i%%25 == 0){
        cat("Done with permutation: ", i, "\n")
      }
    }
  }

  # get emperical p-value
  empPval <- calcPval(empVar, nullVars)

  return(list(pval = empPval,
              empVar = empVar,
              nullVars = nullVars))
}

###########################################################################################################
