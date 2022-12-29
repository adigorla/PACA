#' @title PACA Null Case Evaluation (nullEvalPACA)
#'
#' @name nullEvalPACA
#'
#' @description
#' This method applies a simple case/control label permutation approach to qunatify the
#' statistical significance of the presence of subphenotypic variation the cases, a givne fixed \eqn{k}.
#' The procedure should be able to reject the null (no subphenotypic variation) when there is sufficently
#' strong variation unique to the cases.
#'
#' @usage nullEvalPACA(X, Y, k, nPerm = 100)
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
#' @param nPerm Positive integer, optional (default \eqn{100}); \cr
#'          Number of random permutations to build the emperical null distribution.
#'
#' @return \code{nullEvalPACA} returns a list containing the following components:
#' \describe{
#'    \item{pval}{     non-negative real value; \cr
#'                     the significance of rejecting the null hypothesis that there is no subphenotypic structure in the case data \code{X}.
#'    }
#'    \item{empVar}{   non-negative real value; \cr
#'                     the variance of the top \emph{PACA} PC of the case data (\code{xtil}).
#'    }
#'
#'    \item{nullVars}{ list of size \eqn{nPerm}; \cr
#'                     the variances of the top \emph{PACA} PC of each of the permuted null data.
#'    }
#'}
#'
#' @export
nullEvalPACA <- function(X, Y, k, nPerm = 100){

  # perp data for perm
  xy <- cbind(t(X), t(Y))
  cat("\nStarting permutations...") # \nComb matrix dim : ", dim(xy))
  ids <- c(1:dim(xy)[2])
  colnames(xy) <- ids
  sz <- ceiling(length(ids)/2)

  # Scale input for CCA
  stdDat <- transformCCAinput(X, Y, .center = TRUE, .scale = TRUE)
  rm(X, Y)

  # get point stat for selected k
  empVar <- PACAvarPC1(stdDat$x, stdDat$y, k)
  rm(stdDat)

  # get dist of permuted null
  nullVars <- c()
  for (i in 1:nPerm){
    inCase <- sample(ids, size=sz, replace=F)
    Xs <- xy[,inCase]
    Ys <- xy[,setdiff(ids, inCase)]
    Xs <- scale(Xs, center = TRUE, scale = T)
    Ys <- scale(Ys, center = TRUE, scale = T)

    nv <- PACAvarPC1(Xs, Ys, k)

    nullVars <- c(nullVars, nv)
    if(i%%25 == 0){
      cat("\n\tDone with perm ", i)
    }
  }

  # get emperical p-value
  empPval <- calcPval(empVar, nullVars)

  return(list(pval = empPval,
              empVar = empVar,
              nullVars = nullVars))
}

#' @title Calculate PACA PC1 Variance
#'
#' @name PACAvarPC1
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
#' @param pcRank Positive integer, optional (default \eqn{2}); \cr
#'               Number of dominant principle components to be computed for the corrected case data.
#'
#'
#' @return non-negative real value; the variance of the top \emph{PACA} PC of the input data (\code{X}).
#'
#' @importFrom stats var
PACAvarPC1 <- function(X, Y, k, pcRank = 2){

  res <- PACA(X, Y, k, pcRank = pcRank, residOnly = FALSE)

  return(var(res$x[,1]))
}

#' @title Calculate Emperical P-values
#'
#' @name calcPval
#'
#' @noRd
#'
#' @param accs non-negative real value; \cr
#'             the variance of the top \emph{PACA} PC of the case data (\code{xtil}).
#' @param nulls list of size \eqn{nPerm}; \cr
#'              the variances of the top \emph{PACA} PC of each of the permuted null data.
#'
#' @return non-negative real value; the probability of getting the test statistic at least as extreme as one that was actually observed given the emperical null distrubution
calcPval <- function(accs, nulls){
  if(accs < min(nulls)){
    pval <- 1
  } else if (accs > max(nulls)){
    pval <- 0
  } else {
    n <- length(nulls)
    r <- sum(nulls > accs)
    pval <- (r+1)/(n+1)
  }

  return(pval)
}


