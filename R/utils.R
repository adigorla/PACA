###########################################################################################################

#' @title Transform Input for CCA
#'
#' @description
#' This method transposes the input matrices from the standard \eqn{dim(n,m)} to
#' the form \eqn{dim(m,n)} which PACA uses. Next, it also centers/scales the
#' input data along the sample axis.
#'
#' @param X \eqn{n_1} by \eqn{m} matrix; \cr
#'          Case (foreground) input data matrix. \cr
#'          It is recommended to normailize the feature scales as appropriate for the data modality.
#'          E.g. quantile normalization (or other comparable approaches) for RNAseq data.
#' @param Y \eqn{n_0} by \eqn{m} matrix; \cr
#'          Control (foreground) input data matrix. \cr
#'          It is recommended to normailize the feature scales as appropriate for the data modality.
#'          E.g. quantile normalization (or other comparable approaches) for RNAseq data.
#' @param .center bool, optional (default \eqn{TRUE}); \cr
#'                  If \eqn{TRUE}, data is mean centered.
#' @param .scale bool, optional (default \eqn{FALSE}); \cr
#'                  If \eqn{TRUE}, data is unit scaled.
#'
#' @return \code{transformCCAinput} returns a list containing the following components:
#' \describe{
#'    \item{x}{   \eqn{m} by \eqn{n_1} matrix; \cr
#'                the transposed and sample centered/scaled case data \code{X}.
#'    }
#'    \item{Y}{   \eqn{m} by \eqn{n_0} matrix; \cr
#'                the transposed and sample centered/scaled case data \code{Y}.
#'    }
#'}
#'
#' @keywords internal
#' @noRd
transformCCAinput <- function(X, Y, center = TRUE, scale = FALSE){

  # Transpose and scale (by samples) the input matrix
  x <- scale(t(X), center = .center, scale = .scale)
  y <- scale(t(Y), center = .center, scale = .scale)

  if ( ((sum(sum(is.na(x))) ) > 0) || ((sum(sum(is.na(y))) ) > 0)){
    stop("Division by zero due to constant features in either X or Y.")
  }
  return(list(x = x, y = y))
}

###########################################################################################################

###########################################################################################################

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
#'
#' @keywords internal
#' @noRd
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

###########################################################################################################

###########################################################################################################

#' @title Get names from data object
#'
#' @description
#' Internal function to get names from inputted matrix so we can re-add them to the return output from functions.
#' Need this since Eigen/our cpp code can destroy this infor during processing
#'
#' @param X \eqn{n_1} by \eqn{m} matrix;
#' @param Y \eqn{n_0} by \eqn{m} matrix;
#'
#' @return Returns a list with row and colnames
#' @keywords internal
#' @noRd
getNames <- function(X, Y){

  return(list(
    X = list(
    col = colnames(X),
    row = row.names(X)
  ),
  Y = list(
    col = colnames(Y),
    row = row.names(Y)
  )
                ))
}

###########################################################################################################
