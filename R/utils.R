###########################################################################################################

#' @title Principal Components Analysis using RcppEigen
#'
#' @description
#' This function performs a principal components analysis on the given data matrix
#' using an efficient RcppEigen implementation. It returns the results as a list containing
#' various components of the PCA.
#'
#' @param x A numeric matrix or data frame which provides the data for the
#'          principal components analysis.
#' @param center A logical value indicating whether the variables should be
#'               shifted to be zero centered. Alternately, a vector of length
#'               equal the number of columns of x can be supplied. The value
#'               is passed to scale.
#' @param scale A logical value indicating whether the variables should be
#'              scaled to have unit variance before the analysis takes place.
#'              Alternatively, a vector of length equal the number of columns
#'              of x can be supplied. The value is passed to scale.
#' @param rank A positive integer specifying the desired rank (i.e., number
#'             of principal components to be computed).
#' @param tol A threshold for omitting principal components (see Details).
#' @param info Integer, optional (default: 0); \cr
#'          Verbosity level for the log generated. \cr
#'          0: Errors and warnings only \cr
#'          1: Basic informational messages \cr
#'          2: More detailed informational messages \cr
#'          3: Debug mode, all informational log is dumped
#'
#' @return A list with class "prcomp" containing the following components:
#' \describe{
#'   \item{sdev}{The standard deviations of the principal components.}
#'   \item{rotation}{The matrix of variable loadings (eigenvectors).}
#'   \item{x}{The rotated data (the centred, scaled data multiplied by the rotation matrix).}
#'   \item{center}{The centering used, or FALSE.}
#'   \item{scale}{The scaling used, or FALSE.}
#' }
#'
#' @details
#' The calculation is done using singular value decomposition of the (centered
#' and scaled) data matrix. This is generally the preferred method for numerical
#' accuracy. The `tol` parameter can be used to omit components whose standard
#' deviations are less than or equal to `tol` times the standard deviation of
#' the first component.
#'
#' @export
#'
#' @rdname Eigenprcomp
eigenprcomp <- function(x, center = TRUE, scale = FALSE, rank = NULL, tol = NULL, info = 1) {
  # Check if x is a data frame and convert to matrix if necessary
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }

  # Save original row and column names
  orig_row_names <- rownames(x)
  orig_col_names <- colnames(x)

  # Call the Rcpp function
  result <- cpp_prcomp(x, center, scale,
                       rank, tol, info)
  # Add back original row and column names
  if (!is.null(orig_row_names)) {
    row.names(result$x) <- orig_row_names
  }
  if (!is.null(orig_col_names)) {
    row.names(result$rotation) <- orig_col_names
    names(result$center) <- orig_col_names
    names(result$scale) <- orig_col_names
  }

#   # Name the PC components as PAC1, PAC2, ...
#   colnames(result$x) <- paste0("PAC", seq_len(ncol(result$x)))
#   colnames(result$rotation) <- paste0("PAC", seq_len(ncol(result$rotation)))
  # Name the PC components as PC1, PC2, ...
  colnames(result$x) <- paste0("PC", seq_len(ncol(result$x)))
  colnames(result$rotation) <- paste0("PC", seq_len(ncol(result$rotation)))

  return(result)
}

###########################################################################################################

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
transformCCAinput <- function(X, Y, center = TRUE, scale = FALSE){

  # Transpose and scale (by samples) the input matrix
  x <- scale(t(X), center = center, scale = scale)
  y <- scale(t(Y), center = center, scale = scale)

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
