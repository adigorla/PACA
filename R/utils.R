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
#' @export
transformCCAinput <- function(X, Y, .center = TRUE, .scale = FALSE){

  # Transpose and scale (by samples) the input matrix
  x <- scale(t(X), center = .center, scale = .scale)
  y <- scale(t(Y), center = .center, scale = .scale)

  if ( ((sum(sum(is.na(x))) ) > 0) || ((sum(sum(is.na(y))) ) > 0)){
    stop("Division by zero due to constant features in either X or Y.")
  }
  return(list(x = x, y = y))
}
