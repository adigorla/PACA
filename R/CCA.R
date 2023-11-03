###########################################################################################################

#' @title Highly efficient implementation of canonical correlation analysis (CCA)
#'
#' @description
#' This function performs CCA on two matrices. As input, it takes
#' two matrices, \eqn{X} and \eqn{Y},  of size \eqn{m} by \eqn{n_1} and
#' \eqn{m} by \eqn{n_2} respectively, where \eqn{m} > \eqn{n_1},\eqn{n_2} (i.e., both
#' matrices have the same number of rows, but not necessarily the same
#' number of columns and the number of rows is greater than the number of columns).
#' The canonical variables are returned in decreasing order of correlation.
#' We solve the CCA without using an explicit inverse which should improve numerical
#' stability and speed relative to than other implementations that invert.
#'
#' @param X matrix; \eqn{m} by \eqn{n_1} cases or foreground
#' @param Y matrix; \eqn{m} by \eqn{n_2} controls or backgorund
#' @param scale bool, optional (default: \eqn{TRUE}); normalize (center+scale) each matrix columnwise
#' @param info int, optional (default: 1); verbosity for the log generated
#' \itemize{
#' \item 0 : Errors and warnings only
#' \item 1 : Basic Informational messages
#' \item 2 : More detailed Informational messages
#' \item 3 : Debug mode, all informational log is dumped
#' }
#'@return \code{cca} returns a list containing the following components:
#' \describe{
#'    \item{A}{ the loadings for \eqn{X}
#'    }
#'    \item{B}{ the loadings for \eqn{Y}
#'    }
#'    \item{U}{ canonical variables of \eqn{X}, calculated by column centering \eqn{X} and projecting it on \eqn{A}
#'    }
#'    \item{V}{ canonical variables of \eqn{Y}, calculated by column centering \eqn{Y} and projecting it on \eqn{B}
#'    }
#'}
#' @export
cca = function (X, Y, scale = TRUE, info = 1){

  names_list <- getNames(X, Y)

  tmp<-doCCA(X, Y, scale, info)

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

  return(tmp)
}

###########################################################################################################
