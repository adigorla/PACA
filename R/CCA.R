#' @title General function for quickly performing canonical correlation analysis (CCA)
#'
#' @description
#' This function performs CCA on two matrices. As input, it takes
#' two matrices, \eqn{X} and \eqn{Y},  of size \eqn{m} by \eqn{n_1} and
#' \eqn{m} by \eqn{n_2} respectively, where \eqn{m} > \eqn{n_1},\eqn{n_2} (i.e., both
#' matrices have the same number of rows, but not necessarily the same
#' number of columns and the number of rows is greater than the number of columns).
#' The canonical variables are returned in decreasing order of correlation.
#' This code is  based on the '\code{CONFINED}' package by Mike Thompson.
#'  See \code{https://github.com/cozygene/CONFINED}.
#'
#' @param X \eqn{m} by \eqn{n_1} matrix
#' @param Y \eqn{m} by \eqn{n_2} matrix
#'
#' @return A  -  the loadings for \eqn{X}
#' @return B  -  the loadings for \eqn{Y}
#' @return U  -  canonical variables of \eqn{X}, calculated by column centering \eqn{X} and projecting it on \eqn{A}
#' @return V  -  canonical variables of \eqn{Y}, calculated by column centering \eqn{Y} and projecting it on \eqn{B}
#' @export
CCA = function (X, Y){
  tmp<-doRCCA(X, Y)

  A<-tmp[[1]]
  B<-tmp[[2]]

  d = min(c(dim(A)[2], dim(B)[2]))
  print(d)

  A = A[,1:d]
  B = B[,1:d]

  U = getUV(X, A)
  V = getUV(Y, B)

  return(list(A= A, B= B, U = U, V = V))
}

#' @useDynLib PACA, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp evalCpp
NULL
