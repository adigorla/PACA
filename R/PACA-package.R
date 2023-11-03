################################################################################

#' @import RcppEigen
#' @useDynLib PACA, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp evalCpp
#'
#'
#' @param X \eqn{m} by \eqn{n_1} matrix of cases or a foregroung dataset\cr
#' **No missing values accepted**
#' @param Y \eqn{m} by \eqn{n_2} matrix of controls or a background dataset\cr
#' **No missing values accepted**
#'
#' @keywords internal
"_PACKAGE"

################################################################################
