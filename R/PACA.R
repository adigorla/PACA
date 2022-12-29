#' @title  Phenotype Aware Components Analysis (PACA)
#'
#' @description
#' Phenotype Aware Components Analysis (PACA) is a
#' contrastive learning approach leveraging canonical correlation analysis to robustly capture weak sources of
#' subphenotypic variation. Given case-control data of any modality, PACA highlights the dominant variation in a
#' subspace that is not affected by background variation as a putative representation of phenotypic heterogeneity. We do so by
#' removing the top \code{k} components of shared variation from the cases \code{X}.
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
#' @param k Positive integer, \eqn{k > 1}; \cr
#'          Number of, \eqn{k}, dimensions of shared variation to be removed from case data \code{X}.
#' @param pcRank Positive integer, optional (default \eqn{2}); \cr
#'               Number of dominant principle components to be computed for the corrected case data.
#' @param residOnly bool, optional (default \eqn{FALSE}); \cr
#'                  If \eqn{TRUE}, return the \emph{PACA} corrected case data (\code{xtil}) ONLY.
#'
#'@return By default, \code{PACA} returns a list containing the following components:
#' \describe{
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
PACA <- function(X, Y, k, pcRank = 2, residOnly = FALSE){

  # calculate the CCA components
  pacaFull <- CCA(X, Y)
  rm(Y)

  # k cannot be less than 2
  if (k > 1){
    k0_hat <- k
  } else{
    k0_hat <- 2
  }

  # get the corrected case data matrix
  U_1 <- pacaFull$U[,1:k0_hat] / t(kronecker(matrix(1,1,dim(pacaFull$U)[1]),sqrt(colSums(pacaFull$U[,1:k0_hat]^2))))
  means_matrix <- t(kronecker(matrix(1,1,dim(X)[1]),colMeans(X)))
  X_centered <- X - means_matrix
  X_tilde <- X_centered - (U_1 %*% (t(U_1) %*% X_centered)) + means_matrix
  rm(means_matrix, X_centered, U_1)

  # END algo here and return Xtilde only
  if (residOnly) {
    return(X_tilde)
  }

  res <- rpca(t(X_tilde), k = pcRank, center = TRUE, scale = FALSE, q = 2)

  return(list(x = res$x,
              rotation = res$rotation,
              xtil = X_tilde))
}
