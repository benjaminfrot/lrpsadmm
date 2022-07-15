#' Estimate the correlation matrix of a semiparametric elliptical
#' copula using a modified Kendall rank correlation matrix
#'
#' @title Estimate the scale matrix of a semiparametric elliptical copula (a.k.a.
#' transelliptical distribution).
#' @param X An n x p data matrix.
#' @return A p x p correlation matrix.
#' @details Given X, an n x p data matrix, let K be its Kendall correlation matrix.
#' This function first computes \code{Chat <- sin(0.5*pi*K)} and projects it onto the
#' space of correlation matrices (in Frobenius norm) using the \code{nearPD} function
#' of the \code{Matrix} package
#' @examples
#' set.seed(0)
#' # Data from a mixture semiparametric elliptical copula (Cauchy) / multivariate normal
#' sim.data <- generate.latent.ggm.data(n=2000, p=100, h=5, outlier.fraction = 0.05,
#'                                     sparsity = 0.02, sparsity.latent = 0.7)
#' X <- sim.data$obs.data;
#' Sigma.Kendall <- Kendall.correlation.estimator(X)

#' @import Matrix
#' @export
Kendall.correlation.estimator <- function(X) {
  C <- cor(X, method='kendall')

  # Make sure it is a correlation matrix
  as.matrix(Matrix::nearPD(sin(0.5 * pi * C), corr = T)$mat)
}
