#' Generate data from Gaussian graphical model with hidden variables.
#'
#' @title Generate synthetic data
#' @param n Number of samples.
#' @param p Number of observed variables.
#' @param h Number of hidden variables.
#' @param sparsity Real between 0 and 1. The density of the sparse graph.
#' @param sparsity.latent Real between 0 and 1. Probability of connection between any
#' observed variable and any observed hidden variable.
#' @param outlier.fraction Fraction of the samples that should be drawn from a Cauchy
#' distribution. The remaining ones are drawn from a multivariate normal with the same
#' scale matrix.
#' @return A list with keys:
#'  - obs.data: n x p data matrix of observed variables
#'  - full.data: (n+h) x p data matrix which also contains the hidden variables
#'  - true.precision.matrix: (n+h)x(n+h) matrix from which the data was sampled.
#' @import mvtnorm
#' @export
generate.data <- function(n, p, h, sparsity=0.02, sparsity.latent=0.7,
                          outlier.fraction=0) {
  library(mvtnorm)

  if( h <= 0 ) stop('The number of latent variables h must be >= 0.')

  tot.var <- p + h
  S <- matrix(runif(n=p**2, min=-1) * rbinom(n=p**2, size = 1,
                                                   prob = sparsity),
              ncol=p)
  L <- (diag(h))
  SLX <- matrix(runif(n=p*h, min=-1) * rbinom(n=p*h, size = 1,
                                             prob = sparsity.latent),
              ncol=h, nrow=p)
  S <- 0.5 * (S + t(S))
  true.prec.mat <- rbind(cbind(S, SLX), cbind(t(SLX), L))
  true.prec.mat <- true.prec.mat - min(eigen(true.prec.mat)$val) * diag(tot.var) + diag(tot.var)
  true.prec.mat <- cov2cor(true.prec.mat)
  true.cov.mat <- solve(true.prec.mat)

  #generate data
  if (outlier.fraction==0) {
    dat <- rmvnorm(n,mean=rep(0,tot.var),sigma=true.cov.mat)
  } else {
    f = 1 - outlier.fraction
    dat1 <- rmvnorm(ceiling(f * n), mean=rep(0,tot.var),sigma=true.cov.mat)
    dat2 <- rmvt(n - ceiling(f * n), sigma=true.cov.mat, df=1)
    dat <- rbind(dat1, dat2)
  }
  #remove the columns that corresponds to the latent variables
  obs.dat <- dat[,-((p+1):tot.var)]
  obs.dat <- scale(obs.dat)

  data <- list()
  data$obs.data <- obs.dat
  data$full.data <- dat
  data$true.precision.matrix <- true.prec.mat

  data
}
