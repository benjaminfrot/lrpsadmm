
#' @title Simulate data from a Gaussian graphical model with hidden variables
#' @description
#'   A partitioned inverse covariance (precision) matrix K is constructed as \eqn{K := [K_X, K_{XL}; K_{XL}^T, K_L]}, where \eqn{K_X} is a p x p sparse matrix,
#'   \eqn{K_{LX}} is a p x h matrix connecting the h hidden variables to observed ones and \eqn{K_L} is a h x h diagonal matrix.
#'
#'   n samples are then drawn from a multivariate normal distribution: \eqn{N(0, K^{-1})} and only the first p variables are observed.
#'
#'   This function is here to demonstrate the features of the package.
#' @param n Number of samples.
#' @param p Number of observed variables.
#' @param h Number of hidden variables.
#' @param sparsity Real between 0 and 1. The density of the sparse graph (encoded by \eqn{K_X}).
#' @param sparsity.latent Real between 0 and 1. Probability of connection between any
#' observed variable and any hidden variable (\eqn{K_{LX}}).
#' @param outlier.fraction Fraction of the samples that should be drawn from a Cauchy
#' distribution. The remaining ones are drawn from a multivariate normal with the same
#' scale matrix.
#' @return A list with keys:
#'   \describe{
#'     \item{obs.data}{n x p data matrix of observed variables.}
#'     \item{full.data}{(n+h) x p data matrix which also contains the hidden variables.}
#'     \item{precision.matrix}{ (n+h)x(n+h) matrix from which the data was sampled. Rows/Columns indexed by 1:p correspond to observed variables. The remaining h rows/cols are the hidden ones.}
#'     \item{sparsity}{Realised sparsity of the precision matrix restricted to observed variables.}
#'  }
#' @examples
#'  n <- 2000 # Number of samples
#'  p <- 100 # Number of variables
#'  h <- 5 # Number of hidden variables
#'  sim.data <- generate.latent.ggm.data(n=n, p=p, h=h, outlier.fraction = 0.0,
#'                                sparsity = 0.02, sparsity.latent = 0.7)
#'  true.S <- sim.data$precision.matrix[-((p+1):(p+h)),-((p+1):(p+h))] # The sparse matrix
#'  observed.data <- sim.data$obs.data
#'
#'  # Generate data with 10 of samples drawn from a Cauchy
#'  sim.data <- generate.latent.ggm.data(n=n, p=p, h=h, outlier.fraction = 0.1,
#'                                sparsity = 0.02, sparsity.latent = 0.7)
#'  true.S <- sim.data$precision.matrix[-((p+1):(p+h)),-((p+1):(p+h))] # The sparse matrix
#'  observed.data <- sim.data$obs.data
#'
#' @import mvtnorm
#' @importFrom stats runif rbinom
#' @export
simulate.latent.ggm.data <- function(n, p, h, sparsity=0.02, sparsity.latent=0.7,
                          outlier.fraction=0) {
  sparsity <- sparsity * 0.5
  if( h <= 0 ) stop('The number of latent variables h must be > 0.')

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
  true.prec.mat <- true.prec.mat -
    min(eigen(true.prec.mat)$val) * diag(tot.var) + diag(tot.var)
  true.prec.mat <- cov2cor(true.prec.mat)
  true.cov.mat <- solve(true.prec.mat)

  #generate data
  if (outlier.fraction==0) {
    dat <- rmvnorm(n,mean=rep(0,tot.var),sigma=true.cov.mat)
  } else {
    f <- 1 - outlier.fraction
    dat1 <- rmvnorm(ceiling(f * n), mean=rep(0,tot.var),sigma=true.cov.mat)
    dat2 <- rmvt(n - ceiling(f * n), sigma=true.cov.mat, df=1)
    dat <- rbind(dat1, dat2)
  }
  #remove the columns that corresponds to the latent variables
  obs.dat <- dat[,-((p+1):tot.var)]
  obs.dat <- scale(obs.dat)

  true.S <- true.prec.mat[-((p+1):(p+h)),-((p+1):(p+h))]
  true.S <- (true.S!=0) - diag(diag(true.S!=0))
  realised.sparsity <- 0.5 * sum(true.S) / (choose(p,2))
  data <- list()
  data$obs.data <- obs.dat
  data$full.data <- dat
  data$precision.matrix <- true.prec.mat
  data$sparsity <- realised.sparsity
  data
}