generate.latent.conditional.ggm.data <- function(n, m, p, rankLX, rankLZX, sparsity=0.02,
                                                 sparsity.latent=0.7) {

  if( rankLX <= 0 ) stop('The rank of LX must be >= 0.')
  if( rankLZX <= 0 ) stop('The rank of LZX must be >= 0.')

  # Start by generating a matrix of the form SX - LX
  tot.var <- p + rankLX
  S <- matrix(runif(n=p**2, min=-1) * rbinom(n=p**2, size = 1,
                                             prob = sparsity),
              ncol=p)
  L <- (diag(rankLX))
  SLX <- matrix(runif(n=p*rankLX, min=-1) * rbinom(n=p*rankLX, size = 1,
                                              prob = sparsity.latent),
                ncol=rankLX, nrow=p)
  S <- 0.5 * (S + t(S))
  true.prec.mat <- rbind(cbind(S, SLX), cbind(t(SLX), L))
  true.prec.mat <- true.prec.mat -
    min(eigen(true.prec.mat)$val) * diag(tot.var) + diag(tot.var)
  true.prec.mat <- cov2cor(true.prec.mat)
  true.cov.mat <- solve(true.prec.mat)

  SX <- true.prec.mat[1:p, 1:p]
  LZX <- true.prec.mat[(p+1):(p+rankLX), 1:p]
  LX <- true.prec.mat[(p+1):(p+rankLX), (p+1):(p+rankLX)]
  print(dim(LX))
  print(dim(LZX))
  L <- t(LZX) %*% solve(LX) %*% LZX

  SZX <- matrix(runif(n=p*m, min=-1) * rbinom(n=p*m, size = 1,
                                               prob = sparsity),
                ncol=p)
  # Some p * m rank rankLZX matrix
  U <- qr.Q(qr(rmvnorm(n=m, sigma = diag(m))))
  V <- qr.Q(qr(rmvnorm(n=p, sigma = diag(p))))
  D <- matrix(0, nrow=m, ncol=p)
  D[1, 1] <- runif(n=1, min = 0.1, max=1)
  LZX <- U %*% D %*% t(V)
  if (rankLZX <= 2) {
    for (i in 2:rankLZX) {
      U <- qr.Q(qr(rmvnorm(n=m, sigma = diag(m))))
      V <- qr.Q(qr(rmvnorm(n=p, sigma = diag(p))))
      D <- matrix(0, nrow=m, ncol=p)
      D[1, 1] <- runif(n=1, min = 0.1, max=1)
      LZX <- LZX + U %*% D %*% t(V)
    }
  }

  Z <- rmvnorm(n=n, sigma=diag(m)) # All the Zs are independent
  X <- matrix(0, nrow=n, ncol=p)
  iS <- solve(SX - L)
  for (i in 1:n) {
    z <- matrix(Z[i,], ncol=1)
    mu = -iS %*% t(SZX - LZX) %*% z
    X[i,] <- rmvnorm(n=1, mean=mu, sigma=iS)
  }

  data <- list()
  data$true.SX <- SX
  data$true.LX <- L
  data$true.LZX <- LZX
  data$true.SZX <- SZX
  data$Z <- Z
  data$X <- X
  data
}
