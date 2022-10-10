

lscggm_loglikelihood.R <- function(SigZ, SigZX, SigX, AX, AZX) {
  eigAX <- eigen(AX, symmetric = TRUE)
  s <- eigAX$values
  s[abs(s) < 1e-09] <- 1e-09
  AXi <- eigAX$vectors %*% (diag(1 / s)) %*% t(eigAX$vectors)
  ll <- -sum(log(s))
  ll <-
    ll + sum(diag(SigX %*% AX)) + 2 * sum(diag(SigZX %*% t(AZX)))
  ll <- ll +  sum(diag(AXi %*% t(AZX) %*% SigZ %*% AZX))
  
  return(ll)
}

lsscgm_obj_func.R <- function(SigZ,
                              SigZX,
                              SigX,
                              AX,
                              SX,
                              LX,
                              AZX,
                              SZX,
                              LZX,
                              l1,
                              l2) {
  ll <- lscggm_loglikelihood.R(SigZ, SigZX, SigX,
                               AX, AZX)
  p1 <- l1 * (sum(abs(SX)) + sum(abs(SZX)))
  L <- rbind(LX, LZX)
  p2 <- l2 * sum(svd(L)$d)
  
  return(ll + p1 + p2)
}

solve.sylvester <- function(Q, evalsA, evalsB, evecsA, evecsB) {
  
  U <- evecsA
  R <- diag(evalsA)
  V <- evecsB
  Sinv <- diag(1 / evalsB)
  m <- ncol(R)
  p <- nrow(Sinv)
  
  F <- t(U) %*% Q %*% V
  ones <- matrix(1, nrow=m, ncol=p)
  factor <- R %*% ones %*% Sinv + ones
  Ysol <- (F %*% Sinv) / factor
  Xsol <- U %*% Ysol %*% t(V)
  
  Xsol
}

lscggm_updateA.R <- function(SigZ,
                             SigZX,
                             SigX,
                             inv_sqrt_SigZ,
                             AZX,
                             SX,
                             LX,
                             UX,
                             SZX,
                             LZX,
                             UZX,
                             mu,
                             tol,
                             eigen2SigZ) {
  A = 2 * SigZ # Matrix A in the Sylvester equation.
  AX <- SX * 0
  p = dim(SigX)[1]
  m = dim(SigZ)[1]
  TX = SigX - mu * (SX - LX) + UX
  TZX = 2 * SigZX - mu * (SZX - LZX) + UZX
  ll_prev <- 0.0
  for (i in 1:20) {
    # Solve AX given AZX. That's a 2nd order polynomial
    GammaZX = 0.5 * inv_sqrt_SigZ %*% (mu * AZX + TZX)
    GTG = t(GammaZX) %*% GammaZX
    Delta = (GTG + TX) %*% (GTG + TX) + 4 * mu * diag(p)
    eigDelta = eigen(Delta, symmetric = TRUE)
    sqrtDelta = eigDelta$vectors %*% diag(sqrt(eigDelta$values)) %*% t(eigDelta$vectors)
    AXold <- AX
    AX = 0.5 * mu * (sqrtDelta - (GTG + TX))
    
    # Solve AZX given AX
    B = mu * AX
    eigB <- eigen(B, symmetric = TRUE)
    Q = -TZX %*% AX
    AZXold <- AZX
    # Pure R Solution
    AZX <- solve.sylvester(Q, eigen2SigZ$values, eigB$values, 
                           eigen2SigZ$vectors, eigB$vectors)
    
    d_normX <- norm(AX - AXold, 'F') ** 2
    d_normZX <- norm(AZX - AZXold, 'F') ** 2
    d_norm <- sqrt(d_normX + d_normZX)
    if ((i > 5) && (d_norm < tol)) {
      break()
    }
  }
  
  return(list(AX = AX, AZX = AZX))
}

lscggm_updateS.R <- function(AX, LX, UX, AZX, LZX, UZX, l1, mu) {
  X1_X = AX + LX + UX / mu
  X1_ZX = AZX + LZX + UZX / mu
  X2_X = abs(X1_X) - l1 / mu
  X2_ZX = abs(X1_ZX) - l1 / mu
  X2_X[X2_X < 0] = 0.0
  X2_ZX[X2_ZX < 0] = 0.0
  SX = sign(X1_X) * X2_X
  SZX = sign(X1_ZX) * X2_ZX
  
  return(list(SX = SX, SZX = SZX))
}

lscggm_update_L_proxg.R <- function(BX, BZX, rho) {
  p <- dim(BX)[1]
  m <- dim(BZX)[1]
  B <- rbind(BX, BZX)
  svdB <- svd(B)
  d <- svdB$d - rho
  d[d < 0] <- 0
  B <- svdB$u %*% diag(d) %*% t(svdB$v)
  BX <- B[1:p,]
  BZX <- B[(p + 1):(p + m),]
  
  return(list(BX = BX, BZX = BZX))
}

lscggm_update_L_proxf.R <- function(X) {
  eigX <- eigen(X, symmetric = TRUE)
  s <- eigX$values
  s[s < 0] <- 0
  u <- eigX$vectors
  return(u %*% diag(s) %*% t(u))
}

lscggm_updateL.R <- function(AX, SX, UX,
                             AZX, SZX, UZX,
                             Lambda2, mu,
                             tol) {
  XX = SX - AX - (UX / mu)
  XZX = SZX - AZX - (UZX / mu)
  
  lp2 <- Lambda2 / mu
  PX <- QX <- XX * 0
  PZX <- QZX <- XZX * 0
  ucond <- NULL
  for (i in 1:10) {
    out = lscggm_update_L_proxg.R(XX + PX, XZX + PZX, lp2)
    YX = out$BX
    YZX = out$BZX
    PX = XX + PX - YX
    PZX = XZX + PZX - YZX
    
    old.XX <- XX
    old.XZX <- XZX
    XX = lscggm_update_L_proxf.R(YX + QX)
    XZX <- YZX + QZX
    QX = YX + QX - XX
    QZX = YZX + QZX - XZX
    
    diffX <- norm(old.XX - XX, 'F') ** 2
    diffZX <- norm(old.XZX - XZX, 'F') ** 2
    diff <- sqrt(diffX + diffZX)
    if ((i > 2) && (diff < tol)) {
      break()
    }
  }
  
  return(list(LX = XX, LZX = XZX))
}

lscggm_updateU.R <- function(AX, SX, LX, UX,
                             AZX, SZX, LZX, UZX,
                             mu) {
  UX <- UX + mu * (AX - (SX - LX))
  UZX <- UZX + mu * (AZX - (SZX - LZX))
  
  return(list(UX = UX, UZX = UZX))
}


#'
#' Fit the conditional low-rank plus sparse estimator using the alternating method direction of multipliers
#' @description
#' The estimator fits a low-rank plus sparse estimator to a set of observed variables X
#' conditionally on a set of variables Z. 
#' Given a n x p data matrix X with empirical correlation matrix
#' \eqn{\Sigma_X}, and an n x m data matrix Z with empirical correlation matrix \eqn{\Sigma_Z}, 
#' an alternating direction method of multipliers (ADMM) algorithm is fitted to obtain the solutions \eqn{S, L} to
#' \deqn{(S, L) = argmin_{A, B}  -loglikelihood(A, B, ; \Sigma_Z, \Sigma_X) + \lambda_1 ||A||_1 + \lambda_2 ||B||_\ast),}
#' subject to \eqn{A_X - B_X > 0} and \eqn{B_X \ge 0}.
#' loglikelihood is the log-likelihood of a conditional Gaussian graphical model (see reference
#' for notations and further details).
#' @param SigmaX An estimate of the correlation matrix of X. p x p matrix.
#' @param SigmaZX An estimate of the cross-correlation matrix between Z and X. p x m matrix.
#' @param SigmaZ An estimate of the correlation matrix of Z. m x m matrix.
#' @param Lambda1 Penalty on the l1 norm of S
#' @param Lambda2 Penalty on the sum of the singular values of L
#' @param init The output of a previous run of the algorithm. For warm starts. Defaults to NULL.
#' @param maxiter Maximal number of iterations
#' @param mu Step size of the ADMM algorithm.
#' @param rel_tol Relative tolerance required to stop the algorithm. The algorithm
#' stops when both the change in parameters is below tolerance and the constraints
#' are satisfied. Default 1e-02.
#' @param abs_tol Absolute tolerance required to stop the algorithm. Default 1e-04.
#' @param print_progress Whether the algorithm should report on its progress.
#' @param print_every How often should the algorithm report on its progress (in terms
#' of #iterations).
#' @details
#' Given a n x p data matrix X with empirical correlation matrix
#' \eqn{\Sigma_X}, and an n x m data matrix Z with empirical correlation matrix \eqn{\Sigma_Z}, 
#' an alternating direction method of multipliers (ADMM) algorithm is fitted to obtain the solutions \eqn{S, L} to
#' \deqn{(S, L) = argmin_{A, B}  -loglikelihood(A, B, ; \Sigma_Z, \Sigma_X) + \lambda_1 ||A||_1 + \lambda_2 ||B||_\ast),}
#' subject to \eqn{A_X - B_X > 0} and \eqn{B_X \ge 0}.
#' loglikelihood is the log-likelihood of a conditional Gaussian graphical model (see reference
#' for notations and further details).
#' This is the estimator suggested in Frot, Jostins and McVean.
#'
#' The optimisation problem is decomposed as a three-block ADMM optimisation problem, as described in Ye et al.
#' Because it is a so-called consensus problem, the ADMM is guaranteed to converge.
#'
#' The tuning parameters \eqn{\lambda_1} and \eqn{\lambda_2} are typically reparametrised as
#' \eqn{\lambda_1 = \lambda \gamma} and \eqn{\lambda_2 = \lambda (1 - \gamma)}, for \eqn{\gamma \in (0,1)}.
#' Here, for a fixed \eqn{\gamma}, \eqn{\lambda} controls the overall shrinkage along the path defined by \eqn{\gamma}.
#' \eqn{\gamma} controls the trade off on the penalties between sparse and low-rank components.
#'
#' For numerical stability, a smaller value of \eqn{\gamma} is preferable when the number of samples n is close to p. 
#' See examples.
#'
#' @return
#' An S3 object of class lscggm It is essentially a list with keys:
#' \describe{
#'  \item{SX}{A p x p matrix. The sparse estimate SX.}
#'  \item{LX}{A p x p matrix. The low-rank estimate LX.}
#'  \item{SZX}{A p x p matrix. The sparse estimate SZX.}
#'  \item{LZX}{A p x p matrix. The low-rank estimate LZX.}
#'  \item{UX}{A p x p matrix. Augmented Lagrangian multiplier. It is stored in order to allow warm starts.}
#'  \item{AX}{A p x p matrix. Used interall. Stored to allow warm starts.}
#'  \item{UZX}{A p x p matrix. Augmented Lagrangian multiplier. It is stored in order to allow warm starts.}
#'  \item{AZX}{A p x p matrix. Used internally. Stored to allow warm starts}

#'  \item{termcode}{An integer. Its value determines whether the algorithm terminated normally or with an error.
#'  0: Convergence reached. -1: Maxiter reached. -2: Shrinkage too strong.}
#'  \item{termmsg}{A character vector. The message corresponding to the value \code{termcode}.}
#'  \item{history}{A numerical dataframe with the objective function at each
#'  iterations, the norm and dual norm as well the primal and dual tolerance
#'  criteria for convergence. The algorithm exits when r_norm < eps_pri and s_norm < eps_dual.}
#' }
#'
#' @references
#' B Frot, L Jostins, G McVean. Graphical model selection for Gaussian conditional random fields in the presence of latent variables
#' Journal of the American Statistical Association 114 (526), 723-734
#' 
#' Chandrasekaran, Venkat; Parrilo, Pablo A.; Willsky, Alan S. Latent variable graphical model selection via convex optimization.
#' Ann. Statist. 40 (2012), no. 4, 1935--1967. doi:10.1214/11-AOS949. \url{https://projecteuclid.org/euclid.aos/1351602527}
#'
#' Gui-Bo Ye, Yuanfeng Wang, Yifei Chen, Xiaohui Xie;
#' Efficient Latent Variable Graphical Model Selection via Split Bregman Method.
#' \url{https://arxiv.org/abs/1110.3076}
#'
#' @examples
#' set.seed(1)
#' # Generate data for a well-powered dataset
#' sim.data <- generate.latent.ggm.data(n=2000, p=100, h=5, outlier.fraction = 0.0,
#'                                      sparsity = 0.02, sparsity.latent = 0.7)
#' X <- sim.data$obs.data; Sigma <- cor(X) # Sample correlation matrix
#'
#' @export
#' @seealso lscggm.cv lscggm.path
#' @import MASS
lscggmadmm <- function(SigmaZ,
                       SigmaZX,
                       SigmaX,
                       Lambda1,
                       Lambda2,
                       init = NULL,
                       maxiter = 1000,
                       mu = 1.0,
                       abs_tol = 1e-04,
                       rel_tol = 1e-02,
                       print_progress = TRUE,
                       print_every = 10) {
  # Compute a "pseudo inverse square root" of SigmaZ
  # According to the theory, SigmaZ should not be singular
  # but in practice it might be different
  Zeig = eigen(SigmaZ, symmetric = TRUE)
  s <- Zeig$values
  s[abs(s) < 1e-09] <- 1e-09
  s = sqrt(s)
  s[s != 0] = 1.0 / (s[s != 0])
  inv_sqrt_SigZ = Zeig$vectors %*% diag(s) %*% t(Zeig$vectors)
  
  eigen2SigZ <- eigen(2 * SigmaZ, symmetric = T)
  
  p <- dim(SigmaX)[1]
  m <- dim(SigmaZ)[1]
  
  if (is.null(init)) {
    SX <- diag(p)
    LX <- SX * 0.0
    UX <- SX * 0.0
    SZX <- matrix(0, nrow = m, ncol = p)
    LZX <- SZX * 0.0
    UZX <- SZX * 0.0
    AZX <- SZX * 0.0
    AX <- SX
  } else {
    parameters <- init
    SX <- init$SX
    LX <- init$LX
    UX <- init$UX
    SZX <- init$SZX
    LZX <- init$LZX
    UZX <- init$UZX
    AZX <- init$AZX
    AX <- init$AX
  }
  
  history <- matrix(NA, nrow = 0, ncol = 6)
  colnames(history) <-
    c('Iteration',
      'Objval',
      's_norm',
      'r_norm',
      'eps_pri',
      'eps_dual')
  parameters <- list()
  parameters$termcode <- -1
  parameters$termmsg <- "Maximum number of iterations reached."
  eps_dual <- rel_tol
  learning_rate <- 1.0
  for (i in 1:maxiter) {
    # Update A
    out <-
      lscggm_updateA.R(
        SigZ = SigmaZ,
        SigZX = SigmaZX,
        SigX = SigmaX,
        AZX = AZX,
        SX = SX,
        SZX = SZX,
        LX = LX,
        LZX = LZX,
        UX = UX,
        UZX = UZX,
        mu = mu,
        inv_sqrt_SigZ = inv_sqrt_SigZ,
        tol = eps_dual,
        eigen2SigZ = eigen2SigZ
      )
    AX = out$AX
    AZX = out$AZX
    learning_rate = out$L * 10
    
    # Update S
    SX_old <- SX
    SZX_old <- SZX
    out <- lscggm_updateS.R(AX, LX, UX,
                            AZX, LZX, UZX,
                            Lambda1, mu)
    SX <- out$SX
    SZX <- out$SZX
    if (any(diag(SX) == 0)) {
      SX <- diag(p)
      LX <- UX <- SX * 0
      SZX <- SZX * 0
      AZX <- UZX <- LZX <- SZX
      parameters$termcode <- -2
      parameters$termmsg <-
        "Shrinkage too strong: sparse component is empty."
      break()
    }
    
    # Update L
    LX_old <- LX
    LZX_old <- LZX
    out <- lscggm_updateL.R(AX, SX, UX, AZX, SZX, UZX, Lambda2, mu,
                            tol = eps_dual)
    LX <- out$LX
    LZX <- out$LZX
    
    # Update U
    out <- lscggm_updateU.R(AX, SX, LX, UX,
                            AZX, SZX, LZX, UZX,
                            mu)
    UX <- out$UX
    UZX <- out$UZX
    
    # Diagnostics
    objval <- lsscgm_obj_func.R(SigmaZ,
                                SigmaZX,
                                SigmaX,
                                AX,
                                SX,
                                LX,
                                AZX,
                                SZX,
                                LZX,
                                Lambda1,
                                Lambda2)
    
    r_normX <- norm(AX - (SX - LX), 'F') ** 2
    s_normX <- norm(mu * ((SX - LX) - (SX_old - LX_old)), 'F') ** 2
    r_normZX <- norm(AZX - (SZX - LZX), 'F') ** 2
    s_normZX <- norm(mu * ((SZX - LZX) - (SZX_old - LZX_old)), 'F') ** 2
    r_norm <- sqrt(r_normX + r_normZX)
    s_norm <- sqrt(s_normX + s_normZX)
    
    tmp1 <- max(sqrt(norm(AX, 'F') ** 2 + norm(AZX, 'F') ** 2),
                sqrt(norm(SX - LX, 'F') ** 2 + norm(SZX - LZX, 'F') ** 2))
    eps_pri <- (p + sqrt(p * m)) * abs_tol + rel_tol * tmp1
    tmp2 <- sqrt(norm(mu * UX, 'F') ** 2 + norm(mu * UZX, 'F') ** 2)
    eps_dual <- (p + sqrt(p * m)) * abs_tol + rel_tol * tmp2
    history <- rbind(history, c(i, objval, s_norm,
                                r_norm, eps_pri, eps_dual))
    
    if ((s_norm < eps_dual) && (r_norm < eps_pri)) {
      parameters$termcode <- 0
      parameters$termmsg <- 'Convergence Reached.'
      break()
    }
    if ((print_progress) & (i > print_every)) {
      if ((i %% print_every) == 0) {
        print(paste(
          c(
            "Iteration:",
            "Obj. fun.:",
            "s_norm:",
            "r_norm:",
            "eps_pri:",
            "eps_dual:"
          ),
          history[i, ]
        ))
      }
    }
    
  }
  
  parameters$SX <- SX
  parameters$LX <- LX
  parameters$UX <- UX
  parameters$SZX <- SZX
  parameters$LZX <- LZX
  parameters$UZX <- UZX
  parameters$AZX <- AZX
  parameters$AX <- AX
    
  parameters$history <- as.data.frame(history)
  
  attr(parameters, "class") <- "lscggmadmm"
  
  parameters
}