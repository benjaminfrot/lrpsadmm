

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
                             SylR,
                             SylU) {
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
    out = .Call(
      '_lrpsadmm_rcppeigen_solve_sylvester',
      PACKAGE = 'lrpsadmm',
      A,
      SylR,
      SylU,
      B,
      diag(eigB$values),
      eigB$vectors,
      Q
    )
    AZXold <- AZX
    AZX = out$AZX
    # For Xest a solution of Sylvester's equation,
    # we compute Qest = A Xest + XestB and return
    # ||Q - Qest||_F.
    diff = out$diff
    d_normX <- norm(AX - AXold, 'F') ** 2
    d_normZX <- norm(AZX - AZXold, 'F') ** 2
    d_norm <- sqrt(d_normX + d_normZX)
    if ((i > 1) && (d_norm < tol)) {
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
    if (diff < tol) {
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

#' @title Compute the LRpS estimator along a path (for a fixed value of \eqn{gamma})
#' @description
#' The penalty for the LRpS estimator is written as \eqn{\lambda_1 ||S||_1 + \lambda_2 Trace(L)} in the
#' objective function of \code{lrpsadmm}.
#' This can be equivalently rewritten in terms of the regularisation parameters \eqn{\lambda} and \eqn{\gamma} as follows
#' \deqn{\lambda \gamma ||S||_1 + \lambda (1 - \gamma) Trace(L),}
#' for \eqn{\gamma \in (0, 1)}.
#' This function estimates the path of the estimator for a fixed value of \eqn{\gamma} (which controls the trade-off between
#' the two penalties) by varying the value of \eqn{\lambda}. See the documentation of \code{lrpsadmm} and references therein for
#' more details.
#' @details
#' The function \code{lrpsadmm} is fitted for successive values of \eqn{\lambda} using warm starts.
#' The sequence of values of \eqn{\lambda} can be provided directly by the user.
#' It is automatically sorted in decreasing order.
#' By default, a decreasing sequence of 20 values within a reasonable range is selected as follows.
#' We set \eqn{\lambda_{max} = \max_{ij, i \neq j} |\Sigma_{ij}|/\gamma} and
#' \eqn{\lambda_{min} = \lambda_{max}} * \code{lambda.ratio}; then
#' 20 values between \eqn{\lambda_{max}} and \eqn{\lambda_{min}} are taken following
#' a geometric progression.
#'
#' Because it does not make much sense to fit this estimator when the sparse estimate S becomes too dense
#' or if the rank of the low-rank estimate L becomes too high, the computation of the path is aborted
#' when the sparsity of S reaches \code{max.sparsity} or when the rank of L reaches \code{max.rank}.
#'
#' @param Sigma A p x p matrix. An estimate of the correlation matrix
#' @param gamma A real between 0 and 1. The value of the tuning parameter gamma in the
#' parametrisation of the penalty described avove. This is the trade-off between the sparse and trace penalties.
#' @param lambdas A decreasing sequence of values of lambda. See Details for the default value.
#' @param lambda.max A positive real. Maximum value of lambda. See Details.
#' @param lambda.ratio A real between 0 and 1. The smallest value of lambda is given by lambda.max * lambda.ratio. See Details.
#' @param n.lambdas A positive integer. The number of values of lambda to
#' generate according a geometric sequence between lambda.max and lambda.max * lambda.ratio. See Details.
#' @param max.sparsity A real between 0 and 1. Abort the computation of the path if S becomes denser than this value.
#' @param max.rank A real between 0 and 1. Abort the computuation of the path if the rank of L becomes higher than this value.
#' @param rel_tol \code{rel_tol} parameter of the \code{lrpsadmm} function.
#' @param abs_tol \code{rel_tol} parameter of the \code{lrpsadmm} function.
#' @param max.iter \code{max.iter} parameter of the \code{lrpsadmm} function.
#' @param mu \code{mu} parameter of the \code{lrpsadmm} function.
#' @param verbose A boolean. Whether to print the value of lambda, gamma, sparsity of S and rank of L after each fit.
#' @param zeros A p x p matrix with entries set to 0 or 1. Whereever its entries are
#' 0, the entries of the estimated S will be forced to 0.
#
#' @return
#'   An object of class lrpsaddmpath. This is essentially a list (see examples). Each element is itself a list with keys:
#'   \describe{
#'      \item{lambda}{Value of lambda used for that fit. Recall that the value of the tuning parameters \eqn{\lambda_1, \lambda_2}
#'      is given by \eqn{\lambda_1 = lambda * gamma} and \eqn{\lambda_2 = lambda * (1 - gamma)}}
#'      \item{gamma}{Value of gamma used for that fit. }
#'      \item{lambda1}{Corresponds to the parameter \code{l1} given as argument to the \code{lrpsadmm} function.
#'      lambda1 = lambda * gamma}
#'      \item{lambda2}{Corresponds to the parameter \code{l2} given as argument to the \code{lrpsadmm} function.
#'      lambda2 = lambda * (1 - gamma)}
#'      \item{number.of.edges}{Number of edges in the estimated sparse graphical model.}
#'      \item{rank.L}{Rank of the estimated low-rank matrix L.}
#'      \item{sparsity}{Sparsity of the estimated sparse matrix. This is fraction of entries that are non-zero.}
#'      \item{fit}{An object of class lrpsadmm.
#'      This is the outcome of calling \code{lrpsadmm} with tuning parameters l1 = lambda * gamma and l2 = lambda (1 - gamma).
#'      See the documentation of the function \code{lrpsadmm} for more information.}
#'   }
#' @examples
#'
#' set.seed(0)
#' # Generate data with a well-powered dataset
#' sim.data <- generate.latent.ggm.data(n=2000, p=100, h=5, outlier.fraction = 0.0,
#'                                      sparsity = 0.02, sparsity.latent = 0.7)
#' X <- sim.data$obs.data; Sigma <- cor(X) # Sample correlation matrix
#'
#'
#' gamma <- 0.1 # Some reasonble value for gamma
#' # We ask for 30 lambdas, but the sparse graph becomes too dense so the
#' # computation is stopped.
#' my.path <- lrpsadmm.path(Sigma = Sigma, gamma = gamma,
#'                          lambda.ratio = 1e-03, n.lambdas = 30,
#'                          verbose = TRUE, rel_tol = 1e-04, abs_tol=1e-06)
#'
#' # This time let us ask for 30 values,
#' # but let us narrow down the range by using a
#' # a smaller ratio
#' my.path <- lrpsadmm.path(Sigma = Sigma, gamma = gamma,
#'                          lambda.max = 0.96, lambda.ratio = 0.1, n.lambdas = 30,
#'                          verbose = TRUE, rel_tol = 1e-04, abs_tol=1e-06)
#'
#' # Plot some basic information about the path
#' plot(my.path)
#' # Look at the first graph in the path
#' plot(my.path[[1]]$fit)
#' # Because this is simulated data, we know the ground truth
#' # Let us use it to compute the precsion and recall metrics
#' # along the path
#' ground.truth <- sim.data$precision.matrix[1:100, 1:100]
#' # Remove the elements along the diagonal. Keep a matrix of 0s and 1s
#' ground.truth <- 1 * (( ground.truth - diag(diag(ground.truth)) ) !=0)
#' # There is a new plot with the precision / recall curve
#' plot(my.path, ground.truth = ground.truth)
#'
#' ### Let us use a robust estimator of the correlation matrix
#' # Generate data with 5% of outliers
#' set.seed(0)
#' sim.data <- generate.latent.ggm.data(n=2000, p=100, h=5, outlier.fraction = 0.05,
#'                                      sparsity = 0.02, sparsity.latent = 0.7)
#' ground.truth <- sim.data$precision.matrix[1:100, 1:100]
#' # Remove the elements along the diagonal. Keep a matrix of 0s and 1s
#' ground.truth <- 1 * (( ground.truth - diag(diag(ground.truth)) ) !=0)
#' X <- sim.data$obs.data;
#' Sigma <- cor(X) # Sample correlation matrix
#' Sigma.Kendall <- Kendall.correlation.estimator(X) # The robust estimator
#'
#' # With that many strong outliers, using the sample corr. mat.
#' # is not going to work well
#' gamma <- 0.2
#' my.path <- lrpsadmm.path(Sigma = Sigma, gamma = gamma,
#'                          lambda.ratio = 1e-02, n.lambdas = 30, verbose = TRUE)
#' # Use another estimator for the correlation matrix:
#' my.robust.path <- lrpsadmm.path(Sigma = Sigma.Kendall, gamma = gamma,
#'                                 lambda.ratio = 1e-01, n.lambdas = 30, verbose = TRUE)
#' # The output of the sample correlation path is poor (in terms of prec/recall)
#' # This is pretty much noise
#' plot(my.path, ground.truth)
#' # The Kendall estimator produces far better results.
#' # It is not affected by the 5% of outliers
#' plot(my.robust.path, ground.truth)
#' @useDynLib lrpsadmm
#' @importFrom Rcpp evalCpp
#' @import MASS RcppEigen
#' @export
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
  
  out = .Call('_lrpsadmm_rcppeigen_Schur',
              PACKAGE = 'lrpsadmm', 2 * SigmaZ)
  SylU <- out$U
  SylR <- out$R
  
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
        SylR = SylR,
        SylU = SylU
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
  
  parameters$history <- as.data.frame(history)
  
  attr(parameters, "class") <- "lscggmadmm"
  
  parameters
}