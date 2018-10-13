
.obj_func <- function(ps, opts) {
  p <- opts$p
  n <- opts$n
  C <- opts$Sigma
  l1 <- opts$lp1 * opts$mu
  l2 <- opts$lp2 * opts$mu
  X <- ps$S - ps$L
  X <- 0.5 * (X + t(X))
  evals <- eigen(X, symmetric = T)$val
  evals[abs(evals) < 1e-08] <- 1e-16
  if (any(evals < 0)) {
    return(NaN)
  }
  # Log-Likelihood
  ll <- sum(diag(C %*% X)) - sum(log(evals))
  # + Penalty
  ll <- ll + l1 * sum(abs(ps$S)) + l2 * sum(diag(ps$L))

  ll
}

.updateA <- function(ps, opts) {
  lp1 <- opts$lp1
  lp2 <- opts$lp2
  mu <- opts$mu
  C <- opts$Sigma
  Shat <- ps$Shat
  Lhat <- ps$Lhat
  Uhat <- ps$Uhat
  n <- opts$n
  p <- opts$p

  X1 <- mu * (Shat - Lhat) - C - Uhat
  X2 <- X1 %*% X1 + 4 * mu * diag(dim(X1)[1])
  X2 <- 0.5 * (X2 + t(X2))
  eig <- eigen(X2, symmetric = T)
  sqrtX2 <- eig$vectors %*% diag(sqrt(eig$values)) %*% t(eig$vectors)
  A <- (X1 + sqrtX2) / (2 * mu)
  A <- 0.5 * (A + t(A))

  ps$A <- A
  ps$exit <- F
  ps
}

.updateAlpha <- function(ps, opts) {
  alpha <- ps$alpha[length(ps$alpha)]
  alpha <- 0.5 * (1 + sqrt(1 + 4 * alpha**2))
  ps$alpha <- c(ps$alpha, alpha)

  ps
}

.updateL <- function(ps, opts) {
  lp1 <- opts$lp1
  lp2 <- opts$lp2
  mu <- opts$mu
  C <- opts$Sigma
  A <- ps$A
  S <- ps$S
  Uhat <- ps$Uhat
  ps$prevL <- ps$L

  X1 <- S - A - (Uhat / mu)
  X1 <- 0.5 * (X1 + t(X1))
  eig <- tryCatch({
    RSpectra::eigs_sym(X1, opts$max.rank)
    }, 
    error = function(e) {
      print(e)
      print('hi')
      ps$exit = T
      return(ps)
  }
  )
  eigVal <- eig$values - lp2
  eigVal[eigVal < 0] <- 0
  L <- eig$vectors %*% diag(eigVal) %*% t(eig$vectors)
  L <- 0.5 * (L + t(L))

  ps$L <- L
  ps$exit <- F
  ps
}

.updateS <- function(ps, opts) {
  lp1 <- opts$lp1
  lp2 <- opts$lp2
  mu <- opts$mu
  C <- opts$Sigma
  A <- ps$A
  L <- ps$L
  Uhat <- ps$Uhat
  ps$prevS <- ps$S

  X1 <- A + L + (Uhat / mu)
  X2 <- abs(X1) - lp1
  X2[X2 <= 0] <- 0
  S = X2 * sign(X1)
  S <- S * opts$zeros
  S <- 0.5 * (S + t(S))

  ps$S <- S
  ps$exit <- F
  ps
}

.updateShatLhatUhat <- function(ps, opts) {
  S <- ps$S
  pS <- ps$prevS
  L <- ps$L
  pL <- ps$prevL
  U <- ps$U
  pU <- ps$prevU
  alpha_k_min_1 <- ps$alpha[length(ps$alpha)-2]
  alpha_k_plus_1 <- ps$alpha[length(ps$alpha)]

  Shat <- S + (alpha_k_min_1 / alpha_k_plus_1) * (S - pS)
  Lhat <- L + (alpha_k_min_1 / alpha_k_plus_1) * (L - pL)
  Uhat <- U + (alpha_k_min_1 / alpha_k_plus_1) * (U - pU)

  ps$Shat <- Shat
  ps$Lhat <- Lhat
  ps$Uhat <- Uhat

  ps
}

.updateU <- function(ps, opts) {
  lp1 <- opts$lp1
  lp2 <- opts$lp2
  mu <- opts$mu
  C <- opts$Sigma
  A <- ps$A
  S <- ps$S
  L <- ps$L
  Uhat <- ps$Uhat
  ps$prevU <- ps$U

  U <- Uhat + mu * (A - S + L)
  U <- 0.5 * (U + t(U))

  ps$U <- U
  ps$exit <- F
  ps
}


.split.bregman.low.rank.plus.sparse.update.parameters <- function(ps, opts) {

  pL <- ps$L
  pS <- ps$S
  pU <- ps$U

  fs <- c(.updateA, .updateS, .updateL, .updateU)
  for (f in fs) {
    ps <- f(ps, opts)
    if(ps$exit) {
      break()
    }
  }

  pck <- ps$ck
  ck <- (1.0 / opts$mu) * sum((ps$Uhat - ps$U)**2) +
    opts$mu * (sum((ps$Shat - ps$Lhat - ps$S + ps$L)**2))
  if (ck < (opts$eta * pck)) { # Do not restart
    ps$restarts <- c(ps$restarts, 0)
    ps$ck <- ck
    ps <- .updateAlpha(ps, opts)
    ps <- .updateShatLhatUhat(ps, opts)
  } else { # Then we restart
    ps$restarts <- c(ps$restarts, 1)
    ps$ck <- pck / (opts$eta)
    ps$Shat <- pS
    ps$Uhat <- pU
    ps$Lhat <- pL
    ps$alpha <- c(ps$alpha, 1.0)
  }

  # One problem with restarting when the value of mu is not
  # good is that it "flip-flops": alternating restart and non-restart
  # endlessly. When that happens 100 times in a row, we reduce the value of mu.
  if ((length(ps$restarts) > 200)) {
    if (all(ps$restarts[(length(ps$restarts)-99):length(ps$restarts)] ==
            rep(c(0,1), 50))) {
      # then we know it flip-flops.
      msg1 = "The accelerated ADDM algorithm restarted 100 times in a row."
      msg2 = "It seems like the required convergence tolerance cannot be achieved."
      msg3 = "This might due to the problem being too ill-posed. For example, the value of gamma might be too large given the ratio p/n."
      msg4 = "You might also want to consider trying a smaller value of mu."
      warning(paste(msg1, msg2, msg3, msg4))
      ps$exit <- TRUE
    }
  }

  list(ps=ps, opts=opts)
}

#'
#' Fit the the low-rank plus sparse estimator
#' @details Given a n x p data matrix X and its empirical correlation matrix
#' \eqn{\Sigma}, an alternative direction method of multipliers (ADMM) algorithm
#' is used to obtain the solutions \eqn{S, L} to
#' \deqn{(S, L) = argmin_{A, B}  -logdet(A - B) + Tr((A-B)\Sigma) + \lambda_1 ||A||_1 + \lambda_2Tr(B),}
#' subject to \eqn{A - B > 0} and \eqn{B \ge 0}.
#' @param Sigma An estimate of the correlation matrix.
#' @param Lambda1 Penalty on the l1 norm of S
#' @param Lambda2 Penalty on the sum of the eigenvalues of L
#' @param n Number of samples
#' @param init The output of a previous run of the algorithm. For warm starts.
#' @param maxiter Maximal number of iterations
#' @param mu Stepsize of the ADMM algorithm
#' @param tol Relative tolerance required to stop the algorithm. Algorithm is stopped
#' when either the relative change in log-likelihood is below this threshold. or the
#' relative change in the Frobenius norm of the parameters is below this threshold.
#' @param eta Restart parameter in the accelerated ADMM
#' @param print_progress Whether the algorithm should report on its progress
#' @param print_every How often should the algorithm report on its progress (in terms
#' of #iterations).
#' @param zeros A p x p matrix with entries set to 0 or 1. Whereever its entries are
#' 0, the entries of the estimated S will be set to 0.
#' @param max.rank The maximum rank of the estimated L. Especially useful if n > p.
#' By default it is min(p-1, n-1).
#' @export
#' @import matrixcalc RSpectra MASS
fit.low.rank.plus.sparse <- function(Sigma, Lambda1, Lambda2, n, init=NULL,
                                     maxiter=2000, mu=0.1, tol=1e-05, eta=0.999,
                                     print_progress=T, print_every=20,
                                     zeros=NULL, max.rank=NA) {
  library(matrixcalc)
  library(RSpectra)
  library(MASS)

  if (is.null(zeros)) {
    zeros <- 1
  }
  p <- dim(Sigma)[1]
  if(is.na(max.rank)) {
    max.rank <- ceiling(min(n-1, p-1))
  }

  options <- list()
  options$mu <- mu
  options$Sigma <- Sigma
  options$lp1 <- Lambda1 / mu
  options$lp2 <- Lambda2 / mu
  options$eta <- eta
  options$zeros <- zeros
  options$max.rank <- max.rank
  options$n <- n
  options$p <- p

  if (is.null(init)) {
    S <- MASS::ginv(Sigma)
    L <- S * 0.01
    A <- S - L #+ diag(rep(1, p))
    U <- mu * (A - S + L)
    parameters <- list(S=S, L=L, A=A, U=U)
  } else {
    parameters = init
  }

  parameters$Shat <- parameters$S
  parameters$Lhat <- parameters$L
  parameters$alpha <- c(1.0, 1.0)
  parameters$restarts <- c()
  parameters$Uhat <- parameters$U
  parameters$ck <- p**2
  parameters$exit <- FALSE
  parameters$S <- parameters$S * options$zeros

  diffs <- c()
  lls <- c(NaN)
  parameters$termcode <- -1
  parameters$termmsg <- "Maximum number of iterations reached."
  for (i in 1:maxiter) {
    L <- .split.bregman.low.rank.plus.sparse.update.parameters(parameters, options)
    new_parameters <- L$ps
    options <- L$opts

    if(parameters$exit) {
      parameters$termcode <- -3
      parameters$termmsg <- "Algorithm was restarted 100 times without improvement."
      break()
    }

    if(sum(new_parameters$S) == 0) {
      parameters$S <- diag(p)
      parameters$termcode <- -2
      parameters$termmsg <- "Shrinkage too strong: sparse component is empty."
      break()
    }
    # Compute the relative change in parameters with the previous iteration
    if ((matrixcalc::frobenius.norm(parameters$L) > 0)) {
      diff <- matrixcalc::frobenius.norm(new_parameters$S - parameters$S) / matrixcalc::frobenius.norm(parameters$S) +
        matrixcalc::frobenius.norm(new_parameters$L - parameters$L) / matrixcalc::frobenius.norm(parameters$L)
    } else {
      diff <- matrixcalc::frobenius.norm(new_parameters$S - parameters$S) / matrixcalc::frobenius.norm(parameters$S) +
        matrixcalc::frobenius.norm(new_parameters$L - parameters$L)
    }
    diffs <- c(diffs, diff)

    parameters <- new_parameters
    if (diff < tol) {
      parameters$termcode <- 0
      parameters$termmsg <- "Convergence reached."
      break()
    }
    # Compute the objective function and its relative change
    lls <- c(lls, .obj_func(parameters, options))
    if (!is.nan(lls[length(lls)]) & !is.nan(lls[length(lls)-1])) {
      ll_m_1 <- lls[length(lls)-1]
      ll <- lls[length(lls)]
      if ((abs(ll - ll_m_1) / abs(ll_m_1)) < tol) {
        parameters$termcode <- 0
        parameters$termmsg <- "Convergence reached."
        break()
      }
    }
    if ((print_progress) & (i > print_every)) {
      if((i %% print_every) == 0) {
        last_lls <- lls[(i-print_every):i]
        avg_chg <- mean(abs(diff(last_lls)) / abs(last_lls[1:(print_every)]),na.rm=T)
        print(paste("Iteration:", i,
                    "Log-Likelihood:", ll,
                    "Avg relative log-likelihood change:", avg_chg,
                    "Relative change in parameters:", diff))
      }
    }
  }

  parameters$iter <-i
  parameters$diffs <- diffs
  parameters$lls <- lls

  # Clean up the variables that are specific to the algo
  parameters$Shat <- NULL
  parameters$Lhat <- NULL
  parameters$alpha <- NULL
  parameters$restarts <- NULL
  parameters$Uhat <- NULL
  parameters$ck <- NULL
  parameters$exit <- NULL
  parameters$prevS <- NULL
  parameters$prevL <- NULL
  parameters$prevU <- NULL

  parameters
}
