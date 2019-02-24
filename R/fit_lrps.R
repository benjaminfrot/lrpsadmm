
.obj_func <- function(ps, opts) {
  p <- opts$p
  n <- opts$n
  C <- opts$Sigma
  l1 <- opts$lp1 * opts$mu
  l2 <- opts$lp2 * opts$mu
  X <- ps$S - ps$L
  X <- 0.5 * (X + t(X))
  evals <- eigen(X, symmetric = TRUE)$val
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
  eig <- eigen(X2, symmetric = TRUE)
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
    ps$exit <- TRUE
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
  S <- X2 * sign(X1)
  S <- S * opts$zeros
  S <- 0.5 * (S + t(S))

  ps$S <- S
  ps$exit <- FALSE
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
  ps$exit <- FALSE
  ps
}


.split.bregman.low.rank.plus.sparse.update.parameters <- function(ps, opts) {

  pL <- ps$L
  pS <- ps$S
  pU <- ps$U
  ps$exit <- FALSE
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
      msg1 <- "The accelerated ADMM algorithm restarted 100 times in a row."
      msg2 <- "It seems like the required convergence tolerance cannot be achieved."
      msg3 <- "The problem might be too ill-posed. For example, the value of gamma might be too large given the ratio p/n."
      msg4 <- "You can also consider rying a smaller value of mu.\n"
      warning(paste(msg1, msg2, msg3, msg4))
      ps$exit <- TRUE
    }
  }

  list(ps=ps, opts=opts)
}

#'
#' Fit the the low-rank plus sparse estimator using an accelerated alternating method direction of multipliers
#' @description
#' Given a n x p data matrix X and its empirical correlation matrix
#' \eqn{\Sigma}, an alternating direction method of multipliers (ADMM) algorithm
#' is used to obtain the solutions \eqn{S, L} to
#' \deqn{(S, L) = argmin_{A, B}  -logdet(A - B) + Tr((A-B)\Sigma) + \lambda_1 ||A||_1 + \lambda_2Tr(B),}
#' subject to \eqn{A - B > 0} and \eqn{B \ge 0}.
#' @param Sigma An estimate of the correlation matrix.
#' @param Lambda1 Penalty on the l1 norm of S
#' @param Lambda2 Penalty on the sum of the eigenvalues of L
#' @param n Number of samples. Giving its value is useful when n < p.
#' @param init The output of a previous run of the algorithm. For warm starts.
#' @param maxiter Maximal number of iterations
#' @param mu Stepsize of the ADMM algorithm.
#' @param tol Relative tolerance required to stop the algorithm. Algorithm is stopped
#' when either the relative change in log-likelihood is below this threshold, or the
#' relative change in the Frobenius norm of the parameters is below this threshold.
#' @param eta Restart parameter in the accelerated ADMM.
#' @param print_progress Whether the algorithm should report on its progress.
#' @param print_every How often should the algorithm report on its progress (in terms
#' of #iterations).
#' @param zeros A p x p matrix with entries set to 0 or 1. Whereever its entries are
#' 0, the entries of the estimated S will be forced to 0.
#' @details
##' Given a n x p data matrix X and its empirical correlation matrix
#' \eqn{\Sigma}, an alternating direction method of multipliers (ADMM) algorithm
#' is used to obtain the solutions \eqn{S, L} to
#' \deqn{(S, L) = argmin_{A, B}  -logdet(A - B) + Tr((A-B)\Sigma) + \lambda_1 ||A||_1 + \lambda_2Tr(B),}
#' subject to \eqn{A - B > 0} and \eqn{B \ge 0}.
#' This is the estimator suggested in Chandrasekaran et al.
#'
#' The optimisation problem is decomposed as a three-block ADMM optimisation problem, as described in Ye et al.
#' Because it is a so-called consensus problem, the ADMM is guaranteed to converge.
#' Given this decomposition of the problem, we use the "Fast ADMM with Restart" (Alg. 8) of Goldstein et al.
#' in order to fit the problem (see references).
#'
#' The tuning parameters \eqn{\lambda_1} and \eqn{\lambda_2} are typically reparametrised as
#' \eqn{\lambda_1 = \lambda \gamma} and \eqn{\lambda_2 = \lambda (1 - \gamma)}, for \eqn{\gamma \in (0,1)}.
#' Here, for a fixed \eqn{\gamma}, \eqn{\lambda} controls the overall shrinkage along the path defined by \eqn{\gamma}.
#' \eqn{\gamma} controls the tradeoff on the penalties between sparse and low-rank components.
#'
#' For numerical stability, a smaller value of \eqn{\gamma} is preferable when n is close to p. See examples.
#'
#' @return
#' An S3 object of class lrpsadmm. It is essentially a list with keys:
#' \describe{
#'  \item{S}{A p x p matrix. The sparse estimate S}
#'  \item{L}{A p x p matrix. The low-rank estimate L}
#'  \item{termcode}{An integer. Its value determines whether the algorithm terminated normally or with an error.
#'  0: Convergence reached. -1: Maxiter reached. -2: Shrinkage too strong. -3: Too many restarts in a row.}
#'  \item{termmsg}{A character vector. The message corresponding to the value \code{termcode}.}
#'  \item{iter}{An integer. Number of iterations until convergence.}
#'  \item{diffs}{A vector. Relative difference in change of parameters (in Frobenius norm) at each iteration.}
#'  \item{lls}{A vector. Values taken by the objective function at each iteration.}
#'  \item{A}{A p x p matrix. A variable used by the algorithm. Its value should be close to S - L.
#'  It is stored and returned by the algorithm in order to allow warm starts.}
#'  \item{U}{A p x p matrix. Augmented Lagrangian multiplier. It is stored in order to allow warm starts.}
#' }
#'
#' @references
#' Chandrasekaran, Venkat; Parrilo, Pablo A.; Willsky, Alan S. Latent variable graphical model selection via convex optimization.
#' Ann. Statist. 40 (2012), no. 4, 1935--1967. doi:10.1214/11-AOS949. \url{https://projecteuclid.org/euclid.aos/1351602527}
#'
#' Tom Goldstein, Brendan O'Donoghue, Simon Setzer, and Richard Baraniuk;
#' Fast Alternating Direction Optimization Methods.
#' SIAM Journal on Imaging Sciences, 2014, Vol. 7, No. 3 : pp. 1588-1623
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
#' ### Fit the estimator for some value of the tuning parameters
#' lambda <- 0.7; gamma <- 0.1 # The tuning parameters.
#' l1 <- lambda * gamma; l2 <- lambda * (1 - gamma)
#' fit <- lrpsadmm(Sigma = Sigma, Lambda1 = l1, Lambda2 = l2, n = dim(X)[1])
#' plot(fit) # Use the S3 method plot
#' estS <- fit$S # Sparse estimate
#' image(estS!=0) # Visualise its non-zero pattern
#' estL <- fit$L # Low-rank estimate
#' plot(eigen(estL)$values) # Visualise the spectrum of the low-rank estimate
#' plot(fit$lls[15:50], type='l') # The log-likelihood from iteration 15 onwards
#'
#' ### Fit for another value of the tuning parameters and compare cold/warm starts
#' lambda <- 0.4; gamma <- 0.1 # Change the tuning parameters
#' l1 <- lambda * gamma; l2 <- lambda * (1 - gamma)
#' # Reuse the previous fit as warm start:
#' warm.fit <- lrpsadmm(Sigma = Sigma, Lambda1 = l1, Lambda2 = l2, n = dim(X)[1], init=fit, tol = 1e-09)
#' warm.fit$iter # Number of itereations of the algorithm
#' # Fit without warm start
#' cold.fit <- lrpsadmm(Sigma = Sigma, Lambda1 = l1, Lambda2 = l2, n = dim(X)[1], tol=1e-09)
#' cold.fit$iter # Number of iterations
#' plot(cold.fit$lls[10:50], type='l', col='red', ylab = "Log-Likelihood", xlab = "#Iteration")
#' lines(warm.fit$lls[10:50], col='blue')
#' legend(20, 103, c("Cold Start", "Warm Start"), col=c("red", "blue"), lty=c(1,1))
#'
#' ### Force the sparsity pattern of the sparse component
#' zeros = 1 * (sim.data$precision.matrix != 0) # A mtrix of 0 and 1.
#' zeros = zeros[1:100,1:100] # Keep only the observed part
#' # Whereever zeros[i,j] = 0, the estimated S will be 0.
#' fit.zeros <- lrpsadmm(Sigma, l1, l2, n=dim(X)[1], tol=1e-09, zeros = zeros)
#' fit.no.zeros <- lrpsadmm(Sigma, l1, l2, n=dim(X)[1], tol=1e-09)
#' image(fit.zeros$S!=0) # Comparing the sparsity patterns
#' image(fit.no.zeros$S!=0)
#'
#' ### Fit the estimator when the problem is not so well-posed (n close to p)
#' set.seed(0)
#' # n = 80, p = 100 with 5 latent variables.
#' sim.data <- generate.latent.ggm.data(n=80, p=100, h=5, outlier.fraction = 0.0,
#'                                      sparsity = 0.02, sparsity.latent = 0.7)
#' X <- sim.data$obs.data; Sigma <- cor(X) # Sample correlation matrix
#'
#' lambda <- 2; gamma <- 0.1 # Here gamma is fairly small
#' l1 <- lambda * gamma; l2 <- lambda * (1 - gamma)
#' fit.small.gamma <- lrpsadmm(Sigma = Sigma, Lambda1 = l1, Lambda2 = l2, n = dim(X)[1])
#' plot(eigen(fit.small.gamma$L)$value) # Spectrum of L
#' lambda <- 0.28 ; gamma <- 0.7 # A large gamma, favourising a non low-rank component.
#' # This is too high for this ratio n/p
#' l1 <- lambda * gamma; l2 <- lambda * (1 - gamma)
#' fit.large.gamma <- lrpsadmm(Sigma = Sigma, Lambda1 = l1, Lambda2 = l2, n = dim(X)[1])
#' plot(eigen(fit.large.gamma$L)$value) # Spectrum of L
#' # Numerical stability and convergence are not guaranteed for such an ill-posed proble.
#' # Gamma is too high.
#' plot(fit.large.gamma$lls[350:500], type = 'l', xlab= "#Iterations", ylab = "Log-Likelihood")
#'
#' ### Fit the estimator with a robust estimator of the correlation matrix
#' # Generate data with 5% of outliers
#' set.seed(0)
#' sim.data <- generate.latent.ggm.data(n=2000, p=100, h=5, outlier.fraction = 0.05,
#'                                      sparsity = 0.02, sparsity.latent = 0.7)
#' X <- sim.data$obs.data;
#' Sigma <- cor(X) # Sample correlation matrix
#' Sigma.Kendall <- Kendall.correlation.estimator(X) # The robust estimator
#'
#' lambda <- 0.7; gamma <- 0.1 # The tuning parameters.
#' l1 <- lambda * gamma; l2 <- lambda * (1 - gamma)
#' fit <- lrpsadmm(Sigma = Sigma, Lambda1 = l1, Lambda2 = l2,
#' n = dim(X)[1], tol=1e-08, print_every = 200) # Outliers make the problem very ill-posed
#' Kendall.fit <- lrpsadmm(Sigma = Sigma.Kendall, Lambda1 = l1,
#' Lambda2 = l2, n = dim(X)[1], tol=1e-08) # Use the Kendall based estimator
#' image(fit$S!=0)
#' image(Kendall.fit$S!=0)
#' plot(fit$lls[800:1200], xlab="#Iterations", ylab="Log-Likelihood")
#' plot(Kendall.fit$lls, xlab="#Iterations", ylab="Log-Likelihood")
#'
#' @export
#' @seealso lrpsadmm.cv lrpsadmm.path
#' @import matrixcalc RSpectra MASS
lrpsadmm <- function(Sigma, Lambda1, Lambda2, n=NA, init=NULL,
                                     maxiter=2000, mu=0.1, tol=1e-05, eta=0.999,
                                     print_progress=T, print_every=20,
                                     zeros=NULL) {
  max.rank <- NULL
  if (is.null(zeros)) {
    zeros <- 1
  }
  p <- dim(Sigma)[1]
  if(is.na(n)) {
    n <- p
  }
  if(is.null(max.rank)) {
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
    parameters <- init
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
      parameters$termmsg <-
        "Algorithm was restarted 100 times without improvement."
      break()
    }

    if(any(diag(new_parameters$S) == 0)) {
      parameters$S <- diag(p)
      parameters$termcode <- -2
      parameters$termmsg <- "Shrinkage too strong: sparse component is empty."
      break()
    }
    # Compute the relative change in parameters with the previous iteration
    if ((matrixcalc::frobenius.norm(parameters$L) > 0)) {
      diff <- matrixcalc::frobenius.norm(new_parameters$S - parameters$S) /
        matrixcalc::frobenius.norm(parameters$S) +
        matrixcalc::frobenius.norm(new_parameters$L - parameters$L) /
        matrixcalc::frobenius.norm(parameters$L)
    } else {
      diff <- matrixcalc::frobenius.norm(new_parameters$S - parameters$S) /
        matrixcalc::frobenius.norm(parameters$S) +
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
    } else {
      ll <- NaN
    }
    if ((print_progress) & (i > print_every)) {
      if((i %% print_every) == 0) {
        last_lls <- lls[(i-print_every):i]
        avg_chg <- mean(abs(diff(last_lls)) / abs(last_lls[1:(print_every)]),
                        na.rm=T)
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
  parameters$exit <- NULL

  attr(parameters, "class") <- "lrpsadmm"

  parameters
}

#' @title Plotting function for 'lrpsadmm' Objects
#' @description Plots the sparsity pattern of S, the eigenvalues of L and
#' the values taken by the objective function at each iteration.
#' @param x An object of class lrpsadmm output by the function \code{lrpsadmm}
#' @examples
#' set.seed(1)
#' sim.data <- generate.latent.ggm.data(n=2000, p=100, h=5, outlier.fraction = 0.0,
#'                                     sparsity = 0.02, sparsity.latent = 0.7)
#' X <- sim.data$obs.data; Sigma <- cor(X) # Sample correlation matrix

#' lambda <- 0.7; gamma <- 0.1 # The tuning parameters.
#' l1 <- lambda * gamma; l2 <- lambda * (1 - gamma)
#' fit <- lrpsadmm(Sigma = Sigma, Lambda1 = l1, Lambda2 = l2, n = dim(X)[1])
#' plot(fit)
#' @importFrom graphics par plot title
#' @export
plot.lrpsadmm <- function(x) {
  fit <- x
  par(mfrow=c(2,2))
  image(fit$S!=0)
  title('Non-zero pattern of estimated sparse matrix S')
  plot(eigen(fit$L)$values, ylab = "EValues of L")
  title('Eigenvalues of estimated L')
  plot(fit$lls, xlab="#Iteration")
  title("Objective Function")
  par(mfrow=c(1,1))
}
