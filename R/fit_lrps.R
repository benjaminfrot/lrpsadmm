

.obj_func.R <- function(Sigma, A, S, L, l1, l2) {
  evals <- eigen(A, symmetric = TRUE, only.values = TRUE)$values
  evals[abs(evals) < 1e-09] <- 1e-8
  if (any(evals < 0)) {
    return(NaN)
  }
  
  # Log-Likelihood
  ll <- sum(diag(Sigma %*% A)) - sum(log(evals))
  # + Penalty
  ll <- ll + l1 * sum(abs(S)) + l2 * sum(diag(L))
  
  ll
}

.updateA.R <- function(Sigma, S, L, U, mu) {
  X1 <- mu * (S - L) - Sigma - U
  X2 <- X1 %*% X1 + 4 * mu * diag(dim(X1)[1])
  eig <- eigen(X2, symmetric = TRUE)
  sqrtX2 <-
    eig$vectors %*% diag(sqrt(eig$values)) %*% t(eig$vectors)
  A <- (X1 + sqrtX2) / (2 * mu)
  
  return(A)
}

.updateL.R <- function(A, S, U, l2, mu) {
  lp2 <- l2 / mu
  
  X1 <- S - A - (U / mu)
  eig <- eigen(X1)
  eigVal <- eig$values - lp2
  eigVal[eigVal < 0] <- 0
  L <- eig$vectors %*% diag(eigVal) %*% t(eig$vectors)
  
  return(L)
}

.updateS.R <- function(A, L, U, l1, mu, zeros) {
  lp1 <- l1 / mu
  
  X1 <- A + L + (U / mu)
  X2 <- abs(X1) - lp1
  X2[X2 <= 0] <- 0
  S <- X2 * sign(X1)
  S <- S * zeros
  
  return(S)
}

.updateU.R <- function(A, S, L, U, mu) {
  U <- U + mu * (A - S + L)
  
  return(U)
}


#' @import MASS
lrpsadmm.R <- function(Sigma,
                       Lambda1,
                       Lambda2,
                       init,
                       maxiter,
                       mu,
                       abs_tol,
                       rel_tol,
                       print_progress,
                       print_every,
                       zeros) {
  if (is.null(zeros)) {
    zeros <- 1
  }
  
  p <- dim(Sigma)[1]
  
  if (is.null(init)) {
    S <- diag(p)
    L <- S * 0.0
    U <- S * 0.0
  } else {
    parameters <- init
    S <- init$S
    L <- init$L
    U <- init$U
  }
  
  # Enforce the 0 pattern
  S <- S * zeros
  
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
  for (i in 1:maxiter) {
    # Update A
    A <- .updateA.R(Sigma, S, L, U, mu)
    
    # Update S
    S_old <- S
    S <- .updateS.R(A, L, U, Lambda1, mu, zeros)
    if (any(diag(S) == 0)) {
      S <- diag(p)
      L <- S * 0
      U <- S * 0
      parameters$termcode <- -2
      parameters$termmsg <-
        "Shrinkage too strong: sparse component is empty."
      break()
    }
    
    # Update L
    L_old <- L
    L <- .updateL.R(A, S, U, Lambda2, mu)
    
    # Update U
    U <- .updateU.R(A, S, L, U, mu)
    
    # Diagnostics
    objval <- .obj_func.R(Sigma, A, S, L, Lambda1, Lambda2)
    
    r_norm <- norm(A - (S - L), 'F')
    s_norm <- norm(mu * ((S - L) - (S_old - L_old)), 'F')
    eps_pri <- p * abs_tol + rel_tol *
      max(norm(A, 'F'), norm(S - L, 'F'))
    eps_dual <- p * abs_tol + rel_tol * norm(mu * U, 'F')
    history <- rbind(history, c(i, objval, s_norm,
                                r_norm, eps_pri, eps_dual))
    
    if ((s_norm < eps_dual) && (r_norm < eps_pri)) {
      parameters$termcode <- 0
      parameters$termmsg <- 'Convergence Reached.'
      break()
    }
    if ((print_progress) & (i >= print_every)) {
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
  
  parameters$S <- S
  parameters$L <- L
  parameters$U <- U
  parameters$history <- as.data.frame(history)
  
  attr(parameters, "class") <- "lrpsadmm"
  
  parameters
}


#'
#' Fit the low-rank plus sparse estimator using the alternating method direction of multipliers
#' @description
#' Given a n x p data matrix X and its empirical correlation matrix
#' \eqn{\Sigma}, an alternating direction method of multipliers (ADMM) algorithm
#' is used to obtain the solutions \eqn{S, L} to
#' \deqn{(S, L) = argmin_{A, B}  -logdet(A - B) + Tr((A-B)\Sigma) + \lambda_1 ||A||_1 + \lambda_2Tr(B),}
#' subject to \eqn{A - B > 0} and \eqn{B \ge 0}.
#' @param Sigma An estimate of the correlation matrix.
#' @param Lambda1 Penalty on the l1 norm of S
#' @param Lambda2 Penalty on the sum of the eigenvalues of L
#' @param init The output of a previous run of the algorithm. For warm starts.
#' @param maxiter Maximal number of iterations
#' @param mu Stepsize of the ADMM algorithm.
#' @param rel_tol Relative tolerance required to stop the algorithm. The algorithm
#' stops when both the change in parameters is below tolerance and the constraints
#' are satisfied. Default 1e-02.
#' @param abs_tol Absolute tolerance required to stop the algorithm. Default 1e-04.
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
#'  0: Convergence reached. -1: Maxiter reached. -2: Shrinkage too strong.}
#'  \item{termmsg}{A character vector. The message corresponding to the value \code{termcode}.}
#'  \item{history}{A numerical dataframe with the objective function at each
#'  iterations, the norm and dual norm as well the primal and dual tolerance
#'  criteria for convergence. The algorithm exits when r_norm < eps_pri and s_norm < eps_dual.}
#'  \item{U}{A p x p matrix. Augmented Lagrangian multiplier. It is stored in order to allow warm starts.}
#' }
#'
#' @references
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
#' ### Fit the estimator for some value of the tuning parameters
#' lambda <- 0.7; gamma <- 0.1 # The tuning parameters.
#' l1 <- lambda * gamma; l2 <- lambda * (1 - gamma)
#' fit <- lrpsadmm(Sigma = Sigma, Lambda1 = l1, Lambda2 = l2,
#'                 abs_tol=1e-06, rel_tol=1e-04)
#' plot(fit) # Use the S3 method plot
#' estS <- fit$S # Sparse estimate
#' image(estS!=0) # Visualise its non-zero pattern
#' estL <- fit$L # Low-rank estimate
#' plot(eigen(estL)$values) # Visualise the spectrum of the low-rank estimate
#' plot(fit$history$Objval, type='l') # The log-likelihood from iteration 15 onwards
#'
#' ### Fit for another value of the tuning parameters and compare cold/warm starts
#' lambda <- 0.4; gamma <- 0.1 # Change the tuning parameters
#' l1 <- lambda * gamma; l2 <- lambda * (1 - gamma)
#' # Reuse the previous fit as warm start:
#' warm.fit <- lrpsadmm(Sigma = Sigma, Lambda1 = l1, Lambda2 = l2,
#'                      init=fit, rel_tol = 1e-04, abs_tol=1e-06)
#' # Fit without warm start
#' cold.fit <- lrpsadmm(Sigma = Sigma, Lambda1 = l1, Lambda2 = l2, rel_tol=1e-04,
#'                      abs_tol=1e-06)
#' plot(cold.fit$history$Objval, type='l', col='red', ylab = "Log-Likelihood", xlab = "#Iteration")
#' lines(warm.fit$history$Objval, col='blue')
#' xleg = 0.5 * nrow(cold.fit$history)
#' yleg = 0.5 * (max(cold.fit$history$Objval) + min(cold.fit$history$Objval))
#' legend(x = xleg, y=yleg, legend=c("Cold Start", "Warm Start"), col=c("red", "blue"), lty=c(1,1))
#'
#' ### Force the sparsity pattern of the sparse component
#' zeros = 1 * (sim.data$precision.matrix != 0) # A mtrix of 0 and 1.
#' zeros = zeros[1:100,1:100] # Keep only the observed part
#' # Whereever zeros[i,j] = 0, the estimated S will be 0.
#' fit.zeros <- lrpsadmm(Sigma, l1, l2, abs_tol=1e-06, rel_tol=1e-04, zeros = zeros)
#' fit.no.zeros <- lrpsadmm(Sigma, l1, l2, abs_tol=1e-06, rel_tol=1e-04)
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
#' fit.small.gamma <- lrpsadmm(Sigma = Sigma, Lambda1 = l1, Lambda2 = l2)
#' plot(eigen(fit.small.gamma$L)$value) # Spectrum of L
# lambda <- 0.28 ; gamma <- 0.7 # A large gamma, favourising a non low-rank component.
#' # This is too high for this ratio n/p
#' l1 <- lambda * gamma; l2 <- lambda * (1 - gamma)
#' fit.large.gamma <- lrpsadmm(Sigma = Sigma, Lambda1 = l1, Lambda2 = l2)
#' plot(eigen(fit.large.gamma$L)$value) # Spectrum of L
#' # Numerical stability and convergence are not guaranteed for such an ill-posed proble.
#' # Gamma is too high.
#' plot(fit.large.gamma$history$Objval, type = 'l', xlab= "#Iterations", ylab = "Log-Likelihood")
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
#' # Outliers make the problem very ill-posed
#' fit <- lrpsadmm(Sigma = Sigma, Lambda1 = l1, Lambda2 = l2,
#'                 abs_tol=1e-06, rel_tol=1e-04, print_every = 200)
#' # Use the Kendall based estimator
#' Kendall.fit <- lrpsadmm(Sigma = Sigma.Kendall, Lambda1 = l1,
#'                         Lambda2 = l2, abs_tol=1e-06, rel_tol=1e-04)
#' image(fit$S!=0)
#' image(Kendall.fit$S!=0)
#' plot(fit$history$Objval, xlab="#Iterations", ylab="Log-Likelihood")
#' plot(Kendall.fit$history$Objval, xlab="#Iterations", ylab="Log-Likelihood")
#' @export
#' @seealso lrpsadmm.cv lrpsadmm.path
#' @import MASS
lrpsadmm <- function(Sigma,
                     Lambda1,
                     Lambda2,
                     init = NULL,
                     maxiter = 1000,
                     mu = 1.0,
                     abs_tol = 1e-04,
                     rel_tol = 1e-02,
                     print_progress = TRUE,
                     print_every = 10,
                     zeros = NULL) {
  res <- lrpsadmm.R(
    Sigma,
    Lambda1,
    Lambda2,
    init,
    maxiter,
    mu,
    abs_tol,
    rel_tol,
    print_progress,
    print_every,
    zeros
  )
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
#' fit <- lrpsadmm(Sigma = Sigma, Lambda1 = l1, Lambda2 = l2)
#' plot(fit)
#' @importFrom graphics par plot title
#' @export
plot.lrpsadmm <- function(x) {
  fit <- x
  par(mfrow = c(2, 2))
  image(fit$S != 0)
  title('Non-zero pattern of estimated sparse matrix S')
  plot(eigen(fit$L)$values, ylab = "EValues of L")
  title('Eigenvalues of estimated L')
  plot(fit$history$Objval, xlab = "#Iteration")
  title("Objective Function")
  par(mfrow = c(1, 1))
}