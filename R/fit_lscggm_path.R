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
#'
#' @import Matrix
#' @export
lscggmadmm.path <- function(SigmaZ,
                          SigmZX,
                          SigmaX,
                          gamma,
                          lambdas = NULL,
                          lambda.max = NULL,
                          lambda.ratio = 1e-4,
                          n.lambdas = 20,
                          max.sparsity = 0.5,
                          max.rank = NA,
                          rel_tol = 1e-02,
                          abs_tol = 1e-04,
                          max.iter = 2000,
                          mu = 1.0,
                          verbose = FALSE) {
  p <- dim(SigmaX)[1]
  m <- dim(SigmaZ)[1]
  
  if (is.null(lambdas)) {
    if (is.null(lambda.max)) {
      lambda.max1 <- max(abs(SigmaX - diag(diag(SigmaX)))) / gamma
      lambda.max2 <- max(abs(SigmaZX)) / gamma
      lambda.max <- max(lambda.max2, lambda.max2)
    }
    lambda.min <- lambda.max * lambda.ratio
    reason <- lambda.min / lambda.max
    lambdas <- lambda.max * reason ** ((0:n.lambdas) / n.lambdas)
  } else {
    lambdas <- sort(lambdas, decreasing = TRUE)
  }
  fit <- NULL
  path <- list()
  counter <- 1
  for (lambda in lambdas) {
    l1 <- gamma * lambda
    l2 <- (1 - gamma) * lambda
    fit <- lscggmadmm(
      SigmaZ, SigmaZX, SigmaX,
      l1,
      l2,
      init = fit,
      maxiter = max.iter,
      mu = mu,
      rel_tol = rel_tol,
      abs_tol = abs_tol,
      print_progress = FALSE
    )
    if (fit$termcode == -2) {
      next()
    }
    path[[counter]] <- list()
    path[[counter]]$lambda <- lambda
    path[[counter]]$gamma <- gamma
    path[[counter]]$lambda1 <- lambda * gamma
    path[[counter]]$lambda2 <- lambda * (1 - gamma)
    path[[counter]]$fit <- fit
    
    rank.LX <- as.numeric(Matrix::rankMatrix(fit$LX))
    SX <- (fit$SX != 0) - diag(diag(fit$SX != 0))
    rank.LZX <- as.numeric(Matrix::rankMatrix(fit$LZX))
    SZX <- (fit$SZX != 0)
    sparsity <- 0.5 * sum(SX) + sum(SZX)
    n.edges <- sparsity
    sparsity <- sparsity / (choose(p, 2) + p * m)
    
    path[[counter]]$number.of.edges <- n.edges
    path[[counter]]$sparsity <- sparsity
    path[[counter]]$rank.LX <- rank.LX
    path[[counter]]$rank.LZX <- rank.LZX
    counter <- counter + 1
    
    if (!is.na(max.rank)) {
      if (rank.LX > max.rank)
        break()
    }
    if (sparsity > max.sparsity)
      break()
    
    if (verbose) {
      print(
        paste(
          "Fitting with gamma=",
          gamma,
          " and lambda=",
          lambda,
          "Sparsity:",
          sparsity,
          "Rank of LX:",
          rank.LX,
          "Rank of LZX:",
          rank.LZX
        )
      )
    }
  }
  
  attr(path, "class") <- "lscggmadmmpath"
  path
}

#' @title Plotting for 'lscggmadmmpath' Objects
#' @description
#' Plots the sparsity and number of edges of S as a function of
#' the tuning parameter lambda (see documentation of lscggmadmmpath).
#' The rank of L as a function of lambda.
#'
#' The user can also provide a matrix of 0s and 1s to be used as
#' "ground truth". If the location of non-zero entries of S is known
#' (for example because this is simulated data), then the precision/recall
#' metrics are computed and plotted. See examples below.
#' @param x An object of class lrpsadmmpath output by the function \code{lrpsadmmpath}
#' @param ground.truth A binary matrix representing the adjacency matrix of the "true" graph
#' one seeks to recover. Useful mostly on simulated data where the true parameter is known.
#' @examples
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
#'                          lambda.ratio = 1e-03, n.lambdas = 30, verbose = TRUE)
#'
#' # This time let us ask for 30 values, but let us narrow down the range by using a
#' # a smaller ratio
#' my.path <- lrpsadmm.path(Sigma = Sigma, gamma = gamma,
#'                          lambda.max = 0.96, lambda.ratio = 0.1, n.lambdas = 30, verbose = TRUE)
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
#' @importFrom graphics par plot title
#' @seealso lrpsadmm lrpsadmm.cv
#' @export
plot.lrpsadmmpath <-
  function(x, ground.truth = NULL) {
    lrps.path <- x
    gamma <- lrps.path[[1]]$gamma
    lambdas <- ranks <- sparsities <- edges <- c()
    if (!is.null(ground.truth)) {
      ground.truth <- (ground.truth != 0) - diag(diag(ground.truth != 0))
    }
    prs <- rs <- c()
    for (i in 1:length(lrps.path)) {
      lambdas <- c(lambdas, lrps.path[[i]]$lambda)
      ranks <- c(ranks, lrps.path[[i]]$rank.L)
      edges <- c(edges, lrps.path[[i]]$number.of.edges)
      sparsities <- c(sparsities, lrps.path[[i]]$sparsity)
      if (!is.null(ground.truth)) {
        S <- lrps.path[[i]]$fit$S
        S <- (S != 0) - diag(diag(S != 0))
        if (sum(S) > 0) {
          pr <- sum(S * ground.truth) / sum(S)
          r <- sum(S * ground.truth) / sum(ground.truth)
          prs <- c(prs, pr)
          rs <- c(rs, r)
        }
      }
    }
    
    par(mfrow = c(2, 2))
    plot(-log(lambdas), sparsities, xlab = "-Log10(Lambda)", ylab = "Sparsity of S")
    title("Sparsity of the estimated S\n as a function Lambda")
    plot(-log(lambdas), edges, xlab = "-Log10(Lambda)", ylab = "Number of Edges")
    title("Number of edges in the estimated\n graphical model as a function of Lambda")
    plot(-log(lambdas), ranks, xlab = "-Log10(Lambda)", ylab = "Rank of L")
    title("Rank of estimated low-rank\n component L as a function of Lambda")
    if (!is.null(ground.truth)) {
      plot(rs, prs, xlab = "Recall", ylab = "Precision")
      title("Precision / Recall curve for\n the recovery of the support of S.")
    }
    par(mfrow = c(1, 1))
  }
