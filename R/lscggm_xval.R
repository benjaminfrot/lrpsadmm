#' @title Perform K-fold cross-validation for the Conditional Low-Rank plus Sparse estimator
#' @description
#'   Performs K-fold cross-validation in order to select the tuning parameters
#'   lambda and gamma.
#'
#'   Recall that the penalty for the conditional low-rank plus sparse estimator is written as \eqn{\lambda_1 ||S||_1 + \lambda_2 Trace(L)} in the
#'   objective function of \code{lscggmadmm}.
#'   This can be equivalently rewritten in terms of the regularisation parameters \eqn{\lambda} and \eqn{\gamma} as follows
#'   \deqn{\lambda \gamma ||S||_1 + \lambda (1 - \gamma) Trace(L),}
#'   for \eqn{\gamma \in (0, 1)}.
#'
#'   For a given value of \eqn{\gamma}, one can perform cross-validation along the regularisation path in order to choose
#'   \eqn{\lambda}. This function computes the regularisation paths for each value of \eqn{\gamma} supplied as arguments
#'   and performs cross-validation. The pair (\eqn{\lambda, \gamma}) that produces the smallest cross-validated log-likelihood
#'   is returned.
#' @details
#'   Recall that the penalty for the conditional estimator is written as \eqn{\lambda_1 ||S||_1 + \lambda_2 Trace(L)} in the
#'   objective function of \code{lscggmadmm}.
#'   This can be equivalently rewritten in terms of the regularisation parameters \eqn{\lambda} and \eqn{\gamma} as follows
#'   \deqn{\lambda \gamma ||S||_1 + \lambda (1 - \gamma) Trace(L),}
#'   for \eqn{\gamma \in (0, 1)}.
#'
#'   For a given value of \eqn{\gamma}, one can perform cross-validation along the regularisation path in order to choose
#'   \eqn{\lambda}. This function computes the regularisation paths for each value of \eqn{\gamma} supplied as arguments
#'   and performs cross-validation. One can then select the pair (\eqn{\lambda, \gamma}) that produces the smallest
#'   cross-validated log-likelihood.
#'
#' The function \code{lscggmadmm} is fitted for successive values of \eqn{\lambda} using warm starts.
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
#' Recall that cross-validation overselects. It might be good for prediction purposes but methods such
#' statibility selection are probably better if the support of S is what is of interest.
#' @param X n x p data matrix
#' @param Z n x m data matrix of variables to condition on.
#' @param gammas A real or a vector of reals between 0 and 1 (non inclusive). For each value of gamma
#' the regularisation path is computed and cross-validation is performed. See examples for
#' guidance on how to choose reasonable values for gamma. Too high a value might result in the
#' problem being under-identified and therefore numerical instabilities.
#' @param n.folds Number of folds for cross-validation. Default 5.
#' @param covariance.estimator A function that takes a data matrix and outputs an estimate of its correlation matrix.
#' Default: the sample correlation matrix output by the \code{cor} function.
#' @param lambdas A decreasing sequence of values of lambda. See Details for the default values.
#' @param lambda.max A positive real. Maximum value of lambda. See Details.
#' @param lambda.ratio A real between 0 and 1. The smallest value of lambda is given by lambda.max * lambda.ratio. See Details.
#' @param n.lambdas A positive integer. The number of values of lambda to
#' generate according a geometric sequence between lambda.max and lambda.max * lambda.ratio. See Details.
#' @param max.sparsity A real between 0 and 1. Abort the computation of the path if S becomes denser than this value.
#' @param max.rank A real between 0 and 1. Abort the computation of the path if the rank of L becomes higher than this value.
#' @param rel_tol \code{rel_tol} parameter of the \code{lscggmadmm} function.
#' @param abs_tol \code{abs_tol} parameter of the \code{lscggmadmm} function.
#' @param max.iter \code{max.iter} parameter of the \code{lscggmadmm} function.
#' @param mu \code{mu} parameter of the \code{lscggmadmm} function.
#' @param verbose A boolean. Whether to print the value of lambda, gamma, sparsity of S, etc... after each fit
#' @param seed Set the seed of the random number generator used for the K folds.
#' 
#' @return
#'   An object of class lscggmadmmcv.
#'   It contains the values of the mean cross-validated log-likelihood, its standard deviation for each
#'   pair (lambda, gamma) and an object of class lscggmadmm.path for each value of gamma.
#'   See the examples for how to access the selected tuning parameters, best fit etc...
#' @import cvTools
#' @seealso lscggm lscggm.path
#' @export
lscggm.cv <- function(X, Z,
                        gammas = c(0.05, 0.1, 0.15),
                        covariance.estimator = cor,
                        n.folds = 5,
                        lambdas = NULL,
                        lambda.max = NULL,
                        lambda.ratio = 1e-4,
                        n.lambdas = 20,
                        max.sparsity = 0.5,
                        max.rank = NA,
                        abs_tol = 1e-05,
                        rel_tol = 1e-03,
                        max.iter = 2000,
                        mu = 1.0,
                        verbose = FALSE,
                        seed = NA) {
  n <- dim(X)[1]
  if (!is.na(seed)) {
    set.seed(seed)
  }
  folds <- cvTools::cvFolds(n, K = n.folds)
  
  all.paths <- list()
  counter <- 1
  best.ll <- Inf
  for (gamma in gammas) {
    all.paths[[counter]] <- list()
    all.paths[[counter]]$cross.validated.path <- .cv.lscggm.one.gamma(
      X,Z,
      gamma,
      folds,
      covariance.estimator,
      lambdas,
      lambda.max,
      lambda.ratio,
      n.lambdas,
      max.sparsity,
      max.rank,
      rel_tol,
      abs_tol,
      max.iter,
      mu,
      verbose,
      seed
    )
    all.paths[[counter]]$gamma <- gamma
    all.paths[[counter]]$best.fit <-
      .choose.cross.validate.lscggm(all.paths[[counter]]$cross.validated.path)
    mxll <- all.paths[[counter]]$best.fit$mean_xval_ll
    if (mxll < best.ll) {
      best.ll <- mxll
      best.fit <- all.paths[[counter]]$best.fit
    }
    counter <- counter + 1
  }
  
  toReturn <- list()
  toReturn$best.fit <- best.fit
  toReturn$best.gamma <- best.fit$gamma
  toReturn$best.lambda <- best.fit$lambda
  toReturn$cross.validated.paths <- all.paths
  
  attr(toReturn, "class") <- "lscggmadmmcv"
  toReturn
}

#' @import Matrix
#' @importFrom stats cor sd
.cv.lscggm.one.gamma <- function(X, Z,
                               gamma,
                               folds,
                               covariance.estimator = cor,
                               lambdas,
                               lambda.max,
                               lambda.ratio,
                               n.lambdas,
                               max.sparsity,
                               max.rank,
                               rel_tol,
                               abs_tol,
                               max.iter,
                               mu,
                               verbose,
                               seed) {
  SigmaX <- covariance.estimator(X)
  SigmaZX <- covariance.estimator(Z, X)
  SigmaZ <- covariance.estimator(Z)
  
  p <- dim(SigmaX)[1]
  m <- dim(SigmaZ)[1]
  n <- dim(X)[1]
  n.folds <- folds$K
  
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
  
  # Start by computing the whole path on the full dataset.
  if (verbose) {
    print ("### Computing the path on the full dataset first ###")
  }
  path <-
    lscggmadmm.path(
      SigmaZ, SigmaZX, SigmaZ,
      gamma,
      lambdas = lambdas,
      max.sparsity = max.sparsity,
      max.rank = max.rank,
      rel_tol = rel_tol,
      abs_tol = abs_tol,
      max.iter = max.iter,
      mu = mu,
      verbose = verbose
    )
  if (verbose) {
    print(paste("### Now performing ", n.folds, " fold cross validation. ###"))
  }
  
  valid.lambdas <- c()
  for (i in 1:length(path)) {
    valid.lambdas <- c(valid.lambdas, path[[i]]$lambda)
  }
  X <- X[folds$subsets[, 1], ]
  Z <- Z[folds$subsets[, 1], ]
  
  log.liks <- matrix(NA, ncol = 2)
  for (lambda in valid.lambdas) {
    l1 <- gamma * lambda
    l2 <- (1 - gamma) * lambda
    fit <- NULL
    for (i in 1:length(path)) {
      if (path[[i]]$lambda == lambda) {
        fit <- path[[i]]$fit
        index <- i
        break()
      }
    }
    lls <- c()
    for (k in 1:n.folds) {
      Strain <- covariance.estimator(X[folds$which != k,])
      Stest <- covariance.estimator(X[folds$which == k,])
      SZtrain <- covariance.estimator(Z[folds$which != k,])
      SZtest <- covariance.estimator(Z[folds$which == k,])
      SZXtrain <- covariance.estimator(Z[folds$which != k,], X[folds$which != k,])
      SZXtest <- covariance.estimator(Z[folds$which != k,], X[folds$which == k,])
      fitll <- lscggmadmm(
        Strain, SZXtrain, SZtrain,
        l1,
        l2,
        init = fit,
        print_progress = F,
        rel_tol = rel_tol,
        abs_tol = abs_tol,
        maxiter = max.iter,
        mu = mu
      )
      if (fitll$termcode == -2) {
        ll <- NaN
        lls <- c(lls, ll)
        break()
      }
      
      # Compute the log likelihood on the testing set
      AX <- fitll$SX - fitll$LX
      AZX <- fitll$SZX - fitll$LZX
      eigen.out <- tryCatch({
        # In some cases this does not converge.
        eigen(AX, symmetric = TRUE)$values
      },
      error = function(e) {
        print(e)
        c(-10)
      })
      evecs <- eigen.out@vectors
      evals <- eigen.out$values
      evals[abs(evals) < 1e-05] <- 1e-16
      if (any(evals < 0)) {
        ll <- NaN
        lls <- c(lls, ll)
        break()
      }
      AXi <- evecs %*% (diag(1 / evals)) %*% t(evecs)
      ll <- -sum(log(evals))
      ll <-
        ll + sum(diag(Stest %*% AX)) + 2 * sum(diag(SZXtest %*% t(AZX)))
      ll <- ll +  sum(diag(AXi %*% t(AZX) %*% SZtest %*% AZX))
      lls <- c(lls, ll)
    }
    if (any(is.nan(lls))) {
      path[[index]] <- NULL
      next()
    }
    mean.ll <- mean(lls)
    sd.ll <- sd(lls)
    path[[index]]$mean_xval_ll <- mean.ll
    path[[index]]$sd_xval_ll <- sd.ll
    if (verbose) {
      print(
        paste(
          "Lambda:",
          lambda,
          "X-Val Log-lik:",
          mean.ll,
          "#Edges:",
          path[[index]]$number.of.edges
        )
      )
    }
  }
  
  path
}

.choose.cross.validate.lscggm <-
  function(xval.path, method = "min") {
    min.ll <- Inf
    min.sd <- NULL
    min.index <- NULL
    min.lambda <- NULL
    for (i in 1:length(xval.path)) {
      mean.ll <- xval.path[[i]]$mean_xval_ll
      if (mean.ll < min.ll) {
        min.ll <- mean.ll
        min.sd <- xval.path[[i]]$sd_xval_ll
        min.index <- i
        min.lambda <- xval.path[[i]]$lambda
      }
    }
    if (method == "min") {
      return(xval.path[[min.index]])
    }
    
    if (method == "hastie") {
      for (i in 1:length(xval.path)) {
        if (xval.path[[i]]$lambda < min.lambda) {
          next()
        }
        if (xval.path[[i]]$mean_xval_ll < (min.ll + min.sd)) {
          min.ll <- xval.path[[i]]$mean_xval_ll
          min.sd <- xval.path[[i]]$sd_xval_ll
          min.index <- i
          min.lambda <- xval.path[[i]]$lambda
        }
      }
      
      return(xval.path[[min.index]])
      
    } else {
      stop("Method has to be 'min' or 'hastie'.")
    }
  }