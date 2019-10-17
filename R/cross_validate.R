#' @title Perform K-fold cross-validation for the Low-Rank plus Sparse estimator
#' @description
#'   Performs K-fold cross-validation in order to select the tuning parameters
#'   lambda and gamma.
#'
#'   Recall that the penalty for the LRpS estimator is written as \eqn{\lambda_1 ||S||_1 + \lambda_2 Trace(L)} in the
#'   objective function of \code{lrpsadmm}.
#'   This can be equivalently rewritten in terms of the regularisation parameters \eqn{\lambda} and \eqn{\gamma} as follows
#'   \deqn{\lambda \gamma ||S||_1 + \lambda (1 - \gamma) Trace(L),}
#'   for \eqn{\gamma \in (0, 1)}.
#'
#'   For a given value of \eqn{\gamma}, one can perform cross-validation along the regularisation path in order to choose
#'   \eqn{\lambda}. This function computes the regularisation paths for each value of \eqn{\gamma} supplied as arguments
#'   and performs cross-validation. The pair (\eqn{\lambda, \gamma}) that produces the smallest cross-validated log-likelihood
#'   is returned.
#' @details
#'   Recall that the penalty for the LRpS estimator is written as \eqn{\lambda_1 ||S||_1 + \lambda_2 Trace(L)} in the
#'   objective function of \code{lrpsadmm}.
#'   This can be equivalently rewritten in terms of the regularisation parameters \eqn{\lambda} and \eqn{\gamma} as follows
#'   \deqn{\lambda \gamma ||S||_1 + \lambda (1 - \gamma) Trace(L),}
#'   for \eqn{\gamma \in (0, 1)}.
#'
#'   For a given value of \eqn{\gamma}, one can perform cross-validation along the regularisation path in order to choose
#'   \eqn{\lambda}. This function computes the regularisation paths for each value of \eqn{\gamma} supplied as arguments
#'   and performs cross-validation. One can then select the pair (\eqn{\lambda, \gamma}) that produces the smallest
#'   cross-validated log-likelihood.
#'
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
#' Recall that cross-validation overselects. It might be good for prediction purposes but methods such
#' statibility selection are probably better if the support of S is what is of interest.
#' @param X n x p data matrix
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
#' @param max.rank A real between 0 and 1. Abort the computuation of the path if the rank of L becomes higher than this value.
#' @param rel_tol \code{rel_tol} parameter of the \code{lrpsadmm} function.
#' @param abs_tol \code{abs_tol} parameter of the \code{lrpsadmm} function.
#' @param max.iter \code{max.iter} parameter of the \code{lrpsadmm} function.
#' @param mu \code{mu} parameter of the \code{lrpsadmm} function.
#' @param verbose A boolean. Whether to print the value of lambda, gamma, sparsity of S, etc... after each fit
#' @param seed Set the seed of the random number generator used for the K folds.
#' @param zeros A p x p matrix with entries set to 0 or 1. Whereever its entries are
#' 0, the entries of the estimated S will be forced to 0.
#' @param backend The \code{backend} parameter of lrpsadmm. It is one of 'R' or 'RcppEigen'. 
#' 
#' @return
#'   An object of class lrpsadmmcv.
#'   It contains the values of the mean cross-validated log-likelihood, its standard deviation for each
#'   pair (lambda, gamma) and an object of class lrpsadmm.path for each value of gamma.
#'   See the examples for how to access the selected tuning paramters, best fit etc...
#' @examples
#' set.seed(0)
#' ## 1 - A simple example: Well-powered dataset
#' sim.data <- generate.latent.ggm.data(n=2000, p=100, h=5, outlier.fraction = 0.0,
#'                                     sparsity = 0.02, sparsity.latent = 0.7)
#' ground.truth <- sim.data$precision.matrix[1:100, 1:100]
# Remove the elements along the diagonal. Keep a matrix of 0s and 1s
#' ground.truth <- 1 * (( ground.truth - diag(diag(ground.truth)) ) !=0)
#' X <- sim.data$obs.data;
#' gammas <- c(0.1, 0.15, 0.2)
#' # Let the function decide the range of value of Lambda
#' cvpath <- lrpsadmm.cv(X, gammas = gammas,
#'                          lambda.ratio = 1e-02, n.lambdas = 30, verbose = TRUE)
#' best.gamma <- cvpath$best.gamma
#' best.lambda <- cvpath$best.lambda
#' best.fit <- cvpath$best.fit$fit # Object of class lrpsadmm
#' plot(best.fit)
#' # The value Gamma = 0.15 is selected by X-validation
#' plot(cvpath) # Plot the outcome. Object of class lrpsadmmcv
#' # We can look at the path corresponding to this value
#' best.path <- cvpath$cross.validated.paths[[which(gammas == best.gamma)]]$cross.validated.path
#' # We know the ground truth, so let's use it to see how well we did:
#' plot(best.path, ground.truth = ground.truth)
#'
#' ## 2 - Data with outliers: use a robust estimator
#' set.seed(0)
#' sim.data <- generate.latent.ggm.data(n=2000, p=50, h=5, outlier.fraction = 0.05,
#'                                     sparsity = 0.02, sparsity.latent = 0.7)
#' ground.truth <- sim.data$precision.matrix[1:50, 1:50]
#' # Remove the elements along the diagonal. Keep a matrix of 0s and 1s
#' ground.truth <- 1 * (( ground.truth - diag(diag(ground.truth)) ) !=0)
#' X <- sim.data$obs.data;
#' # We can afford high values of gamma because n >> p
#' gammas <- c(0.2, 0.3, 0.4)
#' # Use the Kendall based estimator:
#' cvpath <- lrpsadmm.cv(X, gammas = gammas, covariance.estimator = Kendall.correlation.estimator,
#'                          lambda.ratio = 1e-03, n.lambdas = 30, verbose = TRUE)
#' plot(cvpath)
#' best.gamma <- cvpath$best.gamma
#' best.lambda <- cvpath$best.lambda
#' best.path <- cvpath$cross.validated.paths[[which(gammas == best.gamma)]]$cross.validated.path
#' plot(best.path, ground.truth = ground.truth)
#'
#' ## 3 - A tougher problem n is close to p
#' set.seed(0)
#' sim.data <- generate.latent.ggm.data(n=150, p=100, h=5, outlier.fraction = 0.0,
#'                                      sparsity = 0.02, sparsity.latent = 0.7)
#' ground.truth <- sim.data$precision.matrix[1:100, 1:100]
#' # Remove the elements along the diagonal. Keep a matrix of 0s and 1s
#' ground.truth <- 1 * (( ground.truth - diag(diag(ground.truth)) ) !=0)
#' X <- sim.data$obs.data;
#' # Since n < p, do not try too high values of gamma. Stay close to 0.
#' gammas <- c(0.07, 0.1, 0.12)
#' cvpath <- lrpsadmm.cv(X, gammas = gammas,
#'                          lambda.ratio = 0.1, n.lambdas = 20, verbose = TRUE)
#' plot(cvpath) # Plot the outcome
#'
#' # Clearly the range selected by the function is not good enough
#' # We need better values for lambda. In that case it is better
#' # to see which value of lambda yields a very sparse graph
#' # for a given gamma:
#' gammas <- c(0.1)
#' # We set the seed so we can compre the x-validated log-likelihoods between
#' # two runs of the function
#' cvpath <- lrpsadmm.cv(X, gammas = gammas, lambda.max = 2.2,
#'                         lambda.ratio = 0.1, n.lambdas = 20, verbose = TRUE, seed = 0)
#' plot(cvpath)
#' gammas <- c(0.12)
#' cvpath <- lrpsadmm.cv(X, gammas = gammas, lambda.max = 1.7,
#'                          lambda.ratio = 0.1, n.lambdas = 20, verbose = TRUE, seed = 0)
#' plot(cvpath)
#'
#' @import cvTools
#' @seealso lrpsadmm lrpadmm.path
#' @export
lrpsadmm.cv <- function(X,
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
                        seed = NA,
                        zeros = NULL,
                        backend='RcppEigen') {
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
    all.paths[[counter]]$cross.validated.path <- .cv.lrps.one.gamma(
      X,
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
      seed,
      zeros,
      backend
    )
    all.paths[[counter]]$gamma <- gamma
    all.paths[[counter]]$best.fit <-
      .choose.cross.validate.low.rank.plus.sparse(all.paths[[counter]]$cross.validated.path)
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

  attr(toReturn, "class") <- "lrpsadmmcv"
  toReturn
}

#' @import RSpectra Matrix
#' @importFrom stats cor sd
.cv.lrps.one.gamma <- function(X,
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
                               seed,
                               zeros,
                               backend) {
  Sigma <- covariance.estimator(X)
  p <- dim(Sigma)[1]
  n <- dim(X)[1]
  n.folds <- folds$K

  if (is.null(lambdas)) {
    if (is.null(lambda.max)) {
      max.cor <- max(abs(Sigma - diag(diag(Sigma)))) * 2
      lambda.max <- max.cor / gamma
    }
    lambda.min <- lambda.max * lambda.ratio
    reason <- lambda.min / lambda.max
    lambdas <- lambda.max * reason ** ((0:n.lambdas) / n.lambdas)
  }

  # Start by computing the whole path on the full dataset.
  if (verbose) {
    print ("### Computing the path on the full dataset first ###")
  }
  path <-
    lrpsadmm.path(
      Sigma,
      gamma,
      lambdas = lambdas,
      max.sparsity = max.sparsity,
      max.rank = max.rank,
      rel_tol = rel_tol,
      abs_tol = abs_tol,
      max.iter = max.iter,
      mu = mu,
      verbose = verbose,
      zeros=zeros,
      backend=backend
    )
  if (verbose) {
    print(paste("### Now performing ", n.folds, " fold cross validation. ###"))
  }

  valid.lambdas <- c()
  for (i in 1:length(path)) {
    valid.lambdas <- c(valid.lambdas, path[[i]]$lambda)
  }
  X <- X[folds$subsets[, 1], ]

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
      fitll <- lrpsadmm(
        Strain,
        l1,
        l2,
        init = fit,
        print_progress = F,
        rel_tol = rel_tol,
        abs_tol = abs_tol,
        maxiter = max.iter,
        mu = mu,
        zeros=zeros,
        backend=backend
      )
      if (fitll$termcode == -2) {
        ll <- NaN
        lls <- c(lls, ll)
        break()
      }

      # Compute the log likelihood on the testing set
      A <- fitll$S - fitll$L
      evals <- tryCatch({
        # In some cases this does not converge.
        RSpectra::eigs_sym(A, min(n - 1, p - 1))$values
      },
      error = function(e) {
        print(e)
        c(-10)
      })
      evals[abs(evals) < 1e-05] <- 1e-16
      if (any(evals < 0)) {
        ll <- NaN
        lls <- c(lls, ll)
        break()
      }

      ll <- sum(diag(Stest %*% A)) - sum(log(evals))
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

.choose.cross.validate.low.rank.plus.sparse <-
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

.prepare.plot.one.cv.path <- function(cvpath) {
  lambdas <- c()
  sds <- c()
  lls <- c()
  bl <- cvpath$best.fit$lambda
  for (i in 1:length(cvpath$cross.validated.path)) {
    lambdas <- c(lambdas, -log10(cvpath$cross.validated.path[[i]]$lambda))
    sds <- c(sds, (cvpath$cross.validated.path[[i]]$sd_xval_ll))
    lls <- c(lls, (cvpath$cross.validated.path[[i]]$mean_xval_ll))
  }
  ymin <- min(lls - sds)
  ymin <- ymin - 2 * mean(sds)
  ymax <- max(lls + sds)
  ymax <- ymax + 2 * mean(sds)

  a <- list()
  a$sds <- sds
  a$lambdas <- lambdas
  a$lls <- lls
  a$bl <- bl
  a$ymin <- ymin
  a$ymax <- ymax

  a
}

#' @title Plotting for Objects of class "lrpsadmmcv"
#' @description  Plots an object of class lrpsadmmcv (output by the \code{lrpsadmm.cv} function).
#' The function plots the support and spectrum of the fit that achieves the best cross-validated negative log-likelihood,
#' along with the cross-validated log-likelihood (and its standard deviation) for each value of gamma
#' @param x An object of class lrpsadmmcv output by the \code{lrpsadmm.cv} function.
#' @importFrom RColorBrewer brewer.pal
#' @importFrom graphics abline legend lines
#' @seealso plot.lrpsadmm plot.lrpsadmmpath
#' @export
plot.lrpsadmmcv <- function(x) {
    my.path <- x
    par(mfrow = c(2, 2))
    image(x$best.fit$fit$S!=0)
    title(paste("Support of selected sparse matrix S.\nLambda=", round(x$best.lambda, 6), "\n Gamma=", x$best.gamma))
    plot(eigen(x$best.fit$fit$L)$values, ylab='Eigenvalue')
    title("Spectrum of selected\n low-rank matrix L")
    A <- lapply(my.path$cross.validated.paths, .prepare.plot.one.cv.path)
    gammas <- unlist(lapply(my.path$cross.validated.paths, function(x){x$gamma}))
    ymin <- min(unlist(lapply(A, function(x) {x$ymin})))
    ymax <- max(unlist(lapply(A, function(x) {x$ymax})))
    xmin <- min(unlist(lapply(A, function(x) {min(x$lambdas)})))
    xmax <- max(unlist(lapply(A, function(x) {max(x$lambdas)})))
    colours <- brewer.pal(n=length(A), name="Set1")
    i <- 1
    plot(A[[i]]$lambdas, A[[i]]$lls, xlim = c(xmin, xmax),
         ylim = c(ymin, ymax), col=colours[i], type = 'l', lwd='2',
         xlab = "-Log10(Lambda)", ylab= "Cross-Validated Log-Lik")
    lines(A[[i]]$lambdas, A[[i]]$lls + A[[i]]$sds, type='l', lty=3, col=colours[i])
    lines(A[[i]]$lambdas, A[[i]]$lls - A[[i]]$sds, type='l', lty=3, col=colours[i])
    for (i in 1:length(A)) {
      lines(A[[i]]$lambdas, A[[i]]$lls, type='l', col=colours[i], lwd='2')
      lines(A[[i]]$lambdas, A[[i]]$lls + A[[i]]$sds, type='l', lty=3, col=colours[i])
      lines(A[[i]]$lambdas, A[[i]]$lls - A[[i]]$sds, type='l', lty=3, col=colours[i])
    }
    abline(v=-log10(my.path$best.lambda), lwd='3')
    y <- (ymax + ymin) / 2
    bidx <- which(gammas == my.path$best.gamma)
    gammas <- as.character(gammas)
    gammas[bidx] <- paste(gammas[bidx],"*", sep='')
    legend("topright", legend = paste("Gamma = ", gammas), col=colours, lty=rep(1, length(colours)))
    title("Cross-Validated log-likelihood\n as a function of Lambda and Gamma")

    par(mfrow = c(1, 1))
}
