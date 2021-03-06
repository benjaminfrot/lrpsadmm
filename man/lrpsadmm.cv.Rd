% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cross_validate.R
\name{lrpsadmm.cv}
\alias{lrpsadmm.cv}
\title{Perform K-fold cross-validation for the Low-Rank plus Sparse estimator}
\usage{
lrpsadmm.cv(X, gammas = c(0.05, 0.1, 0.15), covariance.estimator = cor,
  n.folds = 5, lambdas = NULL, lambda.max = NULL,
  lambda.ratio = 1e-04, n.lambdas = 20, max.sparsity = 0.5,
  max.rank = NA, abs_tol = 1e-05, rel_tol = 0.001, max.iter = 2000,
  mu = 1, verbose = FALSE, seed = NA, zeros = NULL,
  backend = "RcppEigen")
}
\arguments{
\item{X}{n x p data matrix}

\item{gammas}{A real or a vector of reals between 0 and 1 (non inclusive). For each value of gamma
the regularisation path is computed and cross-validation is performed. See examples for
guidance on how to choose reasonable values for gamma. Too high a value might result in the
problem being under-identified and therefore numerical instabilities.}

\item{covariance.estimator}{A function that takes a data matrix and outputs an estimate of its correlation matrix.
Default: the sample correlation matrix output by the \code{cor} function.}

\item{n.folds}{Number of folds for cross-validation. Default 5.}

\item{lambdas}{A decreasing sequence of values of lambda. See Details for the default values.}

\item{lambda.max}{A positive real. Maximum value of lambda. See Details.}

\item{lambda.ratio}{A real between 0 and 1. The smallest value of lambda is given by lambda.max * lambda.ratio. See Details.}

\item{n.lambdas}{A positive integer. The number of values of lambda to
generate according a geometric sequence between lambda.max and lambda.max * lambda.ratio. See Details.}

\item{max.sparsity}{A real between 0 and 1. Abort the computation of the path if S becomes denser than this value.}

\item{max.rank}{A real between 0 and 1. Abort the computuation of the path if the rank of L becomes higher than this value.}

\item{abs_tol}{\code{abs_tol} parameter of the \code{lrpsadmm} function.}

\item{rel_tol}{\code{rel_tol} parameter of the \code{lrpsadmm} function.}

\item{max.iter}{\code{max.iter} parameter of the \code{lrpsadmm} function.}

\item{mu}{\code{mu} parameter of the \code{lrpsadmm} function.}

\item{verbose}{A boolean. Whether to print the value of lambda, gamma, sparsity of S, etc... after each fit}

\item{seed}{Set the seed of the random number generator used for the K folds.}

\item{zeros}{A p x p matrix with entries set to 0 or 1. Whereever its entries are
0, the entries of the estimated S will be forced to 0.}

\item{backend}{The \code{backend} parameter of lrpsadmm. It is one of 'R' or 'RcppEigen'.}
}
\value{
An object of class lrpsadmmcv.
  It contains the values of the mean cross-validated log-likelihood, its standard deviation for each
  pair (lambda, gamma) and an object of class lrpsadmm.path for each value of gamma.
  See the examples for how to access the selected tuning paramters, best fit etc...
}
\description{
Performs K-fold cross-validation in order to select the tuning parameters
  lambda and gamma.

  Recall that the penalty for the LRpS estimator is written as \eqn{\lambda_1 ||S||_1 + \lambda_2 Trace(L)} in the
  objective function of \code{lrpsadmm}.
  This can be equivalently rewritten in terms of the regularisation parameters \eqn{\lambda} and \eqn{\gamma} as follows
  \deqn{\lambda \gamma ||S||_1 + \lambda (1 - \gamma) Trace(L),}
  for \eqn{\gamma \in (0, 1)}.

  For a given value of \eqn{\gamma}, one can perform cross-validation along the regularisation path in order to choose
  \eqn{\lambda}. This function computes the regularisation paths for each value of \eqn{\gamma} supplied as arguments
  and performs cross-validation. The pair (\eqn{\lambda, \gamma}) that produces the smallest cross-validated log-likelihood
  is returned.
}
\details{
Recall that the penalty for the LRpS estimator is written as \eqn{\lambda_1 ||S||_1 + \lambda_2 Trace(L)} in the
  objective function of \code{lrpsadmm}.
  This can be equivalently rewritten in terms of the regularisation parameters \eqn{\lambda} and \eqn{\gamma} as follows
  \deqn{\lambda \gamma ||S||_1 + \lambda (1 - \gamma) Trace(L),}
  for \eqn{\gamma \in (0, 1)}.

  For a given value of \eqn{\gamma}, one can perform cross-validation along the regularisation path in order to choose
  \eqn{\lambda}. This function computes the regularisation paths for each value of \eqn{\gamma} supplied as arguments
  and performs cross-validation. One can then select the pair (\eqn{\lambda, \gamma}) that produces the smallest
  cross-validated log-likelihood.

The function \code{lrpsadmm} is fitted for successive values of \eqn{\lambda} using warm starts.
The sequence of values of \eqn{\lambda} can be provided directly by the user.
It is automatically sorted in decreasing order.
By default, a decreasing sequence of 20 values within a reasonable range is selected as follows.
We set \eqn{\lambda_{max} = \max_{ij, i \neq j} |\Sigma_{ij}|/\gamma} and
\eqn{\lambda_{min} = \lambda_{max}} * \code{lambda.ratio}; then
20 values between \eqn{\lambda_{max}} and \eqn{\lambda_{min}} are taken following
a geometric progression.

Because it does not make much sense to fit this estimator when the sparse estimate S becomes too dense
or if the rank of the low-rank estimate L becomes too high, the computation of the path is aborted
when the sparsity of S reaches \code{max.sparsity} or when the rank of L reaches \code{max.rank}.

Recall that cross-validation overselects. It might be good for prediction purposes but methods such
statibility selection are probably better if the support of S is what is of interest.
}
\examples{
set.seed(0)
## 1 - A simple example: Well-powered dataset
sim.data <- generate.latent.ggm.data(n=2000, p=100, h=5, outlier.fraction = 0.0,
                                    sparsity = 0.02, sparsity.latent = 0.7)
ground.truth <- sim.data$precision.matrix[1:100, 1:100]
ground.truth <- 1 * (( ground.truth - diag(diag(ground.truth)) ) !=0)
X <- sim.data$obs.data;
gammas <- c(0.1, 0.15, 0.2)
# Let the function decide the range of value of Lambda
cvpath <- lrpsadmm.cv(X, gammas = gammas,
                         lambda.ratio = 1e-02, n.lambdas = 30, verbose = TRUE)
best.gamma <- cvpath$best.gamma
best.lambda <- cvpath$best.lambda
best.fit <- cvpath$best.fit$fit # Object of class lrpsadmm
plot(best.fit)
# The value Gamma = 0.15 is selected by X-validation
plot(cvpath) # Plot the outcome. Object of class lrpsadmmcv
# We can look at the path corresponding to this value
best.path <- cvpath$cross.validated.paths[[which(gammas == best.gamma)]]$cross.validated.path
# We know the ground truth, so let's use it to see how well we did:
plot(best.path, ground.truth = ground.truth)

## 2 - Data with outliers: use a robust estimator
set.seed(0)
sim.data <- generate.latent.ggm.data(n=2000, p=50, h=5, outlier.fraction = 0.05,
                                    sparsity = 0.02, sparsity.latent = 0.7)
ground.truth <- sim.data$precision.matrix[1:50, 1:50]
# Remove the elements along the diagonal. Keep a matrix of 0s and 1s
ground.truth <- 1 * (( ground.truth - diag(diag(ground.truth)) ) !=0)
X <- sim.data$obs.data;
# We can afford high values of gamma because n >> p
gammas <- c(0.2, 0.3, 0.4)
# Use the Kendall based estimator:
cvpath <- lrpsadmm.cv(X, gammas = gammas, covariance.estimator = Kendall.correlation.estimator,
                         lambda.ratio = 1e-03, n.lambdas = 30, verbose = TRUE)
plot(cvpath)
best.gamma <- cvpath$best.gamma
best.lambda <- cvpath$best.lambda
best.path <- cvpath$cross.validated.paths[[which(gammas == best.gamma)]]$cross.validated.path
plot(best.path, ground.truth = ground.truth)

## 3 - A tougher problem n is close to p
set.seed(0)
sim.data <- generate.latent.ggm.data(n=150, p=100, h=5, outlier.fraction = 0.0,
                                     sparsity = 0.02, sparsity.latent = 0.7)
ground.truth <- sim.data$precision.matrix[1:100, 1:100]
# Remove the elements along the diagonal. Keep a matrix of 0s and 1s
ground.truth <- 1 * (( ground.truth - diag(diag(ground.truth)) ) !=0)
X <- sim.data$obs.data;
# Since n < p, do not try too high values of gamma. Stay close to 0.
gammas <- c(0.07, 0.1, 0.12)
cvpath <- lrpsadmm.cv(X, gammas = gammas,
                         lambda.ratio = 0.1, n.lambdas = 20, verbose = TRUE)
plot(cvpath) # Plot the outcome

# Clearly the range selected by the function is not good enough
# We need better values for lambda. In that case it is better
# to see which value of lambda yields a very sparse graph
# for a given gamma:
gammas <- c(0.1)
# We set the seed so we can compre the x-validated log-likelihoods between
# two runs of the function
cvpath <- lrpsadmm.cv(X, gammas = gammas, lambda.max = 2.2,
                        lambda.ratio = 0.1, n.lambdas = 20, verbose = TRUE, seed = 0)
plot(cvpath)
gammas <- c(0.12)
cvpath <- lrpsadmm.cv(X, gammas = gammas, lambda.max = 1.7,
                         lambda.ratio = 0.1, n.lambdas = 20, verbose = TRUE, seed = 0)
plot(cvpath)

}
\seealso{
lrpsadmm lrpadmm.path
}
