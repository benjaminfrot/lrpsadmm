#' Fit a path of the LRpS
#' @import Matrix
#' @export
fit.low.rank.plus.sparse.path <- function(Sigma,
                                          gamma,
                                          n,
                                          lambdas = NULL,
                                          lambda.max=NULL,
                                          lambda.ratio=1e-4,
                                          n.lambdas=20,
                                          max.sparsity=0.5,
                                          max.rank=NA,
                                          tol=1e-05,
                                          max.iter=2000,
                                          mu=0.1,
                                          verbose=FALSE) {
  p <- dim(Sigma)[1]

  if(is.null(lambdas)) {
    if (is.null(lambda.max)) {
      lambda.max <- max(abs(Sigma - diag(diag(Sigma)))) * 10
    }
    lambda.min <- lambda.max * lambda.ratio
    reason <- lambda.min / lambda.max
    lambdas <- lambda.max * reason **((0:n.lambdas)/n.lambdas)
  }
  fit <- NULL
  path <- list()
  counter <- 1
  for (lambda in lambdas) {
    l1 <- gamma * lambda
    l2 <- (1 - gamma) * lambda
    fit <- fit.low.rank.plus.sparse(Sigma, l1, l2, n, init=fit, 
                                    maxiter = max.iter,
                                    mu=mu, tol=tol, print_progress = FALSE)
    if (fit$termcode == -2) {
      next()
    }
    path[[counter]] <- list()
    path[[counter]]$lambda <- lambda
    path[[counter]]$gamma <- gamma
    path[[counter]]$fit <- fit

    rank.L <- as.numeric(Matrix::rankMatrix(fit$L))
    S <- (fit$S!=0) - diag(diag(fit$S!=0))
    sparsity <- 0.5 * sum(S)
    n.edges <- sparsity
    sparsity <- sparsity / choose(p, 2)

    path[[counter]]$number.of.edges <- n.edges
    path[[counter]]$sparsity <- sparsity
    path[[counter]]$rank.L <- rank.L
    counter <- counter + 1

    if (!is.na(max.rank)) {
      if (rank.L > max.rank)
        break()
    }
    if (sparsity > max.sparsity)
      break()

    if (verbose) {
      print(paste("Fitting with gamma=", gamma, " and lambda=", lambda,
                  "Sparsity:", sparsity, "Rank of L:", rank.L))
    }
  }

  path
}

#' Plot the output fit.low.rank.plus.sparse.path
#' @export
show.low.rank.plus.sparse.path <- function(lrps.path, ground.truth=NULL) {
  gamma <- lrps.path[[1]]$gamma
  lambdas <- ranks <- sparsities <- edges <- c()
  if (!is.null(ground.truth)) {
    ground.truth <- (ground.truth!=0) - diag(diag(ground.truth!=0))
  }
  prs <- rs <- c()
  for (i in 1:length(lrps.path)) {
    lambdas <- c(lambdas, lrps.path[[i]]$lambda)
    ranks <- c(ranks, lrps.path[[i]]$rank.L)
    edges <- c(edges, lrps.path[[i]]$number.of.edges)
    sparsities <- c(sparsities, lrps.path[[i]]$sparsity)
    if(!is.null(ground.truth)) {
      S <- lrps.path[[i]]$fit$S
      S <- (S!=0) - diag(diag(S!=0))
      if (sum(S) > 0) {
        pr <- sum(S * ground.truth) / sum(S)
        r <- sum(S * ground.truth) / sum(ground.truth)
        prs <- c(prs, pr)
        rs <- c(rs, r)
      }
    }
  }

  par(mfrow=c(2,2))
  plot(-log(lambdas), sparsities, xlab="-Log10(Lambda)", ylab="Sparsity of S")
  plot(-log(lambdas), edges, xlab="-Log10(Lambda)", ylab="Number of Edges")
  plot(-log(lambdas), ranks, xlab="-Log10(Lambda)", ylab="Rank of L")
  if (!is.null(ground.truth)) {
    plot(rs, prs, xlab="Recall", ylab="Precision")
  }
  par(mfrow=c(1,1))
}
