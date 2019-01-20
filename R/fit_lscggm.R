
.lscggm_obj_func <- function(ps, opts) {
  p <- opts$p
  q <- opts$q
  n <- opts$n
  CX <- cor(opts$X)
  CZX <- cor(opts$X, opts$Z)
  CZ <- cor(opts$Z)
  l1 <- opts$lp1 * opts$mu
  l2 <- opts$lp2 * opts$mu
  X <- ps$S - ps$L
  XX <- X[1:p,]
  XXZ <- X[(p+1):(p+q),]
  eig <- eigen(XX, symmetric = TRUE)
  evals <- eig$val
  evals[abs(evals) < 1e-08] <- 1e-9
  XX <- eig$vec %*% diag(evals) %*% t(eig$vec)
  if (any(evals < 0)) {
    return(NaN)
  }

  # Log-Likelihood
  ll1 <- sum(diag(CX %*% XX)) - sum(log(evals))
  ll2 <- 2 * sum(diag(CZX %*% XXZ))
  ll3 <- sum(diag(solve(XX) %*% t(XXZ) %*% CZ %*% XXZ))
  ll <- ll1 + ll2 + ll3
  # + Penalty
  ll <- ll + l1 * sum(abs(ps$S)) + l2 * sum(svd(ps$L)$d)

  ll
}

###### Proximal operator implementation ######
.prox_g <- function(X, rho) {
  svd.res <- svd(X)
  ss <- svd.res$d - rho
  ss[ss < 0] <- 0

  svd.res$u %*% diag(ss) %*% t(svd.res$v)
}

.prox_f <- function(X, p, q) {
  A = X[1:p,]
  eig.res = eigen(0.5 * (A + t(A)), symmetric = TRUE);
  ss <- eig.res$values
  ss[ss < 0] <- 0

  rbind(eig.res$vectors %*% diag(ss) %*% t(eig.res$vectors), X[(p+1):(p+q),])
}

.compute_proximal_operator <- function (X, p, q, rho, tol, maxiter) {
  P <- Q <- X * 0
  for (i in 1:maxiter) {
    Y <- .prox_g(X + P, rho)
    P <- X + P - Y
    oX <- X
    X <- .prox_f(Y + Q, p, q)
    Q <- Y + Q - X
    if(norm(oX, "F") < 1e-06) {
      diff = norm(X - oX, "F")
    } else {
      diff = norm(X - oX, "F") / norm(oX, "F");
    }
    if (diff < tol)
      break()
  }

  X
}



###### Update of A #####
step1_compute_gradient <- function(theta, Sx, Sxy, Sy, N, comp_grad, ps) {
  mu <- ps$mu
  flag        = 0;
  # Check that SX is positive semi-definite
  is.pos.sem.def <- is.positive.definite(theta$yy)
  if (!is.pos.sem.def) {
    return(list(flag=1, obj=Inf, grad=theta))
  }

  cyy  = chol(theta$yy);
  logdetyy <- 2 * sum(log(diag(cyy)));
  if((is.nan(logdetyy) | is.infinite(logdetyy))) {
    return(list(flag=1, obj=Inf, grad=theta))
  }

  #icyy	 = chol2inv(cyy);
  ithetayy = chol2inv(cyy);
  txyityy  = theta$xy %*% ithetayy;
  XtXth    = Sx %*% txyityy;
  txyXtXth = t(theta$xy) %*% Sx %*% txyityy;

  l1 = matrix.trace( theta$yy%*%Sy );
  l2 = matrix.trace( Sxy%*%t(theta$xy) );
  l3 = matrix.trace( txyXtXth );
  value = 0.5*l1 + l2 + 0.5*l3 - 0.5*N*logdetyy ;
  value = value / N;
  A <- rbind(theta$yy, theta$xy)
  value = value + 0.5 * ps$mu * frobenius.norm(A - ps$Shat + ps$Lhat + ps$Uhat/mu)**2


  if(comp_grad) {
    grad <- list()
    grad$xy <- (Sxy + XtXth)/N + ps$Uhat[(ps$p+1):(ps$p + ps$m),] + mu * (theta$xy - ps$Shat[(ps$p+1):(ps$p + ps$m),] + ps$Lhat[(ps$p+1):(ps$p + ps$m),]);
    grad$yy <- 0.5*(Sy - N*ithetayy - ithetayy%*%txyXtXth)/N + ps$Uhat[1:ps$p,] + mu * (theta$yy - ps$Shat[(1):(ps$p),] + ps$Lhat[(1):(ps$p),]);
    grad$yy <- 0.5 * (grad$yy + t(grad$yy))
  } else{
    grad <- NULL
  }
  result <- list()
  result$obj <- value
  result$grad <- grad
  result$flag <- flag

  result
}

latent.scggm.update.A <- function(Z, X, ps, A, maxiter=200, tol=1e-7) {

  Sx	= t(Z) %*% Z
  Sy 	= t(X) %*% X
  Sxy	= t(Z) %*% X
  N 	= dim(Z)[1]
  p = dim(X)[2]
  m = dim(Z)[2]

  nobj	= 10
  bconv	= 0
  obj	= rep(0, maxiter)
  L	= 1
  thk_0 	= 2/3
  ls_maxiter = 300
  eta <- 1.5

  # Initialise
  theta <- list()
  theta$yy <- A[1:p, ]
  theta$yy <- 0.5 * (theta$yy + t(theta$yy))
  theta$xy <- A[(p+1):(p+m),]
  grad  = step1_compute_gradient(theta, Sx, Sxy, Sy, N, F, ps);
  obj1 <- grad$obj
  flag <- grad$flag
  if(flag == 1) {
    warning("Initial estimate of SX is not positive-definite. We use a diagonal matrix instead.")
    theta$yy <- diag(p)
    grad  = step1_compute_gradient(theta, Sx, Sxy, Sy, N, F, ps);
    obj1 <- grad$obj
    flag <- grad$flag
  }
  obj[1] = obj1;

  xk      = theta;
  zk   	= theta;
  thk     = thk_0;

  for (iter in 2:maxiter) {
    thk  = (sqrt( thk^4 + 4 * thk^2 ) - thk^2) / 2;
    y <- list()
    y$xy = (1 - thk) * xk$xy + thk * zk$xy;
    y$yy = (1 - thk) * xk$yy + thk * zk$yy;
    grady = step1_compute_gradient( y, Sx, Sxy, Sy, N, TRUE, ps);
    fyk = grady$obj
    flagy = grady$flag
    grady = grady$grad
    ik = 0;
    while(T) {
      zk_grady <- list()
      zk_grady$xy = zk$xy - 1/(L*thk) * grady$xy;
      zk_grady$yy = zk$yy - 1/(L*thk) * grady$yy;
      zk1 		= zk_grady;

      y_grady <- list()
      y_grady$xy	= y$xy - 1/L * grady$xy;
      y_grady$yy	= y$yy - 1/L * grady$yy;
      xk1         = y_grady;

      gradxk1 <- step1_compute_gradient(xk1, Sx, Sxy, Sy, N, FALSE, ps)
      fxk1 <- gradxk1$obj

      flagxk1 <- gradxk1$flag
      if(is.positive.definite(zk1$yy)) {
        flagzk1 <- 0
      } else{
        flagzk1 <- 1
      }
      if ( flagzk1 == 0 & flagy ==0 & flagxk1 ==0 ) {
        xk1_y <- list()
        xk1_y$xy    = xk1$xy - y$xy;
        xk1_y$yy    = xk1$yy - y$yy;
        lfxk1_y = fyk + sum(grady$xy * xk1_y$xy) + sum(grady$yy * xk1_y$yy)
        diffxk1y <- list()
        diffxk1y$xy = xk1$xy - y$xy;
        diffxk1y$yy = xk1$yy - y$yy;
        RHS         = lfxk1_y + L/2 *(sum(diffxk1y$xy**2) + sum(diffxk1y$yy**2));
        if(fxk1 <= RHS + tol) {
          xk = xk1;
          zk = zk1;
          bconv = 1;
          break()
        }
      }

      ik = ik + 1;

      if ( ik > ls_maxiter ) {
        bconv = 0;
        iter  = max(1, iter - 1);
        Theta = xk;
        break()
      }
      L = L * eta;
    }
    obj[iter]  = fxk1
    if(bconv == 0) {
      break()
    }

    if ( iter > nobj + 1) {
      value           = obj[iter];
      prevVals        = obj[iter - nobj];
      avgimprovement  = abs(prevVals - value)/nobj;
      relAvgImpr      = avgimprovement / abs(value)

      if ( relAvgImpr < tol ) {
        bconv = 1;
        break()
      }
    }
  }
  Theta = xk;
  A <- rbind(Theta$yy, Theta$xy)
  obj   = obj[1:iter];
  termcode <- 0
  if(bconv != 1) {
    termcode <- -1
  }
  list(A=A, objs=obj, iter=length(obj), termcode=termcode)

}


######
.lscggm.updateA <- function(ps, opts) {
  lp1 <- opts$lp1
  lp2 <- opts$lp2
  mu <- opts$mu
  Shat <- ps$Shat
  Lhat <- ps$Lhat
  Uhat <- ps$Uhat
  ps$m <- opts$q
  ps$p <- opts$p
  ps$mu <- mu

  out <- latent.scggm.update.A(opts$Z, opts$X, ps, ps$A, tol = opts$tol*0.01,
                               maxiter = opts$prox_maxiter)
  A <- out$A
  AX <- A[1:opts$p,]
  AX <- 0.5 * (AX + t(AX))
  A[1:opts$p,] <- AX
  ps$A <- A

  ps$exit <- FALSE
  ps
}

.lscggm.updateAlpha <- function(ps, opts) {
  alpha <- ps$alpha[length(ps$alpha)]
  alpha <- 0.5 * (1 + sqrt(1 + 4 * alpha**2))
  ps$alpha <- c(ps$alpha, alpha)

  ps
}

.lscggm.updateL <- function(ps, opts) {
  lp1 <- opts$lp1
  lp2 <- opts$lp2
  mu <- opts$mu
  C <- opts$Sigma
  A <- ps$A
  S <- ps$S
  Uhat <- ps$Uhat
  ps$prevL <- ps$L

  X1 <- S - A - (Uhat / mu)
  L <- .compute_proximal_operator(X1, opts$p, opts$q, tol=0.01*opts$tol, rho = lp2,
                                  maxiter = opts$prox_maxiter)
  LX <- L[1:opts$p,]
  LX <- 0.5 * (LX + t(LX))
  L[1:opts$p,] <- LX

  ps$L <- L
  ps$exit <- FALSE
  ps
}

.lscggm.updateS <- function(ps, opts) {
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
  S <- S

  ps$S <- S
  ps$exit <- FALSE
  ps
}

.lscggm.updateShatLhatUhat <- function(ps, opts) {
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

.lscggm.updateU <- function(ps, opts) {
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

  ps$U <- U
  ps$exit <- FALSE
  ps
}

.lscggm.update.parameters <- function(ps, opts) {

  pL <- ps$L
  pS <- ps$S
  pU <- ps$U
  ps$exit <- FALSE
  fs <- c(.lscggm.updateA, .lscggm.updateS, .lscggm.updateL, .lscggm.updateU)
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
    ps <- .lscggm.updateAlpha(ps, opts)
    ps <- .lscggm.updateShatLhatUhat(ps, opts)
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
      msg3 <- "This might due to the problem being too ill-posed. For example, the value of gamma might be too large given the ratio p/n."
      msg4 <- "You could consider trying a smaller value of mu.\n"
      warning(paste(msg1, msg2, msg3, msg4))
      ps$exit <- TRUE
    }
  }

  list(ps=ps, opts=opts)
}

#######################
fit.lscggm <- function(Z, X, Lambda1, Lambda2,
                       init=NULL, maxiter=2000,
                       mu=0.1, tol=1e-05, eta=0.999,
                       prox_maxiter=50,
                       print_progress=TRUE, print_every=20) {

  n <- dim(X)[1]
  p <- dim(X)[2] # Number of variables we model
  q <- dim(Z)[2] # Number of variables we condition on

  options <- list()
  options$mu <- mu
  options$X <- scale(X)
  options$Z <- scale(Z)
  options$lp1 <- Lambda1 / mu
  options$lp2 <- Lambda2 / mu
  options$eta <- eta
  options$n <- n
  options$p <- p
  options$q <- q
  options$tol <- tol
  options$prox_maxiter <- prox_maxiter

  if (is.null(init)) {
    Sigma <- cor(cbind(X, Z))
    S <- MASS::ginv(Sigma)
    L <- S * 0.01
    A <- S - L
    U <- mu * (A - S + L)
    S <- S[,1:p]
    L <- L[,1:p]
    A <- A[,1:p]
    U <- U[,1:p]
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

  diffs <- c()
  lls <- c(NaN)
  parameters$termcode <- -1
  parameters$termmsg <- "Maximum number of iterations reached."
  for (i in 1:maxiter) {
    L <- .lscggm.update.parameters(parameters, options)
    new_parameters <- L$ps
    options <- L$opts

    if(parameters$exit) {
      parameters$termcode <- -3
      parameters$termmsg <-
        "Algorithm was restarted 100 times without improvement."
      break()
    }

    if(any(diag(new_parameters$S) == 0)) {
      parameters$S[,] <- NA
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
    lls <- c(lls, .lscggm_obj_func(parameters, options))
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
        avg_chg <- mean(abs(diff(last_lls)) / abs(last_lls[1:(print_every)]),
                        na.rm=T)
        print(paste("Iteration:", i,
                    "Log-Likelihood:", ll,
                    "Avg relative log-likelihood change:", avg_chg,
                    "Relative change in parameters:", diff))
      }
    }
  }

  parameters$iter <- i
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

  SX <- parameters$S[1:p,]
  SXZ <- parameters$S[(p+1):(p+q),]
  LX <- parameters$L[1:p,]
  LXZ <- parameters$L[(p+1):(p+q),]
  parameters$SX <- SX
  parameters$LX <- LX
  parameters$SXZ <- SXZ
  parameters$LXZ <- LXZ
  parameters
}
