## EM for MSMM (Multi-variate Student t Mixture Model) with Nc components
## if init is not specified, then knn is used to initialize means,
## cluster weights and covariance matrices (taken from the knn
## determined clusters)

## this function follows the details given in ref. "Robust mixture
## modelling using the t distribution". D. Peel and G.J. McLachlan,
## Statistics and Computing (2000), 10, 339-348,

EM_msmm <- function(X, Nc, init, Ninit = 50, verbose = TRUE, Niter.max = 500, tol = 1e-1) {
  ## in case X is no matrix, interpret it as uni-variate case
  if (!is.matrix(X)) {
    X <- matrix(X, ncol = 1)
  }

  N <- dim(X)[1]
  Nd <- dim(X)[2]

  ## initialize randomly
  if (missing(init)) {
    ## assume that the sample is ordered randomly
    ind <- seq(1, N - Nc, length = Ninit)
    knnInit <- list(mu = matrix(0, nrow = Nc, ncol = Nd), p = (1 / seq(1, Nc)) / sum(seq(1, Nc)))
    for (k in seq(Nc)) {
      knnInit$mu[k, ] <- colMeans(X[ind + k - 1, , drop = FALSE])
    }
    ## use k means clustering with K=Nc as init; ignore warnings
    ## as we may hit the maximal number of iterations
    suppressWarnings(KNN <- knn(X, K = Nc, init = knnInit, verbose = verbose, Niter.max = 50))
    pEst <- KNN$p
    cmin <- which.min(pEst)
    muEst <- KNN$center
    ## nuEst <- rlnorm(Nc, log(20), log(5)/1.96) ## avoid randomness during initialization
    nuEst <- rep(12, times = Nc)
    covEst <- array(0, dim = c(Nc, Nd, Nd))
    Xtau <- sqrt(colVars(X))
    for (i in seq(Nc)) {
      if (i == cmin) next
      ind <- KNN$cluster == i
      if (sum(ind) > 10) {
        covKNN <- as.matrix(cov(X[ind, , drop = FALSE]))
        R <- cov2cor(covKNN)
        tau <- sqrt(diag(covKNN))
        ## set variances below or equal to 0 to the sample variance
        tau[tau <= 0] <- Xtau[tau <= 0]
      } else {
        R <- diag(Nd)
        tau <- Xtau
      }
      ## tauR <- rlnorm(Nd, log(tau), log(3)/1.96)
      ## covEst[i,,] <- diag(tauR, Nd, Nd) %*% R %*% diag(tauR, Nd, Nd)
      ## ensure that the smallest variance is not less than the global
      ## variance divided by 100... which is to stabilize things
      tau <- pmax(tau, Xtau / 100)
      covEst[i, , ] <- diag(tau, Nd, Nd) %*% R %*% diag(tau, Nd, Nd)
    }
    ## tauR <- rlnorm(Nd, log(Xtau), log(3)/1.96)
    ## covEst[cmin,,] <- diag(tauR, Nd, Nd) %*% diag(Nd) %*% diag(tauR, Nd, Nd)
    covEst[cmin, , ] <- diag(Xtau, Nd, Nd) %*% diag(Nd) %*% diag(Xtau, Nd, Nd)
    muEst[cmin] <- sum(pEst * muEst)
  } else {
    pEst <- init$p
    muEst <- init$center
    nuEst <- init$nu
    covEst <- init$cov
  }

  iter <- 0
  logN <- log(N)
  traceLli <- c()
  Dlli <- Inf

  ## initialize component and element wise log-likelihood matrix
  lli <- array(-20, dim = c(N, Nc))

  ## and the log(U) matrix which is the inverse mahalanobis distances
  ## decorated with nuEst
  lU <- array(0, dim = c(N, Nc))

  if (verbose) {
    message("EM multi-variate student t with Nc =", Nc, ":\n")
  }

  nu_ml <- function(c1) {
    function(nu) {
      (log(nu / 2) - digamma(nu / 2) + c1 + digamma((nu + Nd) / 2) - log((nu + Nd) / 2))^2
    }
  }

  while (iter < Niter.max) {
    ## calculate responsabilities from the likelihood terms;
    ## calculations are done in log-space to avoid numerical
    ## difficulties if some points are far away from some
    ## component and hence recieve very low density
    for (i in seq(Nc)) {
      lli[, i] <- log(pEst[i]) + mvtnorm::dmvt(X, muEst[i, ], as.matrix(covEst[i, , ]), nuEst[i], log = TRUE)
    }
    ## ensure that the log-likelihood does not go out of numerical
    ## reasonable bounds
    lli <- apply(lli, 2, pmax, -30)
    ## lnresp <- apply(lli, 1, log_sum_exp)
    lnresp <- matrixStats::rowLogSumExps(lli)
    ## the log-likelihood is then given by the sum of lresp
    lliCur <- sum(lnresp)
    traceLli <- c(traceLli, lliCur)
    if (iter > 1) {
      ## Dlli is the slope of the log-likelihood evaulated with
      ## a second order method
      Dlli <- (traceLli[iter + 1] - traceLli[iter - 1]) / 2
    }
    if (verbose) {
      message("Iteration", iter, ": log-likelihood", lliCur, "; Dlli =", Dlli, "\n")
    }
    if (Dlli < tol) {
      break
    }
    ## ... and the responisbility matrix follows from this by
    ## appropiate normalization.
    lresp <- sweep(lli, 1, lnresp, "-") ## Eq. 16
    resp <- exp(lresp)

    ## calculate additional weights of the U matrix aka latent
    ## tail mass of a point
    for (i in seq(Nc)) {
      Xc_i <- sweep(X, 2L, muEst[i, ])
      Sigma_i <- as.matrix(covEst[i, , ])
      ## in rare cases the covariance matrix becomes (almost)
      ## singular in which case the alternative cholesky
      ## factorization gives more stable results
      maha_dist <- tryCatch(
        mahalanobis(Xc_i, FALSE, Sigma_i),
        error = function(e) {
          ## see https://stats.stackexchange.com/questions/147210/efficient-fast-mahalanobis-distance-computation
          ## also adding eps to the diagonal to further stabilize the computation
          L_i <- t(chol(Sigma_i + diag(5 * .Machine$double.eps, Nd, Nd)))
          y_i <- forwardsolve(L_i, t(Xc_i))
          colSums(y_i^2)
        }
      )
      lU[, i] <- log(nuEst[i] + Nd) - log(nuEst[i] + maha_dist) ## Eq. 20
    }

    ## mean probability to be in a specific mixture component ->
    ## updates pEst
    ## lzSum <- apply(lresp, 2, log_sum_exp)
    lzSum <- colLogSumExps(lresp)
    zSum <- exp(lzSum)
    ## zSum <- colSums(resp)
    ## pEst <- zSum/N ## Eq. 29
    pEst <- exp(lzSum - logN)

    ## make sure it is scale to exactly 1 which may not happen due
    ## to small rounding issues
    pEst <- pEst / sum(pEst)

    ## product of u weight and responsabilities
    lW <- lresp + lU ## intermediate formed, i.e. tau_ij * u_ij
    ## wSum <- exp(apply(lW, 2, log_sum_exp))
    wSum <- exp(colLogSumExps(lW))
    ## wSum <- colSums(W)

    ## now obtain new estimates for each component of the mixtures
    ## of their mu vector and covariance matrices
    for (i in seq(Nc)) {
      ## xc <- sqrt(W[,i]) * sweep(X, 2, muEst[i,], check.margin = FALSE)
      ## covEst[i,,] <- crossprod(xc) / zSum[i]    ## Eq. 31
      ## muEst[i,] <- colSums(W[,i] * X) / wSum[i] ## Eq. 30
      ## xc <- exp(0.5 * lW[,i]) * sweep(X, 2, muEst[i,], check.margin = FALSE)
      ## covEst[i,,] <- crossprod(xc) / zSum[i]          ## Eq. 31
      xc <- exp(0.5 * (lW[, i] - lzSum[i])) * sweep(X, 2, muEst[i, ], check.margin = FALSE)
      covEst[i, , ] <- crossprod(xc) ## Eq. 31 (divisor moved to xc)
      muEst[i, ] <- colSums(exp(lW[, i]) * X) / wSum[i] ## Eq. 30
      ## ensure that diagonal stays non-zero
      for (j in 1:Nd) {
        covEst[i, j, j] <- max(covEst[i, j, j], .Machine$double.eps)
      }
    }

    ## finally get the new nu estimates via numerical solution
    ## first calculate necessary constants which don't involve v_i
    ## c1 <- 1 + colSums(resp * (log(U) - U)) / zSum
    c1 <- 1 + colSums(resp * (lU - exp(lU))) / zSum
    for (i in seq(Nc)) {
      nuEstML <- optimize(nu_ml(c1[i]), c(0, 150)) # Eq. 32
      if (is.na(nuEstML$objective)) {
        warning("Component ", i, " in iteration ", iter, " failed convergence.")
      } else if (nuEstML$objective > 1e-3 & nuEstML$minimum < 50) {
        ## only warn if we had trouble finding the minimum
        ## when below 50, larger values are anyway normals
        warning("Component ", i, " in iteration ", iter, " had convergence problems.\nObjective function = ", nuEstML$objective, "\n")
      }
      nuEst[i] <- nuEstML$minimum
    }

    iter <- iter + 1
  }
  if (iter == Niter.max) {
    warning("Maximum number of iterations reached.")
  }

  ## degrees of freedom
  ## covariance matrix df per component
  cov.df <- (Nd - 1) * Nd / 2 + Nd
  df <- Nc * (Nd + cov.df) + Nc + Nc - 1

  ## sort by largest weight
  o <- order(pEst, decreasing = TRUE)
  pEst <- pEst[o]
  muEst <- muEst[o, , drop = FALSE]
  covEst <- covEst[o, , , drop = FALSE]
  nuEst <- nuEst[o, drop = FALSE]

  if (Nd != 1) {
    rhoEst <- array(apply(covEst, 1, cov2cor), c(Nd, Nd, Nc))
    rhoEst <- apply(rhoEst, 3, function(x) x[lower.tri(x)])
    tauEst <- sqrt(t(apply(covEst, 1, diag)))
  } else {
    rhoEst <- NULL
    tauEst <- sqrt(as.vector(covEst))
  }

  invisible(list(cov = covEst, center = muEst, nu = nuEst, p = pEst, rho = rhoEst, tau = tauEst, lli = lliCur, df = df, Dlli = Dlli, niter = iter))
}
