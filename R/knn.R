#' k nearest neighbor algorithm for multi-variate data
#'
#' @param X data matrix, i.e. observations X dimensions
#' @param K number of clusters to use
#' @param init list of p and mu used for initialization
#' @param Ninit number of samples used per cluster if no init argument is given
#' @param verbose allows print out of progress information; in verbose mode the cluster memberships are added to the output
#' @param tol smaller changes than tol in the objective function indicate convergence, if missing chosen automatically to be 1/5 of the smallest sample variance per dimension
#' @param Niter.max maximum number of admissible iterations
#'
#' @keywords internal
knn <- function(X, K = 2, init, Ninit = 50, verbose = FALSE, tol, Niter.max = 500) {
  ## in case X is no matrix, interpret it as uni-variate case
  if (!is.matrix(X)) {
    X <- matrix(X, ncol = 1)
  }

  if (missing(tol)) {
    ## if tol is missing, then set it well below the minimial
    ## length-scale in the data set
    sdSq <- colVars(X)
    ## use k means clustering with K=Nc as init
    tol <- min(sdSq) / 5
  }

  N <- dim(X)[1]
  Nd <- dim(X)[2]

  if (missing(init)) {
    ## initialize randomly
    pEst <- runif(K) / K
    pEst <- pEst / sum(pEst)
    muEst <- matrix(0, K, Nd)
    ## sample for each component from the base data
    for (i in seq(K)) {
      muEst[i, ] <- colMeans(X[sample.int(N, min(Ninit, N), replace = FALSE), , drop = FALSE])
    }
  } else {
    pEst <- init$p
    muEst <- init$mu
  }

  ## init 1-of-K coding matrix indicating cluster membership
  Kresp <- matrix(1:K, nrow = N, ncol = K, byrow = TRUE)

  ## distance matrix which gets updated in each iteration
  DM <- matrix(0, N, K)

  iter <- 1
  J <- Inf

  if (verbose) {
    message("K nearest neighbors clustering with K =", K, ":\n")
  }

  while (iter < Niter.max) {
    Jprev <- J

    ## "E" step, i.e. find for each data point the cluster with the
    ## smallest euclidean distance

    for (i in seq(K)) {
      DM[, i] <- rowSums(scale(X, muEst[i, ], FALSE)^2)
    }

    ## resp <- 1*(Kresp == apply(DM, 1, which.min))
    resp <- 1 * (1 == matrixStats::rowRanks(DM, ties.method = "first"))
    respM <- matrixStats::colSums2(resp)
    if (any(respM == 0)) {
      warning("Some components are assigned the empty set! Try reducing K.")
      respM[respM == 0] <- 1
    }

    ## "M" step, i.e. given cluster membership, calculate new means
    ## muEst <- sweep(t(resp) %*% X, 1, respM, "/")
    muEst <- sweep(crossprod(resp, X), 1, respM, "/", FALSE)

    ## functional to be minimized
    J <- sum((X - resp %*% muEst)^2)
    delta <- Jprev - J

    if (verbose) {
      message("Iteration", iter, ": J =", J, "; delta =", delta, "\n")
    }
    if (delta < tol) {
      break
    }
    iter <- iter + 1
  }
  if (iter == Niter.max) {
    warning("Maximum number of iterations reached.")
  }

  res <- list(center = muEst, p = colMeans(resp), J = J, delta = delta, niter = iter)
  ## res$cluster <- apply(resp==1, 1, which)
  ## 10x faster
  res$cluster <- ((which(t(resp == 1)) - 1) %% K) + 1

  invisible(res)
}
