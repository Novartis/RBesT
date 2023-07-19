## EM for GMM with Nc components
EM_gmm <- function(x, Nc, mix_init, Ninit=50, verbose=FALSE, Niter.max=500, tol, Neps, eps=c(weight=0.005,alpha=0.005,beta=0.005)) {
    N <- length(x)
    assert_that(N+Nc >= Ninit)

    ## check data for 0  values which are problematic, but may be
    ## valid. Moving these to eps ensures proper handling during fit.
    x0 <- x==0
    if(any(x0)) {
        message("Detected ", sum(x0), " value(s) which are exactly 0.\nTo avoid numerical issues during EM such values are moved to smallest eps on machine.")
        x[x0] <- .Machine$double.eps
    }

    ## temporaries needed during EM
    Lx <- matrix(log(x), ncol=Nc, nrow=N)

    xRep <- matrix(x, ncol=Nc, nrow=N)

    ## initialize randomly using KNN
    if(missing(mix_init)) {
        ## assume that the sample is ordered randomly
        ind <- seq(1,N-Nc,length=Ninit)
        knnInit <- list(mu=matrix(0,nrow=Nc,ncol=1), p=rep(1/Nc, times=Nc))
        for(k in seq(Nc))
            knnInit$mu[k,1] <- mean(x[ind+k-1])
        KNN <- suppressWarnings(knn(x, K=Nc, init=knnInit, Niter.max=50))
        muInit <- rep(mean(x), times=Nc)
        varInit <- rep(1.5*var(x), times=Nc)
        for(k in 1:Nc) {
            kind <- KNN$cluster == k
            if(sum(kind) > 10) {
                muInit[k] <- KNN$center[k]
                varInit[k] <- var(x[kind])
            }
        }
        ## relocate the component with the least weight in the center
        ## and assign the sample variance to it; the idea is that we
        ## expect an informative component and a heavy tailed
        ## background which is best "absorbed" if a wide component is
        ## place initially at the center of the data
        cmin <- which.min(KNN$p)
        varInit[cmin] <- var(x)
        muInit[cmin] <- sum(KNN$center * KNN$p)
        ##varInit <- rlnorm(Nc, log(varInit), log(2.5)/1.96)
        bInit <- muInit / varInit
        aInit <- muInit * bInit
        mixEst <- rbind(KNN$p, aInit, bInit)
        dlink(mixEst) <- identity_dlink
        rownames(mixEst) <- c("w", "a", "b")
    } else {
        mixEst <- mix_init
    }

    if(verbose) {
        message("EM for gamma mixture model.\n")
        message("Initial estimates:\n")
        print(mixEst)
    }

    ## mixEst parametrization during fitting
    mixEstPar <- mixEst
    mixEstPar[1,] <- logit(mixEst[1,,drop=FALSE])
    mixEstPar[2,] <- log(mixEst[2,])
    mixEstPar[3,] <- log(mixEst[3,])
    rownames(mixEstPar) <-  c("w", "la", "lb")

    ## the optimizer needs a fixed range where search log-alpha
    MLrange <- c(min(mixEstPar[2,]) - log(1e4), max(mixEstPar[2,]) + log(1e4))

    ## in case tolerance is not specified, then this criteria is
    ## ignored
    if(missing(tol)) {
        checkTol <- FALSE
        tol <- -1
    } else
        checkTol <- TRUE

    if(missing(Neps)) {
        ## in case tolerance has been declared, but Neps not, we flag
        ## to disable checking of running mean convergence check
        checkEps <- FALSE
        Neps <- 5
    } else
        checkEps <- TRUE

    ## if nothing is specified, we declare convergence based on a
    ## running mean of differences in parameter estimates
    if(!checkTol & !checkEps) {
        checkEps <- TRUE
    }

    assert_that(Neps > 1)
    assert_that(ceiling(Neps) == floor(Neps))


    ## eps can also be given as a single integer which is interpreted
    ## as number of digits
    if(length(eps) == 1) eps <- rep(10^(-eps), 3)

    iter <- 0
    logN <- log(N)
    traceMix <- list()
    traceLli <- c()
    Dlli <- Inf
    runMixPar <- array(-Inf, dim=c(Neps,3,Nc), dimnames=list(NULL, rownames(mixEstPar), NULL ))
    runOrder <- 0:(Neps-1)
    Npar <- Nc + 2*Nc
    if(Nc == 1) Npar <- Npar - 1

    ## find alpha and beta for a given component in log-space
    gmm_ml <- function(c1) {
        function(la) {
            (c1 - digamma(exp(la)) + la)^2
        }
    }

    gmm_ml_grad <- function(c1) {
        function(la) {
            a <- exp(la)
            val <- (c1 - digamma(a) + la)
            grad <- 2 * val * (1 - trigamma(a) * a)
            grad
        }
    }

    while(iter < Niter.max) {
        ## calculate responsabilities from the likelihood terms;
        ## calculations are done in log-space to avoid numerical difficulties if some points are far away from some component and hence recieve very low density
        ##li <- t(matrix(abmEst[,1] * dgamma(xRep, abmEst[,2], abmEst[,3]), nrow=Nc))

        ##lli <- t(matrix(log(mixEst[1,]) + dgamma(xRep, mixEst[2,], mixEst[3,], log=TRUE), nrow=Nc))
        w <- mixEst[1,]
        a <- mixEst[2,]
        b <- mixEst[3,]
        ## Gamma density: x^(a-1) * exp(-b * x) * b^a / Gamma(a)
        lli  <- sweep(sweep(Lx, 2, a-1, "*") - sweep(xRep, 2, b, "*"), 2, a * log(b) - lgamma(a) + log(w), "+")

        ## ensure that the log-likelihood does not go out of numerical
        ## reasonable bounds
        lli <- apply(lli, 2, pmax, -30)

        ##lnresp <- apply(lli, 1, log_sum_exp)
        lnresp <- matrixStats::rowLogSumExps(lli)
        ## the log-likelihood is then given by the sum of lresp norms
        lliCur <- sum(lnresp)
        ## record current state
        traceMix <- c(traceMix, list(mixEst))
        traceLli <- c(traceLli, lliCur)
        if(iter > 1) {
            ## Dlli is the slope of the log-likelihood evaulated with
            ## a second order method
            Dlli <- (traceLli[iter+1] - traceLli[iter - 1])/2
        }
        if(Nc > 1) {
            smean <- apply(runMixPar[order(runOrder),,,drop=FALSE], c(2,3), function(x) mean(abs(diff(x))))
            eps.converged <- sum(sweep(smean, 1, eps, "-") < 0)
        } else {
            smean <- apply(runMixPar[order(runOrder),-1,,drop=FALSE], c(2,3), function(x) mean(abs(diff(x))))
            eps.converged <- sum(sweep(smean, 1, eps[-1], "-") < 0)
        }
        if(is.na(eps.converged)) eps.converged <- 0
        if(verbose) {
            message("Iteration ", iter, ": log-likelihood = ", lliCur, "; Dlli = ", Dlli, "; converged = ", eps.converged, " / ", Npar, "\n", sep="")
        }
        if(checkTol & Dlli < tol) {
            break
        }
        if(iter >= Neps & checkEps & eps.converged == Npar) {
            break
        }
        ## ... and the (log) responseability matrix follows from this by
        ## appropiate normalization.
        lresp <- sweep(lli, 1, lnresp, "-")
        resp <- exp(lresp)

        ## mean probability to be in a specific mixture component -> updates
        ## mixEst first row
        ##lzSum <- apply(lresp, 2, log_sum_exp)
        lzSum <- colLogSumExps(lresp)
        zSum <- exp(lzSum)
        mixEst[1,] <- exp(lzSum - logN)

        ## make sure it is scale to exactly 1 which may not happen due
        ## to small rounding issues
        mixEst[1,] <- mixEst[1,] / sum(mixEst[1,])

        ##lrx <- apply(Lx + lresp, 2, log_sum_exp)
        lrx <- colLogSumExps(Lx + lresp)
        resp_zscaled  <- exp(sweep(lresp, 2, lzSum, "-"))
        c1 <- colSums(Lx * resp_zscaled) + lzSum - lrx
        c2 <- lzSum - lrx

        ## now solve for new alpha and beta estimates
        for(i in 1:Nc) {
            Lest <- optimize(gmm_ml(c1[i]), MLrange)
            ##theta <- c(log(mixEst[2,i]))
            ##Lest <- optim(theta, gmm_ml(c1[i]), gr=gmm_ml_grad(c1[i]), method="BFGS", control=list(maxit=500))
            if(abs(Lest$objective) > 1E-4) {
                warning("Warning: Component", i, "in iteration", iter, "had convergence problems!")
            }
            mixEstPar[2,i] <- Lest$minimum
            mixEstPar[3,i] <- Lest$minimum + c2[i]
            mixEst[c(2,3),i] <- exp(mixEstPar[c(2,3),i])
        }

        mixEstPar[1,] <- logit(mixEst[1,,drop=FALSE])
        ind <- 1 + iter %% Neps
        runMixPar[ind,,] <- mixEstPar
        runOrder[ind] <- iter

        iter <- iter + 1
    }
    if(iter == Niter.max)
        warning("Maximum number of iterations reached.")

    mixEst <- mixEst[,order(mixEst[1,], decreasing=TRUE),drop=FALSE]
    colnames(mixEst) <- paste("comp", seq(Nc), sep="")
    dlink(mixEst) <- identity_dlink
    class(mixEst) <- c("EM", "EMgmm", "gammaMix", "mix")

    ## give further details
    attr(mixEst, "df") <- Nc-1 + 2*Nc
    attr(mixEst, "nobs") <- N
    attr(mixEst, "lli") <- lliCur

    attr(mixEst, "Nc") <- Nc

    attr(mixEst, "tol") <- tol
    attr(mixEst, "traceLli") <- traceLli
    attr(mixEst, "traceMix") <- lapply(traceMix, function(x) {class(x) <- c("gammaMix", "mix"); x})
    attr(mixEst, "x") <- x

    mixEst
}


#' @export
print.EMgmm <- function(x, ...) {
    cat("EM for Gamma Mixture Model\nLog-Likelihood = ", logLik(x), "\n\n",sep="")
    NextMethod()
}

