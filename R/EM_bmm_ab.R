
## EM for Beta Mixture Models (BMM) with Nc components

EM_bmm_ab <- function(x, Nc, mix_init, Ninit=50, verbose=FALSE, Niter.max=500, tol, Neps, eps=c(w=0.005,a=0.005,b=0.005), constrain_gt1=TRUE) {
    N <- length(x)
    assert_that(N+Nc >= Ninit)

    assert_logical(constrain_gt1, any.missing=FALSE, len=1)

    ## check data for 0 and 1 values which are problematic, but may be
    ## valid, depending on a and b. Moving these to eps or 1-eps
    ## ensures proper handling during fit.
    x0 <- x==0
    if(any(x0)) {
        message("Detected ", sum(x0), " value(s) which are exactly 0.\nTo avoid numerical issues during EM such values are moved to smallest eps on machine.")
        x[x0] <- .Machine$double.eps
    }
    x1 <- x==1
    if(any(x1)) {
        message("Detected ", sum(x1), " value(s) which are exactly 1.\nTo avoid numerical issues during EM such values are moved to one minus smallest eps on machine.")
        x[x1] <- 1-.Machine$double.eps
    }

    ## temporaries needed during EM
    Lx <- matrix(log(x), ncol=Nc, nrow=N)
    LxC <- matrix(log1p(-x), ncol=Nc, nrow=N)

    xRep <- rep(x, each=Nc)

    ## initialize randomly using KNN
    if(missing(mix_init)) {
        ##abmEst <- matrix(1+rlnorm(Nc*3, 0, log(5)/1.96), nrow=Nc)
        ##abmEst[,1] <- 1/Nc
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
        nInit <- pmax(muInit*(1-muInit)/varInit - 1, 1, na.rm=TRUE)
        ## place the component which recieved the least weight at the
        ## data center with roughly the variance of the sample
        cmin <- which.min(KNN$p)
        muInit[cmin] <- sum(KNN$p * KNN$center)
        ## muInit[cmin] <- mean(x) ## could be considered here
        nInit[cmin] <- pmax(muInit[cmin]*(1-muInit[cmin])/var(x) - 1, 1, na.rm=TRUE)
        ##Nmax <- max(2, max(nInit))
        ## ensure n is positive for each cluster; if this is not the
        ## case, sample uniformly from the range of n we have
        ##Nneg <- nInit <= .Machine$double.eps
        ##Nsmall <- nInit <= 0.5
        ##if(any(Nsmall))
        ##    nInit[Nsmall] <- runif(sum(Nsmall), 0.5, Nmax)
        ##nInitR <- 0.5 + rlnorm(Nc, log(nInit), log(5)/1.96)
        mixEst <- rbind(KNN$p, nInit*muInit, nInit*(1-muInit))
        dlink(mixEst) <- identity_dlink
        rownames(mixEst) <- c("w", "a", "b")
    } else {
        mixEst <- mix_init
    }

    ## mixEst parametrization during fitting
    mixEstPar <- mixEst
    mixEstPar[1,] <- logit(mixEst[1,,drop=FALSE])
    mixEstPar[2,] <- log(mixEst[2,])
    mixEstPar[3,] <- log(mixEst[3,])
    rownames(mixEstPar) <-  c("w", "la", "lb")

    ## constrain to a>=1 & b>=1 if requested... thus we subtract 1
    if(constrain_gt1) {
        mixEstPar[2,]  <- log(pmax(mixEst[2,] - 1, rep(1E-8, times=Nc)))
        mixEstPar[3,]  <- log(pmax(mixEst[3,] - 1, rep(1E-8, times=Nc)))

        mixEst[2,]  <- 1 + exp(mixEstPar[2,])
        mixEst[3,]  <- 1 + exp(mixEstPar[3,])
    }

    if(verbose) {
        message("EM for beta mixture model.\n")
        message("Initial estimates:\n")
        print(mixEst)
    }

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
    bmm_ml <- function(c1,c2) {
        function(par) {
            ab <- exp(par)
            if(constrain_gt1)
                ab <- 1 + ab
            s <- digamma(sum(ab))
            eq1 <- digamma(ab[1]) - s
            eq2 <- digamma(ab[2]) - s
            (eq1 - c1)^2 + (eq2 - c2)^2
        }
    }

    bmm_ml_grad  <- function(c1,c2) {
        function(par) {
            ab <- exp(par)
            if(constrain_gt1)
                ab <- 1 + ab
            n <- sum(ab)
            s <- digamma(n)
            eq1 <- digamma(ab[1]) - s
            eq2 <- digamma(ab[2]) - s
            sqTerm1 <- (eq1 - c1)
            sqTerm2 <- (eq2 - c2)

            trig_n  <- trigamma(n)

            grad1 <- 2 * sqTerm1 * (trigamma(ab[1]) * ab[1] - trig_n * ab[1] ) -
                2 * sqTerm2 * (trig_n * ab[1])

            grad2 <- -2 * sqTerm1 * (trig_n * ab[2]) +
                2 * sqTerm2 * (trigamma(ab[2]) * ab[2] - trig_n * ab[2])

            c(grad1, grad2)
        }
    }

    while(iter < Niter.max) {
        ## calculate responsabilities from the likelihood terms;
        ## calculations are done in log-space to avoid numerical
        ## difficulties if some points are far away from some
        ## component and hence recieve very low density
        ##lli <- t(matrix(log(mixEst[1,]) + dbeta(xRep, mixEst[2,], mixEst[3,], log=TRUE), nrow=Nc))

        ## Beta: Gamma(a + b) / (Gamma(a) * Gamma(b)) * x^(a-1) * (1-x)^(b-1)
        w <- mixEst[1,]
        a <- mixEst[2,]
        b <- mixEst[3,]
        ##lli <- sweep( sweep(Lx, 2, a - 1, "*", check.margin=FALSE) + sweep(LxC, 2, b - 1, "*", check.margin=FALSE), 2, log(w) + lgamma(a + b) - lgamma(a) - lgamma(b), "+", check.margin=FALSE)
        lli <- sweep( sweep(Lx, 2, a - 1, "*", check.margin=FALSE) + sweep(LxC, 2, b - 1, "*", check.margin=FALSE), 2, log(w) - lbeta(a, b), "+", check.margin=FALSE)
        ##lli <- t(matrix(log(mixEst[1,]) + dbeta(xRep, mixEst[2,], mixEst[3,], log=TRUE), nrow=Nc))

        ## ensure that the log-likelihood does not go out of numerical
        ## reasonable bounds
        lli <- apply(lli, 2, pmax, -30)

        ##lnresp <- apply(lli, 1, log_sum_exp)
        lnresp  <- matrixStats::rowLogSumExps(lli)
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
        lresp <- sweep(lli, 1, lnresp, "-", FALSE)
        ##resp <- exp(lresp)

        ## mean probability to be in a specific mixture component -> updates
        ## abmEst first colum
        ##lzSum <- apply(lresp, 2, log_sum_exp)
        lzSum <- matrixStats::colLogSumExps(lresp)
        ##zSum <- exp(lzSum)
        mixEst[1,] <- exp(lzSum - logN)

        ##c1 <- colSums(Lx * resp)/zSum
        ##c2 <- colSums(LxC * resp)/zSum

        resp_zscaled  <- exp(sweep(lresp, 2, lzSum, "-", FALSE))
        c1 <- matrixStats::colSums2(Lx * resp_zscaled)
        c2 <- matrixStats::colSums2(LxC * resp_zscaled)

        ## now solve for new alpha and beta estimates jointly for each
        ## component
        for(i in 1:Nc) {
            if(constrain_gt1) {
                theta <- c(log(mixEst[2:3,i] - 1))
            } else {
                theta <- c(log(mixEst[2:3,i]))
            }
            ##Lest <- optim(theta, bmm_ml(c1[i], c2[i]), gr=bmm_ml_grad(c1[i], c2[i]), method="BFGS", control=list(maxit=500))
            ## Default would be Nelder-Mead
            Lest <- optim(theta, bmm_ml(c1[i], c2[i]))
            if(Lest$convergence != 0 & Lest$value > 1E-4) {
                warning("Warning: Component", i, "in iteration", iter, "had convergence problems!")
            }
            if(constrain_gt1) {
                mixEst[2:3,i] <- 1+pmax(exp(Lest$par), c(1E-8, 1E-8))
            } else {
                mixEst[2:3,i] <- exp(Lest$par)
            }
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
    class(mixEst) <- c("EM", "EMbmm", "betaMix", "mix")

    ## give further details
    attr(mixEst, "df") <- Nc-1 + 2*Nc
    attr(mixEst, "nobs") <- N
    attr(mixEst, "lli") <- lliCur

    attr(mixEst, "Nc") <- Nc

    attr(mixEst, "tol") <- tol
    attr(mixEst, "traceLli") <- traceLli
    attr(mixEst, "traceMix") <- lapply(traceMix, function(x) {class(x) <- c("betaMix", "mix"); x})
    attr(mixEst, "x") <- x

    mixEst
}


#' @export
print.EMbmm <- function(x, ...) {
    cat("EM for Beta Mixture Model\nLog-Likelihood = ", logLik(x), "\n\n",sep="")
    NextMethod()
}

