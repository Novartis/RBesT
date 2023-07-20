
## EM for MNMM (Multi-variate Normal Mixture Model) with Nc components
## if init is not speciefied, then knn is used to initialize means,
## cluster weights and covariance matrices (taken from the knn
## determined clusters)

EM_mnmm <- function(X, Nc, mix_init, Ninit=50, verbose=FALSE, Niter.max=500, tol, Neps, eps=c(w=0.005,m=0.005,s=0.005)) {
    ## in case X is no matrix, interpret it as uni-variate case
    if(!is.matrix(X))
        X <- matrix(X,ncol=1)

    N <- dim(X)[1]
    Nd <- dim(X)[2]
    assert_that(N+Nc >= Ninit)

    ## initialize normal EM using a Student-t EM which is very robust
    ## against outliers
    if(missing(mix_init)) {
        mix_init <- EM_msmm(X, Nc, Ninit=Ninit, verbose=verbose, Niter.max=round(Niter.max/2), tol=0.1)
    }
    pEst <- mix_init$p
    muEst <- mix_init$center
    covEst <- mix_init$cov

    ## take current estimates and transform to scale on which
    ## convergence is assessed
    est2par <- function(p, mu, cov) {
        est <- rbind(logit(p), matrix(sapply(1:Nc, function(i) mv2vec(mu[i,], cov[i,,])), ncol=Nc))
        est[(1+Nd+1):(1+2*Nd),] <- log(est[(1+Nd+1):(1+2*Nd),])
        est
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

    ## degrees of freedom
    ## covariance matrix df per component
    cov.df <- (Nd-1)*Nd/2 + Nd
    df <- Nc  * (Nd + cov.df) + Nc-1
    df.comp <- cov.df + Nd + 1

    ## expand eps according to the dimensionality
    eps <- c(eps[1], rep(eps[2], Nd), rep(eps[3], Nd), rep(eps[2], (Nd-1)*Nd/2))

    iter <- 0
    logN <- log(N)
    traceMix <- list()
    traceLli <- c()
    Dlli <- Inf
    runMixPar <- array(-Inf, dim=c(Neps,df.comp,Nc))
    runOrder <- 0:(Neps-1)
    if(Nc == 1) Npar <- df else Npar <- df + 1

    ## initialize component and element wise log-likelihood matrix
    lli <- array(-20, dim=c(N,Nc))

    if(verbose) {
        message("EM multi-variate normal with Nc =", Nc, ":\n")
    }

    while(iter < Niter.max) {
        ## calculate responsabilities from the likelihood terms;
        ## calculations are done in log-space to avoid numerical
        ## difficulties if some points are far away from some
        ## component and hence recieve very low density
        for(i in seq(Nc)) {
            lli[,i] <- log(pEst[i]) + dmvnorm(X, muEst[i,], as.matrix(covEst[i,,]), log=TRUE, checkSymmetry=FALSE)
        }
        ## ensure that the log-likelihood does not go out of numerical
        ## reasonable bounds
        lli <- apply(lli, 2, pmax, -30)
        ##lnresp <- apply(lli, 1, log_sum_exp)
        lnresp <- matrixStats::rowLogSumExps(lli)
        ## the log-likelihood is then given by the sum of lresp
        lliCur <- sum(lnresp)
        ## record current state
        traceMix <- c(traceMix, list(list(p=pEst, mean=muEst, sigma=covEst)))
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
        if(verbose)
            message("Iteration ", iter, ": log-likelihood = ", lliCur, "; Dlli = ", Dlli, "; converged = ", eps.converged, " / ", Npar, "\n", sep="")
        if(checkTol & Dlli < tol) {
            break
        }
        if(iter >= Neps & checkEps & eps.converged == Npar) {
            break
        }
        ## ... and the responseability matrix follows from this by
        ## appropiate normalization.
        lresp <- sweep(lli, 1, lnresp, "-")
        ##resp <- exp(lresp)

        ## mean probability to be in a specific mixture component -> updates
        ## pEst
        ##lzSum <- apply(lresp, 2, log_sum_exp)
        lzSum <- colLogSumExps(lresp)
        ##zSum <- exp(lzSum)
        pEst <- exp(lzSum - logN)

        ## make sure it is scale to exactly 1 which may not happen due
        ## to small rounding issues
        pEst <- pEst / sum(pEst)

        ## now obtain new estimates for each component of the mixtures
        ## of their mu vector and covariance matrices
        for(i in seq(Nc)) {
            upd <- cov.wt(X, exp(lresp[,i] - lzSum[i]), method="ML")
            ##upd <- cov.wt(X, resp[,i]/zSum[i], method="ML")
            muEst[i,] <- upd$center
            covEst[i,,] <- upd$cov
            ## ensure that diagonal stays non-zero
            for(j in 1:Nd)
                covEst[i,j,j] <- max(covEst[i,j,j], .Machine$double.eps)
        }

        ind <- 1 + iter %% Neps
        runMixPar[ind,,] <- est2par(pEst, muEst, covEst)
        runOrder[ind] <- iter

        iter <- iter + 1
    }
    if(iter+1 == Niter.max)
        warning("Maximum number of iterations reached.")

    ## sort by largest weight
    o <- order(pEst, decreasing=TRUE)
    pEst <- pEst[o]
    muEst <- muEst[o,,drop=FALSE]
    covEst <- covEst[o,,,drop=FALSE]

    ##if(Nd != 1) {
    ##    rhoEst <- array(apply(covEst, 1, cov2cor), c(Nd,Nd,Nc))
    ##    rhoEst <- apply(rhoEst, 3, function(x) x[lower.tri(x)])
    ##    tauEst <- sqrt(t(apply(covEst, 1, diag)))
    ##} else {
    ##    rhoEst <- NULL
    ##    tauEst <- sqrt(as.vector(covEst))
    ##}

    ##mixEst <- list(p=pEst, mean=muEst, sigma=covEst)

    mixEst <- do.call(mixmvnorm, lapply(1:Nc, function(i) c(pEst[i], muEst[i,,drop=FALSE], matrix(covEst[i,,], Nd, Nd))))
    
    ## give further details
    attr(mixEst, "df") <- df
    attr(mixEst, "nobs") <- N
    attr(mixEst, "lli") <- lliCur

    attr(mixEst, "Nc") <- Nc

    attr(mixEst, "tol") <- tol
    attr(mixEst, "traceLli") <- traceLli
    attr(mixEst, "traceMix") <- traceMix
    attr(mixEst, "x") <- X

    class(mixEst) <- c("EM", "EMmvnmm", "mvnormMix", "mix")

    mixEst
}

## uni-variate case
EM_nmm <- function(X, Nc, mix_init, verbose=FALSE, Niter.max=500, tol, Neps, eps=c(w=0.005,m=0.005,s=0.005)) {
    if(is.matrix(X)) {
        assert_matrix(X, any.missing=FALSE, ncols=1)
    }
    mixEst <- EM_mnmm(X=X, Nc=Nc, mix_init=mix_init, verbose=verbose, Niter.max=Niter.max, tol=tol, Neps=Neps, eps=eps)
    rownames(mixEst) <- c("w", "m", "s")
    class(mixEst) <- c("EM", "EMnmm", "normMix", "mix")
    likelihood(mixEst) <- "normal"
    mixEst
}

#' @export
print.EMnmm <- function(x, ...) {
    cat("EM for Normal Mixture Model\nLog-Likelihood = ", logLik(x), "\n\n",sep="")
    NextMethod()
}

#' @export
print.EMmvnmm <- function(x, ...) {
    cat("EM for Multivariate Normal Mixture Model\nLog-Likelihood = ", logLik(x), "\n\n",sep="")
    NextMethod()
}

## given a vector of means and a covariance matrix for a multi-variate
## normal (and optionally a vector of df), return a vector of
## coefficients with a deterministic mapping
mv2vec <- function(mean, sigma, df, label=TRUE) {
    Nd <- length(mean)
    if(Nd != 1) {
        rho <- cov2cor(sigma)[lower.tri(sigma)]
        tau <- sqrt(diag(sigma))
    } else {
        rho <- NULL
        tau <- sqrt(sigma)
    }
    if(missing(df)) df <- NULL
    res <- c(mean, tau, rho, df)
    if(label) {
        tauN <- paste("sd", 1:Nd, sep="")
        cols <- names(mean)
        if(is.null(cols))
            cols <- paste("mu", 1:Nd, sep="")
        if(Nd == 1)
            corNames <- NULL
        else if(length(cols) == 2)
            corNames <- "cor"
        else {
            corNames <- outer(cols, cols, paste, sep="_")
            corNames <- paste("cor", corNames[lower.tri(corNames)], sep=".")
        }
        if(!is.null(df))
            dfNames <- paste("df", 1:Nd, sep="")
        else
            dfNames <- NULL
        names(res) <- c(cols, tauN, corNames, dfNames)
    }
    return(res)
}

## vec2mv <- function(vec) {
##     ## calculate dimension from the number of parameters
##     Nd <- -3/2 + sqrt( 9/4 + 2*length(vec) )
##     mean <- vec[1:Nd]
##     Tau <- diag(vec[(Nd+1):(2*Nd)])
##     sigma <- diag(Nd)/2
##     L <- lower.tri(sigma)
##     E <- sum(L)
##     sigma[L] <- vec[(2*Nd+1):(2*Nd+E)]
##     sigma <- sigma + t(sigma)
##     sigma <- Tau %*% sigma %*% Tau
##     list(mean=mean, sigma=sigma)
## }
##
## ## extracts results in a flattened form
## extractMVNmix <- function(emFit) {
##     pMix <- emFit$p
##     cols <- colnames(emFit$center)
##     if(is.null(cols))
##         cols <- paste("X", 1:ncol(emFit$center), sep="")
##     if(length(cols) == 1)
##         corNames <- NULL
##     else if(length(cols) == 2)
##         corNames <- "cor"
##     else {
##         corNames <- outer(cols, cols, paste, sep="_")
##         corNames <- paste("cor", corNames[lower.tri(corNames)], sep=".")
##     }
##     MAPmix <- do.call(cbind, list(emFit$center, emFit$tau, emFit$rho))
##     colnames(MAPmix) = c(paste("mean", cols, sep="."), paste("sd", cols, sep="."), corNames)
##     MAPmix <- as.data.frame(t(MAPmix))
##     colnames(MAPmix) <- paste("Mix", 1:length(pMix), sep="")
##     names(pMix) <- names(MAPmix)
##     list(mvn=MAPmix, pMix=pMix)
## }
##
## ## utility functions for multi-variate normal mixtures
## rmvnormMix <- function(n, p, m, sigma){
##   ind <- sample.int(length(p), n, replace = TRUE, prob = p)
##   d <- nrow(m)
##   samp <- matrix(0, nrow=n, ncol=d)
##   for(i in seq(n))
##       samp[i,] <- rmvnorm(1, m[ind[i],], as.matrix(sigma[ind[i],,]))
##   samp
## }
##
## dmvnormMix <- function(x, p, m, sigma) {
##     nc <- length(p)
##     ## logic is to replicate the original data vector such that each
##     ## item appears nc times which allows easy vectorized calls to
##     ## dnorm. Then we cast the result into a matrix with as many rows
##     ## as nc components which we sum together with a fast colSums call.
##     dens <- rep.int(0, nrow(x))
##     for(i in seq_along(p))
##         dens <- dens + p[i] * dmvnorm(x, m[i,], sigma[i,,])
##     dens
## }
##
## ## utility function to plot BVN mixtures
## bicontourNMM <- function(X, bvn, title="", ng=50) {
##     ## some pretty colors
##     ##library(colorspace)
##     k <- 15
##     ##my.cols <- diverge_hcl(k, c = c(100, 0), l = c(50, 90), power = 1.3)
##     my.cols <- c("#4A6FE3", "#6D84E1", "#8898E1", "#A0AAE2", "#B5BBE3", "#C7CBE3", "#D7D9E3",
##                  "#E2E2E2", "#E4D6D8", "#E6C4CA", "#E6AFB9", "#E498A7", "#E07E93", "#DB627F",
##                  "#D33F6A")
##
##     r <- apply(X, 2, range)
##     xg <- seq(r[1,1],r[2,1], length=ng)
##     yg <- seq(r[1,2],r[2,2], length=ng)
##     Z <- outer(xg,yg, function(x,y) dmvnormMix(cbind(x,y), bvn$p, bvn$center, bvn$cov))
##
##     Nc <- length(bvn$p)
##
##     plot(X, pch=19, cex=.2)
##     contour(xg, yg, Z, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE, lwd=2)
##     abline(h=bvn$center[,2], v=bvn$center[,1], lwd=1, col=1:Nc)
##     title(title)
##     legend("bottomleft", legend=paste("Mix", 1:Nc, " ", round(100*bvn$p,1), "%", sep=""), text.col=1:Nc)
## }
##
