#' Effective Sample Size for a Conjugate Prior
#'
#' Calculates the Effective Sample Size (ESS) for a mixture prior. The
#' ESS indicates how many experimental units the prior is roughly
#' equivalent to.
#'
#' @param mix Prior (mixture of conjugate distributions).
#' @param method Selects the used method. Can be either \code{elir}
#'     (default), \code{moment} or \code{morita}.
#' @param s For \code{morita} method large constant to ensure that the
#'     prior scaled by this value is vague (default 100); see Morita
#'     et al. (2008) for details.
#' @param eps Probability mass left out from the numerical integration
#'     of the expected information for the Poisson-Gamma case of
#'     Morita method (defaults to 1E-4).
#' @param sigma reference scale.
#' @param ... Optional arguments applicable to specific methods.
#'
#' @details The ESS is calculated using either the expected local
#'     information ratio (elir) \emph{Neuenschwander et
#'     al. (submitted)}, the moments approach or the method by
#'     \emph{Morita et al. (2008)}.
#'
#' The elir approach is the only ESS which fulfills predictive
#' consistency. The predictive consistency of the ESS requires that
#' the ESS of a prior is the same as averaging the posterior ESS after
#' a fixed amount of events over the prior predictive distribution
#' from which the number of forward simulated events is
#' subtracted. The elir approach results in ESS estimates which are
#' neither conservative nor liberal whereas the moments method yields
#' conservative and the morita method liberal results. See the example
#' section for a demonstration of predictive consistency.
#'
#' For the moments method the mean and standard deviation of the
#' mixture are calculated and then approximated by the conjugate
#' distribution with the same mean and standard deviation. For
#' conjugate distributions, the ESS is well defined. See the examples
#' for a step-wise calculation in the beta mixture case.
#'
#' The Morita method used here evaluates the mixture prior at the mode
#' instead of the mean as proposed originally by Morita. The method
#' may lead to very optimistic ESS values, especially if the mixture
#' contains many components. The calculation of the Morita approach
#' here follows the approach presented in Neuenschwander B. et all
#' (2019) which avoids the need for a minimization and does not
#' restrict the ESS to be an integer.
#'
#' @return Returns the ESS of the prior as floating point number.
#'
#' @template conjugate_pairs
#'
#' @references Morita S, Thall PF, Mueller P.
#' Determining the effective sample size of a parametric prior.
#' \emph{Biometrics} 2008;64(2):595-602.
#'
#' @references Neuenschwander B, Weber S, Schmidli H, O'Hagen A.
#' Predictively Consistent Prior Effective Sample Sizes.
#' \emph{pre-print} 2019; arXiv:1907.04185
#'
#' @examples
#' # Conjugate Beta example
#' a <- 5
#' b <- 15
#' prior <- mixbeta(c(1, a, b))
#'
#' ess(prior)
#' (a+b)
#'
#' # Beta mixture example
#' bmix <- mixbeta(rob=c(0.2, 1, 1), inf=c(0.8, 10, 2))
#'
#' ess(bmix, "elir")
#'
#' ess(bmix, "moment")
#' # moments method is equivalent to
#' # first calculate moments
#' bmix_sum <- summary(bmix)
#' # then calculate a and b of a matching beta
#' ab_matched <- ms2beta(bmix_sum["mean"], bmix_sum["sd"])
#' # finally take the sum of a and b which are equivalent
#' # to number of responders/non-responders respectivley
#' round(sum(ab_matched))
#'
#' ess(bmix, method="morita")
#'
#' # Predictive consistency of elir
#' \donttest{
#' n_forward <- 1E2
#' bmixPred <- preddist(bmix, n=n_forward)
#' pred_samp <- rmix(bmixPred, 1E3)
#' pred_ess <- sapply(pred_samp, function(r) ess(postmix(bmix, r=r, n=n_forward), "elir") )
#' ess(bmix, "elir")
#' mean(pred_ess) - n_forward
#' }
#'
#' # Normal mixture example
#' nmix <- mixnorm(rob=c(0.5, 0, 2), inf=c(0.5, 3, 4), sigma=10)
#'
#' ess(nmix, "elir")
#'
#' ess(nmix, "moment")
#'
#' ## the reference scale determines the ESS
#' sigma(nmix) <- 20
#' ess(nmix)
#'
#' # Gamma mixture example
#' gmix <- mixgamma(rob=c(0.3, 20, 4), inf=c(0.7, 50, 10))
#'
#' ess(gmix) ## interpreted as appropriate for a Poisson likelihood (default)
#'
#' likelihood(gmix) <- "exp"
#' ess(gmix) ## interpreted as appropriate for an exponential likelihood
#'
#'
#' @export
ess <- function(mix, method=c("elir", "moment", "morita"), ...) UseMethod("ess")
#' @export
ess.default <- function(mix, method=c("elir", "moment", "morita"), ...) stop("Unknown density")


calc_loc <- function(mix, loc=c("mode", "median", "mean")) {
    loc <- match.arg(loc)
    if(loc == "mode") {
        tol <- .Machine$double.eps^0.25
        locEst <- mixmode(mix)

        if(length(attr(locEst, "modes")) > 1) {
            warning("Detected multiple modes.\nThe ESS is determined for the largest mode, but ESS concept is ill-defined for multi-modal distributions.")
        } else {
            attr(locEst, "modes") <- NULL
        }
    }
    if(loc == "median") {
        locEst <- qmix(mix, 0.5)
    }
    if(loc == "mean") {
        locEst <- summary(mix, NULL)["mean"]
    }
    names(locEst) <- NULL

    return(unname(locEst))
}
## function to calculate mixture info of arbitrary density; input
## needed is the density function, gradient and hessian of the log
## density with respect to x (data)
mixInfo <- function(mix, x, dens, gradl, hessl) {
    p <- mix[1,]
    a <- mix[2,]
    b <- mix[3,]
    lp <- log(p)
    ldensComp <- dens(x, a, b, log=TRUE)
    ldensMix <- matrixStats::logSumExp(lp + ldensComp)
    lwdensComp <- lp + ldensComp - ldensMix
    dgl <- gradl(x,a,b)
    dhl <- (hessl(x,a,b) + dgl^2)
    ## attempt numerically more stable log calculations if possible,
    ## i.e. if all sings are the same
    if(all(dgl < 0) || all(dgl > 0)) {
        gsum <- exp(2*matrixStats::logSumExp(lwdensComp + log(abs(dgl))))
    } else {
        gsum <- (sum(exp(lwdensComp)*dgl))^2
    }
    if(all(dhl < 0) || all(dhl > 0)) {
        hsum <- sign(dhl[1]) * exp(matrixStats::logSumExp(lwdensComp + log(abs(dhl))))
    } else {
        hsum <- (sum(exp(lwdensComp)*dhl))
    }
    gsum - hsum
}

## local information ratio (which we integrate over the prior)
lir <- function(mix, info, fisher_inverse) {
    fn  <- function(x) {
        info(mix, x) * fisher_inverse(x)
    }
    Vectorize(fn)
}

weighted_lir <- function(mix, info, fisher_inverse) {
    fn  <- function(x) {
        dmix(mix, x) * info(mix, x) * fisher_inverse(x)
    }
    Vectorize(fn)
}

## not used ATM as there have been numerical issues
weighted_lir_link <- function(mix, info, fisher_inverse, link) {
    dlink(mix) <- link_map[[link]]
    fn  <- function(x) {
        x_orig  <- mixinvlink(mix, x)
        dmix(mix, x) * info(mix, x_orig) * fisher_inverse(x_orig)
    }
    Vectorize(fn)
}

## function to calculate the gradient of the log mixture
## mixLogGrad <- function(mix, x, dens, gradl) {
##     p <- mix[1,]
##     a <- mix[2,]
##     b <- mix[3,]
##     densMix <- dmix(mix,x)
##     densComp <- dens(x, a, b)
##     dgl <- gradl(x,a,b)
##     sum(p*densComp*dgl)/densMix
## }


# prior effective sample size ESS for Beta-mixture priors
# based on
# Morita, Thall, Mueller (MTM) 2008 Biometrics
# only difference: evaluated at mode of prior rather than at mean; and the flattened prior are derived with respect to the scale of 1 instead of being relative to the input scale
#' @describeIn ess ESS for beta mixtures.
#' @export
ess.betaMix <- function(mix, method=c("elir", "moment", "morita"), ..., s=100) {

    method <- match.arg(method)

    if(method == "elir") {
        if(!test_numeric(mix[2,], lower=1, finite=TRUE, any.missing=FALSE) ||
           !test_numeric(mix[3,], lower=1, finite=TRUE, any.missing=FALSE)) {
            stop("At least one parameter of the beta mixtures is less than 1.\n",
                 "This leads to an ill-defined elir ess since the defining integral diverges.\n",
                 "Consider constraining all parameters to be greater than 1 (use constrain_gt1=TRUE argument for EM fitting functions).")
        }
        elir <- integrate_density(lir(mix, betaMixInfo, bernoulliFisherInfo_inverse), mix)
        return(elir)
    }

    ## simple and conservative moment matching
    if(method == "moment") {
        smix <- summary(mix)
        res <- sum(ms2beta(smix["mean"], smix["sd"]))
        names(res) <- NULL
        return( res )
    }

    alphaP <- mix[2,]
    betaP <- mix[3,]

    locEst <- calc_loc(mix, "mode")

    deriv2.prior <- betaMixInfo(mix, locEst)

    ESSmax <- ceiling(sum(alphaP+betaP)) * 2

    ## alpha and beta of "flattened" priors
    alphaP0 <- locEst / s
    betaP0 <- (1-locEst) / s

    info.prior0 <- betaInfo(locEst, alphaP0, betaP0)

    ## MTM paper would follow this:
    ##priorN <- sum(mix[1,] * (alphaP + betaP))
    ##alphaP0 <- mode * priorN / s
    ##betaP0 <- (1-mode) * priorN / s

    ## we warn if any of the mixture components has a scale (n) which
    ## is less than 10/s such that the
    if(any(rowSums(mix[2:3,,drop=FALSE]) < 10/s )) {
        warning("Some of the mixture components have a scale which is large compared to the rescaling factor s. Consider increasing s.")
    }

    ## expected 2nd derivative at mode
    ed2p <- function(m) {
        yn <- seq(0,m)
        ## negative 2nd log-derivative at mode
        info <- betaInfo(locEst, alphaP0 + yn, betaP0 + m - yn)
        ## prior predictive
        sum(info * dmix(preddist(mix,n=m), yn) )
    }
    ## function to search for change of sign
    ed2pDiff <- function(m) {
        deriv2.prior - ed2p(m)
    }

    pd0 <- dmix(preddist(mix, n=1), 0)
    Einfo <- binomialInfo(0,locEst,1) * pd0 + binomialInfo(1,locEst,1) * (1-pd0)

    ##return(unname(uniroot_int(ed2pDiff, c(0,ESSmax))))
    ## Eq. 9 of Neuenschwander et al. (2019)
    return( unname((deriv2.prior - info.prior0) / Einfo ) )
}

## derivative of a single log-beta
betaLogGrad <- function(x,a,b) {
    lxm1 <- log1p(-x)
    lx <- log(x)
    - (b-1) * exp(- lxm1) + (a-1)*exp(- lx)
}

## second derivative of a single log-beta
betaLogHess <- function(x,a,b) {
    lxm1 <- log1p(-x)
    lx <- log(x)
    - (b-1) * exp(-2 * lxm1) - (a-1)*exp(-2 * lx)
}

betaMixInfo <- function(mix,x) {
    mixInfo(mix, x, dbeta, betaLogGrad, betaLogHess)
}

## info metric for a single beta, i.e. negative second derivative of log beta
betaInfo <- function(x,a,b) {
    -betaLogHess(x,a,b)
}

## 1/i_F(x): The inverse of the fisher information for a Bernoulli
## experiment (binomial with n=1)
bernoulliFisherInfo_inverse <- function(x) {
    x - x^2
}

## info metric for a binomial, second derivative wrt to theta of the
## log binomial
binomialInfo <- function(r,theta,n) {
    r / theta^2 + (n-r)/(1-theta)^2
}


#' @describeIn ess ESS for gamma mixtures.
#' @export
ess.gammaMix <- function(mix, method=c("elir", "moment", "morita"), ..., s=100, eps=1E-4) {

    method <- match.arg(method)
    lik <- likelihood(mix)

    if(method == "elir") {
        if(lik == "poisson")
            return(integrate_density(lir(mix, gammaMixInfo, poissonFisherInfo_inverse), mix))
        if(lik == "exp")
            return(integrate_density(lir(mix, gammaMixInfo, expFisherInfo_inverse), mix))
    }

    ## simple and conservative moment matching
    if(method == "moment") {
        smix <- summary(mix)
        coef <- ms2gamma(smix["mean"], smix["sd"])
        names(coef) <- NULL
        if(lik == "poisson")
            return(unname(coef[2]))
        if(lik == "exp")
            return(unname(coef[1]))
        stop("Unkown likelihood")
    }

    ## Morita method
    locEst <- calc_loc(mix, "mode")

    deriv2.prior <- gammaMixInfo(mix, locEst)

    if(lik == "poisson") {
        meanPrior <- summary(mix)["mean"]
        names(meanPrior) <- NULL
        priorN <- mix[3,,drop=FALSE]

        ## function to search for change of sign
        ed2pDiff <- function(m) {
            deriv2.prior - ( gammaInfo(locEst, locEst/s + m * meanPrior, 1/s))
        }

        ESSmax <- ceiling(sum(mix[3,])) * 2

        info.prior0  <- gammaInfo(locEst, locEst/s, 1/s)

        ## E_Y1 ( i_F ) using numerical integration
        pred_pmf <- preddist(mix, n=1)
        lim <- qmix(pred_pmf, c(eps/2, 1-eps/2))
        y1 <- seq(lim[1], lim[2])
        Einfo <- sum(dmix(pred_pmf, y1) * poissonInfo(y1, locEst))
    }

    if(lik == "exp") {
        priorN <- mix[2,,drop=FALSE]
        ## function to search for change of sign
        ed2pDiff <- function(m) {
            deriv2.prior - gammaInfo(locEst, 1/s + m, 1/(s*locEst))
        }

        ESSmax <- ceiling(sum(mix[2,])) * 2

        info.prior0  <- gammaInfo(locEst, 1/s, 1/(s*locEst))

        ## E_Y1 ( i_F ) (the i_F does not depend on the data)
        Einfo  <- expInfo(1,locEst)
    }

    if(any(priorN < 10/s )) {
        warning("Some of the mixture components have a scale which is large compared to the rescaling factor s. Consider increasing s.")
    }

    return(unname( (deriv2.prior - info.prior0)/Einfo ) )
}

## respective functions for the gamma distribution
gammaLogGrad <- function(x,a,b) {
    (a-1)/x - b
}
gammaLogHess <- function(x,a,b) {
    -(a-1)/x^2
}
gammaInfo <- function(x,a,b) {
    -gammaLogHess(x,a,b)
}
gammaMixInfo <- function(mix,x) {
    mixInfo(mix, x, dgamma, gammaLogGrad, gammaLogHess)
}
poissonFisherInfo_inverse <- function(x) {
    x
}
expFisherInfo_inverse <- function(x) {
    x^2
}
poissonInfo <- function(y,theta) {
    y/theta^2
}
expInfo <- function(y,theta) {
    1/theta^2
}


#' @describeIn ess ESS for normal mixtures.
#' @export
ess.normMix <- function(mix, method=c("elir", "moment", "morita"), ..., sigma, s=100) {

    method <- match.arg(method)

    if(missing(sigma)) {
        sigma <- RBesT::sigma(mix)
        message("Using default prior reference scale ", sigma)
    }
    assert_number(sigma, lower=0)
    tauSq <- sigma^2

    ## note: sigma is reassigned the sd's of each mixture component
    mu <- mix[2,]
    sigma <- mix[3,]
    sigmaSq <- sigma^2

    if(method == "elir") {
        return(tauSq * integrate_density(lir(mix, normMixInfo, normStdFisherInfo_inverse), mix))
    }

    ## simple and conservative moment matching
    if(method == "moment") {
        smix <- summary(mix)
        res <- tauSq / smix["sd"]^2
        return( unname(res) )
    }

    locEst <- calc_loc(mix, "mode")

    deriv2.prior <- normMixInfo(mix, locEst)

    ESSmax <- ceiling(sum( (1-1/s) * tauSq/sigmaSq )) * 2

    ## "flattened" priors
    muP0 <- locEst
    sigmaP0Sq <- s * max(sigmaSq)

    ## difference of info at locEst and the expected 2nd derivative at
    ## locEst
    ed2pDiff <- function(m) {
        deriv2.prior - normInfo(locEst, muP0, sqrt(1/(m/tauSq + 1/sigmaP0Sq)))
    }

    info.prior0  <- normInfo(locEst, muP0, sqrt(sigmaP0Sq))

    Einfo  <- 1/tauSq

    return(unname( (deriv2.prior - info.prior0)/Einfo ))
}

## derivative of a single log-normal
normLogGrad <- function(x,mu,sigma) {
    -1 * (x-mu)/(sigma)^2
}

## second derivative of a single log-normal
normLogHess <- function(x,mu,sigma) {
    -1/sigma^2
}

normMixInfo <- function(mix,x) {
    mixInfo(mix, x, dnorm, normLogGrad, normLogHess)
}

## info metric for a single norm, i.e. negative second derivative of log norm
normInfo <- function(x,mean,sigma) {
    -normLogHess(x,mean,sigma)
}

## Fisher info for normal sampling sd with known unit variance
normStdFisherInfo_inverse <- function(x) {
    1.0
}
