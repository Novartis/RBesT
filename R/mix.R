#' @rdname mix
#' @name mix
#'
#' @title Mixture Distributions
#'
#' @description Density, cumulative distribution function, quantile
#' function and random number generation for supported mixture
#' distributions.  (d/p/q/r)mix are generic and work with any mixture
#' supported by BesT (see table below).
#'
#' @param mix mixture distribution object
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number
#' required
#' @param log,log.p logical; if \code{TRUE} (not default), probabilities \eqn{p} are given as \eqn{\log(p)}
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are \eqn{P[X\leq x]} otherwise, \eqn{P[X>x]}
#' @param rescale logical; if \code{TRUE}, mixture weights will be rescaled to sum to 1
#' @param ... components to subset given mixture.
#'
#' @details A mixture distribution is defined as a linear
#' superposition of \eqn{K} densities of the same distributional
#' class. The mixture distributions supported have the form
#'
#' \deqn{f(x,\mathbf{w},\mathbf{a},\mathbf{b}) = \sum_{k=1}^K w_k \, f_k(x,a_k,b_k).}{f(x,w,a,b) = \sum_{k=1}^K w_k * f(x,a_k,b_k).}
#'
#' The \eqn{w_k} are the mixing coefficients which must sum to
#' \eqn{1}. Moreover, each density \eqn{f} is assumed to be
#' parametrized by two parameters such that each component \eqn{k} is
#' defined by a triplet, \eqn{(w_k,a_k,b_k)}.
#'
#' Individual mixture components can be extracted using the \code{[[}
#' operator, see examples below.
#'
#' The supported densities are normal, beta and gamma which can be
#' instantiated with \code{\link{mixnorm}}, \code{\link{mixbeta}}, or
#' \code{\link{mixgamma}}, respectively. In addition, the respective
#' predictive distributions are supported. These can be obtained by
#' calling \code{\link{preddist}} which returns appropriate normal,
#' beta-binomial or Poisson-gamma mixtures.
#'
#' For convenience a \code{summary} function is defined for all
#' mixtures. It returns the mean, standard deviation and the requested
#' quantiles which can be specified with the argument \code{probs}.
#'
#' @return
#'
#' \code{dmix} gives the weighted sum of the densities of each
#' component.
#'
#' \code{pmix} calculates the distribution function by
#' evaluating the weighted sum of each components distribution
#' function.
#'
#' \code{qmix} returns the quantile for the given \code{p}
#' by using that the distribution function is monotonous and hence a
#' gradient based minimization scheme can be used to find the matching
#' quantile \code{q}.
#'
#' \code{rmix} generates a random sample of size
#' \code{n} by first sampling a latent component indicator in the
#' range \eqn{1..K} for each draw and then the function samples from
#' each component a random draw using the respective sampling
#' function. The \code{rnorm} function returns the random draws as
#' numerical vector with an additional attribute \code{ind} which
#' gives the sampled component indicator.
#'
#' @template conjugate_pairs
#'
#' @seealso \code{\link{plot.mix}}
#' @family mixdist
#'
#' @examples
#' ## a beta mixture
#' bm <- mixbeta(weak=c(0.2, 2, 10), inf=c(0.4, 10, 100), inf2=c(0.4, 30, 80))
#'
#' ## extract the two most informative components
#' bm[[c(2,3)]]
#' ## rescaling needed in order to plot
#' plot(bm[[c(2,3),rescale=TRUE]])
#'
#' summary(bm)
#'
#'
NULL


### DECLARATION
#' @export
#' @rdname mix
dmix <- function(mix, x, log=FALSE) UseMethod("dmix")
#' @export
#' @rdname mix
pmix <- function(mix, q, lower.tail = TRUE, log.p=FALSE) UseMethod("pmix")
#' @export
#' @rdname mix
qmix <- function(mix, p, lower.tail = TRUE, log.p=FALSE) UseMethod("qmix")
#' @export
#' @rdname mix
rmix <- function(mix, n) UseMethod("rmix")
#' @export
#' @rdname mix
"[[.mix" <- function(mix, ..., rescale=FALSE) {
    ## ensure that the resulting object is a mixture object only
    cl <- grep("mix$", class(mix), ignore.case=TRUE, value=TRUE)
    dl <- dlink(mix)
    if(inherits(mix, "normMix")) s <- sigma(mix)
    mix <- mix[,...,drop=FALSE]
    if(rescale) mix[1,] <- mix[1,] / sum(mix[1,])
    class(mix) <- cl
    dlink(mix) <- dl
    if(inherits(mix, "normMix")) sigma(mix) <- s
    mix
}
#' @export
#' @keywords internal
"[[.betaBinomialMix" <- function(mix, ..., rescale=FALSE) {
    ## ensure that the resulting object has still the N attribute
    n <- attr(mix, "n")
    mix <- NextMethod()
    attr(mix, "n") <- n
    mix
}


## IMPLEMENTATION DETAILS

#' @export
dmix.default <- function(mix, x, log=FALSE) "Unknown mixture"

## default implementation which only needs the density function;
## assumption is that the first argument of the density corresponds to
## the second entry and the third to the last entry in the mix matrix
dmix_impl <- function(dens, mix, x, log) {
    Nc <- ncol(mix)
    ## logic is to replicate the original data vector such that each
    ## item appears nc times which allows easy vectorized calls to
    ## dgamma. Then we cast the result into a matrix with as many rows
    ## as nc components which we sum together with a fast colSums call.
    Nx <- length(x)
    if(is.mixidentity_link(mix)) {
        log_dens  <- matrixStats::colLogSumExps(matrix(log(mix[1,]) + dens(rep(x, each=Nc), rep(mix[2,], times=Nx), rep(mix[3,], times=Nx), log=TRUE), nrow=Nc))
    } else {
        ox <- rep(mixinvlink(mix, x), each=Nc)
        log_dens  <- matrixStats::colLogSumExps(matrix(log(mix[1,]) + rep(mixlJinv_link(mix, x), each=Nc) + dens(ox, rep(mix[2,], times=Nx), rep(mix[3,], times=Nx), log=TRUE), nrow=Nc))
    }
    if(!log)
        return(exp(log_dens))
    return(log_dens)
}

#' @export
dmix.gammaMix <- function(mix, x, log=FALSE) dmix_impl(dgamma, mix, x, log)
#' @export
dmix.betaMix  <- function(mix, x, log=FALSE) dmix_impl(dbeta,  mix, x, log)
#' @export
dmix.normMix  <- function(mix, x, log=FALSE) dmix_impl(dnorm,  mix, x, log)

#' @export
dmix.betaBinomialMix <- function(mix, x, log=FALSE) dmix_impl(Curry(dBetaBinomial, n=attr(mix, "n")),  mix, x, log)

## internal redefinition of negative binomial
##.dnbinomAB <- function(x, a, b, n=1, log=FALSE) dnbinom(x, size=a, prob=(b/n)/((b/n)+1), log=log)
.dnbinomAB <- function(x, a, b, n=1, log=FALSE) dnbinom(x, size=a, prob=b/(b+n), log=log)
#' @export
dmix.gammaPoissonMix <- function(mix, x, log=FALSE) dmix_impl(Curry(.dnbinomAB, n=attr(mix, "n")), mix, x, log)

## DISTRIBUTION FUNCTIONS
#' @export
pmix.default <- function(mix, q, lower.tail = TRUE, log.p=FALSE) "Unknown mixture"

pmix_impl <- function(dist, mix, q, lower.tail = TRUE, log.p=FALSE) {
    Nc <- ncol(mix)
    ## logic is to replicate the original data vector such that each
    ## item appears nc times which allows easy vectorized calls to
    ## dgamma. Then we cast the result into a matrix with as many rows
    ## as nc components which we sum together with a fast colSums call.
    oq <- mixinvlink(mix, q)
    Nq <- length(q)
    if(!log.p)
        return(matrixStats::colSums2(matrix(mix[1,] * dist(rep(oq, each=Nc), rep(mix[2,], times=Nq), rep(mix[3,], times=Nq), lower.tail=lower.tail, log.p=FALSE), nrow=Nc)))
    ## log version is slower, but numerically more stable
    res <- matrixStats::colLogSumExps(matrix(log(mix[1,]) + dist(rep(oq, each=Nc), rep(mix[2,], times=Nq), rep(mix[3,], times=Nq), lower.tail=lower.tail, log.p=TRUE), nrow=Nc))
    return(res)
}

#' @export
pmix.gammaMix <- function(mix, q, lower.tail = TRUE, log.p=FALSE) pmix_impl(pgamma, mix, q, lower.tail, log.p)
#' @export
pmix.betaMix  <- function(mix, q, lower.tail = TRUE, log.p=FALSE) pmix_impl(pbeta,  mix, q, lower.tail, log.p)
#' @export
pmix.normMix  <- function(mix, q, lower.tail = TRUE, log.p=FALSE) pmix_impl(pnorm,  mix, q, lower.tail, log.p)

#' @export
## pmix.betaBinomialMix <- function(mix, q, lower.tail = TRUE, log.p=FALSE) pmix_impl(Curry(pBetaBinomial, n=attr(mix, "n")), mix, q, lower.tail, log.p)
pmix.betaBinomialMix <- function(mix, q, lower.tail = TRUE, log.p=FALSE) {
    assert_that(is.dlink_identity(attr(mix, "link")))
##     ## the analytic solution needs the generalized hypergeometric
##     ## function which is only available in the hyper-geo package which
##     ## is why a numeric integration is performed here
##     ## treat out-of-bounds quantiles special
     out_low  <- q<0
     out_high <- q>attr(mix, "n")
     q[out_low | out_high] <- 0
     dist <- cumsum(dmix.betaBinomialMix(mix, seq(0,max(q))))
     p <- dist[q+1]
     p[out_low] <- 0
     p[out_high] <- 1
     if(!lower.tail) p <- 1-p
     if(log.p) p <- log(p)
     return(p)
}

## internal redefinition of negative binomial
##.pnbinomAB <- function(q, a, b, lower.tail = TRUE, log.p = FALSE ) pnbinom(q, size=a, prob=b/(b+1), lower.tail = lower.tail, log.p = log.p )
.pnbinomAB <- function(q, a, b, n=1, lower.tail = TRUE, log.p = FALSE ) pnbinom(q, size=a, prob=b/(b+n), lower.tail = lower.tail, log.p = log.p )
#' @export
pmix.gammaPoissonMix <- function(mix, q, lower.tail = TRUE, log.p=FALSE) pmix_impl(Curry(.pnbinomAB, n=attr(mix, "n")), mix, q, lower.tail, log.p)

## QUANTILE FUNCTION

#' @export
qmix.default <- function(mix, p, lower.tail = TRUE, log.p=FALSE) "Unknown mixture"

qmix_impl <- function(quant, mix, p, lower.tail = TRUE, log.p=FALSE) {
    Nc <- ncol(mix)
    if(log.p)
        assert_that(all(p <= 0))
    else
        assert_that(all(p >= 0 & p <= 1))
    ## in the simple case of a single component, use R's functions
    if(Nc == 1)
        return(mixlink(mix, quant(p, mix[2,1], mix[3,1], lower.tail=lower.tail, log.p=log.p)))
    assert_that(abs(sum(mix["w",])-1) < sqrt(.Machine$double.eps))
    ## first get the support of the mixture, i.e. the 99.9% CI of each
    ## mixture or lower, if the requested quantile is more in the
    ## tails
    eps <- 1E-1
    plow <- if(log.p) min(c(eps, exp(p), (1-exp(p)))) / 2 else min(c(eps, p, (1-p))) / 2
    phigh <- 1-plow
    qlow <- mixlink(mix, min(quant(rep.int(plow, Nc), mix[2,], mix[3,])))
    qhigh  <- mixlink(mix, max(quant(rep.int(phigh,  Nc), mix[2,], mix[3,])))
    if(is.infinite(qlow )) qlow  <- -sqrt(.Machine$double.xmax)
    if(is.infinite(qhigh)) qhigh <-  sqrt(.Machine$double.xmax)
    res <- rep.int(NA, length(p))
    pboundary <- pmix(mix, c(qlow, qhigh), lower.tail=lower.tail, log.p=log.p)
    for(i in seq_along(p)) {
        ## take advantage of the monotonicity of the CDF function such
        ## that we can use a gradient based method to find the root
        ## 13th Aug 2019: disabled for now; unreliable in rare cases
        ##o <- optimise(function(x) { (pmix(mix,x,lower.tail=lower.tail,log.p=log.p) - p[i])^2 }, c(qlow, qhigh))
        ##res[i] <- o$minimum
        ##if(o$objective > 1E-3) {
        u <- uniroot(function(x) { pmix(mix,x,lower.tail=lower.tail,log.p=log.p) - p[i] },
                     c(qlow, qhigh),
                     f.lower=pboundary[1] - p[i],
                     f.upper=pboundary[2] - p[i],
                     extendInt="upX")
        res[i] <- u$root
        if(u$estim.prec > 1E-3)
            warning("Quantile ", p[i], " possibly imprecise.\nEstimated precision= ", u$estim.prec, ".\nRange = ", qlow, " to ", qhigh, "\n")
        ##}
    }
    res
}

#' @export
qmix.gammaMix <- function(mix, p, lower.tail = TRUE, log.p=FALSE) qmix_impl(qgamma, mix, p, lower.tail, log.p)
#' @export
qmix.betaMix  <- function(mix, p, lower.tail = TRUE, log.p=FALSE) qmix_impl(qbeta,  mix, p, lower.tail, log.p)
#' @export
qmix.normMix  <- function(mix, p, lower.tail = TRUE, log.p=FALSE) qmix_impl(qnorm,  mix, p, lower.tail, log.p)

#' @export
qmix.betaBinomialMix  <- function(mix, p, lower.tail = TRUE, log.p=FALSE) {
    assert_that(is.dlink_identity(attr(mix, "link")))
    ## numerical evaluation
    n <- attr(mix, "n")
    dist <- pmix.betaBinomialMix(mix, seq(0,n-1))
    if(log.p) p <- exp(p)
    ind <- findInterval(p, dist)
    if(!lower.tail) ind <- n - ind
    ind
}

## internal redefinition of negative binomial
##.qnbinomAB <- function(p, a, b, lower.tail = TRUE, log.p = FALSE ) qnbinom(p, size=a, prob=b/(b+1), lower.tail = lower.tail, log.p = log.p )
.qnbinomAB <- function(p, a, b, n=1, lower.tail = TRUE, log.p = FALSE ) qnbinom(p, size=a, prob=b/(b+n), lower.tail = lower.tail, log.p = log.p )
#' @export
##qmix.gammaPoissonMix <- function(mix, p, lower.tail = TRUE, log.p=FALSE) qmix_impl(Curry(.qnbinomAB, n=attr(mix, "n")), mix, p, lower.tail, log.p, discrete=TRUE)

## switched to numeric implementation as discretization seems to cause
## some trouble in the above definitions
qmix.gammaPoissonMix <- function(mix, p, lower.tail = TRUE, log.p=FALSE) {
    assert_that(is.dlink_identity(attr(mix, "link")))
    ## numerical evaulation
    n <- attr(mix, "n")
    eps <- 1e-6
    plow <- if(log.p) min(c(eps, exp(p), (1-exp(p)))) / 2 else min(c(eps, p, (1-p))) / 2
    phigh <- 1-plow
    Nc <- ncol(mix)
    qhigh  <- max(.qnbinomAB(rep.int(phigh,  Nc), mix[2,], mix[3,], n=n))

    dist <- pmix.gammaPoissonMix(mix, seq(0,qhigh-1))
    if(log.p) p <- exp(p)
    ind <- findInterval(p, dist)
    if(!lower.tail) ind <- qhigh - ind
    ind
}


### RANDOM NUMBER GENERATION

#' @export
rmix.default <- function(mix, n) "Unknown mixture"

rmix_impl <- function(rng, mix, n) {
    ind <-  sample.int(ncol(mix), n, replace = TRUE, prob = mix[1,])
    samp <- rng(n, mix[2,ind], mix[3,ind])
    attr(samp, "ind") <- ind
    mixlink(mix, samp)
}

#' @export
rmix.gammaMix <- function(mix, n) rmix_impl(rgamma, mix, n)
#' @export
rmix.betaMix  <- function(mix, n) rmix_impl(rbeta,  mix, n)
#' @export
rmix.normMix  <- function(mix, n) rmix_impl(rnorm,  mix, n)

#' @export
rmix.betaBinomialMix  <- function(mix, n) {
    assert_that(is.dlink_identity(attr(mix, "link")))
    p <- rmix_impl(rbeta,  mix, n)
    samp <- rbinom(n, attr(mix, "n"), p)
    attr(samp, "ind") <- attr(samp, "ind")
    samp
}

## internal redefinition of negative binomial
##.rnbinomAB <- function(n, a, b) rnbinom(n, size=a, prob=b/(b+1))
.rnbinomAB <- function(N, a, b, n=1) rnbinom(N, size=a, prob=b/(b+n))
#' @export
rmix.gammaPoissonMix  <- function(mix, n) rmix_impl(Curry(.rnbinomAB, n=attr(mix, "n")),  mix, n)


#' @export
print.mix <- function(x, digits, ...) {
    tr <- attr(x, "link")
    if(tr$name != "identity") print(tr)
    cat("Mixture Components:\n")
    if(missing(digits)) digits <- NULL
    print.default(format(x, digits=digits), quote=FALSE)
}

#' takes x and transforms it according to the defined link function of
#' the mixture
#' @keywords internal
mixlink <- function(mix, x)
    attr(mix, "link")$link(x)

mixinvlink <- function(mix, x)
    attr(mix, "link")$invlink(x)

mixJinv_orig <- function(mix, x)
    attr(mix, "link")$Jinv_orig(x)

mixlJinv_orig <- function(mix, x)
    attr(mix, "link")$lJinv_orig(x)

mixlJinv_link <- function(mix, l)
    attr(mix, "link")$lJinv_link(l)

is.mixidentity_link <- function(mix, l)
    is.dlink_identity(attr(mix, "link"))
