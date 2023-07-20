#' @name mixgamma
#'
#' @title The Gamma Mixture Distribution
#'
#' @description The gamma mixture density and auxiliary functions.
#'
#' @param ... List of mixture components.
#' @param param Determines how the parameters in the list are
#' interpreted. See details.
#' @param likelihood Defines with what likelihood the Gamma density is used (Poisson or Exp). Defaults to \code{poisson}.
#' @param m Vector of means of the Gamma mixture components
#' @param s Vector of standard deviations of the gamma mixture components,
#' @param n Vector of sample sizes of the Gamma mixture components.
#' @param drop Delete the dimensions of an array which have only one level.
#' @param object Gamma mixture object.
#' @param probs Quantiles reported by the \code{summary} function.
#'
#' @details Each entry in the \code{...} argument list is expected to
#' be a triplet of numbers which defines the weight \eqn{w_k}, first
#' and second parameter of the mixture component \eqn{k}. A triplet
#' can optionally be named which will be used appropriately.
#'
#' The first and second parameter can be given in different
#' parametrizations which is set by the \code{param} option:
#' \describe{
#' \item{ab}{Natural parametrization of Gamma density (\code{a}=shape and \code{b}=rate). Default. }
#' \item{ms}{Mean and standard deviation, \eqn{m=a/b} and \eqn{s=\sqrt{a}/b}.}
#' \item{mn}{Mean and number of observations. Translation to natural
#' parameter depends on the \code{likelihood} argument. For a Poisson
#' likelihood \eqn{n=b} (and \eqn{a=m \cdot n}{a=m n}), for an Exp
#' likelihood \eqn{n=a} (and \eqn{b=n/m}).}
#' }
#'
#' @family mixdist
#'
#' @return \code{mixgamma} returns a gamma mixture with the specified mixture components.
#' \code{ms2gamma} and
#' \code{mn2gamma} return the equivalent natural \code{a} and \code{b} parametrization given
#' parameters \code{m}, \code{s}, or \code{n}.
#'
#' @examples
#' # Gamma mixture with robust and informative component
#' gmix <- mixgamma(rob=c(0.3, 20, 4), inf=c(0.7, 50, 10))
#'
#' # objects can be printed
#' gmix
#' # or explicitly
#' print(gmix)
#'
#' # summaries are defined
#' summary(gmix)
#'
#' # sub-components may be extracted
#' # by component number
#' gmix[[2]]
#' # or component name
#' gmix[["inf"]]
#'
#' # alternative mean and standard deviation parametrization
#' gmsMix <- mixgamma(rob=c(0.5, 8, 0.5), inf=c(0.5, 9, 2), param="ms")
#'
#' # or mean and number of observations parametrization
#' gmnMix <- mixgamma(rob=c(0.2, 2, 1), inf=c(0.8, 2, 5), param="mn")
#'
#' # and mixed parametrizations are also possible
#' gfmix <- mixgamma(rob1=c(0.15, mn2gamma(2, 1)), rob2=c(0.15, ms2gamma(2, 5)), inf=c(0.7, 50, 10))
NULL

#' @rdname mixgamma
#' @export
mixgamma <- function(..., param=c("ab", "ms", "mn"), likelihood=c("poisson", "exp")) {
    mix <- mixdist3(...)
    assert_matrix(mix, nrows=3, any.missing=FALSE)
    param <- match.arg(param)
    likelihood <- match.arg(likelihood)
    mix[c(2,3),] <- switch(param,
                           ab=mix[c(2,3),],
                           ms=t(ms2gamma(mix[2,], mix[3,], FALSE)),
                           mn=t(mn2gamma(mix[2,], mix[3,], likelihood, FALSE))
                           )
    rownames(mix) <- c("w", "a", "b")
    assert_that(all(mix["a",] > 0))
    assert_that(all(mix["b",] > 0))
    class(mix) <- c("gammaMix", "mix")
    likelihood(mix) <- likelihood
    mix
}

#' @rdname mixgamma
#' @export
ms2gamma <- function(m, s, drop=TRUE) {
    b <- m/s^2
    ab <- cbind(a=m*b, b=b)
    if(drop) ab <- drop(ab)
    ab
}

#' @rdname mixgamma
#' @export
mn2gamma <- function(m, n, likelihood=c("poisson", "exp"), drop=TRUE) {
    assert_that(all(n>=0))
    likelihood <- match.arg(likelihood)
    ab <- switch(likelihood, poisson=cbind(a=m*n, b=n), exp=cbind(a=n, b=n/m))
    if(drop) ab <- drop(ab)
    ab
}

#' @rdname mixgamma
#' @method print gammaMix
#' @param x The mixture to print
#' @export
print.gammaMix <- function(x, ...) {
    cat("Univariate Gamma mixture\n")
    NextMethod()
}

#' @rdname mixgamma
#' @method print gammaPoissonMix
#' @param x The mixture to print
#' @export
print.gammaPoissonMix <- function(x, ...) {
    cat("Univariate Gamma-Poisson mixture\n")
    NextMethod()
}

#' @rdname mixgamma
#' @method print gammaExpMix
#' @param x The mixture to print
#' @export
print.gammaExpMix <- function(x, ...) {
    cat("Univariate Gamma-Exponential mixture\n")
    NextMethod()
}

#' @rdname mixgamma
#' @method summary gammaMix
#' @export
summary.gammaMix <- function(object, probs=c(0.025,0.5,0.975), ...) {
    p <- object[1,]
    a <- object[2,]
    b <- object[3,]
    m <- a/b
    v <- a/b^2
    ## calculate mean of the second moment
    m2 <- v + m^2
    ## from this we can get the mean and variance of the mixture
    mmix <- sum(p * m)
    vmix <- sum(p * (m2 - (mmix)^2 ) )
    q <- c()
    if(length(probs) != 0) {
        q <- qmix.gammaMix(object, p=probs)
        names(q) <- paste(format(probs*100,digits=2), "%", sep="")
    }
    c(mean=mmix, sd=sqrt(vmix), q)
}

#' @rdname mixgamma
#' @method summary gammaPoissonMix
#' @export
summary.gammaPoissonMix <- function(object, probs=c(0.025,0.5,0.975), ...) {
    n <- attr(object, "n")
    p <- object[1,]
    a <- object[2,]
    b <- object[3,]/n
    m <- a/b
    v <- a *(b+1)/b^2
    ## calculate mean of the second moment
    m2 <- v + m^2
    ## from this we can get the mean and variance of the mixture
    mmix <- sum(p * m)
    vmix <- sum(p * (m2 - (mmix)^2 ) )
    q <- qmix.gammaPoissonMix(object, p=probs)
    if(length(q) != 0)
        names(q) <- paste(format(probs*100,digits=2), "%", sep="")
    c(mean=mmix, sd=sqrt(vmix), q)
}

