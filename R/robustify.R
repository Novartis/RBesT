#' Robustify Mixture Priors
#'
#' Add a non-informative component to a mixture prior.
#'
#' @param priormix orior (mixture of conjugate distributions).
#' @param weight weight given to the non-informative component (0 < \code{weight} < 1).
#' @param mean mean of the non-informative component. It is recommended to set this parameter explicitly.
#' @param n number of observations the non-informative prior
#' corresponds to, defaults to 1.
#' @param ... optional arguments are ignored.
#'
#' @details It is recommended to robustify informative priors derived
#' with \code{\link{gMAP}} using unit-information priors . This
#' protects against prior-data conflict, see for example
#' \emph{Schmidli et al., 2015}.
#'
#' The procedure can be used with beta, gamma and normal mixture
#' priors. A unit-information prior (see \emph{Kass and Wasserman,
#' 1995}) corresponds to a prior which represents the observation of
#' n=1 at the null hypothesis. As the null is problem dependent we
#' \emph{strongly recommend} to make use of the \code{mean} argument
#' accordingly. See below for the definition of the default mean.
#'
#' The weights of the mixture priors are rescaled to \code{(1-weight)}
#' while the non-informative prior is assigned the \code{weight}
#' given.
#'
#' @return New mixture with an extra non-informative component named
#' \code{robust}.
#' 
#' @references Schmidli H, Gsteiger S, Roychoudhury S, O'Hagan A,
#' Spiegelhalter D, Neuenschwander B.  Robust meta-analytic-predictive
#' priors in clinical trials with historical control information.
#' \emph{Biometrics} 2014;70(4):1023-1032.
#'
#' Kass RE, Wasserman L A Reference Bayesian Test for Nested
#' Hypotheses and its Relationship to the Schwarz Criterion \emph{J
#' Amer Statist Assoc} 1995; 90(431):928-934.
#'
#' @seealso \code{\link{mixcombine}}
#' 
#' @examples
#' bmix <- mixbeta(inf1=c(0.2, 8, 3), inf2=c(0.8, 10, 2))
#' plot(bmix)
#' rbmix <- robustify(bmix, weight=0.1, mean=0.5)
#' rbmix
#' plot(rbmix)
#' 
#' gmnMix <- mixgamma(inf1=c(0.2, 2, 3), inf2=c(0.8, 2, 5), param="mn")
#' plot(gmnMix)
#' rgmnMix <- robustify(gmnMix, weight=0.1, mean=2)
#' rgmnMix
#' plot(rgmnMix)
#' 
#' nm <- mixnorm(inf1=c(0.2, 0.5, 0.7), inf2=c(0.8, 2, 1), sigma=2)
#' plot(nm)
#' rnMix <- robustify(nm, weight=0.1, mean=0, sigma=2)
#' rnMix
#' plot(rnMix)
#' 
#' @export
robustify <- function(priormix, weight, mean, n=1, ...) UseMethod("robustify")

#' @export
robustify.default <- function(priormix, weight, mean, n=1, ...) "Unknown density"

#' @describeIn robustify The default \code{mean} is set to 1/2 which
#' represents no difference between the occurrence rates for one of the
#' two outcomes. As the uniform \code{Beta(1,1)} is more appropriate in
#' practical applications, \code{RBesT} uses \code{n+1} as the sample
#' size such that the default robust prior is the uniform instead of
#' the \code{Beta(1/2,1/2)} which strictly defined would be the unit
#' information prior in this case.
#' @export
robustify.betaMix <- function(priormix, weight, mean, n=1, ...) {
    assert_number(weight, lower=0, upper=1)
    assert_number(n, lower=0, finite=TRUE)
    if(missing(mean)) {
        message("Using default mean for robust component of 1/2.")
        mean <- 1/2
    }
    assert_number(mean, lower=0, upper=1)
    rob <- mixbeta(robust=c(1, mean, n+1), param="mn")
    mixcombine(priormix, rob, weight=c(1-weight, weight))
}

#' @describeIn robustify The default \code{mean} is set to the mean of the
#' prior mixture. It is strongly recommended to explicitly set the
#' mean to the location of the null hypothesis.
#' @export
robustify.gammaMix <- function(priormix, weight, mean, n=1, ...) {
    assert_number(weight, lower=0, upper=1)
    assert_number(n, lower=0, finite=TRUE)
    if(missing(mean)) {
        s <- summary(priormix)
        message(paste("Using default mean for robust component; the mean of the prior which is", s["mean"], "."))
        mean <- s["mean"]
    }
    assert_number(mean, lower=0, finite=TRUE)
    rob <- mixgamma(robust=c(1, mean, n), param="mn", likelihood=likelihood(priormix))
    mixcombine(priormix, rob, weight=c(1-weight, weight))
}

#' @describeIn robustify The default \code{mean} is set to the mean
#' of the prior mixture. It is strongly recommended to explicitly set
#' the mean to the location of the null hypothesis, which is very
#' often equal to 0. It is also recommended to explicitly set the
#' sampling standard deviation using the \code{sigma} argument.
#' @param sigma Sampling standard deviation for the case of Normal
#' mixtures.
#' @export
robustify.normMix <- function(priormix, weight, mean, n=1, ..., sigma) {
    assert_number(weight, lower=0, upper=1)
    assert_number(n, lower=0, finite=TRUE)
    if(missing(mean)) {
        s <- summary(priormix)
        message(paste("Using default mean for robust component; the mean of the prior which is", s["mean"], "."))
        mean <- s["mean"]
    }
    assert_number(mean, finite=TRUE)
    if(missing(sigma)) {
        message("Using default prior reference scale ", RBesT::sigma(priormix))
        sigma <- RBesT::sigma(priormix)
    }
    rob <- mixnorm(robust=c(1, mean, n), param="mn", sigma=sigma)
    mixcombine(priormix, rob, weight=c(1-weight, weight))
}

