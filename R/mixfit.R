#' Fit of Mixture Densities to Samples
#'
#' Expectation-Maximization (EM) based fitting of parametric mixture
#' densities to numerical samples. This provides a convenient approach
#' to approximate MCMC samples with a parametric mixture distribution.
#'
#' @param sample Sample to be fitted.
#' @param type Mixture density to use. Can be either norm, beta or gamma.
#' @param thin Thinning applied to the sample. See description for default behavior.
#' @param ... Parameters passed to the low-level EM fitting functions. Parameter \code{Nc} is mandatory.
#'
#' @details
#' Parameters of EM fitting functions
#' \describe{
#' \item{Nc}{Number of mixture components. Required parameter.}
#' \item{mix_init}{Initial mixture density. If missing (default) then a k-nearest-neighbor algorithm is used to find an initial mixture density.}
#' \item{Ninit}{Number of data points used for initialization. Defaults to 50.}
#' \item{verbose}{If set to \code{TRUE} the function will inform about fitting process}
#' \item{maxIter}{Maximal number of iterations. Defaults to 500.}
#' \item{tol}{Defines a convergence criteria as an upper bound for the change in the log-likelihood, i.e. once the derivative (with respect to iterations) of the log-likelihood falls below \code{tol}, the function declares convergence and stops.}
#' \item{eps}{Must be a triplet of numbers which set the desired accuracy of the inferred parameters per mixture component. See below for a description of the parameters used during EM. EM is stopped once a running mean of the absolute difference between the last successive \code{Neps} estimates is below the given \code{eps} for all parameters. Defaults to 5E-3 for each parameter.}
#' \item{Neps}{Number of iterations used for the running mean of parameter estimates to test for convergence. Defaults to 5.}
#' \item{constrain_gt1}{Logical value controlling if the Beta EM constrains all parameters a & b to be greater than 1. By default constraints are turned on (new since 1.6-0).}
#' }
#'
#' By default the EM convergence is declared when
#' the desired accuracy of the parameters has been reached over the last
#' \code{Neps} estimates. If \code{tol} and \code{Neps} is specified, then
#' whatever criterion is met first will stop the EM.
#'
#' The parameters per component \eqn{k} used internally during fitting
#' are for the different EM procedures:
#' \describe{
#' \item{normal}{\eqn{logit(w_k), \mu_k, \log(\sigma_k)}}
#' \item{beta}{\eqn{logit(w_k), \log(a_k), \log(b_k)}}
#' \item{constrained beta}{\eqn{logit(w_k), \log(a_k-1), \log(b_k-1)}}
#' \item{gamma}{\eqn{logit(w_k), \log(\alpha_k), \log(\beta_k)}}
#' }
#'
#' \emph{Note:} Whenever no \code{mix_init} argument is given,
#' the EM fitting routines assume that the data vector is given in
#' random order. If in the unlikely event that the EM gets caught in a
#' local extremum, then random reordering of the data vector may
#' alleviate the issue.
#'
#' @return
#'
#' A mixture object according the requested \code{type} is
#' returned. The object has additional information attached, i.e. the
#' log-likelihood can be queried and diagnostic plots can be
#' generated. See links below.
#'
#' @family EM
#'
#' @references Dempster A.P., Laird N.M., Rubin D.B. Maximum
#' Likelihood from Incomplete Data via the EM Algorithm. \emph{Journal
#' of the Royal Statistical Society, Series B} 1977; 39 (1): 1-38.
#'
#' @examples
#' bmix <- mixbeta(rob=c(0.2, 1, 1), inf=c(0.8, 10, 2))
#'
#' bsamp <- rmix(bmix, 1000)
#'
#' bfit <- mixfit(bsamp, type="beta", Nc=2)
#'
#' # diagnostic plots can easily by generated from the EM fit with
#' bfit.check <- plot(bfit)
#'
#' names(bfit.check)
#'
#' # check convergence of parameters
#' bfit.check$mix
#' bfit.check$mixdens
#' bfit.check$mixecdf
#'
#' # obtain the log-likelihood
#' logLik(bfit)
#'
#' # or AIC
#' AIC(bfit)
#'

#' @export
mixfit <- function(sample, type=c("norm", "beta", "gamma"), thin, ...) UseMethod("mixfit")

#' @describeIn mixfit Performs an EM fit for the given
#' sample. Thinning is applied only if thin is specified.
#' @export
mixfit.default <- function(sample, type=c("norm", "beta", "gamma"), thin, ...)
{
    type <- match.arg(type)
    assert_that(type %in% c("norm", "beta", "gamma"))
    EM <- switch(type, norm=EM_nmm, beta=EM_bmm_ab, gamma=EM_gmm)
    if(!missing(thin)) {
        assert_that(thin >= 1)
        sample <- sample[seq(1,NROW(sample),by=thin),drop=FALSE]
    }
    EM(sample, ...)
}

#' @describeIn mixfit Fits the default predictive distribution from a
#' gMAP analysis. Automatically obtains the predictive distribution of
#' the intercept only case on the response scale mixture from the
#' \code{\link{gMAP}} object. For the binomial case a beta mixture,
#' for the gaussian case a normal mixture and for the Poisson case a
#' gamma mixture will be used. In the gaussian case, the resulting
#' normal mixture will set the reference scale to the estimated
#' sigma in \code{\link{gMAP}} call.
#' @export
mixfit.gMAP <- function(sample, type, thin, ...)
{
    family <- sample$family$family
    ## automatically thin sample as estimated by gMAP function
    if(missing(thin)) {
        thin <- sample$thin
    }
    assert_that(thin >= 1)
    type <- switch(sample$family$family, binomial = "beta", gaussian = "norm", poisson = "gamma", "unknown")
    sim <- rstan::extract(sample$fit, pars="theta_resp_pred", inc_warmup=FALSE, permuted=FALSE)
    sim <- as.vector(sim[seq(1,dim(sim)[1],by=thin),,])
    mix <- mixfit.default(sim, type, thin=1, ...)
    ## for the case of normal data, read out the estimated reference
    ## scale
    if(family == "gaussian" & !is.null(sample$sigma_ref))
        sigma(mix) <- sample$sigma_ref
    set_likelihood(mix, family)
}

#' @describeIn mixfit Fits a mixture density for each prediction from
#' the \code{\link{gMAP}} prediction.
#' @export
mixfit.gMAPpred <- function(sample, type, thin, ...)
{
    if(attr(sample, "type") == "response") {
        type <- switch(attr(sample, "family")$family, binomial = "beta", gaussian = "norm", poisson = "gamma", "unknown")
        family <- attr(sample, "family")$family
    } else {
        type <- "norm"
        family <- "gaussian"
    }
    res <- list()
    for(i in 1:dim(sample)[2]) {
        ## for the case of normal data, read out the estimated reference
        ## scale
        ## note: gMAPpred objects are already thinned down
        res[[i]] <- set_likelihood(mixfit.default(sample[,i], type=type, thin=1, ...), family)
        if(family == "gaussian" & !is.null(attr(sample, "sigma_ref"))) {
            sigma(res[[i]]) <- attr(sample, "sigma_ref")
        }
    }
    res
}
#' @describeIn mixfit Fits a mixture density for an MCMC sample. It is
#' recommended to provide a thinning argument which roughly yields
#' independent draws (i.e. use \code{\link{acf}} to identify a
#' thinning lag with small auto-correlation). The input array is
#' expected to have 3 dimensions which are nested as iterations,
#' chains, and draws.
#' @export
mixfit.array <- function(sample, type, thin, ...)
{
    Nmcmc <- prod(dim(sample))
    if(dim(sample)[3] != 1)
        stop("Only univariate data is supported.")
    mixfit.default(sample, type, thin, ...)
}

set_likelihood <- function(mix, family) {
    if(family == "binomial") {
        likelihood(mix) <- "binomial"
    } else if(family == "gaussian") {
        likelihood(mix) <- "normal"
    } else if(family == "poisson") {
        likelihood(mix) <- "poisson"
    }
    mix
}
