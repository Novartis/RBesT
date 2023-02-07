#' Predictive Distributions for Mixture Distributions
#'
#' @description
#' Predictive distribution for mixture of conjugate distributions
#' (beta, normal, gamma).
#'
#' @param mix mixture distribution
#' @param n predictive sample size, set by default to 1
#' @param ... includes arguments which depend on the specific prior-likelihood pair, see description below.
#'
#' @details Given a mixture density (either a posterior or a prior)
#'
#' \deqn{p(\theta,\mathbf{w},\mathbf{a},\mathbf{b})}{p(\theta,w,a,b)}
#'
#' and a data likelihood of
#'
#' \deqn{y|\theta \sim f(y|\theta),}{y|\theta ~ f(y|\theta),}
#'
#' the predictive distribution of a one-dimensional summary \eqn{y_n}
#' of $n$ future observations is distributed as
#'
#' \deqn{y_n \sim \int p(\theta,\mathbf{w},\mathbf{a},\mathbf{b}) \, f(y_n|\theta) \, d\theta .}{y_n ~ \int p(u,w,a,b) \, f(y_n|u) du .}
#'
#' This distribution is the marginal distribution of the data under
#' the mixture density. For binary and Poisson data \eqn{y_n =
#' \sum_{i=1}^n y_i} is the sum over future events. For normal data,
#' it is the mean\eqn{\bar{y}_n = 1/n \sum_{i=1}^n y_i}.
#'
#' @return The function returns for a normal, beta or gamma mixture
#' the matching predictive distribution for \eqn{y_n}.
#'
#' @template conjugate_pairs
#'
#' @examples
#'
#' # Example 1: predictive distribution from uniform prior.
#' bm <- mixbeta(c(1,1,1))
#' bmPred <- preddist(bm, n=10)
#' # predictive proabilities and cumulative predictive probabilities
#' x <- 0:10
#' d <- dmix(bmPred, x)
#' names(d) <- x
#' barplot(d)
#' cd <- pmix(bmPred, x)
#' names(cd) <- x
#' barplot(cd)
#' # median
#' mdn <- qmix(bmPred,0.5)
#' mdn
#'
#' # Example 2: 2-comp Beta mixture
#'
#' bm <- mixbeta( inf=c(0.8,15,50),rob=c(0.2,1,1))
#' plot(bm)
#' bmPred <- preddist(bm,n=10)
#' plot(bmPred)
#' mdn <- qmix(bmPred,0.5)
#' mdn
#' d <- dmix(bmPred,x=0:10)
#' \donttest{
#' n.sim <- 100000
#' r <-  rmix(bmPred,n.sim)
#' d
#' table(r)/n.sim
#' }
#'
#' # Example 3: 3-comp Normal mixture
#'
#' m3 <- mixnorm( c(0.50,-0.2,0.1),c(0.25,0,0.2), c(0.25,0,0.5), sigma=10)
#' print(m3)
#' summary(m3)
#' plot(m3)
#' predm3 <- preddist(m3,n=2)
#' plot(predm3)
#' print(predm3)
#' summary(predm3)
#'
#' @export
preddist <- function(mix, ...) UseMethod("preddist")
#' @export
preddist.default <- function(mix, ...) "Unknown distribution"

#' @describeIn preddist Obtain the matching predictive distribution
#' for a beta distribution, the BetaBinomial.
#' @export
preddist.betaMix <- function(mix, n=1, ...) {
    attr(mix, "n") <- n
    class(mix) <- c("betaBinomialMix", "mix")
    mix
}

#' @describeIn preddist Obtain the matching predictive distribution
#' for a Normal with constant standard deviation. Note that the
#' reference scale of the returned Normal mixture is meaningless as the
#' individual components are updated appropriatley.
#' @param sigma The fixed reference scale of a normal mixture. If left
#' unspecified, the default reference scale of the mixture is assumed.
#' @export
preddist.normMix <- function(mix, n=1, sigma, ...) {
    if(missing(sigma)) {
        sigma <- RBesT::sigma(mix)
        message("Using default mixture reference scale ", sigma)
    }
    assert_number(sigma, lower=0)
    sigma_ref <- sigma
    ## note: this is effectivley a hierarchical model as we give the
    ## distribution of the sum of n variables which have exactly the
    ## same mean (which is sampled from the normal)
    ## old: sum over y_i
    ##mix[2,] <- mix[2,] * n
    ##mix[3,] <- sqrt(mix[3,]^2 * n^2 + tau^2 * n )
    ## now: predictive for \bar{y}_n, the mean
    mix[3,] <- sqrt(mix[3,]^2  + sigma_ref^2/n )
    class(mix) <- c("normMix", "mix")
    mix
}

#' @describeIn preddist Obtain the matching predictive distribution
#' for a Gamma. Only Poisson likelihoods are supported.
#' @export
preddist.gammaMix <- function(mix, n=1, ...) {
    assert_set_equal(likelihood(mix), "poisson")
    attr(mix, "n") <- n
    class(mix) <- c(switch(likelihood(mix),
                           poisson="gammaPoissonMix",
                           exp="gammaExpMix"), "mix")
    mix
}

