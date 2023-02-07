#' Probability of Success for a 1 Sample Design
#'
#' The \code{pos1S} function defines a 1 sample design (prior, sample
#' size, decision function) for the calculation of the frequency at
#' which the decision is evaluated to 1 when assuming a distribution
#' for the parameter. A function is returned which performs the
#' actual operating characteristics calculations.
#'
#' @template args-boundary1S
#'
#' @details The \code{pos1S} function defines a 1 sample design and
#' returns a function which calculates its probability of success.
#' The probability of success is the frequency with which the decision
#' function is evaluated to 1 under the assumption of a given true
#' distribution of the data implied by a distirbution of the parameter
#' \eqn{\theta}.
#'
#' Calling the \code{pos1S} function calculates the critical value
#' \eqn{y_c} and returns a function which can be used to evaluate the
#' PoS for different predictive distributions and is evaluated as
#'
#' \deqn{ \int F(y_c|\theta) p(\theta) d\theta, }
#'
#' where \eqn{F} is the distribution function of the sampling
#' distribution and \eqn{p(\theta)} specifies the assumed true
#' distribution of the parameter \eqn{\theta}. The distribution
#' \eqn{p(\theta)} is a mixture distribution and given as the
#' \code{mix} argument to the function.
#'
#' @return Returns a function that takes as single argument
#' \code{mix}, which is the mixture distribution of the control
#' parameter. Calling this function with a mixture distribution then
#' calculates the PoS.
#'
#' @family design1S
#'
#' @examples
#'
#' # non-inferiority example using normal approximation of log-hazard
#' # ratio, see ?decision1S for all details
#' s <- 2
#' flat_prior <- mixnorm(c(1,0,100), sigma=s)
#' nL <- 233
#' theta_ni <- 0.4
#' theta_a <- 0
#' alpha <- 0.05
#' beta  <- 0.2
#' za <- qnorm(1-alpha)
#' zb <- qnorm(1-beta)
#' n1 <- round( (s * (za + zb)/(theta_ni - theta_a))^2 )
#' theta_c <- theta_ni - za * s / sqrt(n1)
#'
#' # assume we would like to conduct at an interim analysis
#' # of PoS after having observed 20 events with a HR of 0.8.
#' # We first need the posterior at the interim ...
#' post_ia <- postmix(flat_prior, m=log(0.8), n=20)
#'
#' # dual criterion
#' decComb <- decision1S(c(1-alpha, 0.5), c(theta_ni, theta_c), lower.tail=TRUE)
#'
#' # ... and we would like to know the PoS for a successful
#' # trial at the end when observing 10 more events
#' pos_ia <- pos1S(post_ia, 10, decComb)
#'
#' # our knowledge at the interim is just the posterior at
#' # interim such that the PoS is
#' pos_ia(post_ia)
#'
#'
#' @export
pos1S <- function(prior, n, decision, ...) UseMethod("pos1S")
#' @export
pos1S.default <- function(prior, n, decision, ...) "Unknown density"

#' @templateVar fun pos1S
#' @template design1S-binomial
#' @export
pos1S.betaMix <- function(prior, n, decision, ...) {

    crit <- decision1S_boundary(prior, n, decision)
    lower.tail <- attr(decision, "lower.tail")

    design_fun <- function(mix) {
        pred_dtheta <- preddist(mix, n=n)
        pmix(pred_dtheta, crit, lower.tail=lower.tail)
    }
    design_fun
}

#' @templateVar fun pos1S
#' @template design1S-normal
#' @export
pos1S.normMix <- function(prior, n, decision, sigma, eps=1e-6, ...) {
    ## distributions of the means of the data generating distributions
    ## for now we assume that the underlying standard deviation
    ## matches the respective reference scales
    if(missing(sigma)) {
        sigma <- RBesT::sigma(prior)
        message("Using default prior reference scale ", sigma)
    }
    assert_number(sigma, lower=0)

    sigma(prior) <- sigma

    crit <- decision1S_boundary(prior, n, decision, sigma, eps)

    ## check where the decision is 1, i.e. left or right
    lower.tail <- attr(decision, "lower.tail")

    design_fun <- function(mix) {
        pred_dtheta_mean <- preddist(mix, n=n, sigma=sigma)
        pmix(pred_dtheta_mean, crit, lower.tail=lower.tail)
    }
    design_fun
}

#' @templateVar fun pos1S
#' @template design1S-poisson
#' @export
pos1S.gammaMix <- function(prior, n, decision, eps=1e-6, ...) {

    crit <- decision1S_boundary(prior, n, decision, eps)
    lower.tail <- attr(decision, "lower.tail")

    design_fun <- function(mix) {
        assert_that(likelihood(prior) == "poisson")
        pred_dtheta_sum <- preddist(mix, n=n)
        pmix(pred_dtheta_sum, crit, lower.tail=lower.tail)
    }
    design_fun
}
