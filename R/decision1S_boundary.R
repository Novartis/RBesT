#' Decision Boundary for 1 Sample Designs
#'
#' Calculates the decision boundary for a 1 sample design. This is the
#' critical value at which the decision function will change from 0
#' (failure) to 1 (success).
#' 
#' @template args-boundary1S
#' 
#' @details The specification of the 1 sample design (prior, sample
#' size and decision function, \eqn{D(y)}), uniquely defines the
#' decision boundary
#' 
#' \deqn{y_c = \max_y\{D(y) = 1\},}{y_c = max_{y}{D(y) = 1},}
#'
#' which is the maximal value of \eqn{y} whenever the decision \eqn{D(y)}
#' function changes its value from 1 to 0 for a decision function
#' with \code{lower.tail=TRUE} (otherwise the definition is \eqn{y_c =
#' \max_{y}\{D(y) = 0\}}{y_c = max_{y}{D(y) = 0}}). The decision
#' function may change at most at a single critical value as only
#' one-sided decision functions are supported. Here,
#' \eqn{y} is defined for binary and Poisson endpoints as the sufficient
#' statistic \eqn{y = \sum_{i=1}^{n} y_i} and for the normal
#' case as the mean \eqn{\bar{y} = 1/n \sum_{i=1}^n y_i}.
#' 
#' The convention for the critical value \eqn{y_c} depends on whether
#' a left (\code{lower.tail=TRUE}) or right-sided decision function
#' (\code{lower.tail=FALSE}) is used. For \code{lower.tail=TRUE} the
#' critical value \eqn{y_c} is the largest value for which the
#' decision is 1, \eqn{D(y \leq y_c) = 1}, while for
#' \code{lower.tail=FALSE} then \eqn{D(y > y_c) = 1} holds. This is
#' aligned with the cumulative density function definition within R
#' (see for example \code{\link{pbinom}}).
#' 
#' @return Returns the critical value \eqn{y_c}.
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
#' # double criterion design
#' # statistical significance (like NI design)
#' dec1 <- decision1S(1-alpha, theta_ni, lower.tail=TRUE)
#' # require mean to be at least as good as theta_c
#' dec2 <- decision1S(0.5, theta_c, lower.tail=TRUE)
#' # combination
#' decComb <- decision1S(c(1-alpha, 0.5), c(theta_ni, theta_c), lower.tail=TRUE)
#' 
#' # critical value of double criterion design
#' decision1S_boundary(flat_prior, nL, decComb)
#'
#' # ... is limited by the statistical significance ...
#' decision1S_boundary(flat_prior, nL, dec1)
#'
#' # ... or the indecision point (whatever is smaller)
#' decision1S_boundary(flat_prior, nL, dec2)
#'
#' @export
decision1S_boundary <- function(prior, n, decision, ...) UseMethod("decision1S_boundary")
#' @export
decision1S_boundary.default <- function(prior, n, decision, ...) "Unknown density"

#' @templateVar fun decision1S_boundary
#' @template design1S-binomial
#' @export
decision1S_boundary.betaMix <- function(prior, n, decision, ...) {

    VdecisionLazy <- Vectorize(function(r) { decision(postmix(prior, r=r, n=n)) - 0.25 } )

    ## find decision boundary
    bounds <- VdecisionLazy(c(0,n))
    lower.tail <- attr(decision, "lower.tail")
    if(prod(bounds) > 0) {
        ## decision is always the same
        if(lower.tail) {
            crit <- ifelse(bounds[1] < 0, 0, n+1)
        } else {
            crit <- ifelse(bounds[1] < 0, n, -1)
        }        
    } else {        
        crit <- uniroot_int(VdecisionLazy, c(0,n),
                            f.lower=bounds[1], f.upper=bounds[2])
    }

    ## crit is always pointing to the 0 just before the decision which
    ## is why we need a discrimination here
    if(lower.tail) {
        crit <- crit - 1
    }        

    crit
}


## returns a function object which is the decision boundary. That is
## the function finds at a regular grid between llim1 and ulim1 the
## roots of the decision function and returns an interpolation
## function object
solve_boundary1S_normMix <- function(decision, mix, n, lim) {

    sigma <- sigma(mix)
    
    cond_decisionStep <- function() {
        fn <- function(m) {
            decision(postmix(mix, m=m, se=sigma/sqrt(n))) - 0.75
        }
        Vectorize(fn)
    }

    ## ensure that at the limiting boundaries the decision function
    ## has a different sign (which must be true)
    ind_fun <- cond_decisionStep()
    dec_bounds <- ind_fun(lim)
    while(prod(dec_bounds) > 0) {
        w <- diff(lim)
        lim <- c(lim[1] - w/2, lim[2] + w/2)
        dec_bounds <- ind_fun(lim)
    }

    uniroot(ind_fun, interval=lim,
            f.lower=dec_bounds[1], f.upper=dec_bounds[2])$root
}

#' @templateVar fun decision1S_boundary
#' @template design1S-normal
#' @export
decision1S_boundary.normMix <- function(prior, n, decision, sigma, eps=1e-6, ...) {
    ## distributions of the means of the data generating distributions
    ## for now we assume that the underlying standard deviation
    ## matches the respective reference scales
    if(missing(sigma)) {
        sigma <- RBesT::sigma(prior)
        message("Using default prior reference scale ", sigma)
    }
    assert_number(sigma, lower=0)
    
    sd_samp <- sigma / sqrt(n)

    sigma(prior) <- sigma

    ## change the reference scale of the prior such that the prior
    ## represents the distribution of the respective means
    ##mean_prior <- prior
    ##sigma(mean_prior) <- sd_samp

    m <- summary(prior, probs=c())["mean"]
    
    lim <- qnorm(p=c(eps/2, 1-eps/2), mean=m, sd=sd_samp)

    ## find the boundary of the decision function within the domain we integrate
    crit <- solve_boundary1S_normMix(decision, prior, n, lim)
    
    crit
}


#' @templateVar fun decision1S_boundary
#' @template design1S-poisson
#' @export
decision1S_boundary.gammaMix <- function(prior, n, decision, eps=1e-6, ...) {
    assert_that(likelihood(prior) == "poisson")

    cond_decisionStep <- function() {
        fn <- function(m) {
            decision(postmix(prior, n=n, m=m/n)) - 0.25
        }
        Vectorize(fn)
    }

    Vdecision <- cond_decisionStep()

    m <- summary(prior, probs=c())["mean"]
    
    lambda_prior <- m * n
    
    lim <- qpois(p=c(eps/2, 1-eps/2), lambda=lambda_prior)
    lim[1] <- 0
    
    bounds <- Vdecision(lim)
    ## make sure there is a decision somewhere
    while(prod(bounds) > 0) {
        lim[2] <- round(lim[2]*2)
        bounds[2] <- Vdecision(lim[2])
        if(diff(lim) > 1e9)
            break;
    }
    
    lower.tail <- attr(decision, "lower.tail")
    
    ## check if the decision is constantly 1 or 0
    if(prod(bounds) > 0) {
        ## decision is always the same
        if(lower.tail) {
            crit <- ifelse(bounds[1] < 0, 0, Inf)
        } else {
            crit <- ifelse(bounds[1] < 0, Inf, -1)
        }        
    } else {        
        ## find decision boundary
        crit <- uniroot_int(Vdecision,
                            c(lim[1],lim[2]),
                            f.lower=bounds[1],
                            f.upper=bounds[2])
    }

    ## crit is always pointing to the 0 just before the decision which
    ## is why we need a discrimination here
    if(lower.tail) {
        crit <- crit - 1
    }        

    crit
}


