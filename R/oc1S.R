#' Operating Characteristics for 1 Sample Design
#'
#' The \code{oc1S} function defines a 1 sample design (prior, sample
#' size, decision function) for the calculation of the frequency at
#' which the decision is evaluated to 1 conditional on assuming
#' known parameters.  A function is returned which performs the actual
#' operating characteristics calculations.
#'
#' @template args-boundary1S
#'
#' @details The \code{oc1S} function defines a 1 sample design and
#' returns a function which calculates its operating
#' characteristics. This is the frequency with which the decision
#' function is evaluated to 1 under the assumption of a given true
#' distribution of the data defined by a known parameter
#' \eqn{\theta}. The 1 sample design is defined by the prior, the
#' sample size and the decision function, \eqn{D(y)}. These uniquely
#' define the decision boundary, see
#' \code{\link{decision1S_boundary}}.
#' 
#' When calling the \code{oc1S} function, then internally the critical
#' value \eqn{y_c} (using \code{\link{decision1S_boundary}}) is
#' calculated and a function is returns which can be used to
#' calculated the desired frequency which is evaluated as
#'
#' \deqn{ F(y_c|\theta). }
#'
#' @return Returns a function with one argument \code{theta} which
#' calculates the frequency at which the decision function is
#' evaluated to 1 for the defined 1 sample design. Note that the
#' returned function takes vectors arguments.
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
#' # standard NI design
#' decA <- decision1S(1 - alpha, theta_ni, lower.tail=TRUE)
#' 
#' # double criterion design
#' # statistical significance (like NI design)
#' dec1 <- decision1S(1-alpha, theta_ni, lower.tail=TRUE)
#' # require mean to be at least as good as theta_c
#' dec2 <- decision1S(0.5, theta_c, lower.tail=TRUE)
#' # combination
#' decComb <- decision1S(c(1-alpha, 0.5), c(theta_ni, theta_c), lower.tail=TRUE)
#' 
#' theta_eval  <- c(theta_a, theta_c, theta_ni)
#' 
#' # evaluate different designs at two sample sizes
#' designA_n1  <- oc1S(flat_prior, n1, decA)
#' designA_nL  <- oc1S(flat_prior, nL, decA)
#' designC_n1  <- oc1S(flat_prior, n1, decComb)
#' designC_nL  <- oc1S(flat_prior, nL, decComb)
#'
#' # evaluate designs at the key log-HR of positive treatment (HR<1),
#' # the indecision point and the NI margin
#' 
#' designA_n1(theta_eval)
#' designA_nL(theta_eval)
#' designC_n1(theta_eval)
#' designC_nL(theta_eval)
#'
#' # to understand further the dual criterion design it is useful to
#' # evaluate the criterions separatley:
#' # statistical significance criterion to warrant NI...
#' designC1_nL <- oc1S(flat_prior, nL, dec1)
#' # ... or the clinically determined indifference point
#' designC2_nL <- oc1S(flat_prior, nL, dec2)
#'
#' designC1_nL(theta_eval)
#' designC2_nL(theta_eval)
#'
#' # see also ?decision1S_boundary to see which of the two criterions
#' # will drive the decision
#' 
#'
#' @export
oc1S <- function(prior, n, decision, ...) UseMethod("oc1S")
#' @export
oc1S.default <- function(prior, n, decision, ...) "Unknown density"

#' @templateVar fun oc1S
#' @template design1S-binomial
#' @export
oc1S.betaMix <- function(prior, n, decision, ...) {

    crit <- decision1S_boundary(prior, n, decision)
    lower.tail <- attr(decision, "lower.tail")
    
    design_fun <- function(theta) {
        if(missing(theta)) {
            deprecated("Use of no argument", "decision1S_boundary")
            return(crit)
        }
        pbinom(crit, n, theta, lower.tail=lower.tail)
    }
    design_fun 
}


#' @templateVar fun oc1S
#' @template design1S-normal
#' @export
oc1S.normMix <- function(prior, n, decision, sigma, eps=1e-6, ...) {
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

    crit <- decision1S_boundary(prior, n, decision, sigma, eps)
    
    ## check where the decision is 1, i.e. left or right
    lower.tail <- attr(decision, "lower.tail")

    design_fun <- function(theta) {
        if(missing(theta)) {
            deprecated("Use of no argument", "decision1S_boundary")
            return(crit)
        }
        pnorm(crit, theta, sd_samp, lower.tail=lower.tail)
    }
    design_fun 
}

#' @templateVar fun oc1S
#' @template design1S-poisson
#' @export
oc1S.gammaMix <- function(prior, n, decision, eps=1e-6, ...) {

    crit <- decision1S_boundary(prior, n, decision, eps)
    lower.tail <- attr(decision, "lower.tail")

    design_fun <- function(theta) {
        if(missing(theta)) {
            deprecated("Use of no argument", "decision1S_boundary")
            return(crit)
        }
        ppois(crit, n*theta, lower.tail=lower.tail)
    }
    design_fun 
}
