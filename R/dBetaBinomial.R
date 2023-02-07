#' Beta-Binomial Probabilities
#'
#' @param r,n number of successes (responders) out of n
#' @param a,b parameters of the Beta distribution for response probability
#'
#' @details
#' r,n,a,b can be scalar or vectors. If vectors are used, they must be of the same length
#'
#' @keywords internal
#'
dBetaBinomial <- function(r, n, a, b, log=FALSE) {
    assert_integerish(n, lower=0L, any.missing=FALSE)
    assert_integerish(r, lower=0L, upper=max(n), any.missing=FALSE)
    p <- lgamma(n+1)-lgamma(r+1)-lgamma(n-r+1)+lgamma(a+b)-lgamma(a)-lgamma(b)+lgamma(a+r)+lgamma(b+n-r)-lgamma(a+b+n);
    if(log)
        return(p)
    exp(p);
}


##pBetaBinomial <- function(r, n, a, b, lower.tail=TRUE, log.p=FALSE) {
##    return(.pBetaBinomial(r, n, a, b, lower.tail, log.p))
##}
