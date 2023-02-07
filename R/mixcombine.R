#' Combine Mixture Distributions
#'
#' Combining mixture distributions of the same class to form a new mixture.
#'
#' @param ... arbitrary number of mixtures with same distributional class.
#' Each component with values of mixture weight and model parameters.
#' @param weight relative weight for each component in new mixture
#' distribution. The vector must be of the same length as input
#' mixtures components. The default value gives equal weight to each
#' component.
#' @param rescale boolean value indicates if the weights are
#' rescaled to sum to 1.
#'
#' @details Combines mixtures of the same class of random
#' variable to form a new mixture distribution.
#'
#' @return
#' A R-object with the new mixture distribution.
#' @family mixdist
#' @seealso \code{\link{robustify}}
#'
#' @examples
#' # beta with two informative components
#' bm <- mixbeta(inf=c(0.5, 10, 100), inf2=c(0.5, 30, 80))
#'
#' # robustified with mixcombine, i.e. a 10% uninformative part added
#' unif <- mixbeta(rob=c(1,1,1))
#' mixcombine(bm, unif, weight=c(9, 1))
#'
#' @export
mixcombine <- function(..., weight, rescale=TRUE) {
    comp <- list(...)
    ## ensure that the resulting object is a mixture object only
    cl <- grep("mix$", class(comp[[1]]), ignore.case=TRUE, value=TRUE)
    dl <- dlink(comp[[1]])
    lik <- likelihood(comp[[1]])
    assert_that(all(sapply(comp, inherits, "mix")), msg="All components must be mixture objects.")
    assert_that(all(sapply(comp, likelihood) == lik), msg="All components must have the same likelihood set.")
    mix <- do.call(cbind, comp)
    if(!missing(weight)) {
        assert_that(length(weight) == length(comp))
        mix[1,] <-  mix[1,] * rep(weight, times=sapply(comp, ncol))
    }
    if(rescale)
        mix[1,] <-  mix[1,] / sum(mix[1,])
    class(mix) <- cl
    dlink(mix) <- dl
    likelihood(mix) <- lik
    if("normMix" %in% cl) sigma(mix) <- sigma(comp[[1]])
    mix
}
