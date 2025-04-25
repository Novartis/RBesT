#' @name likelihood
#'
#' @title Read and Set Likelihood to the Corresponding Conjugate Prior
#'
#' @description Read and set the likelihood distribution corresponding to the conjugate prior distribution.
#'
#' @param mix Prior mixture distribution.
#' @param value New likelihood. **Should** only be changed for Gamma priors as these are supported
#' with either Poisson (`value="poisson"`) or Exponential (`value="exp"`) likelihoods.
#'
#' @details
#' If the prior and posterior distributions are in the same family, then the prior distribution
#' is called a conjugate prior for the likelihood function.
#'
#' @template conjugate_pairs
#'
#' @examples
#'
#' # Gamma mixture
#' gmix <- mixgamma(c(0.3, 20, 4), c(0.7, 50, 10))
#'
#' # read out conjugate partner
#' likelihood(gmix)
#'
#' ess(gmix)
#'
#' # set conjugate partner
#' likelihood(gmix) <- "exp"
#'
#' # ... which changes the interpretation of the mixture
#' ess(gmix)
NULL

#' @rdname likelihood
#' @export
likelihood <- function(mix) {
  likelihood <- attr(mix, "likelihood")
  check_choice(likelihood, c("poisson", "exp", "normal", "binomial"))
  likelihood
}

#' @rdname likelihood
#' @export
"likelihood<-" <- function(mix, value) {
  check_choice(value, c("poisson", "exp", "normal", "binomial"))
  attr(mix, "likelihood") <- value
  mix
}
