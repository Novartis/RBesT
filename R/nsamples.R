#' Return the number of posterior samples
#'
#' @param object fitted model object
#' @template args-dots-ignored
#'
#'
#' @template example-start
#' @examples
#'
#' set.seed(34563)
#' map_AS <- gMAP(cbind(r, n - r) ~ 1 | study,
#'   family = binomial,
#'   data = AS,
#'   tau.dist = "HalfNormal", tau.prior = 1,
#'   beta.prior = 2
#' )
#'
#' nsamples(map_AS)
#'
#' @template example-stop
#'
#' @method nsamples gMAP
#' @aliases nsamples
#' @export
nsamples.gMAP <- function(object, ...) {
  return(
    object$fit@sim$chains *
      (object$fit@sim$iter - object$fit@sim$warmup)
  )
}
