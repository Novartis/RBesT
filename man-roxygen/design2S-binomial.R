#' @describeIn <%= fun %> Applies for binomial model with a mixture
#' beta prior. The calculations use exact expressions.  If the
#' optional argument \code{eps} is defined, then an approximate method
#' is used which limits the search for the decision boundary to the
#' region of \code{1-eps} probability mass. This is useful for designs
#' with large sample sizes where an exact approach is very costly to
#' calculate.
