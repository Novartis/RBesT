#' Support of Distributions
#'
#' Returns the support of a distribution.
#'
#' @param mix Mixture distribution.
#'
#' @keywords internal
support <- function(mix) UseMethod("support")
#' @export
support.default <- function(mix) "Unknown mixture"
#' @export
support.betaMix <- function(mix) mixlink(mix, c(0,1))
#' @export
support.gammaMix <- function(mix) mixlink(mix, c(0,Inf))
#' @export
support.normMix <- function(mix) mixlink(mix, c(-Inf,Inf))
