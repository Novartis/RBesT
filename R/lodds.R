#' @rdname lodds
#' @name lodds
#'
#' @title Logit (log-odds) and inverse-logit function.
#'
#' @description
#' Calculates the logit (log-odds) and inverse-logit.
#'
#' @param mu A numeric object with probabilies, with values in the in
#' the range \eqn{[0,1]}. Missing values (`NA`s) are allowed.
#' @param eta A numeric object with log-odds values, with values in
#' the range \eqn{[-\infty,\infty]}. Missing values (`NA`s) are allowed.
#'
#' @details Values of mu equal to 0 or 1 will return \eqn{-\infty} or
#'   \eqn{\infty} respectively.
#'
#' @return A numeric object of the same type as mu and eta containing
#' the logits or inverse logit of the input values.  The logit and
#' inverse transformation equates to
#'
#' \deqn{\text{logit}(\mu) = \log(\mu/(1-\mu))}{logit(\mu) = log(\mu/(1-\mu))}
#' \deqn{\text{logit}^{-1}(\eta)= \exp(\eta)/(1 + \exp(\eta)).}{logit^-1(\eta) = exp(\eta)/(1 + exp(\eta)).}
#'
#' @examples
#' logit(0.2)
#' inv_logit(-1.386)
#'
NULL
bin <- binomial()

#' @export
#' @rdname lodds
logit <- bin$linkfun

#' @export
#' @rdname lodds
inv_logit <- bin$linkinv
