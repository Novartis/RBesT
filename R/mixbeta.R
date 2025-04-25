#' @name mixbeta
#'
#' @title Beta Mixture Density
#'
#' @description The Beta mixture density and auxilary functions.
#'
#' @param ... List of mixture components.
#' @param param Determines how the parameters in the list
#' are interpreted. See details.
#' @param m Vector of means of beta mixture components.
#' @param s Vector of standard deviations of beta mixture components.
#' @param n Vector of number of observations.
#' @param drop Delete the dimensions of an array which have only one level.
#' @param object Beta mixture object.
#' @param probs Quantiles reported by the `summary` function.
#'
#' @details Each entry in the `...` argument list is expected to
#' be a triplet of numbers which defines the weight \eqn{w_k}, first
#' and second parameter of the mixture component \eqn{k}. A triplet
#' can optionally be named which will be used appropriately.
#'
#' The first and second parameter can be given in different
#' parametrizations which is set by the `param` option:
#' \describe{
#' \item{ab}{Natural parametrization of Beta density (`a`=shape1 and `b`=shape2). Default. }
#' \item{ms}{Mean and standard deviation, \eqn{m=a/(a+b)} and \eqn{s=\sqrt{\frac{m(1-m)}{1+n}}}, where \eqn{n=a+b} is the number of observations. Note that \eqn{s} must be less than \eqn{\sqrt{m(1-m)}}.}
#' \item{mn}{Mean and number of observations, \eqn{n=a+b}.}
#' }
#'
#' @family mixdist
#'
#' @return `mixbeta` returns a beta mixture with the specified mixture components. `ms2beta` and
#' `mn2beta` return the equivalent natural `a` and `b` parametrization given parameters `m`,
#' `s`, or `n`.
#'
#' @examples
#' ## a beta mixture
#' bm <- mixbeta(rob = c(0.2, 2, 10), inf = c(0.4, 10, 100), inf2 = c(0.4, 30, 80))
#'
#' # mean/standard deviation parametrization
#' bm2 <- mixbeta(rob = c(0.2, 0.3, 0.2), inf = c(0.8, 0.4, 0.01), param = "ms")
#'
#' # mean/observations parametrization
#' bm3 <- mixbeta(rob = c(0.2, 0.3, 5), inf = c(0.8, 0.4, 30), param = "mn")
#'
#' # even mixed is possible
#' bm4 <- mixbeta(rob = c(0.2, mn2beta(0.3, 5)), inf = c(0.8, ms2beta(0.4, 0.1)))
#'
#' # print methods are defined
#' bm4
#' print(bm4)
#'
NULL

#' @rdname mixbeta
#' @export
mixbeta <- function(..., param = c("ab", "ms", "mn")) {
  mix <- mixdist3(...)
  assert_matrix(mix, nrows = 3, any.missing = FALSE)
  param <- match.arg(param)
  mix[c(2, 3), ] <- switch(
    param,
    ab = mix[c(2, 3), ],
    ms = t(ms2beta(mix[2, ], mix[3, ], FALSE)),
    mn = t(mn2beta(mix[2, ], mix[3, ], FALSE))
  )
  rownames(mix) <- c("w", "a", "b")
  assert_that(all(mix["a", ] >= 0))
  assert_that(all(mix["b", ] >= 0))
  class(mix) <- c("betaMix", "mix")
  likelihood(mix) <- "binomial"
  mix
}

#' @rdname mixbeta
#' @export
ms2beta <- function(m, s, drop = TRUE) {
  n <- m * (1 - m) / s^2 - 1
  assert_that(all(n >= 0))
  ab <- cbind(a = n * m, b = n * (1 - m))
  if (drop) ab <- drop(ab)
  ab
}

#' @rdname mixbeta
#' @export
mn2beta <- function(m, n, drop = TRUE) {
  assert_that(all(n >= 0))
  ab <- cbind(a = n * m, b = n * (1 - m))
  if (drop) ab <- drop(ab)
  ab
}

#' @rdname mixbeta
#' @method print betaMix
#' @param x The mixture to print
#' @export
print.betaMix <- function(x, ...) {
  cat("Univariate beta mixture\n")
  NextMethod()
}

#' @rdname mixbeta
#' @method print betaBinomialMix
#' @param x The mixture to print
#' @export
print.betaBinomialMix <- function(x, ...) {
  cat("Univariate beta binomial mixture\nn = ", attr(x, "n"), "\n", sep = "")
  NextMethod()
}

#' @rdname mixbeta
#' @method summary betaMix
#' @export
summary.betaMix <- function(object, probs = c(0.025, 0.5, 0.975), ...) {
  p <- object[1, ]
  a <- object[2, ]
  b <- object[3, ]
  m <- a / (a + b)
  v <- m * (1 - m) / (a + b + 1)
  ## calculate mean of the second moment
  m2 <- v + m^2
  ## from this we can get the mean and variance of the mixture
  mmix <- sum(p * m)
  vmix <- sum(p * (m2 - (mmix)^2))
  q <- c()
  if (length(probs) != 0) {
    q <- qmix.betaMix(object, p = probs)
    names(q) <- paste(format(probs * 100, digits = 2), "%", sep = "")
  }
  c(mean = mmix, sd = sqrt(vmix), q)
}

#' @rdname mixbeta
#' @method summary betaBinomialMix
#' @export
summary.betaBinomialMix <- function(object, probs = c(0.025, 0.5, 0.975), ...) {
  n <- attr(object, "n")
  p <- object[1, ]
  a <- object[2, ]
  b <- object[3, ]
  m <- n * a / (a + b)
  v <- n * a * b * (a + b + n) / ((a + b)^2 * (a + b + 1))
  ## calculate mean of the second moment
  m2 <- v + m^2
  ## from this we can get the mean and variance of the mixture
  mmix <- sum(p * m)
  vmix <- sum(p * (m2 - (mmix)^2))
  q <- qmix.betaBinomialMix(object, p = probs)
  if (length(q) != 0) {
    names(q) <- paste(format(probs * 100, digits = 2), "%", sep = "")
  }
  c(mean = mmix, sd = sqrt(vmix), q)
}
