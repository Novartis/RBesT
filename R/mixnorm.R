#' @name mixnorm
#' @aliases sigma
#'
#' @title Normal Mixture Density
#'
#' @description The normal mixture density and auxiliary functions.
#'
#' @param ... List of mixture components.
#' @param sigma Reference scale.
#' @param param Determines how the parameters in the list are
#' interpreted. See details.
#' @param m Vector of means
#' @param n Vector of sample sizes.
#' @param drop Delete the dimensions of an array which have only one level.
#' @param object Normal mixture object.
#' @param probs Quantiles reported by the \code{summary} function.
#' @param value New value of the reference scale \code{sigma}.
#'
#' @details Each entry in the \code{...} argument list is expected to
#' be a triplet of numbers which defines the weight \eqn{w_k}, first
#' and second parameter of the mixture component \eqn{k}. A triplet
#' can optionally be named which will be used appropriately.
#'
#' The first and second parameter can be given in different
#' parametrizations which is set by the \code{param} option:
#' \describe{
#' \item{ms}{Mean and standard deviation. Default.}
#' \item{mn}{Mean and number of observations. \code{n} determines \code{s} via the relation \eqn{s=\sigma/\sqrt{n}} with \eqn{\sigma} being the fixed reference scale.}
#' }
#'
#' The reference scale \eqn{\sigma} is the fixed standard deviation in
#' the one-parameter normal-normal model (observation standard
#' deviation). The function \code{sigma} can be used to query the
#' reference scale and may also be used to assign a new reference
#' scale, see examples below. In case the \code{sigma} is not
#' specified, the user has to supply \code{sigma} as argument to
#' functions which require a reference scale.
#'
#' @family mixdist
#'
#' @return Returns a normal mixture with the specified mixture
#' components. \code{mn2norm} returns the mean and standard deviation
#' given a mean and sample size parametrization.
#'
#' @examples
#'
#' nm <- mixnorm(rob = c(0.2, 0, 2), inf = c(0.8, 2, 2), sigma = 5)
#'
#' print(nm)
#' summary(nm)
#' plot(nm)
#'
#' set.seed(57845)
#' mixSamp <- rmix(nm, 500)
#' plot(nm, samp = mixSamp)
#'
#' # support defined by quantiles
#' qmix(nm, c(0.01, 0.99))
#'
#' # density function
#' dmix(nm, seq(-5, 5, by = 2))
#'
#' # distribution function
#' pmix(nm, seq(-5, 5, by = 2))
#'
#' # the reference scale can be changed (it determines the ESS)
#' ess(nm)
#'
#' sigma(nm) <- 10
#' ess(nm)
NULL

#' @rdname mixnorm
#' @export
mixnorm <- function(..., sigma, param = c("ms", "mn")) {
  mix <- mixdist3(...)
  assert_matrix(mix, nrows = 3, any.missing = FALSE)
  param <- match.arg(param)
  mix[c(2, 3), ] <- switch(param,
    ms = mix[c(2, 3), ],
    mn = t(mn2norm(mix[2, ], mix[3, ], sigma, FALSE))
  )
  rownames(mix) <- c("w", "m", "s")
  assert_that(all(mix["s", ] > 0))
  if (!missing(sigma)) {
    assert_number(sigma, lower = 0)
    attr(mix, "sigma") <- sigma
  }
  class(mix) <- c("normMix", "mix")
  likelihood(mix) <- "normal"
  mix
}

#' @rdname mixnorm
#' @export
mn2norm <- function(m, n, sigma, drop = TRUE) {
  assert_number(sigma, lower = 0)
  sigma_n <- sigma / sqrt(n)
  ms <- cbind(m = m, s = sigma_n)
  if (drop) ms <- drop(ms)
  ms
}

#' @rdname mixnorm
#' @method print normMix
#' @param x The mixture to print
#' @export
print.normMix <- function(x, ...) {
  cat("Univariate normal mixture\n")
  if (!is.null(sigma(x))) {
    cat("Reference scale: ", sigma(x), "\n", sep = "")
  }
  NextMethod()
}

#' @rdname mixnorm
#' @method summary normMix
#' @export
summary.normMix <- function(object, probs = c(0.025, 0.5, 0.975), ...) {
  p <- object[1, ]
  m <- object[2, ]
  v <- object[3, ]^2
  ## calculate mean of the second moment
  m2 <- v + m^2
  ## from this we can get the mean and variance of the mixture
  mmix <- sum(p * m)
  vmix <- sum(p * (m2 - (mmix)^2))
  q <- c()
  if (length(probs) != 0) {
    q <- qmix.normMix(object, p = probs)
    names(q) <- paste(format(probs * 100, digits = 2), "%", sep = "")
  }
  c(mean = mmix, sd = sqrt(vmix), q)
}

#' @rdname mixnorm
#' @method sigma normMix
#' @export
#' @export sigma
#' @rawNamespace importFrom(stats, sigma)
sigma.normMix <- function(object, ...) {
  attr(object, "sigma")
}

#' @describeIn mixnorm Allows to assign a new reference scale \code{sigma}.
#' @export
"sigma<-" <- function(object, value) {
  assert_number(value, lower = 0, null.ok = TRUE)
  attr(object, "sigma") <- unname(value)
  object
}
