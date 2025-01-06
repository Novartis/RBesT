#' Difference of mixture distributions
#'
#' Density, cumulative distribution function, quantile
#' function and random number generation for the difference of two mixture
#' distributions.
#'
#' @param mix1 first mixture density
#' @param mix2 second mixture density
#' @param x vector of values for which density values are computed
#' @param q vector of quantiles for which cumulative probabilities are computed
#' @param p vector of cumulative probabilities for which quantiles are computed
#' @param n size of random sample
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are P[X <= x], otherwise P[X > x].
#'
#' @details If \eqn{x_1 \sim f_1(x_1)}{x_1 ~ f_1(x_1)} and \eqn{x_2 \sim
#' f_2(x_2)}{x_2 ~ f_2(x)}, the density of the difference \eqn{d
#' \equiv x_1 - x_2}{d = x_1 - x_2} is given by
#'
#' \deqn{f_d(d) = \int f_1(u) \, f_2(u - d) \, du.}{f_d(d) = \int f_1(u) f_2(u - d)  du.}
#'
#' The cumulative distribution function equates to
#'
#' \deqn{F_d(d) = \int f_1(u) \, (1-F_2(u-d)) \, du.}{F_d(d) = \int f_1(u) (1-F_2(u-d)) du.}
#'
#' Both integrals are performed over the full support of the
#' densities and use the numerical integration function
#' \code{\link{integrate}}.
#'
# The quantile function is implemented using a gradient based
# minimization of the squared difference of the cumulative
# distribution function with the quantile requested, i.e. making use
# of the fact that the cumulative distribution increases
# monotonically.
#'
#' @return Respective density, quantile, cumulative density or random
#' numbers.
#'
#' @examples
#'
#' # 1. Difference between two beta distributions, i.e. Pr( mix1 - mix2 > 0)
#' mix1 <- mixbeta(c(1, 11, 4))
#' mix2 <- mixbeta(c(1, 8, 7))
#' pmixdiff(mix1, mix2, 0, FALSE)
#'
#' # Interval probability, i.e. Pr( 0.3 > mix1 - mix2 > 0)
#' pmixdiff(mix1, mix2, 0.3) - pmixdiff(mix1, mix2, 0)
#'
#' # 2. two distributions, one of them a mixture
#' m1 <- mixbeta(c(1, 30, 50))
#' m2 <- mixbeta(c(0.75, 20, 50), c(0.25, 1, 1))
#'
#' # random sample of difference
#' set.seed(23434)
#' rM <- rmixdiff(m1, m2, 1E4)
#'
#' # histogram of random numbers and exact density
#' hist(rM, prob = TRUE, new = TRUE, nclass = 40)
#' curve(dmixdiff(m1, m2, x), add = TRUE, n = 51)
#'
#' # threshold probabilities for difference, at 0 and 0.2
#' pmixdiff(m1, m2, 0)
#' mean(rM < 0)
#' pmixdiff(m1, m2, 0.2)
#' mean(rM < 0.2)
#'
#' # median of difference
#' mdn <- qmixdiff(m1, m2, 0.5)
#' mean(rM < mdn)
#'
#' # 95%-interval
#' qmixdiff(m1, m2, c(0.025, 0.975))
#' quantile(rM, c(0.025, 0.975))
#'
#' @name mixdiff
NULL

mixnormdiff <- function(mix1, mix2) {
  w <- outer(mix1[1, ], mix2[1, ], "*")
  mu <- outer(mix1[2, ], mix2[2, ], "-")
  v <- outer(mix1[3, ]^2, mix2[3, ]^2, "+")
  mixdiff <- rbind(w = as.vector(w), m = as.vector(mu), s = sqrt(as.vector(v)))
  class(mixdiff) <- c("normMix", "mix")
  if (!is.null(sigma(mix1))) {
    sigma(mixdiff) <- sigma(mix1)
  }
  likelihood(mixdiff) <- "normal"
  dlink(mixdiff) <- dlink(mix1)
  mixdiff
}

#' @rdname mixdiff
#' @export
dmixdiff <- function(mix1, mix2, x) {
  assert_that(!inherits(mix1, "mvnormMix") & !inherits(mix2, "mvnormMix"), msg = "Multivariate normal mixture density not supported.")
  if (inherits(mix1, "normMix") & inherits(mix2, "normMix")) {
    return(dmix(mixnormdiff(mix1, mix2), x = x))
  }
  interval <- support(mix1)
  if (!all(interval == support(mix2))) {
    warning("Support of variates mix1 and mix2 do not match.")
  }

  Nc1 <- ncol(mix1)
  Nc2 <- ncol(mix2)

  ## rescale integrand to ensure stability of integration when
  ## default 1E-4 tolerances are used
  scale <- 1E3
  lscale <- log(scale)

  ## in case we know that mix2 is only a single-component density,
  ## then we swap the integration order
  if (Nc2 == 1) {
    .dens <- function(sx) {
      integrate_density_log(function(x) lscale + dmix(mix1, x + sx, log = TRUE), mix2) / scale
    }
  } else {
    .dens <- function(sx) {
      integrate_density_log(function(x) lscale + dmix(mix2, x - sx, log = TRUE), mix1) / scale
    }
  }

  vapply(x, .dens, c(d = 0.1))
}


#' @rdname mixdiff
#' @export
pmixdiff <- function(mix1, mix2, q, lower.tail = TRUE) {
  assert_that(!inherits(mix1, "mvnormMix") & !inherits(mix2, "mvnormMix"), msg = "Multivariate normal mixture density not supported.")
  if (inherits(mix1, "normMix") & inherits(mix2, "normMix")) {
    return(pmix(mixnormdiff(mix1, mix2), q = q, lower.tail = lower.tail))
  }
  interval <- support(mix1)
  if (!all(interval == support(mix2))) {
    warning("Support of variates mix1 and mix2 do not match.")
  }

  Nc1 <- ncol(mix1)
  Nc2 <- ncol(mix2)

  ## should the second density only have a single component, then we
  ## take advantage of this
  if (Nc2 == 1) {
    .prob <- function(sx) {
      integrate_density_log(function(x) pmix(mix1, sx + x, lower.tail = TRUE, log.p = TRUE), mix2)
    }
  } else {
    .prob <- function(sx) {
      ## integrate_density_log(function(x) pmix(mix2, x-sx, lower.tail=FALSE, log.p=TRUE), mix1)
      1 - integrate_density_log(function(x) pmix(mix2, x - sx, lower.tail = TRUE, log.p = TRUE), mix1)
    }
  }

  p <- pmin(pmax(vapply(q, .prob, c(p = 0.1)), 0), 1)
  if (!lower.tail) {
    p <- 1 - p
  }
  return(p)
}

#' @rdname mixdiff
#' @export
qmixdiff <- function(mix1, mix2, p, lower.tail = TRUE) {
  assert_that(!inherits(mix1, "mvnormMix") & !inherits(mix2, "mvnormMix"), msg = "Multivariate normal mixture density not supported.")
  if (inherits(mix1, "normMix")) {
    return(qmix(mixnormdiff(mix1, mix2), p = p, lower.tail = lower.tail))
  }
  interval <- support(mix1)
  if (!all(interval == support(mix2))) {
    warning("Support of variates mix1 and mix2 do not match.")
  }
  assert_that(all(p >= 0 & p <= 1))
  assert_that(abs(sum(mix1["w", ]) - 1) < sqrt(.Machine$double.eps))
  assert_that(abs(sum(mix2["w", ]) - 1) < sqrt(.Machine$double.eps))
  ## first get the support of the mixture, i.e. the 99% CI of each
  ## mixture or lower, if the requested quantile is more in the
  ## tails
  plow <- min(c(p, (1 - p))) / 2
  ## plow <- if(log.p) min(c(0.01, exp(p), (1-exp(p)))) / 2 else min(c(0.01, p, (1-p))) / 2
  phigh <- 1 - plow
  qlow <- qmix(mix1, plow) - qmix(mix2, phigh)
  qhigh <- qmix(mix1, phigh) - qmix(mix2, plow)
  res <- rep.int(NA, length(p))
  for (i in seq_along(p)) {
    ## take advantage of the monotonicity of the CDF function such
    ## that we can use a gradient based method to find the root
    o <- optimise(function(x) {
      (pmixdiff(mix1, mix2, x, lower.tail = lower.tail) - p[i])^2
    }, c(qlow, qhigh))
    res[i] <- o$minimum
    if (o$objective > 1e-3) {
      ## in that case fall back to binary search which is more robust
      u <- uniroot(function(x) {
        pmixdiff(mix1, mix2, x, lower.tail = lower.tail) - p[i]
      }, c(qlow, qhigh))
      res[i] <- u$root
      if (u$estim.prec > 1E-3) {
        warning("Quantile ", p[i], " possibly imprecise.\nEstimated precision= ", u$estim.prec, ".\nRange = ", qlow, " to ", qhigh, "\n")
      }
    }
  }
  res
}

#' @rdname mixdiff
#' @export
rmixdiff <- function(mix1, mix2, n) {
  rmix(mix1, n) - rmix(mix2, n)
}
