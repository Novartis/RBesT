#' Effective Sample Size for a Conjugate Prior
#'
#' Calculates the Effective Sample Size (ESS) for a mixture prior. The
#' ESS indicates how many experimental units the prior is roughly
#' equivalent to.
#'
#' @param mix Prior (mixture of conjugate distributions).
#' @param method Selects the used method. Can be either \code{elir}
#'     (default), \code{moment} or \code{morita}.
#' @param s For \code{morita} method large constant to ensure that the
#'     prior scaled by this value is vague (default 100); see Morita
#'     et al. (2008) for details.
#' @param eps Probability mass left out from the numerical integration
#'     of the expected information for the Poisson-Gamma case of
#'     Morita method (defaults to 1E-4).
#' @param sigma reference scale.
#' @param family defines data likelihood and link function
#'     (\code{binomial}, \code{gaussian}, or \code{poisson}).
#' @param ... Optional arguments applicable to specific methods.
#'
#' @details The ESS is calculated using either the expected local
#'     information ratio (elir) \emph{Neuenschwander et
#'     al. (2020)}, the moments approach or the method by
#'     \emph{Morita et al. (2008)}.
#'
#' The elir approach measures effective sample size in terms of the
#' average curvature of the prior in relation to the Fisher
#' information. Informally this corresponds to the average peakiness
#' of the prior in relation to the information content of a single
#' observation. The elir approach is the only ESS which fulfills
#' predictive consistency. The predictive consistency of the ESS
#' requires that the ESS of a prior is consistent when considering an
#' averaged posterior ESS of additional data distributed according to
#' the predictive distribution of the prior. The expectation of the
#' posterior ESS is taken wrt to the prior predictive distribution and
#' the averaged posterior ESS corresponds to the sum of the prior ESS
#' and the number of forward simulated data items. The elir approach
#' results in ESS estimates which are neither conservative nor liberal
#' whereas the moments method yields conservative and the morita
#' method liberal results. See the example section for a demonstration
#' of predictive consistency.
#'
#' For the moments method the mean and standard deviation of the
#' mixture are calculated and then approximated by the conjugate
#' distribution with the same mean and standard deviation. For
#' conjugate distributions, the ESS is well defined. See the examples
#' for a step-wise calculation in the beta mixture case.
#'
#' The Morita method used here evaluates the mixture prior at the mode
#' instead of the mean as proposed originally by Morita. The method
#' may lead to very optimistic ESS values, especially if the mixture
#' contains many components. The calculation of the Morita approach
#' here follows the approach presented in Neuenschwander B. et all
#' (2019) which avoids the need for a minimization and does not
#' restrict the ESS to be an integer.
#'
#' The arguments \code{sigma} and \code{family} are specific for
#' normal mixture densities. These specify the sampling standard
#' deviation for a \code{gaussian} family (the default) while also
#' allowing to consider the ESS of standard one-parameter exponential
#' families, i.e. \code{binomial} or \code{poisson}. The function
#' supports non-gaussian families with unit dispersion only.
#'
#' @return Returns the ESS of the prior as floating point number.
#'
#' @template conjugate_pairs
#'
#' @references Morita S, Thall PF, Mueller P.  Determining the
#'     effective sample size of a parametric prior.  \emph{Biometrics}
#'     2008;64(2):595-602.
#'
#' @references Neuenschwander B., Weber S., Schmidli H., O’Hagan
#'     A. (2020). Predictively consistent prior effective sample
#'     sizes. \emph{Biometrics}, 76(2),
#'     578–587. https://doi.org/10.1111/biom.13252
#'
#' @example inst/examples/ess.R
#'
#' @export
ess <- function(mix, method = c("elir", "moment", "morita"), ...) UseMethod("ess")
#' @export
ess.default <- function(mix, method = c("elir", "moment", "morita"), ...) stop("Unknown density")


calc_loc <- function(mix, loc = c("mode", "median", "mean")) {
  loc <- match.arg(loc)
  if (loc == "mode") {
    tol <- .Machine$double.eps^0.25
    locEst <- mixmode(mix)

    if (length(attr(locEst, "modes")) > 1) {
      warning("Detected multiple modes.\nThe ESS is determined for the largest mode, but ESS concept is ill-defined for multi-modal distributions.")
    } else {
      attr(locEst, "modes") <- NULL
    }
  }
  if (loc == "median") {
    locEst <- qmix(mix, 0.5)
  }
  if (loc == "mean") {
    locEst <- summary(mix, NULL)["mean"]
  }
  names(locEst) <- NULL

  return(unname(locEst))
}
## function to calculate mixture info of arbitrary density; input
## needed is the density function, gradient and hessian of the log
## density with respect to x (data)
mixInfo <- function(mix, x, dens, gradl, hessl) {
  p <- mix[1, ]
  a <- mix[2, ]
  b <- mix[3, ]
  lp <- log(p)
  ldensComp <- dens(x, a, b, log = TRUE)
  ldensMix <- matrixStats::logSumExp(lp + ldensComp)
  lwdensComp <- lp + ldensComp - ldensMix
  dgl <- gradl(x, a, b)
  dhl <- (hessl(x, a, b) + dgl^2)
  ## attempt numerically more stable log calculations if possible,
  ## i.e. if all signs are the same
  if (all(!is.na(dgl)) && (all(dgl < 0) || all(dgl > 0))) {
    gsum <- exp(2 * matrixStats::logSumExp(lwdensComp + log(abs(dgl))))
  } else {
    gsum <- (sum(exp(lwdensComp) * dgl))^2
  }
  if (all(!is.na(dhl)) && (all(dhl < 0) || all(dhl > 0))) {
    hsum <- sign(dhl[1]) * exp(matrixStats::logSumExp(lwdensComp + log(abs(dhl))))
  } else {
    hsum <- (sum(exp(lwdensComp) * dhl))
  }
  gsum - hsum
}

## local information ratio (which we integrate over the prior)
lir <- function(mix, info, fisher_inverse) {
  fn <- function(x) {
    info(mix, x) * fisher_inverse(x)
  }
  Vectorize(fn)
}

weighted_lir <- function(mix, info, fisher_inverse) {
  fn <- function(x) {
    dmix(mix, x) * info(mix, x) * fisher_inverse(x)
  }
  Vectorize(fn)
}

## not used ATM as there have been numerical issues
weighted_lir_link <- function(mix, info, fisher_inverse, link) {
  dlink(mix) <- link_map[[link]]
  fn <- function(x) {
    x_orig <- mixinvlink(mix, x)
    dmix(mix, x) * info(mix, x_orig) * fisher_inverse(x_orig)
  }
  Vectorize(fn)
}

## function to calculate the gradient of the log mixture
## mixLogGrad <- function(mix, x, dens, gradl) {
##     p <- mix[1,]
##     a <- mix[2,]
##     b <- mix[3,]
##     densMix <- dmix(mix,x)
##     densComp <- dens(x, a, b)
##     dgl <- gradl(x,a,b)
##     sum(p*densComp*dgl)/densMix
## }

# prior effective sample size ESS for Beta-mixture priors
# based on
# Morita, Thall, Mueller (MTM) 2008 Biometrics
# only difference: evaluated at mode of prior rather than at mean; and the flattened prior are derived with respect to the scale of 1 instead of being relative to the input scale
#' @describeIn ess ESS for beta mixtures.
#' @export
ess.betaMix <- function(mix, method = c("elir", "moment", "morita"), ..., s = 100) {
  method <- match.arg(method)

  call_arg_names <- names(match.call())
  family_arg_is_set <- "family" %in% call_arg_names
  assert_that(!family_arg_is_set, msg = "Argument family is only supported for normal mixtures.")

  if (method == "elir") {
    if (!test_numeric(mix[2, ], lower = 1, finite = TRUE, any.missing = FALSE) ||
      !test_numeric(mix[3, ], lower = 1, finite = TRUE, any.missing = FALSE)) {
      stop(
        "At least one parameter of the beta mixtures is less than 1.\n",
        "This leads to an ill-defined elir ess since the defining integral diverges.\n",
        "Consider constraining all parameters to be greater than 1 (use constrain_gt1=TRUE argument for EM fitting functions)."
      )
    }
    elir <- integrate_density(lir(mix, betaMixInfo, bernoulliFisherInfo_inverse), mix)
    if (elir < 0) {
      warning("Negative ESS elir found indicating unstable integration of the elir ratio.\nConsider estimating the ESS elir on the logit scale for the respective transformed density and use the family=binomial argument.")
    }
    return(elir)
  }

  ## simple and conservative moment matching
  if (method == "moment") {
    smix <- summary(mix)
    res <- sum(ms2beta(smix["mean"], smix["sd"]))
    names(res) <- NULL
    return(res)
  }

  locEst <- calc_loc(mix, "mode")

  deriv2.prior <- betaMixInfo(mix, locEst)

  ## alpha and beta of "flattened" priors
  alphaP0 <- locEst / s
  betaP0 <- (1 - locEst) / s

  info.prior0 <- betaInfo(locEst, alphaP0, betaP0)

  ## MTM paper would follow this:
  ## priorN <- sum(mix[1,] * (alphaP + betaP))
  ## alphaP0 <- mode * priorN / s
  ## betaP0 <- (1-mode) * priorN / s

  ## we warn if any of the mixture components has a scale (n) which
  ## is less than 10/s such that the
  if (any(rowSums(mix[2:3, , drop = FALSE]) < 10 / s)) {
    warning("Some of the mixture components have a scale which is large compared to the rescaling factor s. Consider increasing s.")
  }

  pd0 <- dmix(preddist(mix, n = 1), 0)
  Einfo <- binomialInfo(0, locEst, 1) * pd0 + binomialInfo(1, locEst, 1) * (1 - pd0)

  ## Eq. 9 of Neuenschwander et al. (2019)
  return(unname((deriv2.prior - info.prior0) / Einfo))
}

## derivative of a single log-beta
betaLogGrad <- function(x, a, b) {
  -(b - 1) / (1 - x) + (a - 1) / x
}

## second derivative of a single log-beta
betaLogHess <- function(x, a, b) {
  -(b - 1) / (x - 1)^2 - (a - 1) / x^2
}

betaMixInfo <- function(mix, x) {
  x <- pmin(pmax(x, .Machine$double.eps), 1.0 - .Machine$double.eps)
  mixInfo(mix, x, dbeta, betaLogGrad, betaLogHess)
}

## info metric for a single beta, i.e. negative second derivative of log beta
betaInfo <- function(x, a, b) {
  x <- pmin(pmax(x, .Machine$double.eps), 1.0 - .Machine$double.eps)
  -betaLogHess(x, a, b)
}

## 1/i_F(x): The inverse of the fisher information for a Bernoulli
## experiment (binomial with n=1)
bernoulliFisherInfo_inverse <- function(x) {
  x - x^2
}

## info metric for a binomial, second derivative wrt to theta of the
## log binomial
binomialInfo <- function(r, theta, n) {
  r / theta^2 + (n - r) / (1 - theta)^2
}


#' @describeIn ess ESS for gamma mixtures.
#' @export
ess.gammaMix <- function(mix, method = c("elir", "moment", "morita"), ..., s = 100, eps = 1E-4) {
  method <- match.arg(method)

  call_arg_names <- names(match.call())
  family_arg_is_set <- "family" %in% call_arg_names
  assert_that(!family_arg_is_set, msg = "Argument family is only supported for normal mixtures.")

  lik <- likelihood(mix)

  if (method == "elir") {
    if (lik == "poisson") {
      return(integrate_density(lir(mix, gammaMixInfo, poissonFisherInfo_inverse), mix))
    }
    if (lik == "exp") {
      return(integrate_density(lir(mix, gammaMixInfo, expFisherInfo_inverse), mix))
    }
  }

  ## simple and conservative moment matching
  if (method == "moment") {
    smix <- summary(mix)
    coef <- ms2gamma(smix["mean"], smix["sd"])
    names(coef) <- NULL
    if (lik == "poisson") {
      return(unname(coef[2]))
    }
    if (lik == "exp") {
      return(unname(coef[1]))
    }
    stop("Unkown likelihood")
  }

  ## Morita method
  locEst <- calc_loc(mix, "mode")

  deriv2.prior <- gammaMixInfo(mix, locEst)

  if (lik == "poisson") {
    meanPrior <- summary(mix)["mean"]
    names(meanPrior) <- NULL
    priorN <- mix[3, , drop = FALSE]

    info.prior0 <- gammaInfo(locEst, locEst / s, 1 / s)

    ## E_Y1 ( i_F ) using numerical integration
    pred_pmf <- preddist(mix, n = 1)
    lim <- qmix(pred_pmf, c(eps / 2, 1 - eps / 2))
    y1 <- seq(lim[1], lim[2])
    Einfo <- sum(dmix(pred_pmf, y1) * poissonInfo(y1, locEst))
  }

  if (lik == "exp") {
    priorN <- mix[2, , drop = FALSE]

    info.prior0 <- gammaInfo(locEst, 1 / s, 1 / (s * locEst))

    ## E_Y1 ( i_F ) (the i_F does not depend on the data)
    Einfo <- expInfo(1, locEst)
  }

  if (any(priorN < 10 / s)) {
    warning("Some of the mixture components have a scale which is large compared to the rescaling factor s. Consider increasing s.")
  }

  return(unname((deriv2.prior - info.prior0) / Einfo))
}

## respective functions for the gamma distribution
gammaLogGrad <- function(x, a, b) {
  (a - 1) / x - b
}
gammaLogHess <- function(x, a, b) {
  -(a - 1) / x^2
}
gammaInfo <- function(x, a, b) {
  -gammaLogHess(x, a, b)
}
gammaMixInfo <- function(mix, x) {
  mixInfo(mix, x, dgamma, gammaLogGrad, gammaLogHess)
}
poissonFisherInfo_inverse <- function(x) {
  x
}
expFisherInfo_inverse <- function(x) {
  x^2
}
poissonInfo <- function(y, theta) {
  y / theta^2
}
expInfo <- function(y, theta) {
  1 / theta^2
}


#' @describeIn ess ESS for normal mixtures.
#' @export
ess.normMix <- function(mix, method = c("elir", "moment", "morita"), ..., family = gaussian, sigma, s = 100) {
  method <- match.arg(method)

  if (is.character(family)) {
    family <- get(family, mode = "function", envir = parent.frame())
  }
  if (is.function(family)) {
    family <- family()
  }
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }

  is_gaussian_family <- family$family == "gaussian"
  if (!is_gaussian_family) {
    assert_that(family$dispersion == 1, msg = "Only dispersion unity is supported for non-gaussian families.")
  }

  normTransformedFisherInfo_inverse <- make_normTransformedFisherInfo_inverse(family)

  if (!is_gaussian_family) {
    sigma <- 1
  } else {
    if (missing(sigma)) {
      sigma <- RBesT::sigma(mix)
      message("Using default prior reference scale ", sigma)
    }
  }
  assert_number(sigma, lower = 0)
  tauSq <- sigma^2

  ## note: sigma is reassigned the sd's of each mixture component
  mu <- mix[2, ]
  sigma <- mix[3, ]
  sigmaSq <- sigma^2

  if (method == "elir") {
    elir <- integrate_density(lir(mix, normMixInfo, normTransformedFisherInfo_inverse), mix)
    if (is_gaussian_family) {
      ## in this case we have to account for the non-unity scale
      ## of the Gaussian likelihood
      elir <- elir * tauSq
    }
    return(elir)
  }


  ## simple and conservative moment matching compared to the
  ## expected fisher information over the prior parameter space
  if (method == "moment") {
    smix <- summary(mix)
    if (is_gaussian_family && family$link == "identity") {
      expected_info <- 1
    } else {
      expected_info <- integrate_density(normTransformedFisherInfo_inverse, mix)
    }
    if (is_gaussian_family) {
      expected_info <- tauSq * expected_info
    }
    res <- expected_info / smix["sd"]^2
    return(unname(res))
  }

  locEst <- calc_loc(mix, "mode")

  deriv2.prior <- normMixInfo(mix, locEst)

  ## "flattened" priors
  muP0 <- locEst
  sigmaP0Sq <- s * max(sigmaSq)

  info.prior0 <- normInfo(locEst, muP0, sqrt(sigmaP0Sq))

  ## info at mode
  mode_info <- 1 / normTransformedFisherInfo_inverse(locEst)

  if (is_gaussian_family) {
    ## in this case we have to account for the non-unity scale
    ## of the Gaussian likelihood
    mode_info <- mode_info / tauSq
  }

  return(unname((deriv2.prior - info.prior0) / mode_info))
}

## derivative of a single log-normal
normLogGrad <- function(x, mu, sigma) {
  -1 * (x - mu) / (sigma)^2
}

## second derivative of a single log-normal
normLogHess <- function(x, mu, sigma) {
  -1 / sigma^2
}

normMixInfo <- function(mix, x) {
  mixInfo(mix, x, dnorm, normLogGrad, normLogHess)
}

## info metric for a single norm, i.e. negative second derivative of log norm
normInfo <- function(x, mean, sigma) {
  -normLogHess(x, mean, sigma)
}

## Fisher info for normal sampling sd with known unit variance
normStdFisherInfo_inverse <- function(x) {
  1.0
}

make_normTransformedFisherInfo_inverse <- function(family) {
  function(x) {
    1 / family$variance(family$linkinv(x))
  }
}
