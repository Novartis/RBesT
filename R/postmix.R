#' Conjugate Posterior Analysis
#'
#' @description
#' Calculates the posterior distribution for data \code{data} given a prior
#' \code{priormix}, where the prior is a mixture of conjugate distributions.
#' The posterior is then also a mixture of conjugate distributions.
#'
#' @param priormix prior (mixture of conjugate distributions).
#' @param data individual data. If the individual data is not given, then
#' summary data has to be provided (see below).
#' @param n sample size.
#' @param ... includes arguments which depend on the specific case, see description below.
#'
#' @details A conjugate prior-likelihood pair has the convenient
#' property that the posterior is in the same distributional class as
#' the prior. This property also applies to mixtures of conjugate
#' priors. Let
#'
#' \deqn{p(\theta;\mathbf{w},\mathbf{a},\mathbf{b})}{p(\theta;w,a,b)}
#'
#' denote a conjugate mixture prior density for data
#'
#' \deqn{y|\theta \sim f(y|\theta),}{y|\theta ~ f(y|\theta),}
#'
#' where \eqn{f(y|\theta)} is the likelihood. Then the posterior is
#' again a mixture with each component \eqn{k} equal to the respective
#' posterior of the \eqn{k}th prior component and updated weights
#' \eqn{w'_k},
#'
#' \deqn{p(\theta;\mathbf{w'},\mathbf{a'},\mathbf{b'}|y) = \sum_{k=1}^K w'_k \, p_k(\theta;a'_k,b'_k|y).}{p(\theta;w',a',b'|y) = \sum_{k=1}^K w'_k * p(\theta;a'_k,b'_k|y).}
#'
#' The weight \eqn{w'_k} for \eqn{k}th component is determined by the
#' marginal likelihood of the new data \eqn{y} under the \eqn{k}th prior
#' distribution which is given by the predictive distribution of the
#' \eqn{k}th component,
#'
#' \deqn{w'_k \propto w_k \, \int p_k(\theta;a_k,b_k) \, f(y|\theta) \, d\theta \equiv w^\ast_k .}{w'_k \propto w_k \int p_k(u;a_k,b_k) f(y|u) du = w^*_k .}
#'
#' The final weight \eqn{w'_k} is then given by appropriate
#' normalization, \eqn{w'_k = w^\ast_k / \sum_{k=1}^K w^\ast_k}{w'_k =
#' w^*_k / \sum_{k=1}^K w^*_k}. In other words, the weight of
#' component \eqn{k} is proportional to the likelihood that data
#' \eqn{y} is generated from the respective component, i.e. the
#' marginal probability; for details, see for example \emph{Schmidli
#' et al., 2015}.
#'
#' \emph{Note:} The prior weights \eqn{w_k} are fixed, but the
#' posterior weights \eqn{w'_k \neq w_k} still change due to the
#' changing normalization.
#'
#' The data \eqn{y} can either be given as individual data or as
#' summary data (sufficient statistics). See below for details for the
#' implemented conjugate mixture prior densities.
#'
#' @template conjugate_pairs
#'
#' @references Schmidli H, Gsteiger S, Roychoudhury S, O'Hagan A, Spiegelhalter D, Neuenschwander B.
#' Robust meta-analytic-predictive priors in clinical trials with historical control information.
#' \emph{Biometrics} 2014;70(4):1023-1032.
#'
#' @examples
#'
#' # binary example with individual data (1=event,0=no event), uniform prior
#' prior.unif <- mixbeta(c(1, 1, 1))
#' data.indiv <- c(1, 0, 1, 1, 0, 1)
#' posterior.indiv <- postmix(prior.unif, data.indiv)
#' print(posterior.indiv)
#' # or with summary data (number of events and number of patients)
#' r <- sum(data.indiv)
#' n <- length(data.indiv)
#' posterior.sum <- postmix(prior.unif, n = n, r = r)
#' print(posterior.sum)
#'
#' # binary example with robust informative prior and conflicting data
#' prior.rob <- mixbeta(c(0.5, 4, 10), c(0.5, 1, 1))
#' posterior.rob <- postmix(prior.rob, n = 20, r = 18)
#' print(posterior.rob)
#'
#' # normal example with individual data
#' sigma <- 88
#' prior.mean <- -49
#' prior.se <- sigma / sqrt(20)
#' prior <- mixnorm(c(1, prior.mean, prior.se), sigma = sigma)
#' data.indiv <- c(-46, -227, 41, -65, -103, -22, 7, -169, -69, 90)
#' posterior.indiv <- postmix(prior, data.indiv)
#' # or with summary data (mean and number of patients)
#' mn <- mean(data.indiv)
#' n <- length(data.indiv)
#' posterior.sum <- postmix(prior, m = mn, n = n)
#' print(posterior.sum)
#'
#' @export
postmix <- function(priormix, data, ...) UseMethod("postmix")

#' @export
postmix.default <- function(priormix, data, ...) "Unknown distribution"

#' @describeIn postmix Calculates the posterior beta mixture
#' distribution. The individual data vector is expected to be a vector
#' of 0 and 1, i.e. a series of Bernoulli experiments. Alternatively,
#' the sufficient statistics \code{n} and \code{r} can be given,
#' i.e. number of trials and successes, respectively.
#' @param r Number of successes.
#' @export
postmix.betaMix <- function(priormix, data, n, r, ...) {
  if (!missing(data)) {
    assert_that(all(data %in% c(0, 1)))
    r <- sum(data)
    n <- length(data)
  }
  w <- log(priormix[1, , drop = FALSE]) + dBetaBinomial(r, n, priormix[2, , drop = FALSE], priormix[3, , drop = FALSE], log = TRUE)
  priormix[1, ] <- exp(w - matrixStats::logSumExp(w))
  priormix[2, ] <- priormix[2, , drop = FALSE] + r
  priormix[3, ] <- priormix[3, , drop = FALSE] + n - r
  class(priormix) <- c("betaMix", "mix")
  priormix
}


#' @describeIn postmix Calculates the posterior normal mixture
#' distribution with the sampling likelihood being a normal with fixed
#' standard deviation. Either an individual data vector \code{data}
#' can be given or the sufficient statistics which are the standard
#' error \code{se} and sample mean \code{m}. If the sample size
#' \code{n} is used instead of the sample standard error, then the
#' reference scale of the prior is used to calculate the standard
#' error. Should standard error \code{se} and sample size \code{n} be
#' given, then the reference scale of the prior is updated; however it
#' is recommended to use the command \code{\link{sigma}} to set the
#' reference standard deviation.
#' @param m Sample mean.
#' @param se Sample standard error.
#' @export
postmix.normMix <- function(priormix, data, n, m, se, ...) {
  if (!missing(data)) {
    m <- mean(data)
    n <- length(data)
    se <- sd(data) / sqrt(n)
  } else {
    if (missing(m) & (missing(n) | missing(se))) {
      stop("Either raw data or summary data (m and se) must be given.")
    }
    if (!missing(se) & !missing(n)) {
      sigma(priormix) <- se * sqrt(n)
      message(paste0("Updating reference scale to ", sigma(priormix), ".\nIt is recommended to use the sigma command instead.\nSee ?sigma or ?mixnorm."))
    }
    if (missing(se) & !missing(n)) {
      message("Using default prior reference scale ", sigma(priormix))
      se <- sigma(priormix) / sqrt(n)
    }
  }
  dataPrec <- 1 / se^2
  ## prior precision
  priorPrec <- 1 / priormix[3, , drop = FALSE]^2
  ## posterior precision is prior + data precision
  postPrec <- priorPrec + dataPrec
  ## old weights times the likelihood under the predictive of each
  ## component in log space
  sigmaPred <- sqrt(priormix[3, , drop = FALSE]^2 + se^2)
  lw <- log(priormix[1, , drop = FALSE]) + dnorm(m, priormix[2, , drop = FALSE], sigmaPred, log = TRUE)
  priormix[1, ] <- exp(lw - matrixStats::logSumExp(lw))
  ## posterior means are precision weighted average of prior mean
  ## and data
  priormix[2, ] <- (priormix[2, , drop = FALSE] * priorPrec + m * dataPrec) / postPrec
  priormix[3, ] <- 1 / sqrt(postPrec)
  class(priormix) <- c("normMix", "mix")
  priormix
}

#' @describeIn postmix Calculates the posterior gamma mixture
#' distribution for Poisson and exponential likelihoods. Only the
#' Poisson case is supported in this version.
#' @export
postmix.gammaMix <- function(priormix, data, n, m, ...) {
  type <- likelihood(priormix)
  if (type != "poisson") {
    stop("NOT YET SUPPORTED: Updating Gamma priors is not yet supported for", type, "data. Sorry.")
  }
  if (!missing(data)) {
    s <- sum(data)
    n <- length(data)
  } else {
    s <- m * n
  }
  ## assert_int(s)
  ## the predictive distribution for n events with sufficient
  ## statistics s is the negative binomial with beta->beta/n
  if (n > 0) {
    w <- log(priormix[1, , drop = FALSE]) + .dnbinomAB(s, priormix[2, ], priormix[3, ] / n, log = TRUE)
  } else { ## case n=0
    w <- log(priormix[1, , drop = FALSE])
  }
  priormix[1, ] <- exp(w - matrixStats::logSumExp(w))
  priormix[2, ] <- priormix[2, , drop = FALSE] + s
  priormix[3, ] <- priormix[3, , drop = FALSE] + n
  class(priormix) <- c("gammaMix", "mix")
  priormix
}
