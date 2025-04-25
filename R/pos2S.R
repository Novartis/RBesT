#' Probability of Success for 2 Sample Design
#'
#' The `pos2S` function defines a 2 sample design (priors, sample
#' sizes & decision function) for the calculation of the probability
#' of success. A function is returned which calculates the calculates
#' the frequency at which the decision function is evaluated to 1 when
#' parameters are distributed according to the given distributions.
#'
#' @template args-boundary2S
#'
#' @details The `pos2S` function defines a 2 sample design and
#' returns a function which calculates its probability of success.
#' The probability of success is the frequency with which the decision
#' function is evaluated to 1 under the assumption of a given true
#' distribution of the data implied by a distirbution of the
#' parameters \eqn{\theta_1} and \eqn{\theta_2}.
#'
#' The calculation is analogous to the operating characeristics
#' [oc2S()] with the difference that instead of assuming
#' known (point-wise) true parameter values a distribution is
#' specified for each parameter.
#'
#' Calling the `pos2S` function calculates the decision boundary
#' \eqn{D_1(y_2)} and returns a function which can be used to evaluate the
#' PoS for different predictive distributions. It is evaluated as
#'
#' \deqn{ \int\int\int f_2(y_2|\theta_2) \, p(\theta_2) \, F_1(D_1(y_2)|\theta_1) \, p(\theta_1) \, dy_2 d\theta_2 d\theta_1. }
#'
#' where \eqn{F} is the distribution function of the sampling
#' distribution and \eqn{p(\theta_1)} and \eqn{p(\theta_2)} specifies
#' the assumed true distribution of the parameters \eqn{\theta_1} and
#' \eqn{\theta_2}, respectively. Each distribution \eqn{p(\theta_1)}
#' and \eqn{p(\theta_2)} is a mixture distribution and given as the
#' `mix1` and `mix2` argument to the function.
#'
#' For example, in the binary case an integration of the predictive
#' distribution, the BetaBinomial, instead of the binomial
#' distribution will be performed over the data space wherever the
#' decision function is evaluated to 1. All other aspects of the
#' calculation are as for the 2-sample operating characteristics, see
#' [oc2S()].
#'
#' @return Returns a function which when called with two arguments
#' `mix1` and `mix2` will return the frequencies at
#' which the decision function is evaluated to 1. Each argument is
#' expected to be a mixture distribution representing the assumed true
#' distribution of the parameter in each group.
#'
#' @family design2S
#'
#' @examples
#'
#' # see ?decision2S for details of example
#' priorT <- mixnorm(c(1, 0, 0.001), sigma = 88, param = "mn")
#' priorP <- mixnorm(c(1, -49, 20), sigma = 88, param = "mn")
#' # the success criteria is for delta which are larger than some
#' # threshold value which is why we set lower.tail=FALSE
#' successCrit <- decision2S(c(0.95, 0.5), c(0, 50), FALSE)
#'
#' # example interim outcome
#' postP_interim <- postmix(priorP, n = 10, m = -50)
#' postT_interim <- postmix(priorT, n = 20, m = -80)
#'
#' # assume that mean -50 / -80 were observed at the interim for
#' # placebo control(n=10) / active treatment(n=20) which gives
#' # the posteriors
#' postP_interim
#' postT_interim
#'
#' # then the PoS to succeed after another 20/30 patients is
#' pos_final <- pos2S(postP_interim, postT_interim, 20, 30, successCrit)
#'
#' pos_final(postP_interim, postT_interim)
#'
#' @export
pos2S <- function(prior1, prior2, n1, n2, decision, ...) UseMethod("pos2S")
#' @export
pos2S.default <- function(prior1, prior2, n1, n2, decision, ...)
  "Unknown density"

#' @templateVar fun pos2S
#' @template design2S-binomial
#' @export
pos2S.betaMix <- function(prior1, prior2, n1, n2, decision, eps, ...) {
  if (missing(eps) & (n1 * n2 > 1e7)) {
    warning("Large sample space. Consider setting eps=1e-6.")
  }

  crit_y1 <- decision2S_boundary(prior1, prior2, n1, n2, decision, eps)

  lower.tail <- attr(decision, "lower.tail")

  approx_method <- !missing(eps)

  design_fun <- function(mix1, mix2) {
    ## for each 0:n1 of the possible outcomes, calculate the
    ## probability mass past the boundary (log space) weighted with
    ## the density as the value for 1 occures (due to theta1)

    pred_mix1 <- preddist(mix1, n = n1)
    pred_mix2 <- preddist(mix2, n = n2)

    assert_that(inherits(pred_mix1, "betaBinomialMix"))
    assert_that(inherits(pred_mix2, "betaBinomialMix"))

    ## now get the decision boundary in the needed range
    lim1 <- c(0, n1)
    lim2 <- c(0, n2)

    ## in case we use the approximate method, we restrict the
    ## evaluation of the decision function range
    if (approx_method) {
      lim1 <- qmix(pred_mix1, c(eps / 2, 1 - eps / 2))
      lim2 <- qmix(pred_mix2, c(eps / 2, 1 - eps / 2))
    }

    boundary <- crit_y1(lim2[1]:lim2[2], lim1 = lim1)
    res <- rep(-Inf, times = length(boundary))

    for (i in lim2[1]:lim2[2]) {
      y2ind <- i - lim2[1] + 1
      if (boundary[y2ind] == -1) {
        ## decision was always 0
        res[y2ind] <- -Inf
      } else if (boundary[y2ind] == n1 + 1) {
        ## decision was always 1
        res[y2ind] <- 0
      } else {
        ## calculate for the predictive for dtheta1 the
        ## probability mass past (or before) the boundary
        res[y2ind] <- pmix(
          pred_mix1,
          boundary[y2ind],
          lower.tail = lower.tail,
          log.p = TRUE
        )
      }
      ## finally weight with the density according to the occurence
      ## of i due to theta2; the pmax avoids -Inf in a case of Prob==0
      res[y2ind] <- res[y2ind] + dmix(pred_mix2, i, log = TRUE)
    }
    exp(matrixStats::logSumExp(res))
  }
  design_fun
}


#' @templateVar fun pos2S
#' @template design2S-normal
#' @export
pos2S.normMix <- function(
  prior1,
  prior2,
  n1,
  n2,
  decision,
  sigma1,
  sigma2,
  eps = 1e-6,
  Ngrid = 10,
  ...
) {
  ## distributions of the means of the data generating distributions
  ## for now we assume that the underlying standard deviation
  ## matches the respective reference scales

  if (missing(sigma1)) {
    sigma1 <- RBesT::sigma(prior1)
    message("Using default prior 1 reference scale ", sigma1)
  }
  assert_number(sigma1, lower = 0)
  sigma(prior1) <- sigma1

  if (missing(sigma2)) {
    sigma2 <- RBesT::sigma(prior2)
    message("Using default prior 2 reference scale ", sigma2)
  }
  assert_number(sigma2, lower = 0)
  sigma(prior2) <- sigma2

  crit_y1 <- decision2S_boundary(
    prior1,
    prior2,
    n1,
    n2,
    decision,
    sigma1,
    sigma2,
    eps,
    Ngrid
  )

  lower.tail <- attr(decision, "lower.tail")

  design_fun <- function(mix1, mix2) {
    ## get the predictive of the mean
    pred_mix1_mean <- preddist(mix1, n = n1, sigma = sigma1)
    if (n2 == 0) {
      ## gets ignored anyway
      pred_mix2_mean <- preddist(mix2, n = 1, sigma = sigma2)
    } else {
      pred_mix2_mean <- preddist(mix2, n = n2, sigma = sigma2)
    }

    assert_that(inherits(pred_mix1_mean, "normMix"))
    assert_that(inherits(pred_mix2_mean, "normMix"))

    lim1 <- qmix(pred_mix1_mean, c(eps / 2, 1 - eps / 2))
    lim2 <- qmix(pred_mix2_mean, c(eps / 2, 1 - eps / 2))
    crit_y1(lim2, lim1)

    ## return(list(crit=crit_y1, m1=pred_dtheta1_mean, m2=pred_dtheta2_mean))

    if (n2 == 0) {
      mean_prior2 <- summary(prior2, probs = c())["mean"]
      return(pmix(
        pred_mix1_mean,
        crit_y1(mean_prior2),
        lower.tail = lower.tail
      ))
    } else {
      return(integrate_density_log(
        function(x)
          pmix(
            pred_mix1_mean,
            crit_y1(x, lim1 = lim1),
            lower.tail = lower.tail,
            log.p = TRUE
          ),
        pred_mix2_mean,
        logit(eps / 2),
        logit(1 - eps / 2)
      ))
    }
  }

  design_fun
}

#' @templateVar fun pos2S
#' @template design2S-poisson
#' @export
pos2S.gammaMix <- function(prior1, prior2, n1, n2, decision, eps = 1e-6, ...) {
  assert_that(likelihood(prior1) == "poisson")
  assert_that(likelihood(prior2) == "poisson")

  crit_y1 <- decision2S_boundary(prior1, prior2, n1, n2, decision, eps)

  lower.tail <- attr(decision, "lower.tail")

  design_fun <- function(mix1, mix2) {
    assert_that(likelihood(mix1) == "poisson")
    assert_that(likelihood(mix2) == "poisson")

    ## get the predictive of the sum
    pred_mix1_sum <- preddist(mix1, n = n1)
    pred_mix2_sum <- preddist(mix2, n = n2)

    assert_that(inherits(pred_mix1_sum, "gammaPoissonMix"))
    assert_that(inherits(pred_mix2_sum, "gammaPoissonMix"))

    lim1 <- qmix(pred_mix1_sum, c(eps / 2, 1 - eps / 2))
    lim2 <- qmix(pred_mix2_sum, c(eps / 2, 1 - eps / 2))

    ## force lower limit of lim1 to be 0 such that we will get and
    ## answer in most cases; performance wise it should be ok as
    ## we run a O(log(N)) search
    lim1[1] <- 0

    ## ensure that the boundaries are cached
    crit_y1(lim2, lim1 = lim1)
    grid <- seq(lim2[1], lim2[2])
    exp(matrixStats::logSumExp(
      dmix(pred_mix2_sum, grid, log = TRUE) +
        pmix(
          pred_mix1_sum,
          crit_y1(grid, lim1 = lim1),
          lower.tail = lower.tail,
          log.p = TRUE
        )
    ))
  }
  design_fun
}
