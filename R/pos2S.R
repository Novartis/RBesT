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
pos2S.default <- function(prior1, prior2, n1, n2, decision, ...) {
  "Unknown density"
}

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

    res <- rep(-Inf, times = length(lim2[1]:lim2[2]))

    if (is(decision, "decision2S_1sided")) {
      boundary <- crit_y1(lim2[1]:lim2[2], lim1 = lim1)

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
      }
    } else {
      boundary_lower_or_equal_than <- crit_y1$lower_or_equal_than(
        lim2[1]:lim2[2],
        lim1 = lim1
      )
      boundary_higher_than <- crit_y1$higher_than(
        lim2[1]:lim2[2],
        lim1 = lim1
      )

      for (i in lim2[1]:lim2[2]) {
        y2ind <- i - lim2[1] + 1
        lower_or_equal <- boundary_lower_or_equal_than[y2ind]
        higher <- boundary_higher_than[y2ind]

        if (
          lower_or_equal <= higher ||
            lower_or_equal == -1 ||
            higher == -1
        ) {
          ## decision was always 0
          res[y2ind] <- -Inf
        } else if (higher == n1 + 1 && lower_or_equal == n1 + 1) {
          ## both criteria are always 1
          res[y2ind] <- 0
        } else if (higher == n1 + 1) {
          ## upper-tail criterion is always 1, hence only lower-tail criterion matters
          res[y2ind] <- pmix(
            pred_mix1,
            lower_or_equal,
            lower.tail = TRUE,
            log.p = TRUE
          )
        } else if (lower_or_equal == n1 + 1) {
          ## lower-tail criterion is always 1, hence only upper-tail criterion matters
          res[y2ind] <- pmix(
            pred_mix1,
            higher,
            lower.tail = FALSE,
            log.p = TRUE
          )
        } else {
          ## calculate for all requested theta1 the probability mass
          ## <= lower_or_equal boundary and > higher_than boundary
          res[y2ind] <- log(
            pmix(
              pred_mix1,
              lower_or_equal,
              lower.tail = TRUE
            ) -
              pmix(
                pred_mix1,
                higher,
                lower.tail = TRUE
              )
          )
        }
      }
    }

    for (i in lim2[1]:lim2[2]) {
      y2ind <- i - lim2[1] + 1

      ## finally weight with the density according to the occurence
      ## of i due to theta2;
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
  Ngrid = 5,
  family = NULL,
  offset1 = 0,
  offset2 = offset1,
  ...
) {
  # Determine data sigma for each arm
  resolved1 <- resolve_sigma_family(family, missing(sigma1), sigma1, offset1)
  resolved2 <- resolve_sigma_family(family, missing(sigma2), sigma2, offset2)
  sigma_fun1 <- resolved1$sigma_fun
  sigma_fun2 <- resolved2$sigma_fun
  family <- resolved1$family

  # Resolve sigma for the fixed-sigma path
  if (is.null(sigma_fun1)) {
    if (missing(sigma1)) {
      sigma1 <- RBesT::sigma(prior1)
      message("Using default prior 1 reference scale ", sigma1)
    }
    if (missing(sigma2)) {
      sigma2 <- RBesT::sigma(prior2)
      message("Using default prior 2 reference scale ", sigma2)
    }
  } else {
    sigma1 <- RBesT::sigma(prior1)
    sigma2 <- RBesT::sigma(prior2)
  }
  assert_number(sigma1, lower = 0)
  assert_number(sigma2, lower = 0)
  sigma(prior1) <- sigma1
  sigma(prior2) <- sigma2

  # Calculate critical boundary
  if (is.null(sigma_fun1)) {
    crit_y1 <- decision2S_boundary(
      prior1, prior2, n1, n2, decision, sigma1, sigma2, eps, Ngrid
    )
  } else {
    crit_y1 <- decision2S_boundary(
      prior1, prior2, n1, n2, decision,
      family = family, offset1 = offset1, offset2 = offset2, eps = eps, Ngrid = Ngrid
    )
  }

  # Family path: PoS via numerical integration over both priors using oc2S
  if (!is.null(sigma_fun1)) {
    # Build OC function (reuses the same boundary)
    oc_fn <- oc2S(prior1, prior2, n1, n2, decision,
         family = family, offset1 = offset1, offset2 = offset2, eps = eps, Ngrid = Ngrid)

    oc_eval <- environment(oc_fn)$design_eval
    oc_freq <- environment(oc_fn)$Vfreq
    oc_misfit <- environment(oc_fn)$Vmisfit

    design_fun <- function(mix1, mix2) {
      # Integrate fun(theta1, theta2) over theta1 ~ mix1, theta2 ~ mix2
      integrate_truths <- function(fun) {
        # Inner: for a SCALAR theta2, integrate over theta1 ~ mix1
        inner_for_theta2 <- function(theta2_scalar) {
          inner_log_integrand <- function(theta1) {
            # theta1 can be a vector from integrate
            log(pmax(vapply(theta1, function(t1) oc_eval(t1, theta2_scalar, fun), numeric(1)),
                     .Machine$double.xmin))
          }
          integrate_density_log(mix1, inner_log_integrand, logit(eps / 2), logit(1 - eps / 2))
        }
        # Outer: integrate inner result over theta2 ~ mix2
        outer_log_integrand <- function(theta2) {
          # theta2 can be a vector from integrate
          log(pmax(vapply(theta2, inner_for_theta2, numeric(1)), .Machine$double.xmin))
        }
        integrate_density_log(
          mix2, outer_log_integrand,
          logit(eps / 2), logit(1 - eps / 2)
        )
      }
      res <- integrate_truths(oc_freq)
      # Probability of the outcomes for which the decision disagrees
      # with the boundary, see oc2S. As this is a coarse check only, the
      # truths are represented by the centers of ten equally likely bins.
      truths <- expand.grid(
        theta1 = qmix(mix1, (1:10 - 0.5) / 10),
        theta2 = qmix(mix2, (1:10 - 0.5) / 10)
      )
      mass <- mean(oc_eval(truths$theta1, truths$theta2, oc_misfit))
      if (mass > boundary2S_misfit_tol) {
        warn_boundary2S_misfit(
          "for outcomes which are likely under the given distributions",
          "The result is",
          mass
        )
      }
      res
    }
    return(design_fun)
  }

  # Fixed-sigma path: use preddist shortcut
  design_fun <- if (is(decision, "decision2S_1sided")) {
    # Simple case of one-sided boundary.
    assert_function(crit_y1)
    lower.tail <- attr(decision, "lower.tail")

    function(mix1, mix2) {
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

      if (n2 == 0) {
        mean_prior2 <- summary(prior2, probs = c())["mean"]
        pmix(
          pred_mix1_mean,
          crit_y1(mean_prior2, lim1 = lim1),
          lower.tail = lower.tail
        )
      } else {
        ## resolve the boundary where the predictive of sample 2 has its
        ## mass; the integrand then only reads the interpolant
        crit_y1_cache <- boundary2S_cache(crit_y1)
        crit_y1_cache$resolve(mass_intervals(pred_mix2_mean, eps), lim1)
        integrate_density_log(
          pred_mix2_mean,
          function(x) {
            pmix(
              pred_mix1_mean,
              crit_y1_cache$eval(x),
              lower.tail = lower.tail,
              log.p = TRUE
            )
          },
          logit(eps / 2),
          logit(1 - eps / 2)
        )
      }
    }
  } else {
    # Mixed boundary case.
    assert_list(crit_y1, len = 2, types = "function")
    crit_y1_lower_or_equal_than <- crit_y1$lower_or_equal_than
    crit_y1_higher_than <- crit_y1$higher_than

    function(mix1, mix2) {
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

      if (n2 == 0) {
        mean_prior2 <- summary(prior2, probs = c())["mean"]
        bound_lower_or_equal_than <- crit_y1_lower_or_equal_than(
          mean_prior2,
          lim1 = lim1
        )
        bound_higher_than <- crit_y1_higher_than(mean_prior2, lim1 = lim1)
        if (bound_lower_or_equal_than <= bound_higher_than) {
          0
        } else {
          pmix(pred_mix1_mean, bound_lower_or_equal_than, lower.tail = TRUE) -
            pmix(pred_mix1_mean, bound_higher_than, lower.tail = TRUE)
        }
      } else {
        ## resolve the boundaries where the predictive of sample 2 has its
        ## mass; the integrand then only reads the interpolants
        mass2 <- mass_intervals(pred_mix2_mean, eps)
        cache_lower_or_equal_than <- boundary2S_cache(
          crit_y1_lower_or_equal_than
        )
        cache_higher_than <- boundary2S_cache(crit_y1_higher_than)
        cache_lower_or_equal_than$resolve(mass2, lim1)
        cache_higher_than$resolve(mass2, lim1)
        integrand <- function(x) {
          bound_lower_or_equal_than <- cache_lower_or_equal_than$eval(x)
          bound_higher_than <- cache_higher_than$eval(x)
          # We need to expect here a vector x.
          ifelse(
            bound_lower_or_equal_than <= bound_higher_than,
            -Inf,
            log(
              pmix(
                pred_mix1_mean,
                bound_lower_or_equal_than,
                lower.tail = TRUE
              ) -
                pmix(pred_mix1_mean, bound_higher_than, lower.tail = TRUE)
            )
          )
        }
        integrate_density_log(
          pred_mix2_mean,
          integrand,
          logit(eps / 2),
          logit(1 - eps / 2)
        )
      }
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

  design_fun <- if (is(decision, "decision2S_1sided")) {
    # Simple case of one-sided boundary.
    assert_function(crit_y1)
    lower.tail <- attr(decision, "lower.tail")

    function(mix1, mix2) {
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
  } else {
    # Mixed boundary case.
    assert_list(crit_y1, len = 2, types = "function")
    crit_y1_lower_or_equal_than <- crit_y1$lower_or_equal_than
    crit_y1_higher_than <- crit_y1$higher_than

    function(mix1, mix2) {
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
      crit_y1_lower_or_equal_than(lim2, lim1)
      crit_y1_higher_than(lim2, lim1)

      grid <- seq(lim2[1], lim2[2])
      bound_lower_or_equal_than <- crit_y1_lower_or_equal_than(
        grid,
        lim1 = lim1
      )
      bound_higher_than <- crit_y1_higher_than(grid, lim1 = lim1)
      log_prob_in_bounds <- numeric(length(grid))
      has_zero_prob <- bound_lower_or_equal_than <= bound_higher_than
      log_prob_in_bounds[has_zero_prob] <- -Inf
      if (!all(has_zero_prob)) {
        log_prob_in_bounds[!has_zero_prob] <-
          log(
            pmix(
              pred_mix1_sum,
              bound_lower_or_equal_than[!has_zero_prob],
              lower.tail = TRUE
            ) -
              pmix(
                pred_mix1_sum,
                bound_higher_than[!has_zero_prob],
                lower.tail = TRUE
              )
          )
      }
      exp(matrixStats::logSumExp(
        dmix(pred_mix2_sum, grid, log = TRUE) +
          log_prob_in_bounds
      ))
    }
  }
  design_fun
}
