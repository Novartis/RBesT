#' Operating Characteristics for 2 Sample Design
#'
#' The `oc2S` function defines a 2 sample design (priors, sample
#' sizes & decision function) for the calculation of operating
#' characeristics. A function is returned which calculates the
#' calculates the frequency at which the decision function is
#' evaluated to 1 when assuming known parameters.
#'
#' @template args-boundary2S
#'
#' @details The `oc2S` function defines a 2 sample design and
#' returns a function which calculates its operating
#' characteristics. This is the frequency with which the decision
#' function is evaluated to 1 under the assumption of a given true
#' distribution of the data defined by the known parameter
#' \eqn{\theta_1} and \eqn{\theta_2}. The 2 sample design is defined
#' by the priors, the sample sizes and the decision function,
#' \eqn{D(y_1,y_2)}. These uniquely define the decision boundary , see
#' [decision2S_boundary()].
#'
#' Calling the `oc2S` function calculates the decision boundary
#' \eqn{D_1(y_2)} (see [decision2S_boundary()]) and returns
#' a function which can be used to calculate the desired frequency
#' which is evaluated as
#'
#' \deqn{ \int f_2(y_2|\theta_2) F_1(D_1(y_2)|\theta_1) dy_2. }
#'
#' See below for examples and specifics for the supported mixture
#' priors.
#'
#' @return Returns a function which when called with two arguments
#' `theta1` and `theta2` will return the frequencies at
#' which the decision function is evaluated to 1 whenever the data is
#' distributed according to the known parameter values in each
#' sample. Note that the returned function takes vector arguments.
#'
#' @references \insertRef{schmidli2014}{RBesT}
#'
#' @family design2S
#'
#' @examples
#'
#' # example from Schmidli et al., 2014
#' dec <- decision2S(0.975, 0, lower.tail = FALSE)
#'
#' prior_inf <- mixbeta(c(1, 4, 16))
#' prior_rob <- robustify(prior_inf, weight = 0.2, mean = 0.5)
#' prior_uni <- mixbeta(c(1, 1, 1))
#'
#' N <- 40
#' N_ctl <- N - 20
#'
#' # compare designs with different priors
#' design_uni <- oc2S(prior_uni, prior_uni, N, N_ctl, dec)
#' design_inf <- oc2S(prior_uni, prior_inf, N, N_ctl, dec)
#' design_rob <- oc2S(prior_uni, prior_rob, N, N_ctl, dec)
#'
#' # type I error
#' curve(design_inf(x, x), 0, 1)
#' curve(design_uni(x, x), lty = 2, add = TRUE)
#' curve(design_rob(x, x), lty = 3, add = TRUE)
#'
#' # power
#' curve(design_inf(0.2 + x, 0.2), 0, 0.5)
#' curve(design_uni(0.2 + x, 0.2), lty = 2, add = TRUE)
#' curve(design_rob(0.2 + x, 0.2), lty = 3, add = TRUE)
#'
#' @export
oc2S <- function(prior1, prior2, n1, n2, decision, ...) UseMethod("oc2S")
#' @export
oc2S.default <- function(prior1, prior2, n1, n2, decision, ...) {
  "Unknown density"
}

#' @templateVar fun oc2S
#' @template design2S-binomial
#' @export
oc2S.betaMix <- function(prior1, prior2, n1, n2, decision, eps, ...) {
  if (missing(eps) & ((n1 + 1) * (n2 + 1) > 1e7)) {
    warning("Large sample space. Consider setting eps=1e-6.")
  }

  crit_y1 <- decision2S_boundary(prior1, prior2, n1, n2, decision, eps)

  lower.tail <- attr(decision, "lower.tail")

  approx_method <- !missing(eps)

  design_fun <- function(theta1, theta2, y2) {
    ## other-wise we calculate the frequencies at which the
    ## decision is 1 (probability mass with decision==1)

    ## in case n2==0, then theta2 is irrelevant
    if (n2 == 0 & missing(theta2)) {
      theta2 <- 0.5
    }

    if (!missing(y2)) {
      lifecycle::deprecate_warn(
        "1.3.0",
        "design_fun(y2)",
        "decision2S_boundary(y2)",
        "Refer to ?decision2S_boundary() for details."
      )
      return(crit_y1(y2, lim1 = c(0, n1)))
    }

    assert_numeric(theta1, lower = 0, upper = 1, finite = TRUE)
    assert_numeric(theta2, lower = 0, upper = 1, finite = TRUE)

    theta_df <- try(data.frame(
      theta1 = theta1,
      theta2 = theta2,
      row.names = NULL
    ))
    if (inherits(theta_df, "try-error")) {
      stop("theta1 and theta2 need to be of same size")
    }

    ## now get the decision boundary in the needed range
    lim1 <- c(0, n1)
    lim2 <- c(0, n2)

    ## in case we use the approximate method, we restrict the
    ## evaluation of the decision function range
    if (approx_method) {
      lim1[1] <- qbinom(eps / 2, n1, min(theta1))
      lim1[2] <- qbinom(1 - eps / 2, n1, max(theta1))
      lim2[1] <- qbinom(eps / 2, n2, min(theta2))
      lim2[2] <- qbinom(1 - eps / 2, n2, max(theta2))
    }

    ## for each 0:n1 of the possible outcomes, calculate the
    ## probability mass past the boundary (log space) weighted with
    ## the density as the value for 1 occures (due to theta1)
    res <- matrix(-Inf, nrow = diff(lim2) + 1, ncol = nrow(theta_df))

    if (is(decision, "decision2S_1sided")) {
      boundary <- crit_y1(lim2[1]:lim2[2], lim1 = c(0, n1))

      for (i in lim2[1]:lim2[2]) {
        y2ind <- i - lim2[1] + 1
        if (boundary[y2ind] == -1) {
          ## decision was always 0
          res[y2ind, ] <- -Inf
        } else if (boundary[y2ind] == n1 + 1) {
          ## decision was always 1
          res[y2ind, ] <- 0
        } else {
          ## calculate for all requested theta1 the probability mass
          ## past (or before) the boundary
          res[y2ind, ] <- pbinom(
            boundary[y2ind],
            n1,
            theta_df$theta1,
            lower.tail = lower.tail,
            log.p = TRUE
          )
        }
      }
    } else {
      boundary_lower_or_equal_than <- crit_y1$lower_or_equal_than(
        lim2[1]:lim2[2],
        lim1 = c(0, n1)
      )
      boundary_higher_than <- crit_y1$higher_than(
        lim2[1]:lim2[2],
        lim1 = c(0, n1)
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
          res[y2ind, ] <- -Inf
        } else if (higher == n1 + 1 && lower_or_equal == n1 + 1) {
          ## both criteria are always 1
          res[y2ind, ] <- 0
        } else if (higher == n1 + 1) {
          ## upper-tail criterion is always 1, hence only lower-tail criterion matters
          res[y2ind, ] <- pbinom(
            lower_or_equal,
            n1,
            theta_df$theta1,
            lower.tail = TRUE,
            log.p = TRUE
          )
        } else if (lower_or_equal == n1 + 1) {
          ## lower-tail criterion is always 1, hence only upper-tail criterion matters
          res[y2ind, ] <- pbinom(
            higher,
            n1,
            theta_df$theta1,
            lower.tail = FALSE,
            log.p = TRUE
          )
        } else {
          ## calculate for all requested theta1 the probability mass
          ## <= lower_or_equal boundary and > higher_than boundary
          res[y2ind, ] <- log(
            pbinom(
              lower_or_equal,
              n1,
              theta_df$theta1,
              lower.tail = TRUE
            ) -
              pbinom(
                higher,
                n1,
                theta_df$theta1,
                lower.tail = TRUE
              )
          )
        }
      }
    }

    for (i in lim2[1]:lim2[2]) {
      y2ind <- i - lim2[1] + 1

      ## finally weight with the density according to the occurence
      ## of i due to theta2; the pmax avoids -Inf in a case of Prob==0
      ## res[y2ind,] <- res[y2ind,] + pmax(dbinom(i, n2, T$theta2, log=TRUE), -700)
      res[y2ind, ] <- res[y2ind, ] + dbinom(i, n2, theta_df$theta2, log = TRUE)
    }
    exp(matrixStats::colLogSumExps(res))
  }
  design_fun
}

#' @templateVar fun oc2S
#' @template design2S-normal
#' @export
oc2S.normMix <- function(
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

  sem1 <- sigma1 / sqrt(n1)
  sem2 <- sigma2 / sqrt(n2)

  if (n2 == 0) {
    sem2 <- sigma(prior2) / sqrt(1E-1)
  }

  ## change the reference scale of the prior such that the prior
  ## represents the distribution of the respective means
  mean_prior1 <- prior1
  sigma(mean_prior1) <- sem1

  freq <- if (is(decision, "decision2S_1sided")) {
    # Simple case of one-sided boundary.
    assert_function(crit_y1)
    lower.tail <- attr(decision, "lower.tail")
    ## design_fun resolves the boundary where sample 2 has its mass; the
    ## integrand then only reads the interpolant.
    crit_y1_cache <- boundary2S_cache(crit_y1)
    function(theta1, theta2) {
      sem1_t <- if (is.null(sigma_fun1)) sem1 else sigma_fun1(theta1) / sqrt(n1)
      sem2_t <- if (is.null(sigma_fun2)) sem2 else sigma_fun2(theta2) / sqrt(n2)
      lim1 <- qnorm(c(eps / 2, 1 - eps / 2), theta1, sem1_t)
      if (n2 == 0) {
        pnorm(
          crit_y1(theta2, lim1 = lim1),
          theta1,
          sem1_t,
          lower.tail = lower.tail
        )
      } else {
        integrate_density_log(
          mixnorm(c(1, theta2, sem2_t), sigma = sem2_t),
          function(x) {
            pnorm(
              crit_y1_cache$eval(x),
              theta1,
              sem1_t,
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
    cache_lower_or_equal_than <- boundary2S_cache(
      crit_y1_lower_or_equal_than
    )
    cache_higher_than <- boundary2S_cache(crit_y1_higher_than)
    function(theta1, theta2) {
      assert_scalar(theta1)
      assert_scalar(theta2)
      sem1_t <- if (is.null(sigma_fun1)) sem1 else sigma_fun1(theta1) / sqrt(n1)
      sem2_t <- if (is.null(sigma_fun2)) sem2 else sigma_fun2(theta2) / sqrt(n2)
      lim1 <- qnorm(c(eps / 2, 1 - eps / 2), theta1, sem1_t)
      if (n2 == 0) {
        bound_lower_or_equal_than <- crit_y1_lower_or_equal_than(
          theta2,
          lim1 = lim1
        )
        bound_higher_than <- crit_y1_higher_than(theta2, lim1 = lim1)
        if (bound_lower_or_equal_than <= bound_higher_than) {
          0
        } else {
          pnorm(bound_lower_or_equal_than, theta1, sem1_t, lower.tail = TRUE) -
            pnorm(bound_higher_than, theta1, sem1_t, lower.tail = TRUE)
        }
      } else {
        integrand <- function(x) {
          bound_lower_or_equal_than <- cache_lower_or_equal_than$eval(x)
          bound_higher_than <- cache_higher_than$eval(x)
          # We need to expect here a vector x.
          ifelse(
            bound_lower_or_equal_than <= bound_higher_than,
            -Inf,
            ## both branches are evaluated, hence avoid a negative argument
            log(pmax(
              pnorm(
                bound_lower_or_equal_than,
                theta1,
                sem1_t,
                lower.tail = TRUE
              ) -
                pnorm(bound_higher_than, theta1, sem1_t, lower.tail = TRUE),
              0
            ))
          )
        }
        integrate_density_log(
          mixnorm(c(1, theta2, sem2_t), sigma = sem2_t),
          integrand,
          logit(eps / 2),
          logit(1 - eps / 2)
        )
      }
    }
  }

  Vfreq <- Vectorize(freq)

  ## With a family-based sigma the sampling standard error grows without
  ## bound for extreme outcomes, where the posterior reverts to the prior
  ## and the decision can change a second time along y1. The boundary
  ## holds a single critical value of y1 per y2 and ignores this. Estimate
  ## the probability mass of the outcomes under a truth for which the
  ## decision itself disagrees with the boundary on a coarse grid.
  Vmisfit <- NULL
  if (!is.null(sigma_fun1)) {
    boundary_success <- if (is(decision, "decision2S_1sided")) {
      lower.tail <- attr(decision, "lower.tail")
      function(y1, crit) {
        if (lower.tail) y1 <= crit[[1]] else y1 > crit[[1]]
      }
    } else {
      function(y1, crit) y1 > crit[[2]] & y1 <= crit[[1]]
    }
    boundary_at <- function(y2, lim1) {
      if (is(decision, "decision2S_1sided")) {
        if (n2 == 0) {
          return(list(crit_y1(y2, lim1 = lim1)))
        }
        return(list(crit_y1_cache$eval(y2)))
      }
      if (n2 == 0) {
        return(list(
          crit_y1_lower_or_equal_than(y2, lim1 = lim1),
          crit_y1_higher_than(y2, lim1 = lim1)
        ))
      }
      list(cache_lower_or_equal_than$eval(y2), cache_higher_than$eval(y2))
    }
    ## stay inside the region of sample 2 in which the boundary is resolved
    z <- seq(-1, 1, length.out = 5) * min(3, qnorm(1 - eps / 2))
    weight <- dnorm(z) / sum(dnorm(z))
    misfit <- function(theta1, theta2) {
      sem1_t <- sigma_fun1(theta1) / sqrt(n1)
      y1 <- theta1 + sem1_t * z
      post1 <- lapply(y1, function(m) {
        postmix(prior1, m = m, se = sigma_fun1(m) / sqrt(n1))
      })
      if (n2 == 0) {
        y2 <- theta2
        weight2 <- 1
        post2 <- list(prior2)
      } else {
        y2 <- theta2 + sigma_fun2(theta2) / sqrt(n2) * z
        weight2 <- weight
        post2 <- lapply(y2, function(m) {
          postmix(prior2, m = m, se = sigma_fun2(m) / sqrt(n2))
        })
      }
      lim1 <- qnorm(c(eps / 2, 1 - eps / 2), theta1, sem1_t)
      crit <- boundary_at(y2, lim1)
      mass <- 0
      for (j in seq_along(y2)) {
        crit_j <- lapply(crit, "[", j)
        for (i in seq_along(y1)) {
          ## an outcome on the boundary disagrees by the accuracy of the
          ## boundary only
          if (any(abs(y1[i] - unlist(crit_j)) < 0.1 * sem1_t)) {
            next
          }
          disagree <- (decision(post1[[i]], post2[[j]]) == 1) !=
            boundary_success(y1[i], crit_j)
          if (disagree) {
            mass <- mass + weight[i] * weight2[j]
          }
        }
      }
      mass
    }
    Vmisfit <- Vectorize(misfit)
  }

  ## Evaluates fun (Vfreq or Vmisfit) for the truths after resolving the
  ## boundary(s) where needed. Used by pos2S as well, which integrates
  ## both over the distributions of the truths.
  design_eval <- function(theta1, theta2, fun) {
    ## Use reference sem for the range in which critical values of
    ## sample 1 are searched
    sem1_ref <- sem1
    lim1 <- c(
      qnorm(eps / 2, min(theta1), sem1_ref),
      qnorm(1 - eps / 2, max(theta1), sem1_ref)
    )

    ## Resolve the boundary(s) where sample 2 has its probability mass
    ## under each requested truth. Regions of nearby truths merge, while
    ## the gap between truths far apart is left alone.
    if (n2 > 0) {
      sem2_t <- if (is.null(sigma_fun2)) sem2 else sigma_fun2(theta2) / sqrt(n2)
      mass2 <- cbind(
        qnorm(eps / 2, theta2, sem2_t),
        qnorm(1 - eps / 2, theta2, sem2_t)
      )
      if (is(decision, "decision2S_1sided")) {
        crit_y1_cache$resolve(mass2, lim1)
      } else {
        ## The caching is in the closures, therefore we don't need to
        ## worry about conflicts between the caches.
        cache_lower_or_equal_than$resolve(mass2, lim1)
        cache_higher_than$resolve(mass2, lim1)
      }
    }

    theta_df <- try(data.frame(
      theta1 = theta1,
      theta2 = theta2,
      row.names = NULL
    ))
    if (inherits(theta_df, "try-error")) {
      stop("theta1 and theta2 need to be of same size")
    }

    do.call(fun, theta_df)
  }

  design_fun <- function(theta1, theta2, y2) {
    if (!missing(y2)) {
      theta2 <- y2
    }

    ## in case n2==0, then theta2 is irrelevant
    if (n2 == 0 && missing(theta2)) {
      theta2 <- theta1
    }

    res <- design_eval(theta1, theta2, Vfreq)

    if (!is.null(Vmisfit)) {
      mass <- design_eval(theta1, theta2, Vmisfit)
      flagged <- mass > boundary2S_misfit_tol
      if (any(flagged)) {
        warn_boundary2S_misfit(
          paste("for", sum(flagged), "of", length(mass), "truths"),
          "Results for these truths are",
          max(mass)
        )
      }
    }

    res
  }

  design_fun
}

#' @templateVar fun oc2S
#' @template design2S-poisson
#' @export
oc2S.gammaMix <- function(prior1, prior2, n1, n2, decision, eps = 1e-6, ...) {
  assert_that(likelihood(prior1) == "poisson")
  assert_that(likelihood(prior2) == "poisson")

  crit_y1 <- decision2S_boundary(prior1, prior2, n1, n2, decision, eps)

  lower.tail <- attr(decision, "lower.tail")

  freq <- if (is(decision, "decision2S_1sided")) {
    # Simple case of one-sided boundary.
    assert_function(crit_y1)
    function(theta1, theta2) {
      lambda1 <- theta1 * n1
      lambda2 <- theta2 * n2
      lim1 <- qpois(c(eps / 2, 1 - eps / 2), lambda1)
      grid <- seq(qpois(eps / 2, lambda2), qpois(1 - eps / 2, lambda2))

      exp(matrixStats::logSumExp(
        dpois(grid, lambda2, log = TRUE) +
          ppois(
            crit_y1(grid, lim1 = lim1),
            lambda1,
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
    function(theta1, theta2) {
      lambda1 <- theta1 * n1
      lambda2 <- theta2 * n2
      lim1 <- qpois(c(eps / 2, 1 - eps / 2), lambda1)
      grid <- seq(qpois(eps / 2, lambda2), qpois(1 - eps / 2, lambda2))

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
            ppois(
              bound_lower_or_equal_than[!has_zero_prob],
              lambda1,
              lower.tail = TRUE
            ) -
              ppois(
                bound_higher_than[!has_zero_prob],
                lambda1,
                lower.tail = TRUE
              )
          )
      }
      exp(matrixStats::logSumExp(
        dpois(grid, lambda2, log = TRUE) + log_prob_in_bounds
      ))
    }
  }

  Vfreq <- Vectorize(freq)

  design_fun <- function(theta1, theta2, y2) {
    ## in case n2==0, then theta2 is irrelevant
    if (n2 == 0 & missing(theta2)) {
      theta2 <- theta1
    }

    if (!missing(y2)) {
      if (missing(theta1)) {
        theta1 <- summary(prior1)["mean"]
      }
      lambda2 <- y2
    } else {
      lambda2 <- theta2 * n2
    }

    lambda1 <- theta1 * n1

    lim1 <- c(0, 0)
    lim1[1] <- qpois(eps / 2, min(lambda1))
    lim1[2] <- qpois(1 - eps / 2, max(lambda1))

    if (n2 == 0) {
      lim2 <- c(0, 0)
    } else {
      lim2 <- c(0, 0)
      lim2[1] <- qpois(eps / 2, min(lambda2))
      lim2[2] <- qpois(1 - eps / 2, max(lambda2))
    }

    ## force lower limit of lim1 to be 0 such that we will get and
    ## answer in most cases; performance wise it should be ok as
    ## we run a O(log(N)) search
    lim1[1] <- 0

    ## Call boundary function(s) to cache all results for all
    ## requested computations.
    if (is(decision, "decision2S_1sided")) {
      crit_y1(lim2, lim1 = lim1)
    } else {
      ## The caching is in the closures, therefore we don't need to
      ## worry about conflicts between the caches.
      crit_y1$lower_or_equal_than(lim2, lim1 = lim1)
      crit_y1$higher_than(lim2, lim1 = lim1)
    }

    if (!missing(y2)) {
      lifecycle::deprecate_warn(
        "1.3.0",
        "design_fun(y2)",
        "decision2S_boundary(y2)",
        "Refer to ?decision2S_boundary() for details."
      )
      return(crit_y1(y2, lim1 = lim1))
    }

    theta_df <- try(data.frame(
      theta1 = theta1,
      theta2 = theta2,
      row.names = NULL
    ))
    if (inherits(theta_df, "try-error")) {
      stop("theta1 and theta2 need to be of same size")
    }

    do.call(Vfreq, theta_df)
  }
  design_fun
}
