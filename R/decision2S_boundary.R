#' Decision Boundary for 2 Sample Designs
#'
#' The \code{decision2S_boundary} function defines a 2 sample design
#' (priors, sample sizes, decision function) for the calculation of
#' the decision boundary. A function is returned which calculates the
#' critical value of the first sample \eqn{y_{1,c}} as a function of
#' the outcome in the second sample \eqn{y_2}. At the decision
#' boundary, the decision function will change between 0 (failure) and
#' 1 (success) for the respective outcomes.
#'
#' @template args-boundary2S
#'
#' @details For a 2 sample design the specification of the priors, the
#' sample sizes and the decision function, \eqn{D(y_1,y_2)}, uniquely
#' defines the decision boundary
#'
#' \deqn{D_1(y_2) = \max_{y_1}\{D(y_1,y_2) = 1\},}{D_1(y_2) = max_{y_1}{D(y_1,y_2) = 1},}
#'
#' which is the critical value of \eqn{y_{1,c}} conditional on the
#' value of \eqn{y_2} whenever the decision \eqn{D(y_1,y_2)} function
#' changes its value from 0 to 1 for a decision function with
#' \code{lower.tail=TRUE} (otherwise the definition is \eqn{D_1(y_2) =
#' \max_{y_1}\{D(y_1,y_2) = 0\}}{D_1(y_2) = max_{y_1}{D(y_1,y_2) =
#' 0}}). The decision function may change at most at a single critical
#' value for given \eqn{y_{2}} as only one-sided decision functions
#' are supported. Here, \eqn{y_2} is defined for binary and Poisson
#' endpoints as the sufficient statistic \eqn{y_2 = \sum_{i=1}^{n_2}
#' y_{2,i}} and for the normal case as the mean \eqn{\bar{y}_2 = 1/n_2
#' \sum_{i=1}^{n_2} y_{2,i}}.
#'
#' @return Returns a function with a single argument. This function
#' calculates in dependence of the outcome \eqn{y_2} in sample 2 the
#' critical value \eqn{y_{1,c}} for which the defined design will
#' change the decision from 0 to 1 (or vice versa, depending on the
#' decision function).
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
#' # the futility criterion acts in the opposite direction
#' futilityCrit <- decision2S(c(0.90), c(40), TRUE)
#'
#' # success criterion boundary
#' successBoundary <- decision2S_boundary(priorP, priorT, 10, 20, successCrit)
#'
#' # futility criterion boundary
#' futilityBoundary <- decision2S_boundary(priorP, priorT, 10, 20, futilityCrit)
#'
#' curve(successBoundary(x), -25:25 - 49, xlab = "y2", ylab = "critical y1")
#' curve(futilityBoundary(x), lty = 2, add = TRUE)
#'
#' # hence, for mean in sample 2 of 10, the critical value for y1 is
#' y1c <- futilityBoundary(-10)
#'
#' # around the critical value the decision for futility changes
#' futilityCrit(postmix(priorP, m = y1c + 1E-3, n = 10), postmix(priorT, m = -10, n = 20))
#' futilityCrit(postmix(priorP, m = y1c - 1E-3, n = 10), postmix(priorT, m = -10, n = 20))
#'
#' @export
decision2S_boundary <- function(prior1, prior2, n1, n2, decision, ...) UseMethod("decision2S_boundary")
#' @export
decision2S_boundary.default <- function(prior1, prior2, n1, n2, decision, ...) "Unknown density"

#' @templateVar fun decision2S_boundary
#' @template design2S-binomial
# If the \code{eps} argument is specificed, then the
# returned function will use the additional \code{lim2}
# argument to limit the search for the critical value.
#' @export
decision2S_boundary.betaMix <- function(prior1, prior2, n1, n2, decision, eps, ...) {
  ## only n2=0 is supported
  assert_number(n1, lower = 1, finite = TRUE)
  assert_number(n2, lower = 0, finite = TRUE)

  if (!missing(eps)) {
    assert_number(eps, lower = 0, upper = 0.1, finite = TRUE)
  }

  cond_decisionDist <- function(post2cond) {
    fn <- function(m1) {
      ## Note: Subtracting from the decision 0.25 leads to
      ## negative decisions being at -0.25 while positives are
      ## at 0.75; since uniroot_int *always* returns the x which
      ## has lowest absolute value we are guaranteed that y2crit
      ## is just before the jump
      ## decision(post1cond, post2[[m2+1]]) - 0.25
      decision(postmix(prior1, r = m1, n = n1), post2cond) - 0.25
      ## decision(postmix(prior2, r=m2, n=n2), post1cond) - 0.25
    }
    Vectorize(fn)
  }

  ## saves the decision boundary conditional on the outcome of the
  ## second variable
  clim1 <- c(Inf, -Inf)
  clim2 <- c(Inf, -Inf)
  boundary <- c()
  full_boundary <- missing(eps)

  lower.tail <- attr(decision, "lower.tail")

  update_boundary <- function(lim1, lim2) {
    boundary <<- rep(NA, diff(lim2) + 1)
    clim2 <<- lim2
    clim1 <<- lim1
    for (y2 in lim2[1]:lim2[2]) {
      ## find decision point
      if (n2 == 0) {
        decFun <- cond_decisionDist(prior2)
      } else {
        decFun <- cond_decisionDist(postmix(prior2, r = y2, n = n2))
      }
      ind_llim <- decFun(lim1[1])
      ind_ulim <- decFun(lim1[2])
      y2ind <- y2 - lim2[1] + 1
      if (ind_llim < 0 & ind_ulim < 0) {
        ## then the decision is never 1
        boundary[y2ind] <<- -1
        next
      }
      if (ind_llim > 0 & ind_ulim > 0) {
        ## then the decision is always 1
        boundary[y2ind] <<- n1 + 1
        next
      }
      ## find boundary
      boundary[y2ind] <<- uniroot_int(decFun, lim1,
        f.lower = ind_llim,
        f.upper = ind_ulim
      )
    }
    if (lower.tail) {
      ## if lower.tail==TRUE, then the condition becomes true when
      ## going from large to small values, hence we need to integrate from
      ## 0 to boundary
      boundary <<- pmax(boundary - 1, -1)
    }
    return()
  }

  if (full_boundary) {
    update_boundary(c(0, n1), c(0, n2))
  }

  decision_boundary <- function(y2, lim1) {
    ## check if we need to recalculate the decision grid for the
    ## case of enabled approximate method
    assert_integerish(y2, lower = 0, upper = n2, any.missing = FALSE)

    if (!full_boundary) {
      if (missing(lim1)) {
        ## if not hint is given we search the full sample
        ## space which should be OK, as the complexity is
        ## log(N)
        lim1 <- c(0, n1)
      } else {
        assert_integerish(lim1, lower = 0, upper = n1, any.missing = FALSE)
      }
      lim2 <- c(min(y2), max(y2))
      ## check if the decision grid needs to be recomputed
      if (lim1[1] < clim1[1] | lim1[2] > clim1[2] |
        lim2[1] < clim2[1] | lim2[2] > clim2[2]) {
        ## ensure that lim1 never shrinks
        lim1[1] <- min(lim1[1], clim1[1])
        lim1[2] <- max(lim1[2], clim1[2])
        update_boundary(lim1, lim2)
      }
    }

    ## make sure y2 is an integer which is the value of
    ## the second read-out for which we return the decision
    ## boundary
    ## TODO: handle case with eps with care
    crit <- boundary[(y2 - clim2[1]) + 1]
    if (!full_boundary) {
      ## in case the lower boundary of the searched grid is not
      ## zero, then we cannot say anything about cases when the
      ## decision is always negative
      if (!lower.tail) {
        ## in this case the decision changes from negative to
        ## positive when going from small to large
        ## values. Hence, if the decision is always negative,
        ## then we can be sure of that we can never be sure,
        ## but should the decision be negative at all values,
        ## it can change at larger values.
        crit[crit == n1 + 1] <- NA
      } else {
        ## now the decision changes from positive to negative
        ## when going from small to large => should the
        ## decision not change in the clim1 domain then we do
        ## not know if it happens later
        if (clim1[1] > 0) {
          crit[crit == -1] <- NA
        }
        ## however, if crit==Inf then we can be sure that the
        ## decision is indeed always positive
      }
    }
    return(crit)
  }
  decision_boundary
}


## returns a function object which is the decision boundary. That is
## the function finds at a regular grid between llim1 and ulim1 the
## roots of the decision function and returns an interpolation
## function object
solve_boundary2S_normMix <- function(decision, mix1, mix2, n1, n2, lim1, lim2, delta2) {
  grid <- seq(lim2[1], lim2[2], length = diff(lim2) / delta2)

  sigma1 <- sigma(mix1)
  sigma2 <- sigma(mix2)

  sem1 <- sigma1 / sqrt(n1)
  scale1 <- sigma1 / (n1^0.25)

  cond_decisionStep <- function(post2) {
    fn <- function(m1) {
      decision(postmix(mix1, m = m1, se = sem1), post2) - 0.75
    }
    Vectorize(fn)
  }

  Neval <- length(grid)
  # cat("Calculating boundary from", lim2[1], "to", lim2[2], "with", Neval, "points\n")
  tol <- min(delta2 / 100, .Machine$double.eps^0.25)
  ## cat("Using tolerance", tol, "\n")
  crit <- rep(NA, times = Neval)
  for (i in 1:Neval) {
    if (n2 == 0) {
      post2 <- mix2
    } else {
      post2 <- postmix(mix2, m = grid[i], se = sigma2 / sqrt(n2))
    }
    ind_fun <- cond_decisionStep(post2)
    dec_bounds <- ind_fun(lim1)
    ## if decision function is not different at boundaries, lim1
    ## is too narrow and we then enlarge
    while (prod(dec_bounds) > 0) {
      w <- diff(lim1)
      lim1 <- c(lim1[1] - w / 2, lim1[2] + w / 2)
      dec_bounds <- ind_fun(lim1)
    }
    y1c <- uniroot(ind_fun, lim1,
      f.lower = dec_bounds[1], f.upper = dec_bounds[2], tol = tol
    )$root
    crit[i] <- y1c
    ## set lim1 tightly around the current critical value and use the
    ## last boundary limits to not shrink too fast
    lim1 <- c(mean(lim1[1], y1c - 2 * scale1), mean(y1c + 2 * scale1, lim1[2]))
  }

  cbind(grid, crit)
}

#' @templateVar fun decision2S_boundary
#' @template design2S-normal
#' @export
decision2S_boundary.normMix <- function(prior1, prior2, n1, n2, decision, sigma1, sigma2, eps = 1e-6, Ngrid = 10, ...) {
  ## distributions of the means of the data generating distributions
  ## for now we assume that the underlying standard deviation
  ## matches the respective reference scales
  if (missing(sigma1)) {
    sigma1 <- RBesT::sigma(prior1)
    message("Using default prior 1 reference scale ", sigma1)
  }
  if (missing(sigma2)) {
    sigma2 <- RBesT::sigma(prior2)
    message("Using default prior 2 reference scale ", sigma2)
  }
  assert_number(sigma1, lower = 0)
  assert_number(sigma2, lower = 0)

  sem1 <- sigma1 / sqrt(n1)
  sem2 <- sigma2 / sqrt(n2)

  sigma(prior1) <- sigma1
  sigma(prior2) <- sigma2

  ## only n2 can be zero
  assert_that(n1 > 0)
  assert_that(n2 >= 0)

  if (n2 == 0) sem2 <- sigma(prior2) / sqrt(1E-1)

  ## change the reference scale of the prior such that the prior
  ## represents the distribution of the respective means
  mean_prior1 <- prior1
  sigma(mean_prior1) <- sem1
  ## mean_prior2 <- prior2
  ## sigma(mean_prior2) <- sem2

  ## discretization step-size
  delta2 <- sem2 / Ngrid

  ## for the case of mix1 and mix2 having just 1 component, then one
  ## can prove that the decision boundary is a linear function.
  ## Hence we only calculate a very rough grid and apply linear
  ## interpolation.

  linear_boundary <- FALSE
  if (ncol(prior1) == 1 && ncol(prior2) == 1) {
    linear_boundary <- TRUE
    ## we could relax this even further
    delta2 <- sigma2 / Ngrid
  }

  ## the boundary function depends only on the samples sizes n1, n2,
  ## the priors and the decision, but not the assumed truths

  clim2 <- c(Inf, -Inf)

  ## the boundary function which gives conditional on the second
  ## variable the critical value where the decision changes
  boundary <- NA
  boundary_discrete <- matrix(NA, nrow = 0, ncol = 2)

  ## check where the decision is 1, i.e. left or right
  lower.tail <- attr(decision, "lower.tail")

  decision_boundary <- function(y2, lim1) {
    lim2 <- range(y2)

    ## check if boundary function must be recomputed
    if (lim2[1] < clim2[1] | lim2[2] > clim2[2]) {
      new_lim2 <- clim2
      ## note: the <<- assignment is needed to set the variable in the enclosure
      if (missing(lim1)) {
        lim1 <- qmix(mean_prior1, c(eps / 2, 1 - eps / 2))
      }
      if (nrow(boundary_discrete) == 0) {
        ## boundary hasn't been calculated before, do it all
        boundary_discrete <<- solve_boundary2S_normMix(decision, prior1, prior2, n1, n2, lim1, lim2, delta2)
        new_lim2 <- lim2
      } else {
        if (lim2[1] < clim2[1]) {
          ## the lower bound is not low enough... only add the region which is missing
          new_left_lim2 <- min(lim2[1], clim2[1] - 2 * delta2)
          boundary_extra <- solve_boundary2S_normMix(decision, prior1, prior2, n1, n2, lim1, c(new_left_lim2, clim2[1] - delta2), delta2)
          new_lim2[1] <- new_left_lim2
          boundary_discrete <<- rbind(boundary_extra, boundary_discrete)
        }
        if (lim2[2] > clim2[2]) {
          ## the upper bound is not large enough.. again only add whats missing
          new_right_lim2 <- max(lim2[2], clim2[2] + 2 * delta2)
          boundary_extra <- solve_boundary2S_normMix(decision, prior1, prior2, n1, n2, lim1, c(clim2[2] + delta2, new_right_lim2), delta2)
          new_lim2[2] <- new_right_lim2
          boundary_discrete <<- rbind(boundary_discrete, boundary_extra)
        }
      }
      ## only for debugging
      ## assert_that(all(order(boundary_discrete[,1]) == 1:nrow(boundary_discrete)), msg="x grid must stay ordered!")
      if (linear_boundary) {
        boundary <<- approxfun(boundary_discrete[, 1], boundary_discrete[, 2], rule = 2)
      } else {
        boundary <<- splinefun(boundary_discrete[, 1], boundary_discrete[, 2])
      }
      clim2 <<- new_lim2
    }

    return(boundary(y2))
  }

  decision_boundary
}


#' @templateVar fun decision2S_boundary
#' @template design2S-poisson
#' @export
decision2S_boundary.gammaMix <- function(prior1, prior2, n1, n2, decision, eps = 1e-6, ...) {
  assert_that(likelihood(prior1) == "poisson")
  assert_that(likelihood(prior2) == "poisson")

  # only the second n2 argument may be 0
  assert_that(n1 > 0)
  assert_that(n2 >= 0)

  cond_decisionStep <- function(post2) {
    fn <- function(m1) {
      decision(postmix(prior1, n = n1, m = m1 / n1), post2) - 0.25
    }
    Vectorize(fn)
  }

  clim1 <- c(Inf, -Inf)
  clim2 <- c(Inf, -Inf)
  boundary <- NA
  grid <- NA
  lower.tail <- attr(decision, "lower.tail")

  decision_boundary <- function(y2, lim1) {
    if (missing(lim1)) {
      lambda1 <- summary(prior1, probs = c())["mean"] * n1
      lim1 <- qpois(c(eps / 2, 1 - eps / 2), lambda1)
    }

    lim2 <- range(y2)

    assert_number(lim1[1], lower = 0, finite = TRUE)
    assert_number(lim1[2], lower = 0, finite = TRUE)
    assert_number(lim2[1], lower = 0, finite = TRUE)
    assert_number(lim2[2], lower = 0, finite = TRUE)

    ## check if the boundary needs to be recomputed
    if (lim1[1] < clim1[1] | lim1[2] > clim1[2] |
      lim2[1] < clim2[1] | lim2[2] > clim2[2]) {
      ## ensure that lim1 never shrinks in size
      lim1[1] <- min(lim1[1], clim1[1])
      lim1[2] <- max(lim1[2], clim1[2])
      grid <<- lim2[1]:lim2[2]
      Neval <- length(grid)
      boundary <<- rep(NA, Neval)
      for (i in 1:Neval) {
        if (n2 == 0) {
          cond_dec <- cond_decisionStep(prior2)
        } else {
          cond_dec <- cond_decisionStep(postmix(prior2, n = n2, m = grid[i] / n2))
        }
        low <- cond_dec(lim1[1])
        high <- cond_dec(lim1[2])
        if (low < 0 & high < 0) {
          boundary[i] <<- -1
          next
        }
        if (low > 0 & high > 0) {
          boundary[i] <<- Inf
          next
        }
        boundary[i] <<- uniroot_int(cond_dec, lim1,
          f.lower = low,
          f.upper = high
        )
      }

      if (lower.tail) {
        ## if lower.tail==TRUE, then the condition becomes
        ## true when going from large to small values, hence
        ## we need to integrate from 0 to the boundary
        boundary <<- pmax(boundary - 1, -1)
      }

      ## save limits of new grid
      clim1 <<- lim1
      clim2 <<- lim2
    }

    assert_numeric(y2, lower = 0, finite = TRUE, any.missing = FALSE)
    crit <- boundary[y2 - clim2[1] + 1]
    ## in case the lower boundary of the searched grid is not
    ## zero, then we cannot say anything about cases when the
    ## decision is always negative
    if (!lower.tail) {
      ## in this case the decision changes from negative to
      ## positive when going from small to large values. Hence,
      ## if the decision is always negative, then we can be sure
      ## that the decision changes past lim1[2]. We set the
      ## boundary to lim1[2]+1 if lim1 has been given; otherwise
      ## to NA.
      ## Should the decision be negative at all values,
      ## it can change at larger values.
      if (missing(lim1)) {
        crit[crit == -1] <- NA
        crit[crit == Inf] <- NA
      } else {
        crit[crit == -1] <- lim1[2] + 1
        crit[crit == Inf] <- lim1[1] - 1
      }
    } else {
      ## now the decision changes from positive to negative
      ## when going from small to large => should the
      ## decision not change in the clim1 domain then we do
      ## not know if it happens later
      if (clim1[1] > 0) {
        crit[crit == -1] <- NA
      }
      ## however, if crit==Inf then we can be sure that the
      ## decision is indeed always positive
    }
    return(crit)
  }
  decision_boundary
}
