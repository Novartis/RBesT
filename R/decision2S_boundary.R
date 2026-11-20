#' Decision Boundary for 2 Sample Designs
#'
#' The `decision2S_boundary` function defines a 2 sample design
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
#' `lower.tail=TRUE` (otherwise the definition is \eqn{D_1(y_2) =
#' \max_{y_1}\{D(y_1,y_2) = 0\}}{D_1(y_2) = max_{y_1}{D(y_1,y_2) =
#' 0}}). The decision function may change at most at a single critical
#' value for given \eqn{y_{2}} as only one-sided decision functions
#' are supported. Here, \eqn{y_2} is defined for binary and Poisson
#' endpoints as the sufficient statistic \eqn{y_2 = \sum_{i=1}^{n_2}
#' y_{2,i}} and for the normal case as the mean \eqn{\bar{y}_2 = 1/n_2
#' \sum_{i=1}^{n_2} y_{2,i}}.
#'
#' @return For one-sided decision functions, returns a function with a
#' single argument. This function calculates in dependence of the
#' outcome \eqn{y_2} in sample 2 the critical value \eqn{y_{1,c}} for
#' which the defined design will change the decision from 0 to 1 (or
#' vice versa, depending on the decision function).
#' For two-sided decision functions, returns a list with components
#' `lower_or_equal_than` and `higher_than`, containing the critical
#' value functions for the lower and upper one-sided decision
#' boundary components.
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
decision2S_boundary <- function(prior1, prior2, n1, n2, decision, ...) {
  UseMethod("decision2S_boundary")
}
#' @export
decision2S_boundary.default <- function(prior1, prior2, n1, n2, decision, ...) {
  "Unknown density"
}

#' @templateVar fun decision2S_boundary
#' @template design2S-binomial
# If the \code{eps} argument is specificed, then the
# returned function will use the additional \code{lim2}
# argument to limit the search for the critical value.
#' @export
decision2S_boundary.betaMix <- function(
  prior1,
  prior2,
  n1,
  n2,
  decision,
  eps,
  ...
) {
  ## only n2=0 is supported
  assert_number(n1, lower = 1, finite = TRUE)
  assert_number(n2, lower = 0, finite = TRUE)

  if (is(decision, "decision2S_2sided")) {
    decision2S_boundary_betaMix_2sided(
      prior1,
      prior2,
      n1,
      n2,
      decision,
      eps,
      ...
    )
  } else {
    decision2S_boundary_betaMix_1sided(
      prior1,
      prior2,
      n1,
      n2,
      decision,
      eps,
      ...
    )
  }
}

decision2S_boundary_betaMix_2sided <- function(
  prior1,
  prior2,
  n1,
  n2,
  decision,
  eps,
  ...
) {
  crit_lower <- decision2S_boundary_betaMix_atomic(
    prior1,
    prior2,
    n1,
    n2,
    lower(decision),
    eps,
    ...
  )
  crit_upper <- decision2S_boundary_betaMix_atomic(
    prior1,
    prior2,
    n1,
    n2,
    upper(decision),
    eps,
    ...
  )
  list(lower_or_equal_than = crit_lower, higher_than = crit_upper)
}

decision2S_boundary_betaMix_1sided <- function(
  prior1,
  prior2,
  n1,
  n2,
  decision,
  eps,
  ...
) {
  decision <- if (has_lower(decision)) {
    lower(decision)
  } else {
    upper(decision)
  }
  decision2S_boundary_betaMix_atomic(
    prior1,
    prior2,
    n1,
    n2,
    decision,
    eps,
    ...
  )
}

decision2S_boundary_betaMix_atomic <- function(
  prior1,
  prior2,
  n1,
  n2,
  decision,
  eps,
  ...
) {
  if (!missing(eps)) {
    assert_number(eps, lower = 0, upper = 0.1, finite = TRUE)
  }

  cond_decisionDist <- function(post2cond) {
    function(m1) {
      ## Note: Subtracting from the decision 0.25 leads to
      ## negative decisions being at -0.25 while positives are
      ## at 0.75; since uniroot_int *always* returns the x which
      ## has lowest absolute value we are guaranteed that y2crit
      ## is just before the jump
      ## decision(post1cond, post2[[m2+1]]) - 0.25
      decision(postmix(prior1, r = m1, n = n1), post2cond) - 0.25
      ## decision(postmix(prior2, r=m2, n=n2), post1cond) - 0.25
    }
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
    ## The step objective is monotone in m1 in the direction set by
    ## lower.tail; the critical y1 itself is monotone in y2 for the
    ## supported one-sided decisions.
    extend_dir <- if (lower.tail) "downX" else "upX"

    decFun_for <- function(y2) {
      if (n2 == 0) {
        cond_decisionDist(prior2)
      } else {
        cond_decisionDist(postmix(prior2, r = y2, n = n2))
      }
    }
    ## raw integer critical value at y2 within a warm bracket;
    ## uniroot_int extends it (bounded by lim1) if the boundary moved
    ## further and reports +/-Inf where the decision is constant.
    solve_raw <- function(y2, bracket) {
      uniroot_int(decFun_for(y2), bracket, extendInt = extend_dir, clamp = lim1)
    }
    ## map the constant-decision sentinels into the ordered integer
    ## domain [-1, n1 + 1] so they participate in the monotone fill:
    ## never true (-Inf -> -1), always true (Inf -> n1 + 1).
    map_root <- function(root) {
      if (is.finite(root)) {
        root
      } else if (root < 0) {
        -1
      } else {
        n1 + 1
      }
    }
    set_val <- function(y2, val) {
      boundary[y2 - lim2[1] + 1] <<- val
    }

    if (getOption("RBesT.decision2S_boundary", "adaptive") == "adaptive") {
      ## Adaptive solve replacement. Because the critical value is
      ## monotone in y2, if the roots at the two ends of a y2 subinterval
      ## are equal then every y2 between them shares that value -- the run
      ## is filled without any (expensive pmixdiff) uniroot_int calls.
      ## Otherwise the midpoint is solved, tightly bracketed by the two
      ## neighbour roots (crit(m) is guaranteed to lie between them), and
      ## the two halves recurse. Steep boundaries end up solving every
      ## point; flat runs cost O(log length).
      bracket_for <- function(ra, rb) {
        lo <- max(min(ra, rb), lim1[1])
        hi <- min(max(ra, rb), lim1[2])
        if (lo >= hi) lim1 else c(lo, hi)
      }
      fill <- function(a, ra, b, rb) {
        if (ra == rb) {
          for (y2 in a:b) {
            set_val(y2, ra)
          }
          return(invisible())
        }
        if (b - a <= 1L) {
          set_val(a, ra)
          set_val(b, rb)
          return(invisible())
        }
        m <- (a + b) %/% 2L
        rm <- map_root(solve_raw(m, bracket_for(ra, rb)))
        fill(a, ra, m, rm)
        fill(m, rm, b, rb)
      }
      ra <- map_root(solve_raw(lim2[1], lim1))
      rb <- map_root(solve_raw(lim2[2], lim1))
      fill(lim2[1], ra, lim2[2], rb)
    } else {
      ## Legacy sequential sweep: warm-start each y2 point's integer
      ## search with a tight bracket around the previous boundary.
      prev <- NA
      for (y2 in lim2[1]:lim2[2]) {
        bracket <- if (is.na(prev)) {
          lim1
        } else {
          c(max(prev - 2, lim1[1]), min(prev + 2, lim1[2]))
        }
        root <- solve_raw(y2, bracket)
        set_val(y2, map_root(root))
        prev <- if (is.finite(root)) root else NA
      }
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
      if (
        lim1[1] < clim1[1] |
          lim1[2] > clim1[2] |
          lim2[1] < clim2[1] |
          lim2[2] > clim2[2]
      ) {
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


boundary_extreme_value <- 1e100

constant_boundary_value <- function(decision_true, lower.tail) {
  if (identical(isTRUE(decision_true), isTRUE(lower.tail))) {
    boundary_extreme_value
  } else {
    -boundary_extreme_value
  }
}

is_constant_boundary <- function(x) {
  is.finite(x) & abs(x) >= boundary_extreme_value
}

## With a family-based sigma the decision can change more than once along
## the outcome of sample 1, while the boundary holds a single critical
## value. Designs warn if the outcomes for which the decision disagrees
## with the boundary have a probability exceeding this tolerance.
boundary2S_misfit_tol <- 0.005

warn_boundary2S_misfit <- function(where, what, mass) {
  warning(
    "The decision changes more than once along the outcome of sample 1 ",
    where,
    ", which the decision boundary does not account for. ",
    what,
    " inaccurate by a probability of roughly ",
    signif(mass, 1),
    " as estimated on a coarse grid of outcomes. ",
    "Consider to use the exact design for the endpoint.",
    call. = FALSE
  )
}

## Adaptive boundary tracer for the normal endpoint.
##
## Recursively bisects the y2 interval instead of sweeping a fixed grid.
## Each new point is seeded from its two bracketing solved neighbours and
## root-found within their two roots as bracket -- guaranteed to contain
## the root because y1_c(y2) is monotone in y2 for the supported one-sided
## decisions; the solver's extendInt remains the safety net. delta2 is the
## resolution floor.
##
## Acceptance uses a nested (quarter-point) error estimate rather than a
## single-midpoint linearity check. For a subinterval [a, b] with solved
## endpoints and midpoint (a, ra), (m, rm), (b, rb) we form the quadratic
## Lagrange interpolant Q through those three points, solve the two quarter
## points q1 = (a+m)/2, q3 = (m+b)/2, and accept the subinterval when
## max(|rq1 - Q(q1)|, |rq3 - Q(q3)|) < lin_tol; otherwise we subdivide.
## The quarter points q1, q3 of a parent are exactly the midpoints of its
## two children, so with memoization no y2 is ever solved twice.
##
## The residual test alone can alias: on a wide interval all five samples
## may fall on a common line while the boundary bends in between (mixture
## priors change the slope between the regime of the informative component
## and the tails). The coverage guard max_gap therefore forbids accepting
## an interval while its knot spacing (b - a) / 4 exceeds max_gap; the
## caller sets it to the scale below which the boundary has no features,
## so that the residual test only ever decides at feature scale.
##
## known is an optional cbind(y2, crit) of already solved roots; those
## inside lim2 seed the memo and are not solved again.
##
## Returns cbind(grid, crit), sorted and de-duplicated for splinefun, or
## NULL if a constant/no-root (NA) region is hit -- the caller then falls
## back to the robust uniform grid sweep.
trace_boundary2S_adaptive <- function(
  solve_at,
  lim1,
  lim2,
  delta2,
  scale1,
  lin_tol,
  max_gap = Inf,
  known = NULL
) {
  ## memoized solver: one uniroot per distinct y2, ever -- a quarter point
  ## of a parent interval is reused as the midpoint of a child interval.
  cache_y <- numeric(0)
  cache_r <- numeric(0)
  if (!is.null(known)) {
    seed <- known[, 1] >= lim2[1] - 1e-9 & known[, 1] <= lim2[2] + 1e-9
    cache_y <- known[seed, 1]
    cache_r <- known[seed, 2]
  }
  ## a constant-decision knot carried over from an earlier request is no
  ## root either and must send the caller to the grid sweep just the same
  bad <- any(is_constant_boundary(cache_r))

  solve_memo <- function(y2, bracket) {
    hit <- which(abs(cache_y - y2) < 1e-9)
    if (length(hit)) {
      return(cache_r[hit[1L]])
    }
    r <- solve_at(y2, bracket)
    if (is.na(r) || is_constant_boundary(r)) {
      bad <<- TRUE
    }
    cache_y[[length(cache_y) + 1L]] <<- y2
    cache_r[[length(cache_r) + 1L]] <<- r
    r
  }

  ## bracket for a point known to lie between two solved neighbours;
  ## monotonicity of y1_c(y2) guarantees the root is inside [r_lo, r_hi].
  ## A degenerate (equal-root) bracket falls back to a symmetric scale1
  ## window so uniroot's extendInt can recover.
  bracket_between <- function(r_lo, r_hi) {
    br <- sort(c(r_lo, r_hi))
    if (diff(br) < 1e-9) {
      mean(br) + c(-2, 2) * scale1
    } else {
      br
    }
  }

  ## quadratic Lagrange interpolant through (a, ra), (m, rm), (b, rb)
  quad <- function(x, a, ra, m, rm, b, rb) {
    ra *
      ((x - m) * (x - b)) /
      ((a - m) * (a - b)) +
      rm * ((x - a) * (x - b)) / ((m - a) * (m - b)) +
      rb * ((x - a) * (x - m)) / ((b - a) * (b - m))
  }

  ## recurse over [a, b]; the midpoint m is already solved and passed in
  ## (it was a quarter point of the parent -- reused, not re-solved).
  recurse <- function(a, ra, m, rm, b, rb) {
    if (bad || (b - a) <= delta2) {
      return(invisible())
    }
    q1 <- (a + m) / 2
    q3 <- (m + b) / 2
    rq1 <- solve_memo(q1, bracket_between(ra, rm))
    rq3 <- solve_memo(q3, bracket_between(rm, rb))
    if (bad) {
      return(invisible())
    }
    ## nested error estimate: solved roots vs the quadratic they did not build
    err <- max(
      abs(rq1 - quad(q1, a, ra, m, rm, b, rb)),
      abs(rq3 - quad(q3, a, ra, m, rm, b, rb))
    )
    if (err > lin_tol || (b - a) / 4 > max_gap) {
      ## q1, q3 become the midpoints of the two children.
      recurse(a, ra, q1, rq1, m, rm)
      recurse(m, rm, q3, rq3, b, rb)
    }
  }

  a <- lim2[1]
  b <- lim2[2]
  ra <- solve_memo(a, lim1)
  rb <- solve_memo(b, lim1)
  if (!bad) {
    m <- (a + b) / 2
    rm <- solve_memo(m, bracket_between(ra, rb))
    if (!bad) {
      recurse(a, ra, m, rm, b, rb)
    }
  }

  if (bad) {
    return(NULL)
  }
  ord <- order(cache_y)
  y <- cache_y[ord]
  r <- cache_r[ord]
  keep <- c(TRUE, diff(y) > 1e-9)
  cbind(y[keep], r[keep])
}

## returns a function object which is the decision boundary. That is
## the function finds at a regular grid between llim1 and ulim1 the
## roots of the decision function and returns an interpolation
## function object
solve_boundary2S_normMix <- function(
  decision,
  mix1,
  mix2,
  n1,
  n2,
  lim1,
  lim2,
  delta2,
  sigma_fun1 = NULL,
  sigma_fun2 = NULL,
  y2_at = NULL,
  max_gap = Inf,
  known = NULL
) {
  assert_class(decision, "decision2S_atomic")

  sigma2 <- sigma(mix2)

  if (is.null(sigma_fun1)) {
    ## sigma is known and fixed
    sigma1 <- sigma(mix1)
    sem1 <- sigma1 / sqrt(n1)
    scale1 <- sigma1 / (n1^0.25)

    cond_decisionStep <-
      function(post2) {
        function(m1) {
          min(decision(postmix(mix1, m = m1, se = sem1), post2, dist = TRUE))
        }
      }
  } else {
    sigma1 <- sigma_fun1(summary(mix1)["mean"])
    sem1 <- sigma1 / sqrt(n1)
    scale1 <- sigma1 / (n1^0.25)

    ## Search window for the critical value. It must not depend on the
    ## call (lim1), since the solved knots are cached for all later calls.
    ## An outcome m1 is uninformative once its standard error exceeds the
    ## prior standard deviation by far; the posterior is the prior from
    ## there on and so is the decision. Such an edge bounds the window. A
    ## side without such an edge (e.g. large counts) is informative
    ## throughout and the objective is monotone beyond the prior means.
    ## Only the non-monotone part between the edge and the prior means is
    ## scanned on such a side.
    prior1_mean <- summary(mix1)[["mean"]]
    uninformative_se <- 100 * max(mix1[3, ])
    uninformative <- function(m1) {
      se_m1 <- sigma_fun1(m1) / sqrt(n1)
      is.na(se_m1) || se_m1 >= uninformative_se
    }
    informative_edge <- function(direction) {
      inner <- prior1_mean
      for (k in 0:60) {
        outer <- prior1_mean + direction * sigma1 * 2^k
        if (uninformative(outer)) {
          for (i in 1:30) {
            mid <- (inner + outer) / 2
            if (uninformative(mid)) outer <- mid else inner <- mid
          }
          ## the standard error approaches the limit continuously at an
          ## edge; a jump is a numerical breakdown of the family functions
          if (sigma_fun1(inner) / sqrt(n1) < uninformative_se / 2) {
            break
          }
          return(inner)
        }
        inner <- outer
      }
      direction * Inf
    }
    if (uninformative(prior1_mean)) {
      window1 <- prior1_mean + c(-8, 8) * sigma1
      edge_found <- c(TRUE, TRUE)
    } else {
      window1 <- c(informative_edge(-1), informative_edge(1))
      edge_found <- is.finite(window1)
    }
    ## An informative side is limited generously with respect to the prior;
    ## the decision counts as constant if it changes only beyond.
    far1 <- range(mix1[2, ]) + c(-8, 8) * max(sigma1, mix1[3, ])
    window1 <- ifelse(edge_found, window1, far1)
    scan1 <- range(mix1[2, ]) + c(-8, 8) * sem1
    scan1 <- ifelse(edge_found, window1, scan1)

    cond_decisionStep <- function(post2) {
      ## Bound uniroot's interval extension to the search window. Returning
      ## NA outside it stops extension; the fallback scans within it.
      function(m1) {
        if (m1 < window1[1] || m1 > window1[2]) {
          return(NA_real_)
        }
        se_m1 <- sigma_fun1(m1) / sqrt(n1)
        min(decision(postmix(mix1, m = m1, se = se_m1), post2, dist = TRUE))
      }
    }
  }

  tol <- min(delta2 / 100, .Machine$double.eps^0.25)

  ## Direction of the requested crossing: decreasing for lower.tail = TRUE,
  ## increasing otherwise. With constant sigma the objective is monotone;
  ## family-based sigma can require the non-monotone fallback below.
  extend_dir <- if (isTRUE(attr(decision, "lower.tail"))) "downX" else "upX"

  ## Root-find the critical y1 for a second-sample outcome y2 within a
  ## warm bracket (extendInt widens it if the boundary moved out).
  ## If there is no sign change, retain whether the decision is always true
  ## or always false as a tail-aware extreme boundary.
  solve_at <- function(y2, bracket) {
    if (n2 == 0) {
      post2 <- mix2
    } else {
      se2 <- if (is.null(sigma_fun2)) {
        sigma2 / sqrt(n2)
      } else {
        sigma_fun2(y2) / sqrt(n2)
      }
      post2 <- postmix(mix2, m = y2, se = se2)
    }
    ind_fun <- cond_decisionStep(post2)
    if (!is.null(sigma_fun1)) {
      bracket <- pmin(pmax(bracket, window1[1]), window1[2])
    }
    root <- tryCatch(
      uniroot(ind_fun, bracket, extendInt = extend_dir, tol = tol)$root,
      error = function(e) NA_real_
    )
    if (!is.na(root)) {
      return(root)
    }

    ## With a family-based standard error the objective can revert toward the
    ## prior in the tails and is not globally monotone. A warm bracket may
    ## therefore miss the crossing in the requested direction. Scan the
    ## window on a uniform grid before classifying the decision as constant.
    ## This is a numerical search, not proof that no narrower excursion
    ## exists.
    if (!is.null(sigma_fun1)) {
      m1 <- seq(scan1[1], scan1[2], length.out = 65L)
      val <- vapply(m1, ind_fun, numeric(1))
      if (anyNA(val)) {
        ## a dropped knot would be replaced by its neighbours silently
        stop(
          "The decision is undefined in the search window of sample 1 for ",
          "the outcome ",
          signif(y2, 6),
          " of sample 2.",
          call. = FALSE
        )
      }
      n_scan <- length(m1)
      up <- extend_dir == "upX"
      cross <- which(
        if (up) {
          val[-n_scan] <= 0 & val[-1L] >= 0
        } else {
          val[-n_scan] >= 0 & val[-1L] <= 0
        }
      )
      if (length(cross)) {
        k <- if (up) cross[length(cross)] else cross[1L]
        return(uniroot(ind_fun, m1[k + 0:1], tol = tol)$root)
      }
      ## An excursion between two grid points shows as the extremum of the
      ## scan which is closest to zero; look for it there. This is done
      ## only if the scan sees no change of sign at all: a further narrow
      ## excursion next to a crossing seen by the scan is not searched for,
      ## as a single critical value cannot represent it anyway.
      if (all(val < 0) || all(val > 0)) {
        bump <- val[1L] < 0
        k <- if (bump) which.max(val) else which.min(val)
        cell <- m1[c(max(k - 1L, 1L), min(k + 1L, n_scan))]
        extremum <- optimize(ind_fun, cell, maximum = bump, tol = tol)
        if (bump == (extremum$objective > 0)) {
          bracket <- if (up == bump) {
            c(cell[1L], extremum[[1L]])
          } else {
            c(extremum[[1L]], cell[2L])
          }
          return(uniroot(ind_fun, bracket, tol = tol)$root)
        }
      }
      ## No crossing inside the scanned range. On a side which is
      ## informative throughout the objective is monotone beyond the scan,
      ## so a crossing there is found between the end of the scan and the
      ## window if the objective has not crossed yet.
      at_end <- val[c(1L, n_scan)]
      not_crossed <- if (extend_dir == "upX") at_end < 0 else at_end > 0
      for (side in 1:2) {
        ## below the scan the crossing lies behind, above it lies ahead
        pending <- if (side == 1) !not_crossed[1] else not_crossed[2]
        if (edge_found[side] || !pending) {
          next
        }
        outward <- sort(c(scan1[side], window1[side]))
        root <- tryCatch(
          uniroot(ind_fun, outward, tol = tol)$root,
          error = function(e) NA_real_
        )
        if (!is.na(root)) {
          return(root)
        }
      }
      ## Without a crossing in the requested direction the decision is
      ## constant if both ends agree; otherwise it changes in the reverse
      ## direction only, which a critical value cannot represent.
      if (sign(at_end[1]) != sign(at_end[2])) {
        return(NA_real_)
      }
      return(constant_boundary_value(
        at_end[1] > 0,
        attr(decision, "lower.tail")
      ))
    }

    decision_at_center <- ind_fun(mean(lim1))
    if (
      length(decision_at_center) != 1L ||
        is.na(decision_at_center)
    ) {
      return(NA_real_)
    }
    constant_boundary_value(
      decision_at_center > 0,
      attr(decision, "lower.tail")
    )
  }

  ## Direct solve at caller-specified y2 points (used by the exact
  ## linear-boundary fast path, which needs only two points to define a
  ## provably-straight boundary). Bypasses both the adaptive tracer and the
  ## uniform grid. Points with no root (NA) are dropped by the caller.
  if (!is.null(y2_at)) {
    crit <- vapply(y2_at, solve_at, numeric(1), bracket = lim1)
    keep <- !is.na(crit)
    return(cbind(y2_at[keep], crit[keep]))
  }

  ## Adaptive tracing is the default; the uniform grid sweep is kept as a
  ## fallback (constant/no-root regions) and as an opt-out via
  ## options(RBesT.decision2S_boundary = "grid"). lin_tol is the accepted
  ## deviation of the boundary from its quadratic fit at the quarter points
  ## of a subinterval; it is set relative to sem1, the posterior scale of
  ## the arm-1 mean that governs the boundary's own resolution ("boundary
  ## located to ~1% of its natural sem1 scale"). The test is only reliable
  ## at feature scale, which the coverage guard max_gap enforces.
  if (getOption("RBesT.decision2S_boundary", "adaptive") == "adaptive") {
    adaptive <- trace_boundary2S_adaptive(
      solve_at,
      lim1,
      lim2,
      delta2,
      scale1,
      lin_tol = sem1 / 100,
      max_gap = max_gap,
      known = known
    )
    if (!is.null(adaptive)) {
      return(adaptive)
    }
  }

  ## Uniform grid sweep: a fixed grid of y2 points, warm-started with a
  ## symmetric +/-2*scale1 bracket around the previous finite root.
  grid <- if (diff(lim2) == 0) {
    lim2[1]
  } else {
    seq(lim2[1], lim2[2], length = max(2, diff(lim2) / delta2))
  }
  crit <- rep(NA_real_, length(grid))
  solved <- rep(FALSE, length(grid))
  ## A first pass visits points which are max_gap apart, the scale below
  ## which the boundary has no features. The decision is taken as constant
  ## in between two such points whenever it is the same constant at both,
  ## which saves the costly search for a crossing at the points in between.
  stride <- 1L
  if (length(grid) > 2 && is.finite(max_gap)) {
    stride <- max(1L, floor(max_gap / (grid[2] - grid[1])))
  }
  coarse <- unique(c(seq(1L, length(grid), by = stride), length(grid)))
  bracket <- lim1
  for (i in c(coarse, setdiff(seq_along(grid), coarse))) {
    if (solved[i]) {
      next
    }
    if (i > 1 && solved[i - 1] && is.finite(crit[i - 1]) &&
          !is_constant_boundary(crit[i - 1])) {
      bracket <- crit[i - 1] + c(-2, 2) * scale1
    } else {
      bracket <- lim1
    }
    crit[i] <- solve_at(grid[i], bracket)
    solved[i] <- TRUE
    if (i == coarse[length(coarse)] && stride > 1L) {
      for (k in seq_len(length(coarse) - 1L)) {
        ends <- crit[coarse[k + 0:1]]
        between <- seq_len(coarse[k + 1L] - coarse[k] - 1L) + coarse[k]
        if (
          length(between) &&
            !anyNA(ends) &&
            all(is_constant_boundary(ends)) &&
            ends[1] == ends[2]
        ) {
          crit[between] <- ends[1]
          solved[between] <- TRUE
        }
      }
    }
  }

  ## Locate transitions between finite-root and constant-decision runs more
  ## accurately than the base grid. This keeps the piecewise consumer from
  ## assigning a whole grid cell to the wrong constant state.
  transition_tol <- delta2 / 100
  repeat {
    status <- ifelse(is_constant_boundary(crit), sign(crit), 0)
    transition <- which(
      !is.na(crit[-length(crit)]) &
        !is.na(crit[-1L]) &
        status[-length(status)] != status[-1L] &
        diff(grid) > transition_tol
    )
    if (!length(transition)) {
      break
    }
    midpoint <- (grid[transition] + grid[transition + 1L]) / 2
    midpoint_crit <- vapply(
      midpoint,
      solve_at,
      numeric(1),
      bracket = lim1
    )
    classified <- !is.na(midpoint_crit)
    if (!any(classified)) {
      break
    }
    grid <- c(grid, midpoint[classified])
    crit <- c(crit, midpoint_crit[classified])
    ord <- order(grid)
    grid <- grid[ord]
    crit <- crit[ord]
  }

  ## Drop only unclassifiable failures. Constant-decision points are retained
  ## as extreme values and handled piecewise by fit_boundary_consumer().
  keep <- !is.na(crit)
  cbind(grid[keep], crit[keep])
}

#' @templateVar fun decision2S_boundary
#' @template design2S-normal
#' @export
decision2S_boundary.normMix <- function(
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

  if (is(decision, "decision2S_2sided")) {
    decision2S_boundary_normMix_2sided(
      prior1,
      prior2,
      n1,
      n2,
      decision,
      sigma1,
      sigma2,
      eps,
      Ngrid,
      sigma_fun1 = sigma_fun1,
      sigma_fun2 = sigma_fun2,
      ...
    )
  } else {
    decision2S_boundary_normMix_1sided(
      prior1,
      prior2,
      n1,
      n2,
      decision,
      sigma1,
      sigma2,
      eps,
      Ngrid,
      sigma_fun1 = sigma_fun1,
      sigma_fun2 = sigma_fun2,
      ...
    )
  }
}

decision2S_boundary_normMix_2sided <- function(
  prior1,
  prior2,
  n1,
  n2,
  decision,
  sigma1,
  sigma2,
  eps,
  Ngrid,
  sigma_fun1 = NULL,
  sigma_fun2 = NULL,
  ...
) {
  crit_lower <- decision2S_boundary_normMix_atomic(
    prior1,
    prior2,
    n1,
    n2,
    lower(decision),
    sigma1,
    sigma2,
    eps,
    Ngrid,
    sigma_fun1 = sigma_fun1,
    sigma_fun2 = sigma_fun2,
    ...
  )
  crit_upper <- decision2S_boundary_normMix_atomic(
    prior1,
    prior2,
    n1,
    n2,
    upper(decision),
    sigma1,
    sigma2,
    eps,
    Ngrid,
    sigma_fun1 = sigma_fun1,
    sigma_fun2 = sigma_fun2,
    ...
  )
  list(lower_or_equal_than = crit_lower, higher_than = crit_upper)
}

decision2S_boundary_normMix_1sided <- function(
  prior1,
  prior2,
  n1,
  n2,
  decision,
  sigma1,
  sigma2,
  eps,
  Ngrid,
  sigma_fun1 = NULL,
  sigma_fun2 = NULL,
  ...
) {
  decision <- if (has_lower(decision)) {
    lower(decision)
  } else {
    upper(decision)
  }
  decision2S_boundary_normMix_atomic(
    prior1,
    prior2,
    n1,
    n2,
    decision,
    sigma1,
    sigma2,
    eps,
    Ngrid,
    sigma_fun1 = sigma_fun1,
    sigma_fun2 = sigma_fun2,
    ...
  )
}

## Union of intervals given as rows (lower, upper) of a 2-column matrix;
## returns sorted, disjoint intervals.
merge_intervals <- function(iv) {
  n <- nrow(iv)
  if (n < 2) {
    return(iv)
  }
  iv <- iv[order(iv[, 1]), , drop = FALSE]
  upper <- cummax(iv[, 2])
  start <- c(TRUE, iv[-1, 1] > upper[-n])
  cbind(iv[start, 1], upper[c(start[-1], TRUE)])
}

## Parts of the interval lim = c(lower, upper) which are not covered by the
## sorted, disjoint intervals in resolved; slivers narrower than tol count
## as covered.
uncovered_intervals <- function(lim, resolved, tol) {
  lower <- pmax(c(-Inf, resolved[, 2]), lim[1])
  upper <- pmin(c(resolved[, 1], Inf), lim[2])
  keep <- upper - lower > tol
  cbind(lower[keep], upper[keep])
}

## Regions of a normal mixture which hold all but eps of the probability
## mass of each component, as rows (lower, upper).
mass_intervals <- function(mix, eps) {
  z <- qnorm(1 - eps / 2)
  cbind(mix[2, ] - z * mix[3, ], mix[2, ] + z * mix[3, ])
}

## The cache (resolve / eval functions) behind a normal boundary function
## returned by decision2S_boundary_normMix_atomic.
boundary2S_cache <- function(boundary_fun) {
  environment(boundary_fun)$cache
}

## Build the interpolating boundary function from the solved (y2, y1c)
## points produced by solve_boundary2S_normMix. This is the single seam
## through which the boundary consumer is chosen:
##   - constant-decision runs -> exact extreme-value function on that run.
##   - `linear`     -> piecewise-linear approxfun for the provably-linear
##     single-component / constant-sigma case.
##   - otherwise    -> smooth spline through the (non-uniformly graded) points.
fit_boundary_consumer <- function(boundary_discrete, linear, decision) {
  if (nrow(boundary_discrete) == 0) {
    stop("No classifiable decision boundary points were found.")
  }

  boundary_discrete <- boundary_discrete[
    order(boundary_discrete[, 1]),
    ,
    drop = FALSE
  ]
  keep <- !duplicated(boundary_discrete[, 1], fromLast = TRUE)
  boundary_discrete <- boundary_discrete[keep, , drop = FALSE]

  status <- ifelse(
    is_constant_boundary(boundary_discrete[, 2]),
    sign(boundary_discrete[, 2]),
    0
  )
  run_id <- cumsum(c(TRUE, diff(status) != 0))
  runs <- split(seq_len(nrow(boundary_discrete)), run_id)

  consumers <- lapply(runs, function(idx) {
    x <- boundary_discrete[idx, 1]
    y <- boundary_discrete[idx, 2]
    if (status[idx[1]] != 0 || length(idx) == 1L) {
      value <- y[1]
      return(function(y2) rep(value, length(y2)))
    }
    if (linear) {
      return(approxfun(x, y, rule = 2))
    }
    ## Shape-preserving monotone Hermite (Fritsch-Carlson): matches the
    ## adaptive tracer's curvature-controlled point placement, avoids the
    ## overshoot a global cubic ("fmm") exhibits over the non-uniformly
    ## graded points, and respects the boundary's monotonicity where it
    ## holds without imposing it where it does not.
    splinefun(x, y, method = "monoH.FC")
  })

  breaks <- if (length(runs) == 1L) {
    numeric(0)
  } else {
    vapply(seq_len(length(runs) - 1L), function(i) {
      left <- runs[[i]][length(runs[[i]])]
      right <- runs[[i + 1L]][1]
      mean(boundary_discrete[c(left, right), 1])
    }, numeric(1))
  }

  function(y2) {
    selected_run <- findInterval(y2, breaks) + 1L
    value <- numeric(length(y2))
    for (i in seq_along(consumers)) {
      use <- selected_run == i
      value[use] <- consumers[[i]](y2[use])
    }
    value
  }
}

decision2S_boundary_normMix_atomic <- function(
  prior1,
  prior2,
  n1,
  n2,
  decision,
  sigma1,
  sigma2,
  eps,
  Ngrid,
  sigma_fun1 = NULL,
  sigma_fun2 = NULL,
  ...
) {
  assert_class(decision, "decision2S_atomic")

  assert_number(sigma1, lower = 0)
  assert_number(sigma2, lower = 0)

  sem1 <- sigma1 / sqrt(n1)
  sem2 <- sigma2 / sqrt(n2)

  sigma(prior1) <- sigma1
  sigma(prior2) <- sigma2

  ## only n2 can be zero
  assert_that(n1 > 0)
  assert_that(n2 >= 0)

  if (n2 == 0) {
    sem2 <- sigma(prior2) / sqrt(1E-1)
  }

  ## change the reference scale of the prior such that the prior
  ## represents the distribution of the respective means
  mean_prior1 <- prior1
  sigma(mean_prior1) <- sem1

  ## discretization step-size
  delta2 <- sem2 / Ngrid

  ## for the case of mix1 and mix2 having just 1 component, then one
  ## can prove that the decision boundary is a linear function.
  ## Hence we only calculate a very rough grid and apply linear
  ## interpolation.
  ## With family path, linearity no longer holds even for single-component priors.

  linear_boundary <- FALSE
  if (ncol(prior1) == 1 && ncol(prior2) == 1 && is.null(sigma_fun1)) {
    linear_boundary <- TRUE
    ## we could relax this even further
    delta2 <- sigma2 / Ngrid
  }

  ## the boundary function depends only on the samples sizes n1, n2,
  ## the priors and the decision, but not the assumed truths

  ## Cache of the boundary:
  ##   knots    -- all solved roots (y2, y1c), sorted in y2
  ##   resolved -- sorted disjoint y2 intervals (lower, upper) in which the
  ##               knots are dense, i.e. traced with the coverage guard
  ##               max_gap. Knots outside are never relied upon.
  ## sem2 is the scale over which the outcome of sample 2 can change the
  ## posterior mixture weights, hence the boundary has no features below it.
  max_gap <- sem2
  knots <- matrix(NA_real_, nrow = 0, ncol = 2)
  resolved <- matrix(NA_real_, nrow = 0, ncol = 2)

  ## the boundary function which gives conditional on the second
  ## variable the critical value where the decision changes
  boundary <- NA
  exact_line_built <- FALSE
  linear_fallback <- FALSE

  ## Ensure that the boundary is resolved on the y2 intervals given as rows
  ## (lower, upper); only parts not resolved before are solved.
  resolve <- function(intervals, lim1) {
    if (exact_line_built) {
      return(invisible())
    }
    if (missing(lim1)) {
      lim1 <- qmix(mean_prior1, c(eps / 2, 1 - eps / 2))
    }
    ## Exact linear fast path (C1/C2/C3): for single-component priors with
    ## constant sigma the boundary y1_c(y2) is provably a straight line, and
    ## for multi-criteria decisions the per-criterion lines are parallel, so
    ## their min/max envelope is still a single line (no kink). We therefore
    ## solve at just two anchor points, fit the exact line, and cache it as a
    ## closed-form function valid for all y2 -- bypassing the tracing and
    ## the resolved intervals entirely. n2 == 0 collapses to a horizontal
    ## line (two equal anchors -> slope 0) automatically.
    if (linear_boundary && !linear_fallback) {
      anchor2 <- range(intervals)
      if (diff(anchor2) == 0) {
        anchor2 <- anchor2 + c(-1, 1) * max(delta2, abs(anchor2[1]), 1)
      }
      pts <- solve_boundary2S_normMix(
        decision,
        prior1,
        prior2,
        n1,
        n2,
        lim1,
        anchor2,
        delta2,
        sigma_fun1 = sigma_fun1,
        sigma_fun2 = sigma_fun2,
        y2_at = anchor2
      )
      if (
        nrow(pts) == 2 &&
          pts[2, 1] != pts[1, 1] &&
          !any(is_constant_boundary(pts[, 2]))
      ) {
        slope <- (pts[2, 2] - pts[1, 2]) / (pts[2, 1] - pts[1, 1])
        intercept <- pts[1, 2] - slope * pts[1, 1]
        boundary <<- local({
          m <- slope
          b <- intercept
          function(y2) b + m * y2
        })
        exact_line_built <<- TRUE
        return(invisible())
      }
      ## Two-anchor solve failed (e.g. a root did not exist at an anchor);
      ## fall through to the general path, which stays linear via the
      ## approxfun consumer. Latch so later calls skip the retry.
      linear_fallback <<- TRUE
    }

    ## Pad by max_gap: the interpolant slope at the outermost knot of a
    ## resolved interval depends on the next knot beyond it, which a later
    ## request may add far away. The padding keeps the affected
    ## sub-interval (at most max_gap wide) out of the requested range.
    todo <- merge_intervals(
      cbind(intervals[, 1] - max_gap, intervals[, 2] + max_gap)
    )
    changed <- FALSE
    for (i in seq_len(nrow(todo))) {
      ## A gap narrower than the refinement floor delta2 is bounded by
      ## knots which are closer than the floor already. Tracing it would
      ## only add knots so close that their secant is solver noise.
      pieces <- uncovered_intervals(
        todo[i, ],
        resolved,
        tol = min(delta2, max_gap)
      )
      for (j in seq_len(nrow(pieces))) {
        extra <- solve_boundary2S_normMix(
          decision,
          prior1,
          prior2,
          n1,
          n2,
          lim1,
          pieces[j, ],
          delta2,
          sigma_fun1 = sigma_fun1,
          sigma_fun2 = sigma_fun2,
          max_gap = max_gap,
          known = knots
        )
        knots <<- rbind(knots, extra)
        changed <- TRUE
      }
    }
    if (changed) {
      ## no knots at all if the decision could not be classified anywhere
      if (nrow(knots) > 1) {
        knots <<- knots[order(knots[, 1]), , drop = FALSE]
        knots <<- knots[c(TRUE, diff(knots[, 1]) > 1e-9), , drop = FALSE]
      }
      resolved <<- merge_intervals(rbind(resolved, todo))
      boundary <<- fit_boundary_consumer(knots, linear_boundary, decision)
    }
    invisible()
  }

  ## resolve decides where the boundary must be accurate, eval only
  ## reads the interpolant. Consumers knowing the distribution of y2
  ## (oc2S, pos2S) resolve the regions holding its mass and then
  ## integrate using eval. This split ensures that the function is
  ## accurate in regions where a high precision is needed and avoids
  ## additional calls to resolve in regions where the boundary is not
  ## needed with high precision due to low density of the integrand.
  cache <- list(resolve = resolve, eval = function(y2) boundary(y2))

  ## Without knowledge about the distribution of y2 the boundary is
  ## resolved at every requested point.
  decision_boundary <- function(y2, lim1) {
    resolve(cbind(y2, y2), lim1)
    boundary(y2)
  }

  decision_boundary
}

#' @templateVar fun decision2S_boundary
#' @template design2S-poisson
#' @export
decision2S_boundary.gammaMix <- function(
  prior1,
  prior2,
  n1,
  n2,
  decision,
  eps = 1e-6,
  ...
) {
  assert_that(likelihood(prior1) == "poisson")
  assert_that(likelihood(prior2) == "poisson")

  # only the second n2 argument may be 0
  assert_that(n1 > 0)
  assert_that(n2 >= 0)

  if (is(decision, "decision2S_2sided")) {
    decision2S_boundary_gammaMix_2sided(
      prior1,
      prior2,
      n1,
      n2,
      decision,
      eps,
      ...
    )
  } else {
    decision2S_boundary_gammaMix_1sided(
      prior1,
      prior2,
      n1,
      n2,
      decision,
      eps,
      ...
    )
  }
}

#' @keywords internal
decision2S_boundary_gammaMix_2sided <- function(
  prior1,
  prior2,
  n1,
  n2,
  decision,
  eps = 1e-6,
  ...
) {
  crit_lower <- decision2S_boundary_gammaMix_atomic(
    prior1,
    prior2,
    n1,
    n2,
    lower(decision),
    eps,
    ...
  )
  crit_upper <- decision2S_boundary_gammaMix_atomic(
    prior1,
    prior2,
    n1,
    n2,
    upper(decision),
    eps,
    ...
  )
  list(lower_or_equal_than = crit_lower, higher_than = crit_upper)
}

#' @keywords internal
decision2S_boundary_gammaMix_1sided <- function(
  prior1,
  prior2,
  n1,
  n2,
  decision,
  eps = 1e-6,
  ...
) {
  decision <- if (has_lower(decision)) {
    lower(decision)
  } else {
    upper(decision)
  }
  decision2S_boundary_gammaMix_atomic(
    prior1,
    prior2,
    n1,
    n2,
    decision,
    eps,
    ...
  )
}

#' @keywords internal
decision2S_boundary_gammaMix_atomic <- function(
  prior1,
  prior2,
  n1,
  n2,
  decision,
  eps = 1e-6,
  ...
) {
  if (!missing(eps)) {
    assert_number(eps, lower = 0, upper = 0.1, finite = TRUE)
  }

  cond_decisionStep <- function(post2) {
    function(m1) {
      decision(postmix(prior1, n = n1, m = m1 / n1), post2) - 0.25
    }
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
    if (
      lim1[1] < clim1[1] |
        lim1[2] > clim1[2] |
        lim2[1] < clim2[1] |
        lim2[2] > clim2[2]
    ) {
      ## ensure that lim1 never shrinks in size
      lim1[1] <- min(lim1[1], clim1[1])
      lim1[2] <- max(lim1[2], clim1[2])
      grid <<- lim2[1]:lim2[2]
      Neval <- length(grid)
      boundary <<- rep(NA, Neval)
      ## The step objective is monotone in m1 in the direction set by
      ## lower.tail; the critical y1 itself is monotone in y2 for the
      ## supported one-sided decisions.
      extend_dir <- if (lower.tail) "downX" else "upX"

      decFun_for <- function(y2) {
        if (n2 == 0) {
          cond_decisionStep(prior2)
        } else {
          cond_decisionStep(postmix(prior2, n = n2, m = y2 / n2))
        }
      }
      solve_raw <- function(y2, bracket) {
        uniroot_int(
          decFun_for(y2),
          bracket,
          extendInt = extend_dir,
          clamp = lim1
        )
      }
      ## map the constant-decision sentinels into the ordered domain
      ## [-1, Inf] so they join the monotone fill: never true
      ## (-Inf -> -1), always true (Inf, no finite upper count).
      map_root <- function(root) {
        if (is.finite(root)) {
          root
        } else if (root < 0) {
          -1
        } else {
          Inf
        }
      }
      set_val <- function(y2, val) {
        boundary[y2 - lim2[1] + 1] <<- val
      }

      if (getOption("RBesT.decision2S_boundary", "adaptive") == "adaptive") {
        ## Adaptive solve replacement. Because the critical value is
        ## monotone in y2, a run of outcomes bracketed by two equal
        ## critical values is provably constant and filled without any
        ## further (expensive pmixdiff) uniroot_int calls; otherwise the
        ## midpoint is solved, bracketed by the two neighbour roots, and
        ## the halves recurse.
        bracket_for <- function(ra, rb) {
          lo <- max(min(ra, rb), lim1[1])
          hi <- min(max(ra, rb), lim1[2])
          if (lo >= hi) lim1 else c(lo, hi)
        }
        fill <- function(a, ra, b, rb) {
          if (ra == rb) {
            for (yy in a:b) {
              set_val(yy, ra)
            }
            return(invisible())
          }
          if (b - a <= 1L) {
            set_val(a, ra)
            set_val(b, rb)
            return(invisible())
          }
          m <- (a + b) %/% 2L
          rm <- map_root(solve_raw(m, bracket_for(ra, rb)))
          fill(a, ra, m, rm)
          fill(m, rm, b, rb)
        }
        ra <- map_root(solve_raw(lim2[1], lim1))
        rb <- map_root(solve_raw(lim2[2], lim1))
        fill(lim2[1], ra, lim2[2], rb)
      } else {
        ## Legacy sequential sweep: warm-start each grid point's integer
        ## search with a tight bracket around the previous boundary.
        ## Note: the loop variable must not shadow the `y2` argument of
        ## decision_boundary, which is read again by the lookup below.
        prev <- NA
        for (yy in lim2[1]:lim2[2]) {
          bracket <- if (is.na(prev)) {
            lim1
          } else {
            c(max(prev - 2, lim1[1]), min(prev + 2, lim1[2]))
          }
          root <- solve_raw(yy, bracket)
          set_val(yy, map_root(root))
          prev <- if (is.finite(root)) root else NA
        }
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
