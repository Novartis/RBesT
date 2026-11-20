test_that("decision2S_boundary works for normal outcome", {
  priorT <- mixnorm(c(1, 0, 0.001), sigma = 88, param = "mn")
  priorP <- mixnorm(c(1, -49, 20), sigma = 88, param = "mn")

  successCrit <- decision2S(c(0.95, 0.5), c(0, 50), FALSE)
  futilityCrit <- decision2S(c(0.90), c(40), TRUE)

  successBoundary <- decision2S_boundary(priorP, priorT, 10, 20, successCrit)
  futilityBoundary <- decision2S_boundary(priorP, priorT, 10, 20, futilityCrit)

  gridVals <- -25:25 - 49

  successBounds <- successBoundary(gridVals)
  futilityBounds <- futilityBoundary(gridVals)
  expect_snapshot_value(successBounds, style = "deparse", tolerance = 1e-4)
  expect_snapshot_value(futilityBounds, style = "deparse", tolerance = 1e-4)

  # Now we define the criterion for the gray zone using mixed lower.tail.
  grayzoneCrit <- decision2S(
    c(0.95, 0.5, 0.9),
    c(0, 50, 40),
    c(TRUE, TRUE, FALSE)
  )
  grayzoneBoundary <- decision2S_boundary(priorP, priorT, 10, 20, grayzoneCrit)
  grayzoneBoundsLower <- grayzoneBoundary$lower_or_equal_than(gridVals)
  grayzoneBoundsHigher <- grayzoneBoundary$higher_than(gridVals)

  expect_snapshot_value(grayzoneBoundsLower, style = "deparse", tolerance = 1e-4)
  expect_snapshot_value(grayzoneBoundsHigher, style = "deparse", tolerance = 1e-4)

  # In this case there is no gray zone:
  expect_false(any(grayzoneBoundsHigher < grayzoneBoundsLower))
})

test_that("Adaptive normal boundary matches the uniform grid fallback", {
  skip_on_cran()

  ## Multi-component priors give a genuinely curved (non-linear) boundary,
  ## which is where the adaptive tracer and the grid sweep could diverge.
  prior1 <- mixnorm(rob = c(0.2, 0, 2), inf = c(0.8, 2, 2), sigma = 5)
  prior2 <- mixnorm(rob = c(0.2, 1, 2), inf = c(0.8, 3, 2), sigma = 5)

  dec <- decision2S(pc = 0.5, qc = 1.5, lower.tail = TRUE)
  gridVals <- -25:25 - 49

  withr::with_options(
    list(RBesT.decision2S_boundary = "adaptive"),
    {
      boundary_adaptive <- decision2S_boundary(prior1, prior2, 50, 50, dec)
      result_adaptive <- boundary_adaptive(gridVals)
    }
  )
  withr::with_options(
    list(RBesT.decision2S_boundary = "grid"),
    {
      boundary_grid <- decision2S_boundary(prior1, prior2, 50, 50, dec)
      result_grid <- boundary_grid(gridVals)
    }
  )

  expect_equal(result_adaptive, result_grid, tolerance = 1e-3)
})

## Mixture priors: the boundary follows one slope where the informative
## component dominates and another one in the tails. On a broad, symmetric
## range the samples of a coarse trace all fall on a common line, which the
## adaptive tracer used to accept.
boundary2S_curved_scenario <- function() {
  list(
    prior1 = mixnorm(c(1, 0, 100), sigma = 1),
    prior2 = mixnorm(c(0.8, 0, 1), c(0.2, 0, 2), sigma = 1),
    decisions = list(
      median = decision2S(0.5, 2.5, FALSE),
      dual = decision2S(c(0.95, 0.5), c(0, 2.5), FALSE),
      dual_lower = decision2S(c(0.95, 0.5), c(0, -2.5), TRUE)
    )
  )
}

test_that("Adaptive normal boundary separates the posterior decisions on a broad range", {
  s <- boundary2S_curved_scenario()
  y2 <- c(-20, -3, 20)
  for (decision in s$decisions) {
    boundary <- decision2S_boundary(
      s$prior1, s$prior2, 1, 1, decision, sigma1 = 1, sigma2 = 1
    )
    crit <- boundary(y2)
    lower.tail <- attr(decision, "lower.tail")
    for (i in seq_along(y2)) {
      post2 <- postmix(s$prior2, m = y2[i], se = 1)
      below <- decision(postmix(s$prior1, m = crit[i] - 0.02, se = 1), post2)
      above <- decision(postmix(s$prior1, m = crit[i] + 0.02, se = 1), post2)
      expect_equal(below, as.numeric(lower.tail))
      expect_equal(above, as.numeric(!lower.tail))
    }
  }
})

test_that("Adaptive normal boundary matches directly solved roots on a broad range", {
  s <- boundary2S_curved_scenario()
  y2 <- seq(-40, 40, length.out = 81)
  for (decision in s$decisions) {
    boundary <- decision2S_boundary(
      s$prior1, s$prior2, 1, 1, decision, sigma1 = 1, sigma2 = 1
    )
    direct <- vapply(
      y2,
      function(y) {
        post2 <- postmix(s$prior2, m = y, se = 1)
        uniroot(
          ## the decision is a step in y1; bisect to the jump
          function(y1) {
            2 * decision(postmix(s$prior1, m = y1, se = 1), post2) - 1
          },
          c(-100, 100),
          tol = 1e-8
        )$root
      },
      numeric(1)
    )
    ## accepted quarter-point residual is sem1 / 100 = 0.01
    expect_lt(max(abs(boundary(y2) - direct)), 0.02)
  }
})

test_that("Normal boundary in a resolved range is not altered by later requests", {
  s <- boundary2S_curved_scenario()
  boundary <- decision2S_boundary(
    s$prior1, s$prior2, 1, 1, s$decisions$dual, sigma1 = 1, sigma2 = 1
  )
  y2 <- seq(-5, 5, by = 0.25)
  first <- boundary(y2)
  boundary(c(-300, 300))
  boundary(seq(20, 60, by = 5))
  expect_equal(boundary(y2), first, tolerance = 1e-10)
})

test_that("Normal boundary is not degraded by a tiny gap between resolved ranges", {
  s <- boundary2S_curved_scenario()
  decision <- s$decisions$dual
  y2 <- seq(0, 2, by = 0.05)
  fresh <- decision2S_boundary(
    s$prior1, s$prior2, 1, 1, decision, sigma1 = 1, sigma2 = 1
  )(y2)
  for (gap in c(3e-9, 1e-6, 1e-4)) {
    boundary <- decision2S_boundary(
      s$prior1, s$prior2, 1, 1, decision, sigma1 = 1, sigma2 = 1
    )
    ## the padded ranges of these two requests miss each other by gap
    boundary(0)
    boundary(2 + gap)
    expect_lt(max(abs(boundary(y2) - fresh)), 0.005)
  }
})

test_that("Normal boundary without any critical value in the requested range is constant", {
  ## family path at low rates: the decision does not change with y1
  prior <- mixnorm(c(0.7, log(5), 0.3), c(0.3, log(5), 1.5), sigma = 1)
  decision <- decision2S(0.95, 0, lower.tail = TRUE)
  boundary <- decision2S_boundary(
    prior, prior, 20, 20, decision, family = poisson()
  )
  expect_equal(boundary(-3), -1e100)
  design <- oc2S(prior, prior, 20, 20, decision, family = poisson())
  expect_true(is.finite(design(log(0.1) - 0.3, log(0.1))))
})

test_that("Constant-sigma normal boundary is exactly linear (single and multi-criteria)", {
  skip_on_cran()

  ## Single-component priors + constant sigma => y1_c(y2) is provably a
  ## straight line. For multi-criteria (min-combined) decisions the
  ## per-criterion lines are parallel, so the envelope is still a single
  ## line with no kink. The exact 2-point line path must reproduce this to
  ## machine precision: second differences of an evenly-spaced evaluation
  ## must vanish.
  priorT <- mixnorm(c(1, 0, 0.001), sigma = 88, param = "mn")
  priorP <- mixnorm(c(1, -49, 20), sigma = 88, param = "mn")
  gridVals <- seq(-80, -20, length.out = 41)

  single <- decision2S(0.95, 0, FALSE)
  multi <- decision2S(c(0.95, 0.5, 0.9), c(0, 50, 40), c(FALSE, FALSE, FALSE))

  b_single <- decision2S_boundary(priorP, priorT, 10, 20, single)(gridVals)
  b_multi <- decision2S_boundary(priorP, priorT, 10, 20, multi)(gridVals)

  ## second difference == 0 for a straight line
  expect_lt(max(abs(diff(diff(b_single)))), 1e-6)
  expect_lt(max(abs(diff(diff(b_multi)))), 1e-6)
})

test_that("Normal boundary with n2 == 0 is exactly constant in y2", {
  skip_on_cran()

  priorT <- mixnorm(c(1, 0, 0.001), sigma = 88, param = "mn")
  priorP <- mixnorm(c(1, -49, 20), sigma = 88, param = "mn")
  dec <- decision2S(0.95, 0, FALSE)

  boundary <- decision2S_boundary(priorP, priorT, 10, 0, dec)
  vals <- boundary(seq(-80, -20, length.out = 41))

  expect_lt(max(abs(vals - vals[1])), 1e-6)
})

test_that("Mixture normal boundary consumer (monoH.FC) does not overshoot", {
  skip_on_cran()

  ## The non-linear consumer must be the shape-preserving monotone Hermite:
  ## on a sharp monotone step a global cubic ("fmm", the old default)
  ## overshoots outside [min, max] of the anchor points, while monoH.FC
  ## stays within them. Exercise fit_boundary_consumer directly so the guard
  ## actually discriminates between the two interpolants (a monotone mixture
  ## boundary alone cannot, since its min/max are its endpoints).
  x <- 1:6
  y <- c(0, 0, 0, 1, 1, 1)
  dec <- decision2S(0.95, 0, FALSE)
  consumer <- fit_boundary_consumer(cbind(x, y), linear = FALSE, decision = dec)

  xs <- seq(1, 6, length.out = 201)
  ys <- consumer(xs)
  expect_lte(max(ys), 1 + 1e-9)
  expect_gte(min(ys), 0 - 1e-9)

  ## discrimination check: the old global cubic really would overshoot here
  fmm <- splinefun(x, y, method = "fmm")
  expect_gt(max(fmm(xs)), 1 + 1e-3)
})

test_that("Normal boundary consumer keeps constant runs out of the spline", {
  dec <- decision2S(0.95, 0, FALSE)
  points <- cbind(
    1:6,
    c(-boundary_extreme_value, -boundary_extreme_value, 0, 1,
      boundary_extreme_value, boundary_extreme_value)
  )
  consumer <- fit_boundary_consumer(points, linear = FALSE, decision = dec)

  expect_equal(
    consumer(c(1, 2.4, 4.6, 6)),
    c(-1e100, -1e100, 1e100, 1e100)
  )
  expect_equal(consumer(c(3, 4)), c(0, 1))
  expect_true(all(abs(consumer(seq(2.51, 4.49, length.out = 21))) < 1e100))
})

test_that("Mixed lower.tail usage works for normal decision boundary calculation", {
  skip_on_cran()

  prior1 <- mixnorm(rob = c(0.2, 0, 2), inf = c(0.8, 2, 2), sigma = 5)
  prior2 <- mixnorm(rob = c(0.2, 1, 2), inf = c(0.8, 3, 2), sigma = 5)

  dec_lower <- decision2S(pc = 0.5, qc = 1.5, lower.tail = TRUE)
  boundary_fn_lower <- decision2S_boundary(
    prior1,
    prior2,
    n1 = 50,
    n2 = 50,
    decision = dec_lower
  )

  gridVals <- -25:25 - 49
  result_lower <- boundary_fn_lower(gridVals)

  dec_upper <- decision2S(pc = 0.6, qc = 0.5, lower.tail = FALSE)
  boundary_fn_upper <- decision2S_boundary(
    prior1,
    prior2,
    n1 = 50,
    n2 = 50,
    decision = dec_upper
  )
  result_upper <- boundary_fn_upper(gridVals)

  decMixed <- decision2S(
    qc = c(1.5, 0.5),
    pc = c(0.5, 0.6),
    lower.tail = c(TRUE, FALSE)
  )
  boundary_fn_mixed <- decision2S_boundary(prior1, prior2, 50, 50, decMixed)

  result_mixed_lower <- boundary_fn_mixed$lower_or_equal_than(gridVals)
  result_mixed_upper <- boundary_fn_mixed$higher_than(gridVals)

  expect_equal(result_mixed_lower, result_lower)
  expect_equal(result_mixed_upper, result_upper)
})

test_that("decision2S_boundary works for binomial outcome", {
  priorT <- mixbeta(c(1, 1, 2))
  priorP <- mixbeta(c(1, 10, 40))

  successCrit <- decision2S(c(0.95, 0.5), c(0, 0.2), FALSE)
  futilityCrit <- decision2S(0.90, 0.2, TRUE)

  successBoundary <- decision2S_boundary(priorP, priorT, 20, 20, successCrit)
  futilityBoundary <- decision2S_boundary(priorP, priorT, 20, 20, futilityCrit)

  gridVals <- seq(0, 20)

  successBounds <- successBoundary(gridVals)
  futilityBounds <- futilityBoundary(gridVals)
  expect_snapshot_value(successBounds, style = "deparse")
  expect_snapshot_value(futilityBounds, style = "deparse")
})

test_that("Adaptive binomial boundary matches the uniform grid fallback", {
  skip_on_cran()

  priorT <- mixbeta(c(0.5, 1, 2), c(0.5, 5, 5))
  priorP <- mixbeta(c(0.5, 10, 40), c(0.5, 3, 6))

  successCrit <- decision2S(c(0.95, 0.5), c(0, 0.2), FALSE)
  futilityCrit <- decision2S(0.90, 0.2, TRUE)
  gridVals <- seq(0, 30)

  for (crit in list(successCrit, futilityCrit)) {
    withr::with_options(
      list(RBesT.decision2S_boundary = "adaptive"),
      {
        b_adaptive <- decision2S_boundary(priorP, priorT, 30, 30, crit)
        r_adaptive <- b_adaptive(gridVals)
      }
    )
    withr::with_options(
      list(RBesT.decision2S_boundary = "grid"),
      {
        b_grid <- decision2S_boundary(priorP, priorT, 30, 30, crit)
        r_grid <- b_grid(gridVals)
      }
    )
    ## discrete critical values are integers -- the monotone fill is exact,
    ## so adaptive and grid must agree identically (not merely to tolerance).
    expect_identical(r_adaptive, r_grid)
  }
})

test_that("Mixed lower.tail usage works for binomial decision boundary calculation", {
  skip_on_cran()

  priorT <- mixbeta(c(1, 1, 2))
  priorP <- mixbeta(c(1, 10, 40))

  dec_lower <- decision2S(pc = 0.5, qc = 0.7, lower.tail = TRUE)
  boundary_fn_lower <- decision2S_boundary(
    priorP,
    priorT,
    n1 = 20,
    n2 = 20,
    decision = dec_lower
  )

  gridVals <- seq(0, 20)
  result_lower <- boundary_fn_lower(gridVals)

  dec_upper <- decision2S(pc = 0.6, qc = 0.5, lower.tail = FALSE)
  boundary_fn_upper <- decision2S_boundary(
    priorP,
    priorT,
    n1 = 20,
    n2 = 20,
    decision = dec_upper
  )
  result_upper <- boundary_fn_upper(gridVals)

  decMixed <- decision2S(
    qc = c(0.7, 0.5),
    pc = c(0.5, 0.6),
    lower.tail = c(TRUE, FALSE)
  )
  boundary_fn_mixed <- decision2S_boundary(priorP, priorT, 20, 20, decMixed)

  result_mixed_lower <- boundary_fn_mixed$lower_or_equal_than(gridVals)
  result_mixed_upper <- boundary_fn_mixed$higher_than(gridVals)

  expect_equal(result_mixed_lower, result_lower)
  expect_equal(result_mixed_upper, result_upper)
})

test_that("decision2S_boundary works for Poisson outcome", {
  priorT <- mixgamma(c(1, 0.5, 2))
  priorP <- mixgamma(c(1, 1, 2))

  successCrit <- decision2S(c(0.95, 0.5), c(0, 1), FALSE)
  futilityCrit <- decision2S(0.90, 1, TRUE)

  successBoundary <- decision2S_boundary(priorP, priorT, 20, 20, successCrit)
  futilityBoundary <- decision2S_boundary(priorP, priorT, 20, 20, futilityCrit)

  gridVals <- seq(0, 20)

  successBounds <- successBoundary(gridVals)
  futilityBounds <- futilityBoundary(gridVals)
  expect_snapshot_value(successBounds, style = "deparse")
  expect_snapshot_value(futilityBounds, style = "deparse")
})

test_that("Adaptive Poisson boundary matches the uniform grid fallback", {
  skip_on_cran()

  priorT <- mixgamma(c(0.6, 0.5, 2), c(0.4, 2, 3))
  priorP <- mixgamma(c(0.6, 1, 2), c(0.4, 3, 4))

  successCrit <- decision2S(c(0.95, 0.5), c(0, 1), FALSE)
  futilityCrit <- decision2S(0.90, 1, TRUE)
  gridVals <- seq(0, 30)

  for (crit in list(successCrit, futilityCrit)) {
    withr::with_options(
      list(RBesT.decision2S_boundary = "adaptive"),
      {
        b_adaptive <- decision2S_boundary(priorP, priorT, 30, 30, crit)
        r_adaptive <- b_adaptive(gridVals)
      }
    )
    withr::with_options(
      list(RBesT.decision2S_boundary = "grid"),
      {
        b_grid <- decision2S_boundary(priorP, priorT, 30, 30, crit)
        r_grid <- b_grid(gridVals)
      }
    )
    ## discrete critical values are integers -- the monotone fill is exact,
    ## so adaptive and grid must agree identically (not merely to tolerance).
    expect_identical(r_adaptive, r_grid)
  }
})

test_that("Mixed lower.tail usage works for Poisson decision boundary calculation", {
  skip_on_cran()

  priorT <- mixgamma(c(1, 0.5, 2))
  priorP <- mixgamma(c(1, 1, 2))

  dec_lower <- decision2S(pc = 0.5, qc = 0.7, lower.tail = TRUE)
  boundary_fn_lower <- decision2S_boundary(
    priorP,
    priorT,
    n1 = 20,
    n2 = 20,
    decision = dec_lower
  )

  gridVals <- seq(0, 20)
  result_lower <- boundary_fn_lower(gridVals)

  dec_upper <- decision2S(pc = 0.6, qc = 0.5, lower.tail = FALSE)
  boundary_fn_upper <- decision2S_boundary(
    priorP,
    priorT,
    n1 = 20,
    n2 = 20,
    decision = dec_upper
  )
  result_upper <- boundary_fn_upper(gridVals)

  decMixed <- decision2S(
    qc = c(0.7, 0.5),
    pc = c(0.5, 0.6),
    lower.tail = c(TRUE, FALSE)
  )
  boundary_fn_mixed <- decision2S_boundary(priorP, priorT, 20, 20, decMixed)

  result_mixed_lower <- boundary_fn_mixed$lower_or_equal_than(gridVals)
  result_mixed_upper <- boundary_fn_mixed$higher_than(gridVals)

  expect_equal(result_mixed_lower, result_lower)
  expect_equal(result_mixed_upper, result_upper)
})
