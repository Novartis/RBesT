# Tests for family-based sigma in 2S normal design functions

# =========================================================================
# Test Group 1: Backward compatibility
# =========================================================================

test_that("decision2S_boundary without family is unchanged", {
  prior1 <- mixnorm(c(1, 0, 1.5), sigma = 2)
  prior2 <- mixnorm(c(1, 0, 1.5), sigma = 2)
  dec <- decision2S(0.95, 0, lower.tail = FALSE)

  crit_fn <- decision2S_boundary(prior1, prior2, 30, 30, dec, sigma1 = 2, sigma2 = 2)
  # Should return a function
  expect_true(is.function(crit_fn))
  # Boundary over a range should give finite values
  y2_grid <- seq(-1, 1, by = 0.5)
  vals <- crit_fn(y2_grid)
  expect_true(all(is.finite(vals)))
})

test_that("oc2S without family is unchanged", {
  prior <- mixnorm(c(1, 0, 100), sigma = 2)
  dec <- decision2S(0.975, 0, lower.tail = FALSE)
  n1 <- 30
  n2 <- 30

  oc_fn <- oc2S(prior, prior, n1, n2, dec)
  # At equal thetas (type I error), should be < 0.05
  expect_true(oc_fn(0, 0) < 0.05)
  # With a true difference, power should be moderate
  expect_true(oc_fn(1, 0) > 0.1)
})

test_that("pos2S without family is unchanged", {
  prior <- mixnorm(c(1, 0, 100), sigma = 2)
  dec <- decision2S(0.975, 0, lower.tail = FALSE)
  n1 <- 30
  n2 <- 30

  pos_fn <- pos2S(prior, prior, n1, n2, dec)
  # PoS with the priors themselves
  result <- pos_fn(prior, prior)
  expect_true(result > 0 && result < 1)
})


# =========================================================================
# Test Group 2: Gaussian(identity) + sigma short-circuits
# =========================================================================

test_that("decision2S_boundary with gaussian(identity) matches fixed sigma", {
  prior1 <- mixnorm(c(1, 0, 1.5), sigma = 2)
  prior2 <- mixnorm(c(1, 0, 1.5), sigma = 2)
  dec <- decision2S(0.95, 0, lower.tail = FALSE)
  n1 <- 30
  n2 <- 30

  crit_fixed <- decision2S_boundary(prior1, prior2, n1, n2, dec, sigma1 = 2, sigma2 = 2)
  crit_family <- decision2S_boundary(prior1, prior2, n1, n2, dec,
                                      sigma1 = 2, sigma2 = 2,
                                      family = gaussian())

  # Evaluate at a few y2 points
  y2_pts <- seq(-2, 2, by = 1)
  expect_equal(crit_fixed(y2_pts), crit_family(y2_pts), tolerance = 1e-4)
})

test_that("oc2S with gaussian(identity) matches fixed sigma", {
  prior1 <- mixnorm(c(1, 0, 1.5), sigma = 2)
  prior2 <- mixnorm(c(1, 0, 1.5), sigma = 2)
  dec <- decision2S(0.95, 0, lower.tail = FALSE)
  n1 <- 30
  n2 <- 30

  oc_fixed <- oc2S(prior1, prior2, n1, n2, dec, sigma1 = 2, sigma2 = 2)
  oc_family <- oc2S(prior1, prior2, n1, n2, dec,
                     sigma1 = 2, sigma2 = 2, family = gaussian())

  expect_equal(oc_fixed(1, 0), oc_family(1, 0), tolerance = 1e-4)
  expect_equal(oc_fixed(0, 0), oc_family(0, 0), tolerance = 1e-4)
})


# =========================================================================
# Test Group 3: Binomial family boundary
# =========================================================================

test_that("decision2S_boundary with binomial() family differs from fixed sigma", {
  prior1 <- mixnorm(c(1, 0, 2), sigma = 2)
  prior2 <- mixnorm(c(1, 0, 2), sigma = 2)
  dec <- decision2S(0.95, qlogis(0.10), lower.tail = FALSE)
  n1 <- 40
  n2 <- 40

  y2_pts <- qlogis(c(0.05, 0.10, 0.15, 0.20))
  crit_family <- decision2S_boundary(prior1, prior2, n1, n2, dec,
                                      family = binomial())
  crit_fixed <- decision2S_boundary(prior1, prior2, n1, n2, dec,
                                     sigma1 = 2, sigma2 = 2)

  # Boundaries should differ at some points
  expect_false(isTRUE(all.equal(crit_family(y2_pts), crit_fixed(y2_pts),
                                tolerance = 1e-2)))
})

test_that("decision2S_boundary with binomial() returns a function", {
  prior1 <- mixnorm(c(1, 0, 2), sigma = 2)
  prior2 <- mixnorm(c(1, 0, 2), sigma = 2)
  dec <- decision2S(0.95, qlogis(0.10), lower.tail = FALSE)
  n1 <- 40
  n2 <- 40

  crit_fn <- decision2S_boundary(prior1, prior2, n1, n2, dec,
                                  family = binomial())
  expect_true(is.function(crit_fn))

  # Should give finite values across a range of y2
  y2_pts <- qlogis(c(0.05, 0.10, 0.15, 0.20))
  vals <- crit_fn(y2_pts)
  expect_true(all(is.finite(vals)))
})


# =========================================================================
# Test Group 4: OC with family
# =========================================================================

test_that("oc2S with binomial() family gives reasonable type I error", {
  prior1 <- mixnorm(c(1, 0, 2), sigma = 2)
  prior2 <- mixnorm(c(1, 0, 2), sigma = 2)
  dec <- decision2S(0.95, 0, lower.tail = FALSE)
  n1 <- 60
  n2 <- 60

  oc_fn <- oc2S(prior1, prior2, n1, n2, dec, family = binomial())

  # Under null (both arms at same rate, logit(0.30)):
  theta_null <- qlogis(0.30)
  type_I <- oc_fn(theta_null, theta_null)
  expect_true(type_I < 0.10)

  # Under alternative (arm 1 higher than arm 2):
  theta_alt1 <- qlogis(0.50)
  theta_alt2 <- qlogis(0.20)
  power <- oc_fn(theta_alt1, theta_alt2)
  expect_true(power > type_I)
})

test_that("oc2S with binomial() family differs from fixed sigma in tails", {
  prior1 <- mixnorm(c(1, 0, 2), sigma = 2)
  prior2 <- mixnorm(c(1, 0, 2), sigma = 2)
  dec <- decision2S(0.95, 0, lower.tail = FALSE)
  n1 <- 40
  n2 <- 40

  oc_family <- oc2S(prior1, prior2, n1, n2, dec, family = binomial())
  oc_fixed <- oc2S(prior1, prior2, n1, n2, dec, sigma1 = 2, sigma2 = 2)

  # At extreme rates where sigma differs most from 2
  theta1 <- qlogis(0.05)
  theta2 <- qlogis(0.05)
  expect_false(isTRUE(all.equal(oc_family(theta1, theta2),
                                oc_fixed(theta1, theta2),
                                tolerance = 1e-2)))
})

test_that("family boundary recovers a narrow binomial decision region", {
  prior1 <- mixnorm(c(1, 0, 1), sigma = 2)
  prior2 <- mixnorm(c(1, 0, 0.2), sigma = 2)
  posterior1 <- function(y) {
    p <- plogis(y)
    postmix(prior1, m = y, se = 1 / sqrt(20 * p * (1 - p)))
  }

  for (lower.tail in c(FALSE, TRUE)) {
    direction <- if (lower.tail) -1 else 1
    qc <- direction * 0.178
    dec <- decision2S(0.975, qc, lower.tail)
    objective <- function(y) {
      post1 <- posterior1(y)
      pnorm(
        qc,
        mean = post1[2, 1],
        sd = sqrt(post1[3, 1]^2 + prior2[3, 1]^2),
        lower.tail = lower.tail
      ) - 0.975
    }
    expected <- uniroot(
      objective, sort(direction * c(2, 2.1)), tol = 1e-10
    )$root
    for (mode in c("adaptive", "grid")) {
      withr::local_options(list(RBesT.decision2S_boundary = mode))
      boundary <- decision2S_boundary(
        prior1, prior2, 20, 0, dec, family = binomial()
      )
      critical <- boundary(0)
      expect_equal(critical, expected, tolerance = 1e-4)
      expect_equal(dec(posterior1(critical + direction * 0.01), prior2), 1)
      expect_equal(dec(posterior1(critical - direction * 0.01), prior2), 0)
    }
  }
})

test_that("family boundary retains arm-1 resolution with a diffuse arm-2 prior", {
  prior1 <- mixnorm(c(1, 0, 1), sigma = 2)

  for (lower.tail in c(FALSE, TRUE)) {
    direction <- if (lower.tail) -1 else 1
    prior2 <- mixnorm(c(1, -direction * 195, 100), sigma = 2)
    dec <- decision2S(0.975, 0, lower.tail)

    for (mode in c("adaptive", "grid")) {
      withr::local_options(list(RBesT.decision2S_boundary = mode))
      boundary <- decision2S_boundary(
        prior1, prior2, 20, 0, dec, family = binomial()
      )
      expect_equal(boundary(0), direction * 1.294373, tolerance = 1e-4)
    }
  }
})

test_that("family boundary preserves an always-success binomial decision", {
  prior1 <- mixnorm(c(1, 0.11, 0.17), sigma = 2)
  prior2 <- mixnorm(c(1, -0.69, 0.17), sigma = 2)
  dec <- decision2S(0.975, 0, lower.tail = FALSE)
  theta1 <- qlogis(c(0.3, 0.5, 0.6))
  theta2 <- rep(qlogis(0.3), 3)

  for (mode in c("adaptive", "grid")) {
    withr::local_options(list(RBesT.decision2S_boundary = mode))
    boundary <- decision2S_boundary(
      prior1,
      prior2,
      20,
      20,
      dec,
      family = binomial()
    )

    expect_equal(boundary(seq(-2, 1, 0.5)), rep(-1e100, 7))
    expect_equal(
      oc2S(prior1, prior2, 20, 20, dec, family = binomial())(theta1, theta2),
      rep(1, 3)
    )
  }
})

test_that("family boundary does not interpolate across a constant region", {
  prior1 <- mixnorm(c(1, 0.11, 0.25), sigma = 2)
  prior2 <- mixnorm(c(1, -0.69, 0.25), sigma = 2)
  dec <- decision2S(0.975, 0, lower.tail = FALSE)
  y2 <- seq(-1.5, 0, by = 0.05)

  ## reference for the frequency of success at a response rate of 0.3 in
  ## both arms: the decision is integrated directly over the outcomes
  sem <- function(y) 1 / sqrt(20 * plogis(y) * (1 - plogis(y)))
  z <- seq(-6, 6, length.out = 401)
  weight <- dnorm(z) / sum(dnorm(z))
  y <- qlogis(0.3) + sem(qlogis(0.3)) * z
  posterior <- function(prior) {
    sapply(y, function(m) postmix(prior, m = m, se = sem(m))[2:3, 1])
  }
  post1 <- posterior(prior1)
  post2 <- posterior(prior2)
  success <- outer(seq_along(y), seq_along(y), function(i, j) {
    pnorm(
      0,
      post1[1, i] - post2[1, j],
      sqrt(post1[2, i]^2 + post2[2, j]^2),
      lower.tail = FALSE
    ) > 0.975
  })
  oc_reference <- sum(outer(weight, weight) * success)

  result <- lapply(c("adaptive", "grid"), function(mode) {
    withr::local_options(list(RBesT.decision2S_boundary = mode))
    boundary <- decision2S_boundary(
      prior1,
      prior2,
      20,
      20,
      dec,
      family = binomial()
    )
    incremental_boundary <- decision2S_boundary(
      prior1,
      prior2,
      20,
      20,
      dec,
      family = binomial()
    )
    incremental_boundary(-1.176)
    list(
      boundary = boundary(y2),
      incremental_boundary = incremental_boundary(y2),
      oc = oc2S(
        prior1,
        prior2,
        20,
        20,
        dec,
        family = binomial()
      )(qlogis(0.3), qlogis(0.3))
    )
  })

  for (x in result) {
    constant <- which(x$boundary == -1e100)
    incremental_constant <- which(x$incremental_boundary == -1e100)
    expect_equal(constant, seq_len(max(constant)))
    expect_true(all(is.finite(x$boundary[-constant])))
    expect_equal(incremental_constant, constant)
    expect_equal(
      x$incremental_boundary[-constant],
      x$boundary[-constant],
      tolerance = 0.005
    )
    expect_equal(x$oc, oc_reference, tolerance = 0.005)
  }
  ## pieces without a constant region are traced adaptively, so the two
  ## schemes agree up to the interpolation accuracy only
  expect_equal(result[[1]], result[[2]], tolerance = 0.005)
})

test_that("family boundary does not depend on the truths of earlier calls", {
  dec <- decision2S(0.975, 0, lower.tail = FALSE)

  ## the critical values lie far from the outcomes of the first call
  prior_weak <- mixnorm(c(1, 0, 2), sigma = 1)
  prior_robust <- mixnorm(c(0.8, 0, 0.3), c(0.2, 0, 2), sigma = 1)
  design <- function() {
    oc2S(
      prior_weak, prior_robust, 20, 20, dec,
      family = poisson(), offset1 = log(5), offset2 = log(2)
    )
  }
  fresh <- design()(3, 0)
  used <- design()
  used(-3, 0)
  expect_equal(fresh, 1, tolerance = 1e-6)
  expect_equal(used(3, 0), fresh, tolerance = 1e-6)

  prior1 <- mixnorm(c(1, log(20), 1), sigma = 1)
  prior2 <- mixnorm(c(0.8, log(20), 0.2), c(0.2, log(20), 1), sigma = 1)
  theta1 <- log(c(20, 25, 30, 40))
  design <- function() oc2S(prior1, prior2, 20, 20, dec, family = poisson())
  fresh <- design()(theta1, log(20))
  used <- design()
  used(log(1), log(20))
  expect_true(all(diff(fresh) > 0))
  expect_equal(used(theta1, log(20)), fresh, tolerance = 1e-3)

  boundary <- decision2S_boundary(prior1, prior2, 20, 20, dec, family = poisson())
  expect_true(abs(boundary(log(20), lim1 = c(-0.5, 0.5))) < 1e100)
})

test_that("family boundary finds critical values far from the prior of sample 1", {
  offset <- log(10)
  posterior <- function(prior, y, n) {
    postmix(prior, m = y, se = exp(-(y + offset) / 2) / sqrt(n))
  }
  expect_separates <- function(prior1, post2, n1, dec, critical) {
    direction <- if (attr(dec, "lower.tail")) -1 else 1
    expect_lt(abs(critical), 1e100)
    expect_equal(dec(posterior(prior1, critical + direction * 0.01, n1), post2), 1)
    expect_equal(dec(posterior(prior1, critical - direction * 0.01, n1), post2), 0)
  }

  ## large counts are informative throughout; the critical value lies
  ## beyond the prior of sample 1 where the decision is monotone
  cases <- list(
    list(c(0.647, 0.121), c(2.804, 0.660), 10, decision2S(0.849, 0, TRUE)),
    list(c(-0.356, 0.342), c(2.080, 0.164), 50, decision2S(0.967, 0, FALSE))
  )
  for (mode in c("adaptive", "grid")) {
    withr::local_options(list(RBesT.decision2S_boundary = mode))
    for (case in cases) {
      prior1 <- mixnorm(c(1, case[[1]]), sigma = 1)
      prior2 <- mixnorm(c(1, case[[2]]), sigma = 1)
      boundary <- decision2S_boundary(
        prior1, prior2, case[[3]], 0, case[[4]],
        family = poisson(), offset1 = offset
      )
      expect_separates(prior1, prior2, case[[3]], case[[4]], boundary(0))
    }
  }

  ## components of the prior which are far apart
  prior1 <- mixnorm(c(0.7, log(5), 0.1), c(0.3, log(50), 0.1), sigma = 1)
  prior2 <- mixnorm(c(1, log(5), 0.1), sigma = 1)
  dec <- decision2S(0.975, 0, lower.tail = FALSE)
  boundary <- decision2S_boundary(
    prior1, prior2, 20, 20, dec,
    family = poisson(), offset1 = offset
  )
  for (y2 in c(3.3, 3.9)) {
    expect_separates(prior1, posterior(prior2, y2, 20), 20, dec, boundary(y2))
  }
})

test_that("oc2S warns if the family decision changes twice along sample 1", {
  prior1 <- mixnorm(c(1, 0.11, 0.25), sigma = 2)
  prior2 <- mixnorm(c(1, -0.69, 0.25), sigma = 2)
  dec <- decision2S(0.975, 0, lower.tail = FALSE)
  dec_2sided <- decision2S(c(0.975, 0.6), c(0, 2), c(FALSE, TRUE))
  pattern <- "changes more than once along the outcome of sample 1 for 1 of 2"

  ## with few expected events sample 1 has mass where the posterior
  ## reverts to the prior, which on its own fulfills the decision
  design <- oc2S(prior1, prior2, 20, 20, dec, family = binomial())
  expect_warning(design(qlogis(c(0.1, 0.5)), qlogis(c(0.3, 0.3))), pattern)
  expect_no_warning(design(qlogis(0.5), qlogis(0.3)))

  design_2sided <- oc2S(prior1, prior2, 20, 20, dec_2sided, family = binomial())
  expect_warning(
    design_2sided(qlogis(c(0.1, 0.5)), qlogis(c(0.3, 0.3))),
    pattern
  )

  design_n2_0 <- oc2S(prior1, prior2, 20, 0, dec, family = binomial())
  expect_warning(design_n2_0(qlogis(c(0.02, 0.5))), pattern)

  ## a truth on the boundary disagrees by the accuracy of the boundary only
  prior_weak <- mixnorm(c(1, 0, 2), sigma = 2)
  boundary <- decision2S_boundary(prior_weak, prior2, 40, 20, dec, family = binomial())
  design_weak <- oc2S(prior_weak, prior2, 40, 20, dec, family = binomial())
  for (shift in c(-1e-6, 1e-6)) {
    expect_no_warning(design_weak(boundary(0) + shift, 0))
  }

  ## the decision changes only once for a constant sigma
  design_fixed <- oc2S(prior1, prior2, 20, 20, dec, sigma1 = 2, sigma2 = 2)
  expect_no_warning(design_fixed(qlogis(c(0.1, 0.5)), qlogis(c(0.3, 0.3))))
})

test_that("constant family boundaries have the correct tail polarity", {
  cases <- list(
    binomial = list(family = binomial()),
    poisson = list(family = poisson()),
    negative_binomial = list(family = MASS::negative.binomial(1)),
    gaussian_log = list(
      family = gaussian("log"),
      sigma1 = 1,
      sigma2 = 1
    )
  )

  for (mode in c("adaptive", "grid")) {
    withr::local_options(list(RBesT.decision2S_boundary = mode))
    for (case_name in names(cases)) {
      case <- cases[[case_name]]
      tail_directions <- if (
        case_name %in% c("binomial", "negative_binomial")
      ) {
        c(FALSE, TRUE)
      } else {
        FALSE
      }
      for (lower.tail in tail_directions) {
        means <- if (lower.tail) c(-1.2, -0.4) else c(-0.4, -1.2)
        prior1 <- mixnorm(c(1, means[1], 0.03), sigma = 2)
        prior2 <- mixnorm(c(1, means[2], 0.03), sigma = 2)
        args <- c(
          list(
            prior1 = prior1,
            prior2 = prior2,
            n1 = 20,
            n2 = 20,
            decision = decision2S(0.975, 0, lower.tail),
            family = case$family
          ),
          case[setdiff(names(case), "family")]
        )
        boundary <- do.call(decision2S_boundary, args)
        expected <- if (lower.tail) 1e100 else -1e100

        expect_equal(
          boundary(seq(-1.1, -0.9, length.out = 3)),
          rep(expected, 3),
          info = paste(case$family$family, case$family$link, mode, lower.tail)
        )
      }
    }
  }
})


# =========================================================================
# Test Group 5: Offset
# =========================================================================

test_that("offset changes boundary for NB family in 2S case", {
  theta_NB <- 0.5
  fam_nb <- MASS::negative.binomial(theta = theta_NB)
  prior1 <- mixnorm(c(1, 0.5, 1), sigma = 1.7)
  prior2 <- mixnorm(c(1, 0.5, 1), sigma = 1.7)
  n1 <- 30
  n2 <- 30
  dec <- decision2S(0.95, 0, lower.tail = FALSE)

  crit_no_offset <- decision2S_boundary(prior1, prior2, n1, n2, dec,
                                         family = fam_nb)
  crit_with_offset <- decision2S_boundary(prior1, prior2, n1, n2, dec,
                                           family = fam_nb,
                                           offset1 = log(2), offset2 = log(2))

  # Boundaries should differ with different exposure
  y2_pts <- seq(0, 1, by = 0.25)
  expect_false(isTRUE(all.equal(crit_no_offset(y2_pts),
                                crit_with_offset(y2_pts),
                                tolerance = 1e-2)))
})

test_that("offset2 defaults to offset1", {
  theta_NB <- 0.5
  fam_nb <- MASS::negative.binomial(theta = theta_NB)
  prior1 <- mixnorm(c(1, 0.5, 1), sigma = 1.7)
  prior2 <- mixnorm(c(1, 0.5, 1), sigma = 1.7)
  n1 <- 30
  n2 <- 30
  dec <- decision2S(0.95, 0, lower.tail = FALSE)

  # offset1 = log(2), offset2 should default to log(2)
  crit_default <- decision2S_boundary(prior1, prior2, n1, n2, dec,
                                       family = fam_nb, offset1 = log(2))
  crit_explicit <- decision2S_boundary(prior1, prior2, n1, n2, dec,
                                        family = fam_nb,
                                        offset1 = log(2), offset2 = log(2))

  y2_pts <- seq(0, 1, by = 0.25)
  expect_equal(crit_default(y2_pts), crit_explicit(y2_pts), tolerance = 1e-6)
})


# =========================================================================
# Test Group 6: Error conditions
# =========================================================================

test_that("specifying sigma with non-gaussian family raises error in 2S", {
  prior1 <- mixnorm(c(1, 0, 2), sigma = 2)
  prior2 <- mixnorm(c(1, 0, 2), sigma = 2)
  dec <- decision2S(0.95, 0, lower.tail = FALSE)

  expect_error(
    decision2S_boundary(prior1, prior2, 30, 30, dec,
                        sigma1 = 2, family = binomial()),
    "sigma.*non-gaussian|Cannot specify"
  )
  expect_error(
    oc2S(prior1, prior2, 30, 30, dec,
         sigma1 = 2, family = binomial()),
    "sigma.*non-gaussian|Cannot specify"
  )
  expect_error(
    pos2S(prior1, prior2, 30, 30, dec,
          sigma1 = 2, family = binomial()),
    "sigma.*non-gaussian|Cannot specify"
  )
})

test_that("gaussian() family without sigma raises error in 2S", {
  prior1 <- mixnorm(c(1, 0, 2), sigma = 2)
  prior2 <- mixnorm(c(1, 0, 2), sigma = 2)
  dec <- decision2S(0.95, 0, lower.tail = FALSE)

  expect_error(
    decision2S_boundary(prior1, prior2, 30, 30, dec, family = gaussian()),
    "sigma.*must be specified|gaussian.*sigma"
  )
})

test_that("invalid family object raises error in 2S", {
  prior1 <- mixnorm(c(1, 0, 2), sigma = 2)
  prior2 <- mixnorm(c(1, 0, 2), sigma = 2)
  dec <- decision2S(0.95, 0, lower.tail = FALSE)

  expect_error(
    decision2S_boundary(prior1, prior2, 30, 30, dec, family = "binomial"),
    "family.*must be"
  )
})


# =========================================================================
# Test Group 7: PoS consistency
# =========================================================================

test_that("pos2S with family at tight priors matches oc2S", {
  prior1 <- mixnorm(c(1, 0, 2), sigma = 2)
  prior2 <- mixnorm(c(1, 0, 2), sigma = 2)
  dec <- decision2S(0.95, qlogis(0.10), lower.tail = FALSE)
  n1 <- 40
  n2 <- 40

  oc_fn <- oc2S(prior1, prior2, n1, n2, dec, family = binomial())
  pos_fn <- pos2S(prior1, prior2, n1, n2, dec, family = binomial())

  # Tight priors centered at specific values
  theta1 <- qlogis(0.20)
  theta2 <- qlogis(0.10)
  tight1 <- mixnorm(c(1, theta1, 0.001), sigma = 2)
  tight2 <- mixnorm(c(1, theta2, 0.001), sigma = 2)

  expect_equal(pos_fn(tight1, tight2), oc_fn(theta1, theta2), tolerance = 0.01)
})

test_that("pos2S with gaussian(identity) matches fixed-sigma pos2S", {
  prior1 <- mixnorm(c(1, 0, 1.5), sigma = 2)
  prior2 <- mixnorm(c(1, 0, 1.5), sigma = 2)
  dec <- decision2S(0.95, 0, lower.tail = FALSE)
  n1 <- 30
  n2 <- 30

  pos_fixed <- pos2S(prior1, prior2, n1, n2, dec, sigma1 = 2, sigma2 = 2)
  pos_family <- pos2S(prior1, prior2, n1, n2, dec,
                       sigma1 = 2, sigma2 = 2, family = gaussian())

  expect_equal(pos_fixed(prior1, prior2), pos_family(prior1, prior2),
               tolerance = 1e-3)
})
