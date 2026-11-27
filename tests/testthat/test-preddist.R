## check that predictive distributions hold what they promise,
## i.e. that they describe the sum of n new data points.

preddist_tolerance <- function() 1e-2

preddist_quantiles <- function() seq(0.1, 0.9, by = 0.1)

preddist_sample_size <- function() 1e5

preddist_trial_count <- function() 25

preddist_mixes <- function() {
  list(
    beta = mixbeta(c(1, 11, 4)),
    betaMix = mixbeta(c(0.8, 11, 4), c(0.2, 1, 1)),
    gamma = mixgamma(c(1, 65, 50)),
    gammaMix = mixgamma(
      rob = c(0.5, 8, 0.5),
      inf = c(0.5, 9, 2),
      param = "ms"
    ),
    norm = mixnorm(c(1, 0, 0.5), sigma = 1),
    normMix = mixnorm(c(0.2, 0, 2), c(0.8, 2, 1), sigma = 1)
  )
}

preddist_cmp <- function(
  mix,
  n,
  n_rng,
  N = preddist_sample_size(),
  qntls = preddist_quantiles(),
  stat = c("sum", "mean"),
  Teps = preddist_tolerance()
) {
  skip_on_cran()
  withr::local_seed(42343)

  ## sample for each draw a single hyper-parameter which is then
  ## used n times in the rng function to return n samples from the
  ## sampling distribution
  stat <- match.arg(stat)
  test <- replicate(N, n_rng(rmix(mix, 1)))
  if (stat == "sum") {
    test_stat <- colSums(test)
  }
  if (stat == "mean") {
    test_stat <- colMeans(test)
  }

  quants_stest <- quantile(test_stat, qntls)
  ## note: in particular for the discrete/counting distributions,
  ## sampling gives better estimates
  quants_sref <- qmix(preddist(mix, n = n), qntls)
  res_sum <- abs(quants_sref - quants_stest)
  ## note: errors are accumulating with n, hence to check the mean,
  ## we scale eps with n
  if (stat == "sum") {
    expect_true(all(res_sum < n * Teps))
  }
  if (stat == "mean") {
    expect_true(all(res_sum < Teps))
  }

  quants_test <- quantile(test[1, ], qntls)
  quants_ref <- qmix(preddist(mix, n = 1), qntls)
  res <- abs(quants_ref - quants_test)
  expect_true(all(res < Teps))

  if (inherits(mix, "betaMix")) {
    ## specifically test BetaBinomial (which is from Stan)
    predmix <- preddist(mix, n = 30)
    dens <- dmix(predmix, 0:30)
    cdens <- cumsum(dens)
    pc <- pmix(predmix, 0:30)
    upc <- pmix(predmix, 0:30, FALSE)
    expect_true(all(abs(cdens - pc) < Teps))
    expect_true(all(abs(1 - cdens - upc) < Teps))
  }

  ## test that output length of the predictive is the same as the
  ## input data vector, even for negative values for values outside
  ## the valid support of the postive only distribtions
  cprob_ltest <- pmix(preddist(mix, n = 1), c(-1, quants_ref))
  expect_true(length(cprob_ltest) == length(c(-1, quants_ref)))
}


test_that("Predictive for a beta evaluates correctly (binary)", {
  n <- preddist_trial_count()
  mix <- preddist_mixes()$beta
  preddist_cmp(mix, n, Curry(rbinom, n = n, size = 1))
})
test_that("Predictive for a beta mixture evaluates correctly (binary)", {
  n <- preddist_trial_count()
  mix <- preddist_mixes()$betaMix
  preddist_cmp(mix, n, Curry(rbinom, n = n, size = 1))
})

test_that("Predictive for a gamma evaluates correctly (poisson)", {
  n <- preddist_trial_count()
  mix <- preddist_mixes()$gamma
  preddist_cmp(mix, n, Curry(rpois, n = n))
})
test_that("Predictive for a gamma mixture evaluates correctly (poisson)", {
  n <- preddist_trial_count()
  mix <- preddist_mixes()$gammaMix
  preddist_cmp(mix, n, Curry(rpois, n = n), Teps = 1E-1)
})

test_that("Predictive for a normal evaluates correctly (normal)", {
  n <- preddist_trial_count()
  mix <- preddist_mixes()$norm
  preddist_cmp(mix, n, Curry(rnorm, n = n, sd = sigma(mix)), stat = "mean")
})
test_that("Predictive for a normal mixture evaluates correctly (normal)", {
  n <- preddist_trial_count()
  mix <- preddist_mixes()$normMix
  preddist_cmp(
    mix,
    n,
    Curry(rnorm, n = n, sd = sigma(mix)),
    stat = "mean",
    Teps = 1E-1
  )
})

test_that("qmix.betaBinomialMix handles lower.tail, log.p, and exact CDF jumps", {
  probs <- c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)
  pd_single <- preddist(mixbeta(c(1, 10, 90)), n = 100)
  pd_mix <- preddist(mixbeta(c(0.7, 10, 90), c(0.3, 2, 8)), n = 100)

  for (pd in list(pd_single, pd_mix)) {
    q_lo <- qmix(pd, 1 - probs, lower.tail = TRUE)
    q_hi <- qmix(pd, probs, lower.tail = FALSE)
    expect_equal(q_hi, q_lo)
    expect_equal(qmix(pd, log(probs), lower.tail = FALSE, log.p = TRUE), q_hi)
    expect_equal(qmix(pd, log(probs), lower.tail = TRUE, log.p = TRUE), qmix(pd, probs, lower.tail = TRUE))
    expect_equal(qmix(pd, c(0, 1), lower.tail = TRUE), c(0, 100))
    expect_equal(qmix(pd, c(0, 1), lower.tail = FALSE), c(100, 0))
  }

  pd_small <- preddist(mixbeta(c(0.6, 2, 5), c(0.4, 4, 2)), n = 10)
  ks <- 0:10
  expect_equal(qmix(pd_small, pmix(pd_small, ks, lower.tail = TRUE), lower.tail = TRUE), ks)
  expect_equal(qmix(pd_small, pmix(pd_small, ks, lower.tail = FALSE), lower.tail = FALSE), ks)
})

test_that("qmix.gammaPoissonMix handles lower.tail, log.p, and exact CDF jumps", {
  probs <- c(0.001, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 0.999)
  gm <- mixgamma(c(1, 20, 10), likelihood = "poisson")
  pd_single <- preddist(gm, n = 10)

  expect_equal(
    qmix(pd_single, probs, lower.tail = TRUE),
    qnbinom(probs, size = 20, prob = 10 / 20, lower.tail = TRUE)
  )
  expect_equal(
    qmix(pd_single, probs, lower.tail = FALSE),
    qnbinom(probs, size = 20, prob = 10 / 20, lower.tail = FALSE)
  )
  expect_equal(
    qmix(pd_single, log(probs), lower.tail = FALSE, log.p = TRUE),
    qnbinom(log(probs), size = 20, prob = 10 / 20, lower.tail = FALSE, log.p = TRUE)
  )
  expect_equal(qmix(pd_single, c(0, 1), lower.tail = TRUE), c(0, Inf))
  expect_equal(qmix(pd_single, c(0, 1), lower.tail = FALSE), c(Inf, 0))

  gm_mix <- mixgamma(c(0.6, 15, 10), c(0.4, 30, 5), likelihood = "poisson")
  pd_mix <- preddist(gm_mix, n = 10)
  q_lo <- qmix(pd_mix, 1 - probs, lower.tail = TRUE)
  q_hi <- qmix(pd_mix, probs, lower.tail = FALSE)
  expect_equal(q_hi, q_lo)
  expect_equal(qmix(pd_mix, log(probs), lower.tail = FALSE, log.p = TRUE), q_hi)
  expect_true(all(pmix(pd_mix, q_hi, lower.tail = FALSE) <= probs))
  expect_true(all(pmix(pd_mix, q_hi - 1, lower.tail = FALSE) > probs))
  expect_equal(qmix(pd_mix, c(0, 1), lower.tail = TRUE), c(0, Inf))
  expect_equal(qmix(pd_mix, c(0, 1), lower.tail = FALSE), c(Inf, 0))

  ks <- 0:20
  expect_equal(qmix(pd_mix, pmix(pd_mix, ks, lower.tail = TRUE), lower.tail = TRUE), ks)
  expect_equal(qmix(pd_mix, pmix(pd_mix, ks, lower.tail = FALSE), lower.tail = FALSE), ks)
})
