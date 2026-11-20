## test the analytical OC function via brute force simulation

oc2S_simulation_size <- function() 1e4

oc2S_normal_scenario <- function() {
  prior1 <- mixnorm(c(0.3, -0.2, 2), c(0.7, 0, 50), sigma = 1)
  prior2 <- mixnorm(c(1.0, 0, 50), sigma = 1)

  ## type I error fairly large to 20% to make it easier to test (less
  ## simulations needed for accurate results).
  pcrit <- 0.80
  qcrit <- 0

  list(
    n1 = 10,
    n2 = 20,
    pcrit = pcrit,
    prior1 = prior1,
    prior2 = prior2,
    qcrit = qcrit,
    theta1 = 0,
    ## theta2 set such that we have about 75% power under this truth.
    theta2 = 0.5
  )
}

run_on_cran <- function() {
  if (identical(Sys.getenv("NOT_CRAN"), "true")) {
    return(FALSE)
  }
  return(TRUE)
}

oc2S_normal_MC <- function(
  prior1,
  prior2,
  N1,
  N2,
  theta1,
  theta2,
  pcrit = 0.975,
  qcrit = 0
) {
  withr::local_seed(491087)

  mean_sd1 <- sigma(prior1) / sqrt(N1)
  mean_sd2 <- sigma(prior2) / sqrt(N2)

  mean_prior1 <- prior1
  sigma(mean_prior1) <- mean_sd1
  mean_prior2 <- prior2
  sigma(mean_prior2) <- mean_sd2

  mean_samp1 <- rnorm(oc2S_simulation_size(), theta1, mean_sd1)
  mean_samp2 <- rnorm(oc2S_simulation_size(), theta2, mean_sd2)

  dec <- rep(NA, oc2S_simulation_size())

  for (i in seq_len(oc2S_simulation_size())) {
    post1 <- postmix(mean_prior1, m = mean_samp1[i], se = mean_sd1)
    post2 <- postmix(mean_prior2, m = mean_samp2[i], se = mean_sd2)
    dec[i] <- as.numeric(pmix(RBesT:::mixnormdiff(post1, post2), qcrit) > pcrit)
  }

  mean(dec)
}

Voc2S_normal_MC <- function(...) {
  Vectorize(oc2S_normal_MC, c("theta1", "theta2"))(...)
}

## first test that the analytic difference distribution for normal
## works as expected

test_that("Analytical convolution of normal mixture matches numerical integration result", {
  skip_on_cran()

  sc <- oc2S_normal_scenario()
  pdiff <- RBesT:::mixnormdiff(sc$prior1, sc$prior2)
  x <- seq(-20, 20, length = 21)
  d1 <- dmix(pdiff, x)
  d2 <- dmixdiff(sc$prior1, sc$prior2, x)
  dres <- abs(d1 - d2)
  expect_equal(sum(dres > 1e-5), 0)
  p1 <- pmix(pdiff, x)
  p2 <- pmixdiff(sc$prior1, sc$prior2, x)
  pres <- 100 * abs(p1 - p2)
  expect_equal(sum(pres > 2), 0)
})

## test that the type I error is matching, i.e. is not off by more than 2%
test_that("Type I error is matching between MC and analytical computations in the normal mixture case", {
  skip_on_cran()

  sc <- oc2S_normal_scenario()
  x <- c(-2, 0)
  alpha <- oc2S(
    sc$prior1,
    sc$prior2,
    sc$n1,
    sc$n2,
    decision2S(sc$pcrit, sc$qcrit),
    sigma1 = sigma(sc$prior1),
    sigma2 = sigma(sc$prior2)
  )(x, x)
  alphaMC <- Voc2S_normal_MC(
    sc$prior1,
    sc$prior2,
    sc$n1,
    sc$n2,
    x,
    x,
    sc$pcrit,
    sc$qcrit
  )
  res <- 100 * abs(alpha - alphaMC)
  expect_equal(sum(res > 2), 0)
})


## test that the power is matching, i.e. is not off by more than 2%
test_that("Power is matching between MC and analytical computations in the normal mixture case", {
  skip_on_cran()

  sc <- oc2S_normal_scenario()
  power <- oc2S(
    sc$prior1,
    sc$prior2,
    sc$n1,
    sc$n2,
    decision2S(sc$pcrit, sc$qcrit),
    sigma1 = sigma(sc$prior1),
    sigma2 = sigma(sc$prior2)
  )(sc$theta1, sc$theta2)
  powerMC <- oc2S_normal_MC(
    sc$prior1,
    sc$prior2,
    sc$n1,
    sc$n2,
    sc$theta1,
    sc$theta2,
    sc$pcrit,
    sc$qcrit
  )
  res <- 100 * abs(power - powerMC)
  expect_equal(sum(res > 2), 0)
})

## further test by cross-checking with Gsponer et. al, "A practical
## guide to Bayesian group sequential designs", Pharmaceut. Statist.
## (2014), 13 71-80, Table 1, Probability at interim

test_that("Gsponer et al. results match (normal end-point)", {
  skip_on_cran()

  ocRef <- data.frame(
    delta = c(0, 40, 50, 60, 70),
    success = c(1.1, 32.2, 50.0, 67.6, 82.2),
    futile = c(63.3, 6.8, 2.5, 0.8, 0.2)
  )
  sigmaFixed <- 88

  priorT <- mixnorm(c(1, 0, 0.001), sigma = sigmaFixed, param = "mn")
  priorP <- mixnorm(c(1, -49, 20), sigma = sigmaFixed, param = "mn")

  ## the success criteria is for delta which are larger than some
  ## threshold value which is why we set lower.tail=FALSE
  successCrit <- decision2S(c(0.95, 0.5), c(0, 50), FALSE)
  ## the futility criterion acts in the opposite direction
  futilityCrit <- decision2S(c(0.90), c(40), TRUE)

  nT1 <- 20
  nP1 <- 10

  oc <- data.frame(delta = c(0, 40, 50, 60, 70))

  ## Note that due to the fact that only a single mixture component is
  ## used, the decision boundary is a linear function such that only few
  ## evaluations of the boundary are needed to estimate reliably the
  ## spline function

  ## Table 1, probability for interim for success
  oc$success <- oc2S(
    priorP,
    priorT,
    nP1,
    nT1,
    successCrit,
    Ngrid = 1,
    sigma1 = sigmaFixed,
    sigma2 = sigmaFixed
  )(-49, -49 - oc$delta)

  ## Table 1, probability for interim for futility
  oc$futile <- oc2S(
    priorP,
    priorT,
    nP1,
    nT1,
    futilityCrit,
    Ngrid = 1,
    sigma1 = sigmaFixed,
    sigma2 = sigmaFixed
  )(-49, -49 - oc$delta)

  ## Table 1, first three columns, page 74
  oc[-1] <- lapply(100 * oc[-1], round, 1)

  resFutility <- abs(ocRef$futile - oc$futile)
  resSuccess <- abs(ocRef$success - oc$success)

  expect_equal(sum(resFutility > 2), 0, info = "futility")
  expect_equal(sum(resSuccess > 2), 0, info = "success")
})


## failure when doing repeated evaluations which came up in consulting
test_that("Ensure that repeated oc2S evaluation works for normal case", {
  skip_on_cran()

  samp_sigma <- 3

  n_ia <- 38
  n_final <- 2 * n_ia
  n_ia_to_final <- n_final - n_ia
  sem_ia <- samp_sigma / sqrt(n_ia)

  theta_ctl <- 0
  delta <- 1.04

  obs_P <- 0.11
  obs_T <- 1.28

  prior <- mixnorm(c(1, 0, 0.001), sigma = samp_sigma, param = "mn")
  postP_interim <- postmix(prior, m = obs_P, se = sem_ia)
  postT_interim <- postmix(prior, m = obs_T, se = sem_ia)

  successCrit <- decision2S(c(0.9), c(0), FALSE)

  interim_CP <- oc2S(
    postT_interim,
    postP_interim,
    n_ia_to_final,
    n_ia_to_final,
    successCrit,
    sigma1 = samp_sigma,
    sigma2 = samp_sigma
  )

  cpd_ia <- interim_CP(obs_T, obs_P)
  cpd_ia2 <- interim_CP(theta_ctl + delta, theta_ctl)

  expect_number(cpd_ia, lower = 0, upper = 1, finite = TRUE)
  expect_number(cpd_ia2, lower = 0, upper = 1, finite = TRUE)

  ## check that when calculating directly that the results
  ## are close enough
  interim_CPalt <- oc2S(
    postT_interim,
    postP_interim,
    n_ia_to_final,
    n_ia_to_final,
    successCrit,
    sigma1 = samp_sigma,
    sigma2 = samp_sigma
  )
  cpd_ia2alt <- interim_CPalt(theta_ctl + delta, theta_ctl)
  expect_number(
    abs(cpd_ia2 - cpd_ia2alt),
    lower = 0,
    upper = 1E-3,
    finite = TRUE
  )
})

## test against Schmidli et. al, "Robust Meta-Analytic-Predictive
## Priors", Table 2, unif and beta case
## Mixture prior in sample 2: the decision boundary is curved between the
## regime where the informative component dominates and the tails. The
## operating characteristics must not depend on which truths were requested
## before.
oc2S_curved_boundary_scenario <- function() {
  prior1 <- mixnorm(c(1, 0, 100), sigma = 1)
  prior2 <- mixnorm(c(0.8, 0, 1), c(0.2, 0, 2), sigma = 1)
  decision <- decision2S(c(0.95, 0.5), c(0, 2.5), FALSE)
  list(
    make_oc = function() oc2S(prior1, prior2, 1, 1, decision, 1, 1),
    ## independent deterministic integration over directly solved roots
    ref_m3 = 0.31863073,
    ref_p3 = 0.8866051
  )
}

test_that("Normal mixture oc2S is accurate when requested with a broad vector of truths", {
  s <- oc2S_curved_boundary_scenario()
  design <- s$make_oc()
  broad <- design(c(-17, 0, 23), c(-20, -3, 20))[2]
  cached <- design(0, -3)
  fresh <- s$make_oc()(0, -3)
  grid <- withr::with_options(
    list(RBesT.decision2S_boundary = "grid"),
    s$make_oc()(0, -3)
  )
  expect_lt(abs(broad - s$ref_m3), 2e-4)
  expect_lt(abs(cached - s$ref_m3), 2e-4)
  expect_lt(abs(fresh - s$ref_m3), 2e-4)
  expect_lt(abs(broad - grid), 5e-4)
})

test_that("Normal mixture oc2S stays accurate when the probability mass moves", {
  s <- oc2S_curved_boundary_scenario()
  design <- s$make_oc()
  ## request far apart truths first, then truths in between
  design(c(-37, 43), c(-40, 40))
  expect_lt(abs(design(0, -3) - s$ref_m3), 2e-4)
  expect_lt(abs(design(6, 3) - s$ref_p3), 2e-4)
  ## the first result is unchanged by the later request
  expect_lt(abs(design(0, -3) - s$ref_m3), 2e-4)
})

test_that("Normal mixture oc2S does not depend on the order of requests", {
  s <- oc2S_curved_boundary_scenario()
  theta2 <- c(-3, 3, 10, -12, 0, 400, -400)
  at_once <- s$make_oc()(theta2 + 3, theta2)
  design <- s$make_oc()
  one_by_one <- rev(vapply(
    rev(theta2),
    function(t2) design(t2 + 3, t2),
    numeric(1)
  ))
  fresh <- vapply(theta2, function(t2) s$make_oc()(t2 + 3, t2), numeric(1))
  expect_lt(max(abs(at_once - one_by_one)), 2e-4)
  expect_lt(max(abs(at_once - fresh)), 2e-4)
})

test_that("Schmidli et al. results (binary end-point)", {
  skip_on_cran()

  ocRef_inf <- expand.grid(pc = seq(0.1, 0.6, by = 0.1), delta = c(0, 0.3))
  ocRef_inf$ref <- c(
    0,
    1.6,
    6.1,
    13.7,
    26.0,
    44.4, ## beta/delta=0
    81.6,
    87.8,
    93.4,
    97.9,
    99.6,
    100.0 ## beta/delta=0.3
  ) /
    100

  ocRef_uni <- expand.grid(pc = seq(0.1, 0.6, by = 0.1), delta = c(0, 0.3))
  ocRef_uni$ref <- c(
    1.8,
    2.3,
    2.4,
    2.6,
    2.8,
    2.6, ## unif/delta=0
    89.7,
    82.1,
    79.5,
    79.5,
    81.9,
    89.8 ## unif/delta=0.3
  ) /
    100
  dec <- decision2S(0.975, 0, lower.tail = FALSE)

  N <- 40

  prior_inf <- mixbeta(c(1, 4, 16))
  prior_uni <- mixbeta(c(1, 1, 1))

  N_ctl_uni <- N - round(ess(prior_uni, method = "morita"))
  N_ctl_inf <- N - round(ess(prior_inf, method = "morita"))

  design_uni <- oc2S(prior_uni, prior_uni, N, N_ctl_uni, dec)
  design_inf <- oc2S(prior_uni, prior_inf, N, N_ctl_inf, dec)

  res_uni <- design_uni(ocRef_uni$pc + ocRef_uni$delta, ocRef_uni$pc)
  res_inf <- design_inf(ocRef_inf$pc + ocRef_inf$delta, ocRef_inf$pc)

  expect_true(all(abs(100 * (res_uni - ocRef_uni$ref)) < 2.5))
  expect_true(all(abs(100 * (res_inf - ocRef_inf$ref)) < 2.5))
})

## some additional, very simple type I error tests and tests for the
## discrete case of correct critical value behavior

test_scenario <- function(oc_res, ref) {
  resA <- oc_res - ref
  expect_true(all(abs(resA) < oc2S_tolerance()))
}

expect_equal_each <- function(test, expected) {
  for (elem in test) {
    expect_equal(elem, expected)
  }
}

## design object, decision function, posterior function must return
## posterior after updatding the prior with the given value; we assume
## that the priors are the same for sample 1 and 2
test_critical_discrete <- function(boundary_design, decision, posterior, y2) {
  lower.tail <- attr(decision, "lower.tail")
  crit <- boundary_design(y2)
  post2 <- posterior(y2)
  if (lower.tail) {
    expect_equal(decision(posterior(crit - 1), post2), 1)
    expect_equal(decision(posterior(crit), post2), 1)
    expect_equal(decision(posterior(crit + 1), post2), 0)
  } else {
    expect_equal(decision(posterior(crit - 1), post2), 0)
    expect_equal(decision(posterior(crit), post2), 0)
    expect_equal(decision(posterior(crit + 1), post2), 1)
  }
}

oc2S_tolerance <- function() 1e-2

oc2S_alpha <- function() 0.05

oc2S_decision <- function(lower.tail = TRUE) {
  decision2S(1 - oc2S_alpha(), 0, lower.tail = lower.tail)
}

oc2S_binary_scenario <- function(eps = NULL) {
  prior <- mixbeta(c(1, 1, 1))
  dec <- oc2S_decision(TRUE)
  decB <- oc2S_decision(FALSE)
  design_args <- list(prior, prior, 100, 100)
  eps_arg <- if (!is.null(eps)) list(eps = eps)
  list(
    alpha = oc2S_alpha(),
    dec = dec,
    decB = decB,
    design = do.call(oc2S, c(design_args, list(dec), eps_arg)),
    designB = do.call(oc2S, c(design_args, list(decB), eps_arg)),
    boundary_design = decision2S_boundary(prior, prior, 100, 100, dec),
    boundary_designB = decision2S_boundary(prior, prior, 100, 100, decB),
    posterior = function(r) postmix(prior, r = r, n = 100),
    theta = 1:9 / 10
  )
}

oc2S_strong_binary_scenario <- function(reversed = FALSE) {
  prior1 <- mixbeta(c(1, 0.9, 1000), param = "mn")
  prior2 <- mixbeta(c(1, 0.1, 1000), param = "mn")
  if (reversed) {
    prior1 <- mixbeta(c(1, 0.1, 1000), param = "mn")
    prior2 <- mixbeta(c(1, 0.9, 1000), param = "mn")
  }
  dec <- oc2S_decision(TRUE)
  decB <- oc2S_decision(FALSE)
  list(
    dec = dec,
    decB = decB,
    design_lower = oc2S(prior1, prior2, 20, 20, dec),
    design_upper = oc2S(prior1, prior2, 20, 20, decB),
    boundary_design_lower = decision2S_boundary(prior1, prior2, 20, 20, dec),
    boundary_design_upper = decision2S_boundary(prior1, prior2, 20, 20, decB),
    theta = 1:9 / 10
  )
}

oc2S_poisson_scenario <- function() {
  prior <- mixgamma(c(1, 2, 2))
  dec <- oc2S_decision(TRUE)
  decB <- oc2S_decision(FALSE)
  list(
    alpha = oc2S_alpha(),
    dec = dec,
    decB = decB,
    design = oc2S(prior, prior, 100, 100, dec),
    designB = oc2S(prior, prior, 100, 100, decB),
    boundary_design = decision2S_boundary(prior, prior, 100, 100, dec),
    boundary_designB = decision2S_boundary(prior, prior, 100, 100, decB),
    posterior = function(m) postmix(prior, m = m / 100, n = 100),
    theta = seq(0.5, 1.3, by = 0.1)
  )
}

test_that("Binary type I error rate", {
  skip_on_cran()
  sc <- oc2S_binary_scenario()
  test_scenario(sc$design(sc$theta, sc$theta), sc$alpha)
})
test_that("Binary critical value, lower.tail=TRUE", {
  skip_on_cran()
  sc <- oc2S_binary_scenario()
  test_critical_discrete(sc$boundary_design, sc$dec, sc$posterior, 30)
})
test_that("Binary critical value, lower.tail=FALSE", {
  skip_on_cran()
  sc <- oc2S_binary_scenario()
  test_critical_discrete(sc$boundary_designB, sc$decB, sc$posterior, 30)
})
test_that("Binary boundary case, lower.tail=TRUE", {
  skip_on_cran()
  sc <- oc2S_binary_scenario()
  expect_numeric(
    sc$design(1, 1),
    lower = 0,
    upper = 1,
    finite = TRUE,
    any.missing = FALSE
  )
})
test_that("Binary boundary case, lower.tail=FALSE", {
  skip_on_cran()
  sc <- oc2S_binary_scenario()
  expect_numeric(
    sc$designB(0, 0),
    lower = 0,
    upper = 1,
    finite = TRUE,
    any.missing = FALSE
  )
})

test_that("Binary case, no decision change, lower.tail=TRUE, critical value", {
  skip_on_cran()
  sc <- oc2S_strong_binary_scenario()
  expect_equal_each(sc$boundary_design_lower(0:20), -1)
})
test_that("Binary case, no decision change, lower.tail=FALSE, critical value", {
  skip_on_cran()
  sc <- oc2S_strong_binary_scenario()
  expect_equal_each(sc$boundary_design_upper(0:20), 21)
})
test_that("Binary case, no decision change, lower.tail=TRUE, frequency=0", {
  skip_on_cran()
  sc <- oc2S_strong_binary_scenario()
  expect_equal_each(sc$design_lower(sc$theta, sc$theta), 0.0)
})
test_that("Binary case, no decision change, lower.tail=FALSE, frequency=1", {
  skip_on_cran()
  sc <- oc2S_strong_binary_scenario()
  expect_equal_each(sc$design_upper(sc$theta, sc$theta), 1.0)
})

test_that("Binary case, no decision change (reversed), lower.tail=TRUE, critical value", {
  skip_on_cran()
  sc <- oc2S_strong_binary_scenario(reversed = TRUE)
  expect_equal_each(sc$boundary_design_lower(0:20), 20)
})
test_that("Binary case, no decision change (reversed), lower.tail=FALSE, critical value", {
  skip_on_cran()
  sc <- oc2S_strong_binary_scenario(reversed = TRUE)
  expect_equal_each(sc$boundary_design_upper(0:20), -1)
})
test_that("Binary case, no decision change (reversed), lower.tail=TRUE, frequency=0", {
  skip_on_cran()
  sc <- oc2S_strong_binary_scenario(reversed = TRUE)
  expect_equal_each(sc$design_lower(sc$theta, sc$theta), 1.0)
})
test_that("Binary case, no decision change (reversed), lower.tail=FALSE, frequency=1", {
  skip_on_cran()
  sc <- oc2S_strong_binary_scenario(reversed = TRUE)
  expect_equal_each(sc$design_upper(sc$theta, sc$theta), 0.0)
})
test_that("Binary case, log-link", {
  skip_on_cran()
  success <- decision2S(
    pc = c(0.90, 0.50),
    qc = c(log(1), log(0.50)),
    lower.tail = TRUE,
    link = "log"
  )
  prior_pbo <- mixbeta(
    inf1 = c(0.60, 19, 29),
    inf2 = c(0.30, 4, 5),
    rob = c(0.10, 1, 1)
  )
  prior_trt <- mixbeta(c(1, 1 / 3, 1 / 3))
  n_trt <- 50
  n_pbo <- 20
  design_suc <- oc2S(prior_trt, prior_pbo, n_trt, n_pbo, success)
  theta <- seq(0, 1, by = 0.1)
  expect_numeric(
    design_suc(theta, theta),
    lower = 0,
    upper = 1,
    finite = TRUE,
    any.missing = FALSE
  )
})
test_that("Binary case, logit-link", {
  skip_on_cran()
  success <- decision2S(
    pc = c(0.90, 0.50),
    qc = c(log(1), log(0.50)),
    lower.tail = TRUE,
    link = "logit"
  )
  prior_pbo <- mixbeta(
    inf1 = c(0.60, 19, 29),
    inf2 = c(0.30, 4, 5),
    rob = c(0.10, 1, 1)
  )
  prior_trt <- mixbeta(c(1, 1 / 3, 1 / 3))
  n_trt <- 50
  n_pbo <- 20
  design_suc <- oc2S(prior_trt, prior_pbo, n_trt, n_pbo, success)
  theta <- seq(0, 1, by = 0.1)
  expect_numeric(
    design_suc(theta, theta),
    lower = 0,
    upper = 1,
    finite = TRUE,
    any.missing = FALSE
  )
})

test_that("Binary type I error rate", {
  skip_on_cran()
  sc <- oc2S_binary_scenario(eps = 1E-3)
  theta <- seq(0.1, 0.9, by = 0.1)
  test_scenario(sc$design(theta, theta), sc$alpha)
})

test_that("Poisson type I error rate", {
  skip_on_cran()
  sc <- oc2S_poisson_scenario()
  test_scenario(sc$design(sc$theta, sc$theta), sc$alpha)
})
test_that("Poisson crticial value, lower.tail=TRUE", {
  skip_on_cran()
  sc <- oc2S_poisson_scenario()
  test_critical_discrete(sc$boundary_design, sc$dec, sc$posterior, 90)
})
test_that("Poisson crticial value, lower.tail=FALSE", {
  skip_on_cran()
  sc <- oc2S_poisson_scenario()
  test_critical_discrete(sc$boundary_designB, sc$decB, sc$posterior, 90)
})

test_that("Normal OC 2-sample case works for n2=0, crohn-1", {
  crohn_sigma <- 88

  map <- mixnorm(c(0.6, -50, 19), c(0.4, -50, 42), sigma = crohn_sigma)

  ## add a 20% non-informative mixture component
  map_robust <- robustify(map, weight = 0.2, mean = -50, sigma = 88)

  poc <- decision2S(pc = c(0.95, 0.5), qc = c(0, -50), lower.tail = TRUE)

  weak_prior <- mixnorm(c(1, -50, 1), sigma = crohn_sigma, param = "mn")
  n_act <- 40
  ## n_pbo <- 20

  design_noprior_b <- oc2S(
    weak_prior,
    map,
    n_act,
    0,
    poc,
    sigma1 = crohn_sigma,
    sigma2 = crohn_sigma
  )

  expect_numeric(
    design_noprior_b(-20, -30),
    lower = 0,
    upper = 1,
    any.missing = FALSE
  )
})

test_that("Normal OC 2-sample case works for n2=0, crohn-2", {
  crohn_sigma <- 88

  map <- mixnorm(c(1.0, -50, 19), sigma = crohn_sigma)

  ## add a 20% non-informative mixture component
  map_robust <- robustify(map, weight = 0.2, mean = -50, sigma = 88)

  poc <- decision2S(pc = c(0.95, 0.5), qc = c(0, -50), lower.tail = TRUE)

  weak_prior <- mixnorm(c(1, -50, 1), sigma = crohn_sigma, param = "mn")
  n_act <- 40
  ## n_pbo <- 20

  design_noprior_b <- oc2S(
    weak_prior,
    map,
    n_act,
    0,
    poc,
    sigma1 = crohn_sigma,
    sigma2 = crohn_sigma
  )

  expect_numeric(
    design_noprior_b(-20, -30),
    lower = 0,
    upper = 1,
    any.missing = FALSE
  )
})

test_that("Normal OC 2-sample avoids undefined behavior, example 1", {
  skip_on_cran()

  sigma_ref <- 3.2
  ## map_ref <- mixnorm(c(0.51, -2.1, 0.39), c(0.42, -2.1, 0.995), c(0.06, -1.99, 2.32), sigma=sigma_ref)
  ## chagned so that weights sum to 1
  map_ref <- mixnorm(
    c(0.52, -2.1, 0.39),
    c(0.42, -2.1, 0.995),
    c(0.06, -1.99, 2.32),
    sigma = sigma_ref
  )
  prior_flat <- mixnorm(c(1, 0, 100), sigma = sigma_ref)
  alpha <- 0.05
  dec <- decision2S(1 - alpha, 0, lower.tail = FALSE)
  n <- 58
  k <- 2
  design_map <- oc2S(
    prior_flat,
    map_ref,
    n,
    n / k,
    dec,
    sigma1 = sigma_ref,
    sigma2 = sigma_ref
  )
  design_map_2 <- oc2S(
    prior_flat,
    map_ref,
    n,
    n / k,
    dec,
    sigma1 = sigma_ref,
    sigma2 = sigma_ref
  )

  x <- seq(-2.6, -1.6, by = 0.1)
  expect_numeric(design_map(x, x), lower = 0, upper = 1, any.missing = FALSE)
  expect_silent(design_map(-3, -4))
  expect_numeric(design_map(-3, -4), lower = 0, upper = 1, any.missing = FALSE)
  expect_numeric(design_map(-3, 4), lower = 0, upper = 1, any.missing = FALSE)
  expect_numeric(
    design_map(-1.6, -1.6),
    lower = 0,
    upper = 1,
    any.missing = FALSE
  )

  expect_numeric(
    design_map_2(-3, -4),
    lower = 0,
    upper = 1,
    any.missing = FALSE
  )
  expect_numeric(design_map_2(-3, 4), lower = 0, upper = 1, any.missing = FALSE)
  expect_numeric(
    design_map_2(-1.6, -1.6),
    lower = 0,
    upper = 1,
    any.missing = FALSE
  )
  expect_numeric(design_map_2(x, x), lower = 0, upper = 1, any.missing = FALSE)
})

test_that("Normal OC 2-sample works with mixed lower.tail decision criterion", {
  skip_on_cran()

  sigma_ref <- 3.2
  map_ref <- mixnorm(
    c(0.52, -2.1, 0.39),
    c(0.42, -2.1, 0.995),
    c(0.06, -1.99, 2.32),
    sigma = sigma_ref
  )
  prior_flat <- mixnorm(c(1, 0, 100), sigma = sigma_ref)
  alpha <- 0.05

  # Here we have 4 decision scenarios: Go, Stop, In-between 1, In-between 2.
  # For any theta combination, the probability sums to 1 across these 4 scenarios.
  # We need to test the cases n2 == 0 and n2 > 0 because they are computed differently.
  # We test two different qc2 value settings to see that both in-between 1 and in-between 2
  # have a case with positive probabilities and everything works as expected.
  n <- 14
  for (n2 in c(0, 7)) {
    for (qc2 in c(0.5, 0.8)) {
      qc1 <- 1
      pc1 <- 0.5
      pc2 <- 0.6

      dec_go <- decision2S(
        qc = c(qc1, qc2),
        pc = c(pc1, pc2),
        lower.tail = c(FALSE, FALSE)
      )
      dec_stop <- decision2S(
        qc = c(qc1, qc2),
        pc = c(1 - pc1, 1 - pc2),
        lower.tail = c(TRUE, TRUE)
      )
      dec_inbetween_1 <- decision2S(
        qc = c(qc1, qc2),
        pc = c(pc1, 1 - pc2),
        lower.tail = c(FALSE, TRUE)
      )
      dec_inbetween_2 <- decision2S(
        qc = c(qc1, qc2),
        pc = c(1 - pc1, pc2),
        lower.tail = c(TRUE, FALSE)
      )

      get_design_map_n2_pos <- function(dec) {
        oc2S(
          prior_flat,
          map_ref,
          n,
          n2,
          dec,
          sigma1 = sigma_ref,
          sigma2 = sigma_ref
        )
      }

      prob_fun_go <- get_design_map_n2_pos(dec_go)
      prob_fun_stop <- get_design_map_n2_pos(dec_stop)
      prob_fun_inbetween_1 <- get_design_map_n2_pos(dec_inbetween_1)
      prob_fun_inbetween_2 <- get_design_map_n2_pos(dec_inbetween_2)

      x <- seq(-3, 3, by = 0.1)

      prob_go <- prob_fun_go(x, x)
      prob_stop <- prob_fun_stop(x, x)
      prob_inbetween_1 <- prob_fun_inbetween_1(x, x)
      prob_inbetween_2 <- prob_fun_inbetween_2(x, x)

      expect_numeric(prob_go, lower = 0, upper = 1, any.missing = FALSE)
      expect_numeric(prob_stop, lower = 0, upper = 1, any.missing = FALSE)
      expect_numeric(
        prob_inbetween_1,
        lower = 0,
        upper = 1,
        any.missing = FALSE
      )
      expect_numeric(
        prob_inbetween_2,
        lower = 0,
        upper = 1,
        any.missing = FALSE
      )

      total_prob <- prob_go + prob_stop + prob_inbetween_1 + prob_inbetween_2
      expect_true(all(abs(total_prob - 1) < 1e-3))
    }
  }
})

test_that("Binomial OC 2-sample works with mixed lower.tail decision criterion", {
  skip_on_cran()

  map_ref <- mixbeta(
    c(0.6, 19, 29),
    c(0.3, 4, 5),
    c(0.1, 1, 1)
  )
  prior_flat <- mixbeta(c(1, 1, 1))
  alpha <- 0.05

  # Again we have 4 decision scenarios: Go, Stop, In-between 1, In-between 2.
  n <- 20
  for (n2 in c(0, 10)) {
    for (qc2 in c(0.5, 0.8)) {
      qc1 <- 0.5
      pc1 <- 0.5
      pc2 <- 0.6

      dec_go <- decision2S(
        qc = c(qc1, qc2),
        pc = c(pc1, pc2),
        lower.tail = c(FALSE, FALSE)
      )
      dec_stop <- decision2S(
        qc = c(qc1, qc2),
        pc = c(1 - pc1, 1 - pc2),
        lower.tail = c(TRUE, TRUE)
      )
      dec_inbetween_1 <- decision2S(
        qc = c(qc1, qc2),
        pc = c(pc1, 1 - pc2),
        lower.tail = c(FALSE, TRUE)
      )
      dec_inbetween_2 <- decision2S(
        qc = c(qc1, qc2),
        pc = c(1 - pc1, pc2),
        lower.tail = c(TRUE, FALSE)
      )

      get_design_map_n2_pos <- function(dec) {
        oc2S(
          prior_flat,
          map_ref,
          n,
          n2,
          dec
        )
      }

      prob_fun_go <- get_design_map_n2_pos(dec_go)
      prob_fun_stop <- get_design_map_n2_pos(dec_stop)
      prob_fun_inbetween_1 <- get_design_map_n2_pos(dec_inbetween_1)
      prob_fun_inbetween_2 <- get_design_map_n2_pos(dec_inbetween_2)

      x <- seq(0.1, 0.9, by = 0.1)

      prob_go <- prob_fun_go(x, x)
      prob_stop <- prob_fun_stop(x, x)
      prob_inbetween_1 <- prob_fun_inbetween_1(x, x)
      prob_inbetween_2 <- prob_fun_inbetween_2(x, x)

      expect_numeric(prob_go, lower = 0, upper = 1, any.missing = FALSE)
      expect_numeric(prob_stop, lower = 0, upper = 1, any.missing = FALSE)
      expect_numeric(
        prob_inbetween_1,
        lower = 0,
        upper = 1,
        any.missing = FALSE
      )
      expect_numeric(
        prob_inbetween_2,
        lower = 0,
        upper = 1,
        any.missing = FALSE
      )

      total_prob <- prob_go + prob_stop + prob_inbetween_1 + prob_inbetween_2
      expect_true(all(abs(total_prob - 1) < 1e-3))
    }
  }
})

test_that("Binomial deprecated y2 argument of oc2S function", {
  skip_on_cran()

  withr::local_options(lifecycle_verbosity = "warning")

  map_ref <- mixbeta(
    c(0.6, 19, 29),
    c(0.3, 4, 5),
    c(0.1, 1, 1)
  )
  prior_flat <- mixbeta(c(1, 1, 1))

  qc1 <- 0.5
  pc1 <- 0.5
  pc2 <- 0.6
  qc2 <- 0.5

  dec_go <- decision2S(
    qc = c(qc1, qc2),
    pc = c(pc1, pc2),
    lower.tail = c(FALSE, FALSE)
  )

  n <- 20
  n2 <- 4
  design <- oc2S(
    prior_flat,
    map_ref,
    n,
    n2,
    dec_go
  )

  expect_warning(
    design(y2 = 2),
    class = "lifecycle_warning_deprecated"
  )
})


test_that("Poisson OC 2-sample works with mixed lower.tail decision criterion", {
  skip_on_cran()

  map_ref <- mixgamma(
    c(0.6, 19, 29),
    c(0.3, 4, 5),
    c(0.1, 1, 1)
  )
  prior_flat <- mixgamma(c(1, 1, 1))
  alpha <- 0.05

  # Again we have 4 decision scenarios: Go, Stop, In-between 1, In-between 2.
  n <- 20
  for (n2 in c(0, 5)) {
    for (qc2 in c(0.5, 0.8)) {
      qc1 <- 0.5
      pc1 <- 0.5
      pc2 <- 0.6

      dec_go <- decision2S(
        qc = c(qc1, qc2),
        pc = c(pc1, pc2),
        lower.tail = c(FALSE, FALSE)
      )
      dec_stop <- decision2S(
        qc = c(qc1, qc2),
        pc = c(1 - pc1, 1 - pc2),
        lower.tail = c(TRUE, TRUE)
      )
      dec_inbetween_1 <- decision2S(
        qc = c(qc1, qc2),
        pc = c(pc1, 1 - pc2),
        lower.tail = c(FALSE, TRUE)
      )
      dec_inbetween_2 <- decision2S(
        qc = c(qc1, qc2),
        pc = c(1 - pc1, pc2),
        lower.tail = c(TRUE, FALSE)
      )

      get_design_map_n2_pos <- function(dec) {
        oc2S(
          prior_flat,
          map_ref,
          n,
          n2,
          dec
        )
      }

      prob_fun_go <- get_design_map_n2_pos(dec_go)
      prob_fun_stop <- get_design_map_n2_pos(dec_stop)
      prob_fun_inbetween_1 <- get_design_map_n2_pos(dec_inbetween_1)
      prob_fun_inbetween_2 <- get_design_map_n2_pos(dec_inbetween_2)

      x <- seq(0.5, 1.5, by = 0.1)

      prob_go <- prob_fun_go(x, x)
      prob_stop <- prob_fun_stop(x, x)
      prob_inbetween_1 <- prob_fun_inbetween_1(x, x)
      prob_inbetween_2 <- prob_fun_inbetween_2(x, x)
      expect_numeric(prob_go, lower = 0, upper = 1, any.missing = FALSE)
      expect_numeric(prob_stop, lower = 0, upper = 1, any.missing = FALSE)
      expect_numeric(
        prob_inbetween_1,
        lower = 0,
        upper = 1,
        any.missing = FALSE
      )
      expect_numeric(
        prob_inbetween_2,
        lower = 0,
        upper = 1,
        any.missing = FALSE
      )
      total_prob <- prob_go + prob_stop + prob_inbetween_1 + prob_inbetween_2
      expect_true(all(abs(total_prob - 1) < 1e-3))
    }
  }
})

test_that("Poisson deprecated y2 argument of oc2S function", {
  skip_on_cran()

  withr::local_options(lifecycle_verbosity = "warning")

  map_ref <- mixgamma(
    c(0.6, 19, 29),
    c(0.3, 4, 5),
    c(0.1, 1, 1)
  )
  prior_flat <- mixgamma(c(1, 1, 1))
  alpha <- 0.05

  qc1 <- 0.5
  pc1 <- 0.5
  pc2 <- 0.6
  qc2 <- 0.5

  dec_go <- decision2S(
    qc = c(qc1, qc2),
    pc = c(pc1, pc2),
    lower.tail = c(FALSE, FALSE)
  )

  n <- 20
  n2 <- 4
  design <- oc2S(
    prior_flat,
    map_ref,
    n,
    n2,
    dec_go
  )

  expect_warning(
    design(y2 = 3),
    class = "lifecycle_warning_deprecated"
  )
})
