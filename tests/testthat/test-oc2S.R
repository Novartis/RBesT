## test the analytical OC function via brute force simulation
set.seed(12354)

prior1 <- mixnorm(c(0.3, -0.2, 2), c(0.7, 0, 50), sigma = 1)
prior2 <- mixnorm(c(1.0, 0, 50), sigma = 1)

N1 <- 10
N2 <- 20

## type I error fairly large to 20% to make it easier to test (less
## simulations needed for accurate results)
pcrit <- 0.80
qcrit <- 0

## theta2 set such that we have about 75% power under this truth
theta1 <- 0
theta2 <- 0.5

Nsim <- 1e4

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
  mean_sd1 <- sigma(prior1) / sqrt(N1)
  mean_sd2 <- sigma(prior2) / sqrt(N2)

  mean_prior1 <- prior1
  sigma(mean_prior1) <- mean_sd1
  mean_prior2 <- prior2
  sigma(mean_prior2) <- mean_sd2

  mean_samp1 <- rnorm(Nsim, theta1, mean_sd1)
  mean_samp2 <- rnorm(Nsim, theta2, mean_sd2)

  dec <- rep(NA, Nsim)

  for (i in 1:Nsim) {
    post1 <- postmix(mean_prior1, m = mean_samp1[i], se = mean_sd1)
    post2 <- postmix(mean_prior2, m = mean_samp2[i], se = mean_sd2)
    dec[i] <- as.numeric(pmix(RBesT:::mixnormdiff(post1, post2), qcrit) > pcrit)
  }

  mean(dec)
}

Voc2S_normal_MC <- Vectorize(oc2S_normal_MC, c("theta1", "theta2"))

## first test that the analytic difference distribution for normal
## works as expected

test_that("Analytical convolution of normal mixture matches numerical integration result", {
  skip_on_cran()

  pdiff <- RBesT:::mixnormdiff(prior1, prior2)
  x <- seq(-20, 20, length = 21)
  d1 <- dmix(pdiff, x)
  d2 <- dmixdiff(prior1, prior2, x)
  dres <- abs(d1 - d2)
  expect_equal(sum(dres > 1e-5), 0)
  p1 <- pmix(pdiff, x)
  p2 <- pmixdiff(prior1, prior2, x)
  pres <- 100 * abs(p1 - p2)
  expect_equal(sum(pres > 2), 0)
})

## test that the type I error is matching, i.e. is not off by more than 2%
test_that("Type I error is matching between MC and analytical computations in the normal mixture case", {
  skip_on_cran()

  x <- c(-2, 0)
  alpha <- oc2S(
    prior1,
    prior2,
    N1,
    N2,
    decision2S(pcrit, qcrit),
    sigma1 = sigma(prior1),
    sigma2 = sigma(prior2)
  )(x, x)
  alphaMC <- Voc2S_normal_MC(prior1, prior2, N1, N2, x, x, pcrit, qcrit)
  res <- 100 * abs(alpha - alphaMC)
  expect_equal(sum(res > 2), 0)
})


## test that the power is matching, i.e. is not off by more than 2%
test_that("Power is matching between MC and analytical computations in the normal mixture case", {
  skip_on_cran()

  power <- oc2S(
    prior1,
    prior2,
    N1,
    N2,
    decision2S(pcrit, qcrit),
    sigma1 = sigma(prior1),
    sigma2 = sigma(prior2)
  )(theta1, theta2)
  powerMC <- oc2S_normal_MC(
    prior1,
    prior2,
    N1,
    N2,
    theta1,
    theta2,
    pcrit,
    qcrit
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
  expect_true(all(abs(resA) < eps))
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

## expect results to be 1% exact
eps <- 1e-2
alpha <- 0.05

dec <- decision2S(1 - alpha, 0, lower.tail = TRUE)
decB <- decision2S(1 - alpha, 0, lower.tail = FALSE)

## test binary case

beta_prior <- mixbeta(c(1, 1, 1))
if (!run_on_cran()) {
  design_binary <- oc2S(beta_prior, beta_prior, 100, 100, dec)
  boundary_design_binary <- decision2S_boundary(
    beta_prior,
    beta_prior,
    100,
    100,
    dec
  )
  design_binaryB <- oc2S(beta_prior, beta_prior, 100, 100, decB)
  boundary_design_binaryB <- decision2S_boundary(
    beta_prior,
    beta_prior,
    100,
    100,
    decB
  )
} else {
  design_binary <- function(...) {
    return(0.1)
  }
  design_binaryB <- function(...) {
    return(0.1)
  }
  boundary_design_binary <- function(...) {
    return(0.1)
  }
  boundary_design_binaryB <- function(...) {
    return(0.1)
  }
}
posterior_binary <- function(r) postmix(beta_prior, r = r, n = 100)
p_test <- 1:9 / 10
test_that("Binary type I error rate", {
  skip_on_cran()
  test_scenario(design_binary(p_test, p_test), alpha)
})
test_that("Binary crticial value, lower.tail=TRUE", {
  skip_on_cran()
  test_critical_discrete(boundary_design_binary, dec, posterior_binary, 30)
})
test_that("Binary crticial value, lower.tail=FALSE", {
  skip_on_cran()
  test_critical_discrete(boundary_design_binaryB, decB, posterior_binary, 30)
})
test_that("Binary boundary case, lower.tail=TRUE", {
  skip_on_cran()
  expect_numeric(
    design_binary(1, 1),
    lower = 0,
    upper = 1,
    finite = TRUE,
    any.missing = FALSE
  )
})
test_that("Binary boundary case, lower.tail=FALSE", {
  skip_on_cran()
  expect_numeric(
    design_binaryB(0, 0),
    lower = 0,
    upper = 1,
    finite = TRUE,
    any.missing = FALSE
  )
})

## check case where decision never changes due to prior being too
## strong

beta_prior1 <- mixbeta(c(1, 0.9, 1000), param = "mn")
beta_prior2 <- mixbeta(c(1, 0.1, 1000), param = "mn")
design_lower <- oc2S(beta_prior1, beta_prior2, 20, 20, dec) ## always 0
design_upper <- oc2S(beta_prior1, beta_prior2, 20, 20, decB) ## always 1
boundary_design_lower <- decision2S_boundary(
  beta_prior1,
  beta_prior2,
  20,
  20,
  dec
) ## always 0
boundary_design_upper <- decision2S_boundary(
  beta_prior1,
  beta_prior2,
  20,
  20,
  decB
) ## always 1

test_that("Binary case, no decision change, lower.tail=TRUE, critical value", {
  skip_on_cran()
  expect_equal_each(boundary_design_lower(0:20), -1)
})
test_that("Binary case, no decision change, lower.tail=FALSE, critical value", {
  skip_on_cran()
  expect_equal_each(boundary_design_upper(0:20), 21)
})
test_that("Binary case, no decision change, lower.tail=TRUE, frequency=0", {
  skip_on_cran()
  expect_equal_each(design_lower(p_test, p_test), 0.0)
})
test_that("Binary case, no decision change, lower.tail=FALSE, frequency=1", {
  skip_on_cran()
  expect_equal_each(design_upper(p_test, p_test), 1.0)
})


if (!run_on_cran()) {
  design_lower_rev <- oc2S(beta_prior2, beta_prior1, 20, 20, dec) ## always 1
  design_upper_rev <- oc2S(beta_prior2, beta_prior1, 20, 20, decB) ## always 0
  boundary_design_lower_rev <- decision2S_boundary(
    beta_prior2,
    beta_prior1,
    20,
    20,
    dec
  ) ## always 1
  boundary_design_upper_rev <- decision2S_boundary(
    beta_prior2,
    beta_prior1,
    20,
    20,
    decB
  ) ## always 0
} else {
  design_lower_rev <- function(...) {
    return(1)
  }
  design_upper_rev <- function(...) {
    return(0)
  }
  boundary_design_lower_rev <- function(...) {
    return(1)
  }
  boundary_design_upper_rev <- function(...) {
    return(0)
  }
}

test_that("Binary case, no decision change (reversed), lower.tail=TRUE, critical value", {
  skip_on_cran()
  expect_equal_each(boundary_design_lower_rev(0:20), 20)
})
test_that("Binary case, no decision change (reversed), lower.tail=FALSE, critical value", {
  skip_on_cran()
  expect_equal_each(boundary_design_upper_rev(0:20), -1)
})
test_that("Binary case, no decision change (reversed), lower.tail=TRUE, frequency=0", {
  skip_on_cran()
  expect_equal_each(design_lower_rev(p_test, p_test), 1.0)
})
test_that("Binary case, no decision change (reversed), lower.tail=FALSE, frequency=1", {
  skip_on_cran()
  expect_equal_each(design_upper_rev(p_test, p_test), 0.0)
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

## check approximate method

beta_prior <- mixbeta(c(1, 1, 1))
design_binary_eps <- oc2S(beta_prior, beta_prior, 100, 100, dec, eps = 1E-3)
p_test <- seq(0.1, 0.9, by = 0.1)
test_that("Binary type I error rate", {
  skip_on_cran()
  test_scenario(design_binary_eps(p_test, p_test), alpha)
})

## 22 Nov 2017: disabled test as we trigger always calculation of the
## boundaries as of now.
## test_that("Binary results cache expands", {
##               design_binary_eps <- oc2S(beta_prior, beta_prior, 100, 100, dec, eps=1E-3)
##               design_binary_eps(theta1=0.99, theta2=0.8)
##               ## in this case the cache boundaries do not cover the
##               ## critical value
##               expect_true(is.na(design_binary_eps(theta1=0.99, y2=80)))
##               ## while now they do as theta1 is set to 0.1 and 0.9
##               ## internally which triggers recalculation of the
##               ## internal boundaries
##               expect_true(!is.na(design_binary_eps(y2=80)))
##           })

## test poisson case

gamma_prior <- mixgamma(c(1, 2, 2))

design_poisson <- oc2S(gamma_prior, gamma_prior, 100, 100, dec)
design_poissonB <- oc2S(gamma_prior, gamma_prior, 100, 100, decB)
boundary_design_poisson <- decision2S_boundary(
  gamma_prior,
  gamma_prior,
  100,
  100,
  dec
)
boundary_design_poissonB <- decision2S_boundary(
  gamma_prior,
  gamma_prior,
  100,
  100,
  decB
)
posterior_poisson <- function(m) postmix(gamma_prior, m = m / 100, n = 100)
lambda_test <- seq(0.5, 1.3, by = 0.1)
test_that("Poisson type I error rate", {
  skip_on_cran()
  test_scenario(design_poisson(lambda_test, lambda_test), alpha)
})
test_that("Poisson crticial value, lower.tail=TRUE", {
  skip_on_cran()
  test_critical_discrete(boundary_design_poisson, dec, posterior_poisson, 90)
})
test_that("Poisson crticial value, lower.tail=FALSE", {
  skip_on_cran()
  test_critical_discrete(boundary_design_poissonB, decB, posterior_poisson, 90)
})
## 22 Nov 2017: disabled test as we trigger always calculation of the
## boundaries as of now.
## test_that("Poisson results cache expands", {
##              design_poisson  <- oc2S(gamma_prior, gamma_prior, 100, 100, dec)
##              design_poisson(theta1=1, theta2=c(0.7,1))
##              expect_true(sum(is.na(design_poisson(y2=70:90)) ) == 4)
##              expect_true(sum(is.na(design_poisson(theta1=c(0.01, 1), y2=70:90)) ) == 0)
##          })

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
