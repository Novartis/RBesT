## test the analytical OC function via brute force simulation

## expect results to be 1% exact
pos2S_tolerance <- function() 1e-2

pos2S_sample_size <- function() 1E4

pos2S_normal_scenario <- function() {
  prior1 <- mixnorm(c(0.3, -0.2, 2), c(0.7, 0, 50), sigma = 1)
  prior2 <- mixnorm(c(1.0, 0, 50), sigma = 1)
  n1 <- 30
  n2 <- 40
  s <- 1
  pcrit <- 0.80
  qcrit <- 0

  list(
    dec = decision2S(pcrit, qcrit),
    decU = decision2S(pcrit, qcrit, lower.tail = FALSE),
    ia1 = postmix(prior1, m = 0.2, se = s / sqrt(15)),
    ia2 = postmix(prior2, m = 0, se = s / sqrt(15)),
    n1 = n1,
    n2 = n2,
    prior1 = prior1,
    prior2 = prior2,
    pcrit = pcrit,
    qcrit = qcrit
  )
}

pos2S_beta_scenario <- function() {
  sc <- pos2S_normal_scenario()
  prior <- mixbeta(c(1, 1, 1))
  list(
    dec = sc$dec,
    decU = sc$decU,
    ia1 = postmix(prior, r = 20, n = 50),
    ia2 = postmix(prior, r = 30, n = 50),
    n1 = sc$n1,
    n2 = sc$n2,
    prior1 = prior,
    prior2 = prior
  )
}

pos2S_gamma_scenario <- function() {
  alpha <- 0.05
  prior <- mixgamma(c(1, 1, 1), param = "mn")
  list(
    dec = decision2S(1 - alpha, 0, lower.tail = TRUE),
    decU = decision2S(1 - alpha, 0, lower.tail = FALSE),
    ia1 = postmix(prior, m = 0.7, n = 60),
    ia2 = postmix(prior, m = 1.2, n = 60),
    n1 = 30,
    n2 = 40,
    prior1 = prior,
    prior2 = prior
  )
}

## we test here that the PoS indeed averages over the predictive of an
## informative prior weighting the conditional power respectively.

test_pos2S <- function(prior1, prior2, ia_dist1, ia_dist2, n1, n2, dec, decU) {
  skip_on_cran()
  withr::local_seed(620339)

  ## the PoS is the expected value of the condition power integrated
  ## over the interim density which is what we check here
  cpo_analytic <- oc2S(prior1, prior2, n1, n2, dec)
  pos_analytic <- pos2S(prior1, prior2, n1, n2, dec)
  samp1 <- rmix(ia_dist1, pos2S_sample_size())
  samp2 <- rmix(ia_dist2, pos2S_sample_size())
  pos_mc <- mean(cpo_analytic(samp1, samp2))
  ## print(pos_mc)
  ## print(pos_analytic(ia_dist1, ia_dist2))
  expect_true(all(
    abs(pos_mc - pos_analytic(ia_dist1, ia_dist2)) < pos2S_tolerance()
  ))
  lower.tail <- attr(dec, "lower.tail")
  if (lower.tail) {
    ## cat("Also testing lower.tail=FALSE\n")
    test_pos2S(prior1, prior2, ia_dist1, ia_dist2, n1, n2, decU)
  }
}

test_that("Normal PoS 2 sample function matches MC integration of CPO", {
  sc <- pos2S_normal_scenario()
  test_pos2S(
    sc$prior1,
    sc$prior2,
    sc$ia1,
    sc$ia2,
    sc$n1,
    sc$n2,
    sc$dec,
    sc$decU
  )
})

## also run a MC comparison
pos2S_normal_MC <- function(
  prior1,
  prior2,
  N1,
  N2,
  dtheta1,
  dtheta2,
  pcrit = 0.975,
  qcrit = 0
) {
  skip_on_cran()
  withr::local_seed(918274)

  mean_sd1 <- sigma(prior1) / sqrt(N1)
  mean_sd2 <- sigma(prior2) / sqrt(N2)

  mean_prior1 <- prior1
  sigma(mean_prior1) <- mean_sd1
  mean_prior2 <- prior2
  sigma(mean_prior2) <- mean_sd2

  pred_dtheta1 <- preddist(dtheta1, n = N1) ## , sigma=mean_sd1)
  pred_dtheta2 <- preddist(dtheta2, n = N2) ## , sigma=mean_sd1)

  ## mean_samp1 <- rnorm(Nsim, theta1, mean_sd1)
  ## mean_samp2 <- rnorm(Nsim, theta2, mean_sd2)
  mean_samp1 <- rmix(pred_dtheta1, pos2S_sample_size())
  mean_samp2 <- rmix(pred_dtheta2, pos2S_sample_size())

  dec <- rep(NA, pos2S_sample_size())

  for (i in seq_len(pos2S_sample_size())) {
    post1 <- postmix(mean_prior1, m = mean_samp1[i], se = mean_sd1)
    post2 <- postmix(mean_prior2, m = mean_samp2[i], se = mean_sd2)
    dec[i] <- as.numeric(pmix(RBesT:::mixnormdiff(post1, post2), qcrit) > pcrit)
  }

  mean(dec)
}

test_that("Normal PoS 2 sample function matches MC integration", {
  sc <- pos2S_normal_scenario()
  pos_mc <- pos2S_normal_MC(
    sc$prior1,
    sc$prior2,
    sc$n1,
    sc$n2,
    sc$ia1,
    sc$ia2,
    pcrit = sc$pcrit,
    qcrit = sc$qcrit
  )
  pos_analytic <- pos2S(sc$prior1, sc$prior2, sc$n1, sc$n2, sc$dec)
  expect_true(all(
    abs(pos_mc - pos_analytic(sc$ia1, sc$ia2)) < pos2S_tolerance()
  ))
})

test_that("Normal mixture PoS 2 sample agrees between adaptive and grid boundary", {
  ## curved decision boundary from a mixture prior in sample 2, broad
  ## predictive distribution for sample 2
  prior1 <- mixnorm(c(1, 0, 100), sigma = 1)
  prior2 <- mixnorm(c(0.8, 0, 1), c(0.2, 0, 2), sigma = 1)
  decision <- decision2S(c(0.95, 0.5), c(0, 2.5), FALSE)
  truth1 <- mixnorm(c(1, 1, 2), sigma = 1)
  truth2 <- mixnorm(c(0.9, -3, 0.5), c(0.1, 0, 8), sigma = 1)
  pos_with <- function(method) {
    withr::with_options(
      list(RBesT.decision2S_boundary = method),
      pos2S(prior1, prior2, 1, 1, decision, sigma1 = 1, sigma2 = 1)(
        truth1,
        truth2
      )
    )
  }
  expect_lt(abs(pos_with("adaptive") - pos_with("grid")), 5e-4)
})

test_that("Binomial PoS 2 sample function matches MC integration of CPO", {
  sc <- pos2S_beta_scenario()
  test_pos2S(
    sc$prior1,
    sc$prior2,
    sc$ia1,
    sc$ia2,
    sc$n1,
    sc$n2,
    sc$dec,
    sc$decU
  )
})


test_that("Poisson PoS 2 sample function matches MC integration of CPO", {
  sc <- pos2S_gamma_scenario()
  test_pos2S(
    sc$prior1,
    sc$prior2,
    sc$ia1,
    sc$ia2,
    sc$n1,
    sc$n2,
    sc$dec,
    sc$decU
  )
})

test_that("Binomial PoS 2 with IA returns results", {
  ## reported by user
  successCrit <- decision2S(c(0.9), c(0), lower.tail = FALSE)
  n0 <- 50
  n <- 100
  n_alt <- 140
  neutr_prior <- mixbeta(c(1, 1 / 3, 1 / 3))
  post_placeboIA <- postmix(neutr_prior, r = 13, n = n0)
  post_treatIA <- postmix(neutr_prior, r = 3, n = n0)
  # Criterion for PPoS at IA
  pos_final <- pos2S(post_treatIA, post_placeboIA, n - n0, n - n0, successCrit)
  pos_final_alt <- pos2S(
    post_treatIA,
    post_placeboIA,
    n_alt - n0,
    n_alt - n0,
    successCrit
  )
  # Predictive proba of success at the end
  expect_number(
    pos_final_alt(post_treatIA, post_placeboIA),
    na.ok = FALSE,
    lower = 0,
    upper = 1,
    finite = TRUE,
    null.ok = FALSE
  )
  expect_number(
    pos_final(post_treatIA, post_placeboIA),
    na.ok = FALSE,
    lower = 0,
    upper = 1,
    finite = TRUE,
    null.ok = FALSE
  )
})

test_that("Mixed lower.tail usage works for normal PoS calculation", {
  prior1 <- mixnorm(rob = c(0.2, 0, 2), inf = c(0.8, 2, 2), sigma = 5)
  prior2 <- mixnorm(rob = c(0.2, 0, 2), inf = c(0.8, 2, 2), sigma = 5)
  post_ia1 <- postmix(prior1, m = -1, n = 15)
  post_ia2 <- postmix(prior2, m = -0.5, n = 15)

  dec_lower <- decision2S(pc = 0.5, qc = 1.5, lower.tail = TRUE)
  pos_lower <- pos2S(
    prior1,
    prior2,
    n1 = 50,
    n2 = 50,
    decision = dec_lower
  )
  result_lower <- pos_lower(post_ia1, post_ia2)

  dec_upper <- decision2S(pc = 0.6, qc = 0.5, lower.tail = FALSE)
  pos_upper <- pos2S(
    prior1,
    prior2,
    n1 = 50,
    n2 = 50,
    decision = dec_upper
  )
  result_upper <- pos_upper(post_ia1, post_ia2)

  dec_mixed <- decision2S(
    qc = c(1.5, 0.5),
    pc = c(0.5, 0.6),
    lower.tail = c(TRUE, FALSE)
  )
  pos_mixed <- pos2S(prior1, prior2, 50, 50, dec_mixed)
  result_mixed <- pos_mixed(post_ia1, post_ia2)

  expected_mixed <- result_lower - (1 - result_upper)
  expect_equal(result_mixed, expected_mixed, tolerance = 1e-5)
})

test_that("Mixed lower.tail usage works for binomial PoS calculation", {
  prior1 <- mixbeta(c(0.3, 0.2, 0.5), c(0.7, 2, 10))
  prior2 <- mixbeta(c(0.3, 0.2, 0.5), c(0.7, 2, 10))
  post_ia1 <- postmix(prior1, r = 10, n = 15)
  post_ia2 <- postmix(prior2, r = 12, n = 15)

  dec_lower <- decision2S(pc = 0.5, qc = 0.5, lower.tail = TRUE)
  pos_lower <- pos2S(
    prior1,
    prior2,
    n1 = 50,
    n2 = 50,
    decision = dec_lower
  )
  result_lower <- pos_lower(post_ia1, post_ia2)

  dec_upper <- decision2S(pc = 0.6, qc = 0.4, lower.tail = FALSE)
  pos_upper <- pos2S(
    prior1,
    prior2,
    n1 = 50,
    n2 = 50,
    decision = dec_upper
  )
  result_upper <- pos_upper(post_ia1, post_ia2)

  dec_mixed <- decision2S(
    qc = c(0.5, 0.4),
    pc = c(0.5, 0.6),
    lower.tail = c(TRUE, FALSE)
  )
  pos_mixed <- pos2S(prior1, prior2, 50, 50, dec_mixed)
  result_mixed <- pos_mixed(post_ia1, post_ia2)

  expected_mixed <- result_lower - (1 - result_upper)
  expect_equal(result_mixed, expected_mixed, tolerance = 1e-5)
})

test_that("Mixed lower.tail usage works for Poisson PoS calculation", {
  prior1 <- mixgamma(c(0.3, 0.2, 0.5), c(0.7, 2, 10), param = "mn")
  prior2 <- mixgamma(c(0.3, 0.2, 0.5), c(0.7, 2, 10), param = "mn")
  post_ia1 <- postmix(prior1, m = 1.0, n = 15)
  post_ia2 <- postmix(prior2, m = 1.2, n = 15)

  dec_lower <- decision2S(pc = 0.5, qc = 0.5, lower.tail = TRUE)
  pos_lower <- pos2S(
    prior1,
    prior2,
    n1 = 50,
    n2 = 50,
    decision = dec_lower
  )
  result_lower <- pos_lower(post_ia1, post_ia2)

  dec_upper <- decision2S(pc = 0.6, qc = 0.4, lower.tail = FALSE)
  pos_upper <- pos2S(
    prior1,
    prior2,
    n1 = 50,
    n2 = 50,
    decision = dec_upper
  )
  result_upper <- pos_upper(post_ia1, post_ia2)

  dec_mixed <- decision2S(
    qc = c(0.5, 0.4),
    pc = c(0.5, 0.6),
    lower.tail = c(TRUE, FALSE)
  )
  pos_mixed <- pos2S(prior1, prior2, 50, 50, dec_mixed)
  result_mixed <- pos_mixed(post_ia1, post_ia2)

  expected_mixed <- result_lower - (1 - result_upper)
  expect_equal(result_mixed, expected_mixed, tolerance = 1e-4)
})

test_that("Negative binomial PoS 2 sample with family-based sigma works", {
  ## Regression test: pos2S with a NB family (sigma determined by family

  ## variance function) previously failed in solve_boundary2S_normMix
  ## because postmix returned zero-row results for extreme eta values,
  ## causing NA in the boundary search loop.
  skip_on_cran()

  nb_family <- MASS::negative.binomial(theta = 1 / 1.9)
  lambda_ctrl <- 1.8
  log_exposure <- log(1) # 1 year follow-up

  ## sigma on the log-rate scale (used only for prior specification)
  sigma_ctrl <- sqrt((1 + 1.9 * lambda_ctrl * 1) / (lambda_ctrl * 1))

  prior_ctrl <- mixnorm(c(1, log(lambda_ctrl), 10), sigma = sigma_ctrl)
  prior_treat <- mixnorm(c(1, log(lambda_ctrl), 10), sigma = sigma_ctrl)

  ## Interim posteriors (as if 45 patients observed per arm)
  se_ia <- sigma_ctrl / sqrt(45)
  post_treat <- postmix(prior_treat, m = log(lambda_ctrl) - 0.3, se = se_ia)
  post_ctrl <- postmix(prior_ctrl, m = log(lambda_ctrl) + 0.1, se = se_ia)

  decision <- decision2S(0.975, 0, lower.tail = TRUE)

  ## pos2S setup should not error
  pos_fn <- pos2S(
    post_treat,
    post_ctrl,
    n1 = 105,
    n2 = 105,
    decision = decision,
    family = nb_family,
    offset1 = log_exposure,
    offset2 = log_exposure
  )

  ## Evaluating PoS should return a valid probability
  pos_val <- pos_fn(post_treat, post_ctrl)
  expect_number(
    pos_val,
    na.ok = FALSE,
    lower = 0,
    upper = 1,
    finite = TRUE,
    null.ok = FALSE
  )
})
