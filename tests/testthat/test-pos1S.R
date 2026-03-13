## test the analytical OC function via brute force simulation
set.seed(12354)

## expect results to be 1% exact
eps <- 1e-2

## we test here that the PoS indeed averages over the predictive of an
## informative prior the conditional power.

## Example for Figure 3
s <- 2
theta_ni <- 0.4

theta_a <- 0
alpha <- 0.05
beta <- 0.2

za <- qnorm(1 - alpha)
n1 <- 155
c1 <- theta_ni - za * s / sqrt(n1)

N_samp <- 1E4

thetaA <- c(theta_a, theta_ni)

## standard NI design, tests only statistical significance to be
## smaller than theta_ni with 1-alpha certainty
decA <- decision1S(1 - alpha, theta_ni, lower.tail = TRUE)
decAU <- decision1S(1 - alpha, theta_ni, lower.tail = FALSE)

prior <- mixnorm(c(1, 0, 100), sigma = s)

## let's say we have 40 events at an interim and a HR of 0.9
ia_dist <- postmix(prior, m = log(0.9), se = s / sqrt(40))

test_pos1S <- function(prior, ia_dist, n, dec, decU) {
  ## the PoS is the expected value of the condition power integrated
  ## over the interim density which is what we check here
  cpo_analytic <- oc1S(prior, n, dec)
  pos_analytic <- pos1S(prior, n, dec)
  samp <- rmix(ia_dist, N_samp)
  pos_mc <- mean(cpo_analytic(samp))
  expect_true(all(abs(pos_mc - pos_analytic(ia_dist)) < eps))
  lower.tail <- attr(dec, "lower.tail")
  if (lower.tail) {
    test_pos1S(prior, ia_dist, n, decU)
  }
}


test_that("Normal PoS 1 sample function matches MC integration of CPO", {
  test_pos1S(prior, ia_dist, n1, decA, decAU)
})

beta_prior <- mixbeta(c(1, 1, 1))
beta_ia <- postmix(beta_prior, r = 20, n = 50)
test_that("Binomial PoS 1 sample function matches MC integration of CPO", {
  test_pos1S(beta_prior, beta_ia, n1, decA, decAU)
})

gamma_prior <- mixgamma(c(1, 1, 1), param = "mn")
dec_count <- decision1S(1 - alpha, 1, lower.tail = TRUE)
dec_countU <- decision1S(1 - alpha, 1, lower.tail = FALSE)
gamma_ia <- postmix(gamma_prior, m = 0.9, n = 40)
test_that("Poisson PoS 1 sample function matches MC integration of CPO", {
  test_pos1S(gamma_prior, gamma_ia, n1, dec_count, dec_countU)
})

prior_unit_inf <- mixnorm(c(1, 0, 1), sigma = s, param = "mn")
post_ia_unit_inf <- postmix(prior_unit_inf, m = -1, n = 162)
test_that("Normal PoS 1 sample function matches MC integration of CPO (more extreme case)", {
  test_pos1S(prior_unit_inf, post_ia_unit_inf, 459 - 162, decA, decAU)
})

test_that("Mixed lower.tail usage works for normal PoS calculation", {
  prior <- mixnorm(rob = c(0.2, 0, 2), inf = c(0.8, 2, 2), sigma = 5)
  post_ia <- postmix(prior, m = -1, n = 15)

  dec_lower <- decision1S(pc = 0.5, qc = 1.5, lower.tail = TRUE)
  pos_lower <- pos1S(
    prior,
    n = 50,
    decision = dec_lower
  )
  result_lower <- pos_lower(post_ia)

  dec_upper <- decision1S(pc = 0.6, qc = 0.5, lower.tail = FALSE)
  pos_upper <- pos1S(
    prior,
    n = 50,
    decision = dec_upper
  )
  result_upper <- pos_upper(post_ia)

  dec_mixed <- decision1S(
    qc = c(1.5, 0.5),
    pc = c(0.5, 0.6),
    lower.tail = c(TRUE, FALSE)
  )
  pos_mixed <- pos1S(prior, 50, dec_mixed)
  result_mixed <- pos_mixed(post_ia)

  expected_mixed <- result_lower - (1 - result_upper)
  expect_equal(result_mixed, expected_mixed)
})

test_that("Mixed lower.tail usage works for binomial PoS calculation", {
  prior <- mixbeta(rob = c(0.2, 0.5, 0.5), inf = c(0.8, 0.5, 0.5))
  post_ia <- postmix(prior, r = 20, n = 50)

  dec_lower <- decision1S(pc = 0.5, qc = 0.8, lower.tail = TRUE)
  pos_lower <- pos1S(
    prior,
    n = 50,
    decision = dec_lower
  )
  result_lower <- pos_lower(post_ia)

  dec_upper <- decision1S(pc = 0.6, qc = 0.5, lower.tail = FALSE)
  pos_upper <- pos1S(
    prior,
    n = 50,
    decision = dec_upper
  )
  result_upper <- pos_upper(post_ia)

  dec_mixed <- decision1S(
    qc = c(0.8, 0.5),
    pc = c(0.5, 0.6),
    lower.tail = c(TRUE, FALSE)
  )
  pos_mixed <- pos1S(prior, 50, dec_mixed)
  result_mixed <- pos_mixed(post_ia)

  expected_mixed <- result_lower - (1 - result_upper)
  expect_equal(result_mixed, expected_mixed)
})

test_that("Mixed lower.tail usage works for Poisson PoS calculation", {
  prior <- mixgamma(
    rob = c(0.2, 0.5, 0.5),
    inf = c(0.8, 0.5, 0.5),
    param = "mn"
  )
  post_ia <- postmix(prior, m = 0.9, n = 40)

  dec_lower <- decision1S(pc = 0.5, qc = 1.5, lower.tail = TRUE)
  pos_lower <- pos1S(
    prior,
    n = 50,
    decision = dec_lower
  )
  result_lower <- pos_lower(post_ia)

  dec_upper <- decision1S(pc = 0.6, qc = 0.5, lower.tail = FALSE)
  pos_upper <- pos1S(
    prior,
    n = 50,
    decision = dec_upper
  )
  result_upper <- pos_upper(post_ia)

  dec_mixed <- decision1S(
    qc = c(1.5, 0.5),
    pc = c(0.5, 0.6),
    lower.tail = c(TRUE, FALSE)
  )
  pos_mixed <- pos1S(prior, 50, dec_mixed)
  result_mixed <- pos_mixed(post_ia)

  expected_mixed <- result_lower - (1 - result_upper)
  expect_equal(result_mixed, expected_mixed)
})
