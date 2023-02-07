context("pos1S: 1-Sample Probability of Success")

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
beta  <- 0.2

za <- qnorm(1-alpha)
n1 <- 155
c1 <- theta_ni - za * s / sqrt(n1)

N_samp <- 1E4

thetaA <- c(theta_a, theta_ni)

## standard NI design, tests only statistical significance to be
## smaller than theta_ni with 1-alpha certainty
decA <- decision1S(1 - alpha, theta_ni, lower.tail=TRUE)
decAU <- decision1S(1 - alpha, theta_ni, lower.tail=FALSE)

prior <- mixnorm(c(1,0,100), sigma=s)

## let's say we have 40 events at an interim and a HR of 0.9
ia_dist <- postmix(prior, m=log(0.9), se=s/sqrt(40))

test_pos1S <- function(prior, ia_dist, n, dec, decU) {
    ## the PoS is the expected value of the condition power integrated
    ## over the interim density which is what we check here
    cpo_analytic <- oc1S(prior, n, dec)
    pos_analytic <- pos1S(prior, n, dec)
    samp <- rmix(ia_dist, N_samp)
    pos_mc <- mean(cpo_analytic(samp))
    expect_true(all(abs(pos_mc - pos_analytic(ia_dist)) < eps))
    lower.tail <- attr(dec,"lower.tail")
    if(lower.tail) {
        test_pos1S(prior, ia_dist, n, decU)
    }
}


test_that("Normal PoS 1 sample function matches MC integration of CPO", test_pos1S(prior, ia_dist, n1, decA, decAU))

beta_prior <- mixbeta(c(1, 1, 1))
beta_ia <- postmix(beta_prior, r=20, n=50)
test_that("Binomial PoS 1 sample function matches MC integration of CPO", test_pos1S(beta_prior, beta_ia, n1, decA, decAU))

gamma_prior <- mixgamma(c(1, 1, 1), param="mn")
dec_count  <- decision1S(1-alpha, 1, lower.tail=TRUE)
dec_countU  <- decision1S(1-alpha, 1, lower.tail=FALSE)
gamma_ia <- postmix(gamma_prior, m=0.9, n=40)
test_that("Poisson PoS 1 sample function matches MC integration of CPO", test_pos1S(gamma_prior, gamma_ia, n1, dec_count, dec_countU))

prior_unit_inf <- mixnorm(c(1, 0, 1), sigma=s, param="mn")
post_ia_unit_inf <- postmix(prior_unit_inf, m=-1, n=162)
test_pos1S(prior_unit_inf, post_ia_unit_inf, 459-162, decA, decAU)
test_that("Normal PoS 1 sample function matches MC integration of CPO (more extreme case)", test_pos1S(prior_unit_inf, post_ia_unit_inf, 459-162, decA, decAU))
