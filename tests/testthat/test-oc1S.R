context("oc1S: 1-Sample Operating Characteristics")

## test the analytical OC function via brute force simulation
set.seed(12354)

## expect results to be 1% exact
eps <- 1e-2

## here we test against the reference Neuenschwander et al.,
## Statist. Med. 2011, 30, 1618 (*the* double criterion paper)

## Example for Figure 3
s <- 2
theta_ni <- 0.4

theta_a <- 0
alpha <- 0.05
beta  <- 0.2

za <- qnorm(1-alpha)
n1 <- 155
c1 <- theta_ni - za * s / sqrt(n1)

thetaA <- c(theta_a, theta_ni)

## standard NI design, tests only statistical significance to be
## smaller than theta_ni with 1-alpha certainty
decA <- decision1S(1 - alpha, theta_ni, lower.tail=TRUE)

prior <- mixnorm(c(1,0,100), sigma=s)

test_scenario <- function(oc_res, ref) {
    resA <- oc_res - ref
    expect_true(all(abs(resA) < eps))
}

test_that("Classical NI design critical value", expect_true( abs(decision1S_boundary(prior, 155, decA) - c1) < eps))

## n set to give power 80% to detect 0 and type I error 5% for no
## better than theta_ni
test_that("Classical NI design at target    sample size for OCs", test_scenario(oc1S(prior, 155, decA)(thetaA), c(1-beta, alpha)) )
test_that("Classical NI design at increased sample size for OCs", test_scenario(oc1S(prior, 233, decA)(thetaA), c(1-0.08, alpha)) )
test_that("Classical NI design at decreased sample size for OCs", test_scenario(oc1S(prior,  77, decA)(thetaA), c(1-0.45, alpha)) )

## now double criterion with indecision point (mean estimate must be
## lower than this)
theta_c <- c1

## statistical significance
dec1 <- decision1S(1-alpha, theta_ni, lower.tail=TRUE)
## require mean to be at least as good as theta_c
dec2 <- decision1S(0.5, theta_c, lower.tail=TRUE)
## combination
decComb <- decision1S(c(1-alpha, 0.5), c(theta_ni, theta_c), lower.tail=TRUE)

thetaD <- c(theta_c, theta_ni)

## since theta_c == c1, both decision criteria are the same for n =
## 155
test_that("Double criterion NI design at target    sample size for OCs, combined      ", test_scenario(oc1S(prior, 155, decComb)(thetaD), c(0.50, alpha)) )
test_that("Double criterion NI design at target    sample size for OCs, stat criterion", test_scenario(oc1S(prior, 155, dec1)(thetaD),    c(0.50, alpha)) )
test_that("Double criterion NI design at target    sample size for OCs, mean criterion", test_scenario(oc1S(prior, 155, dec2)(thetaD),    c(0.50, alpha)) )

## at an increased sample size only the mean criterion is active
test_that("Double criterion NI design at increased sample size for OCs, combined      ", test_scenario(oc1S(prior, 233, decComb)(thetaD), c(0.50, 0.02)) )
test_that("Double criterion NI design at increased sample size for OCs, mean criterion", test_scenario(oc1S(prior, 233, dec2)(thetaD),    c(0.50, 0.02)) )

## at a  decreased sample size only the stat criterion is active
test_that("Double criterion NI design at decreased sample size for OCs, combined      ", test_scenario(oc1S(prior,  78, decComb)(thetaD), c(1-0.68, alpha)) )
test_that("Double criterion NI design at decreased sample size for OCs, stat criterion", test_scenario(oc1S(prior,  78, dec1)(thetaD),    c(1-0.68, alpha)) )

## test type 1 error and correctness of critical values wrt to
## lower.tail=TRUE/FALSE
dec1b <- decision1S(1-alpha, theta_ni, lower.tail=FALSE)

## design object, decision function, posterior function must return
## posterior after updatding the prior with the given value
test_critical_discrete <- function(crit, decision, posterior) {
    lower.tail <- attr(decision, "lower.tail")
    if(lower.tail) {
        expect_equal(decision(posterior(crit-1)), 1)
        expect_equal(decision(posterior(crit  )), 1)
        expect_equal(decision(posterior(crit+1)), 0)
    } else {
        expect_equal(decision(posterior(crit-1)), 0)
        expect_equal(decision(posterior(crit  )), 0)
        expect_equal(decision(posterior(crit+1)), 1)
    }
}

## binary case
beta_prior <- mixbeta(c(1, 1, 1))
design_binary <- oc1S(beta_prior, 1000, dec1)
design_binaryB <- oc1S(beta_prior, 1000, dec1b)
crit1 <- decision1S_boundary(beta_prior, 1000, dec1)
crit1B <- decision1S_boundary(beta_prior, 1000, dec1b)
posterior_binary <- function(r) postmix(beta_prior, r=r, n=1000)
test_that("Binary type I error rate", test_scenario(design_binary(theta_ni), alpha))
test_that("Binary crticial value, lower.tail=TRUE",  test_critical_discrete(crit1,  dec1,  posterior_binary))
test_that("Binary crticial value, lower.tail=FALSE", test_critical_discrete(crit1B, dec1b, posterior_binary))

test_that("Binary boundary case, lower.tail=TRUE",  expect_numeric(design_binary( 1), lower=0, upper=1, finite=TRUE, any.missing=FALSE))
test_that("Binary boundary case, lower.tail=FALSE", expect_numeric(design_binaryB(0), lower=0, upper=1, finite=TRUE, any.missing=FALSE))

## poisson case
gamma_prior <- mixgamma(c(1, 2, 2))
dec_count  <- decision1S(1-alpha, 1, lower.tail=TRUE)
dec_countB <- decision1S(1-alpha, 1, lower.tail=FALSE)
design_poisson  <- oc1S(gamma_prior, 1000, dec_count)
design_poissonB <- oc1S(gamma_prior, 1000, dec_countB)
pcrit1 <- decision1S_boundary(gamma_prior, 1000, dec_count)
pcrit1B <- decision1S_boundary(gamma_prior, 1000, dec_countB)
posterior_poisson <- function(m) postmix(gamma_prior, m=m/1000, n=1000)
test_that("Poisson type I error rate", test_scenario(design_poisson(1), alpha) )
test_that("Poisson critical value, lower.tail=TRUE",  test_critical_discrete(pcrit1,  dec_count,  posterior_poisson))
test_that("Poisson critical value, lower.tail=FALSE", test_critical_discrete(pcrit1B, dec_countB, posterior_poisson))
