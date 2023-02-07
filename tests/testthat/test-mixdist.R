context("mixdist: Mixture Distribution")

## various tests around mixture distributions
set.seed(234534)

## precision at which reference and tested must match
eps_lower <- 1e-5

## sample based match requirements
eps <- 1E-2

## number of samples used for sampling method
Nsamp_quant <- 3

## number of samples used to compare sampled vs estimated quantiles
Nsamp_equant <- 1E6

## define the different test cases
beta <- mixbeta(c(1, 11, 4))
betaMix <- mixbeta(c(0.8, 11, 4), c(0.2, 1, 1))

gamma <- mixgamma(c(1, 5, 10), param="mn")
gammaMix <- mixgamma(rob=c(0.25, 8, 0.5), inf=c(0.75, 8, 10), param="mn")

nm <- mixnorm(rob=c(0.2, 0, 2), inf=c(0.8, 2, 2), sigma=5)

norm <- mixnorm(c(1, 0, sqrt(2)), sigma=1)

normMix <- mixnorm(c(0.2, 0, 2), c(0.8, 2, 2), sigma=1)
normMixWeak <- mixnorm(c(0.2, 0, 2), c(0.8, 2, 2), c(0, 0, 1), sigma=1)

pmix_lower_tail_test <- function(mix, N=Nsamp_quant) {
    ## sample some random quantiles
    do_test <- function(mix) {
        q <- rmix(mix, N)
        pl <- pmix(mix, q, lower.tail=TRUE)
        pu <- pmix(mix, q, lower.tail=FALSE)
        res <- abs(pl  - (1-pu) )
        expect_true(all(res < eps_lower))
    }
    ## now also test the respective predictive
    do_test(mix)
    do_test(preddist(mix, n=100))
}

test_that("Cumulative beta distribution function evaluates lower.tail correctly", pmix_lower_tail_test(beta))
test_that("Cumulative beta mixture distribution function evaluates lower.tail correctly", pmix_lower_tail_test(betaMix))

test_that("Cumulative normal distribution function evaluates lower.tail correctly", pmix_lower_tail_test(norm))
test_that("Cumulative normal mixture distribution function evaluates lower.tail correctly", pmix_lower_tail_test(normMix))

test_that("Cumulative gamma distribution function evaluates lower.tail correctly", pmix_lower_tail_test(gamma))
test_that("Cumulative gamma mixture distribution function evaluates lower.tail correctly", pmix_lower_tail_test(gammaMix))

## tests the quantile and distribution function against simulated samples
mix_simul_test <- function(mix, eps, qtest, ptest = seq(0.1, 0.9, by=0.1), S=Nsamp_equant) {
    samp <- rmix(mix, S)
    qtest_samp <- quantile(samp, ptest)
    qref_qmix <- qmix(mix, ptest)
    res_quants <- abs(qref_qmix - qtest_samp)
    expect_true(all(res_quants < eps))
    ptest_samp <- vapply(qtest, function(q) mean(samp < q), c(0.1))
    pref_pmix <- pmix(mix, qtest)
    res_probs <- abs(pref_pmix - ptest_samp)
    expect_true(all(res_probs < eps))
}

test_that("Beta quantile function is correct", mix_simul_test(beta, eps, c(0.1, 0.9)))
test_that("Beta mixture quantile function is correct", mix_simul_test(betaMix, eps, c(0.1, 0.9)))

test_that("Normal quantile function is correct", mix_simul_test(norm, eps, c(-1, 0)))
test_that("Normal mixture quantile function is correct", mix_simul_test(normMix, eps, c(4, 1)))
test_that("Normal mixture with very weak component quantile function is correct", mix_simul_test(normMixWeak, eps, c(4, 1)))

test_that("Gamma quantile function is correct", mix_simul_test(gamma, eps, c(2, 7)))
test_that("Gamma mixture quantile function is correct", mix_simul_test(gammaMix, eps, c(2, 7), ptest = seq(0.2, 0.8, by=0.1)))

## problematic gamma (triggers internally a fallback to root finding)
gammaMix2 <- mixgamma(c(8.949227e-01, 7.051570e-01, 6.125121e-02),
                      c(1.049106e-01, 3.009986e-01, 5.169626e-04),
                      c(1.666667e-04, 1.836051e+04, 1.044005e-02))

test_that("Singular gamma mixture quantile function is correct", mix_simul_test(gammaMix2, 10*eps, c(1, 1E3), ptest = seq(0.2, 0.8, by=0.1)))


consistent_cdf <- function(mix, values) {
    dens <- dmix(mix, values)
    cdf <- pmix(mix, values)
    expect_true(all(diff(cdf) >= 0))
    expect_numeric(dens, any.missing=FALSE)
    expect_numeric(cdf, any.missing=FALSE)
}

test_that("Beta CDF function is consistent", consistent_cdf(beta, seq(0.1, 0.9, by=0.1)))
test_that("Beta mixture CDF function is consistent", consistent_cdf(betaMix, seq(0.1, 0.9, by=0.1)))

test_that("Normal CDF is consistent", consistent_cdf(norm, seq(-2, 2, by=0.1)))
test_that("Normal mixture CDF is consistent", consistent_cdf(norm, seq(-2, 2, by=0.1)))

test_that("Gamma CDF function is consistent", consistent_cdf(gamma, seq(2, 7, by=0.1)))
test_that("Gamma mixture CDF function is consistent", consistent_cdf(gammaMix, seq(2, 7, by=0.1)))


## problematic beta which triggers that the cumulative of the
## predictive is not monotone (probably fixed with Stan 2.18, check
## again once 2.18 is out)

## problematic Beta density
bm <- mixbeta(c(1.0, 298.30333970, 146.75306521))
test_that("Problematic (1) BetaBinomial CDF function is consistent", consistent_cdf(preddist(bm, n=50), 0:50))
