context("preddist: Mixture Predictive Distribution")

## check that predictive distributions hold what they promise,
## i.e. that they describe the sum of n new data points.

set.seed(42343)

## precision at which reference and tested must match
eps <- 1e-2

## percentiles to test for
p_quants <- seq(0.1, 0.9, by=0.1)

## number of samples used for sampling method
Nsamp <- 1e5

## define the different test cases
beta <- mixbeta(c(1, 11, 4))
betaMix <- mixbeta(c(0.8, 11, 4), c(0.2, 1, 1))

gamma <- mixgamma(c(1, 65, 50))
gammaMix <- mixgamma(rob=c(0.5, 8, 0.5), inf=c(0.5, 9, 2), param="ms")

norm <- mixnorm(c(1, 0, 0.5), sigma=1)
normMix <- mixnorm(c(0.2, 0, 2), c(0.8, 2, 1), sigma=1)

n <- 25

preddist_cmp <- function(mix,  n, n_rng, N=Nsamp, qntls=p_quants, stat=c("sum", "mean"), Teps=eps) {
    skip_on_cran()

    ## sample for each draw a single hyper-parameter which is then
    ## used n times in the rng function to return n samples from the
    ## sampling distribution
    stat <- match.arg(stat)
    test <- replicate(N, n_rng(rmix(mix, 1)))
    if(stat=="sum")  test_stat <- colSums( test)
    if(stat=="mean") test_stat <- colMeans(test)

    quants_stest <- quantile(test_stat, qntls)
    ## note: in particular for the discrete/counting distributions,
    ## sampling gives better estimates
    quants_sref <- qmix(preddist(mix, n=n), qntls)
    res_sum <- abs(quants_sref-quants_stest)
    ## note: errors are accumulating with n, hence to check the mean,
    ## we scale eps with n
    if(stat=="sum")  expect_true(all(res_sum < n * Teps))
    if(stat=="mean") expect_true(all(res_sum < Teps))

    quants_test <- quantile(test[1,], qntls)
    quants_ref <- qmix(preddist(mix, n=1), qntls)
    res <- abs(quants_ref-quants_test)
    expect_true(all(res < Teps))

    if(inherits(mix, "betaMix")) {
        ## specifically test BetaBinomial (which is from Stan)
        predmix <- preddist(mix, n=30)
        dens <- dmix(predmix, 0:30)
        cdens <- cumsum(dens)
        pc <- pmix(predmix, 0:30)
        upc <- pmix(predmix, 0:30, FALSE)
        expect_true(all(abs(cdens - pc) < Teps))
        expect_true(all(abs(1-cdens - upc) < Teps))
    }

    ## test that output length of the predictive is the same as the
    ## input data vector, even for negative values for values outside
    ## the valid support of the postive only distribtions
    cprob_ltest <- pmix(preddist(mix, n=1), c(-1, quants_ref))
    expect_true(length(cprob_ltest) == length(c(-1, quants_ref)))
}


test_that("Predictive for a beta evaluates correctly (binary)", preddist_cmp(beta, n, Curry(rbinom, n=n, size=1)))
test_that("Predictive for a beta mixture evaluates correctly (binary)", preddist_cmp(betaMix, n, Curry(rbinom, n=n, size=1)))

test_that("Predictive for a gamma evaluates correctly (poisson)", preddist_cmp(gamma, n, Curry(rpois, n=n)))
test_that("Predictive for a gamma mixture evaluates correctly (poisson)", preddist_cmp(gammaMix, n, Curry(rpois, n=n), Teps=1E-1))

test_that("Predictive for a normal evaluates correctly (normal)", preddist_cmp(norm, n, Curry(rnorm, n=n, sd=sigma(norm)), stat="mean"))
test_that("Predictive for a normal mixture evaluates correctly (normal)", preddist_cmp(normMix, n, Curry(rnorm, n=n, sd=sigma(normMix)), stat="mean", Teps=1E-1))
