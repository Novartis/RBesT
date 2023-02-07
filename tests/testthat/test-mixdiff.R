context("mixdiff: Mixture Difference Distribution")

## test that calculations for the cumulative distribution function of
## differences in mixture distributions are correct by comparison with
## sampling

set.seed(23534)

## precision at which reference and tested must match
eps <- 1e-2

## number of samples used for sampling method
Nsamp <- 1e6

probs <- seq(0.1, 0.9, by=0.1)

## define the different test cases
beta <- list(mix1=mixbeta(c(1, 11, 4)), mix2=mixbeta(c(1, 8, 7)), q=c(0, 0.3), p=probs)
betaMix <- list(mix1=mixbeta(c(0.8, 11, 4), c(0.2, 1, 1)), mix2=mixbeta(c(0.8, 8, 7), c(0.2, 1, 1)), q=c(0, 0.3), p=probs)

gamma <- list(mix1=mixgamma(c(1, 20, 4)), mix2=mixgamma(c(1, 50, 10)), q=c(0, -2), p=probs)
gammaMix <- list(mix1=mixgamma(rob=c(0.75, 8, 0.5), inf=c(0.25, 9, 2), param="ms"), mix2=mixgamma(c(1, 50, 10)), q=c(0, -2), p=probs)

nm <- mixnorm(rob=c(0.2, 0, 2), inf=c(0.8, 2, 2), sigma=5)

norm <- list(mix1=mixnorm(c(1, 10, sqrt(1/4))), mix2=mixnorm(c(1, 11, sqrt(1/4))), q=c(0, 1.5), p=probs )
norm_ref <- mixnorm(c(1,10-11,sqrt(1/4 + 1/4)))

normMix <- list(mix1=mixnorm(c(0.2, 0, 2), c(0.8, 2, 2)), mix2=mixnorm(c(1, 2, sqrt(4))), q=c(0, 1.5), p=probs )
normMix_ref <- mixnorm(c(0.2, 0-2, sqrt(4 + 4)), c(0.8, 2-2, sqrt(4+4)))

mixdiff_sample <- function(mix1, mix2, q, p, N=Nsamp) {
    samp <- rmix(mix1, N) - rmix(mix2, N)
    list(probs=vapply(q, function(v) mean(samp < v), c(p=0.1)), quants=quantile(samp, p))
}

mixdiff_cmp <- function(case, rev=FALSE) {
    ## skip for speed on CRAN
    skip_on_cran()
    ref_probs <- do.call(pmixdiff, case[c("mix1", "mix2", "q")])
    ref_quants <- do.call(qmixdiff, case[c("mix1", "mix2", "p")])
    test <- do.call(mixdiff_sample, case)
    res_probs <- abs(ref_probs - test$probs)
    res_quants <- abs(ref_quants - test$quants)
    expect_true(all(res_probs < eps))
    expect_true(all(res_quants < eps))
    ## also check the reversed difference case
    if(!rev) {
        case_rev <- case
        case_rev$mix1 <- case$mix2
        case_rev$mix2 <- case$mix1
        mixdiff_cmp(case_rev, TRUE)
    }
}

mixdiff_cmp_norm <- function(case, mixref) {
    test_probs <- do.call(pmixdiff, case[c("mix1", "mix2", "q")])
    test_quants <- do.call(qmixdiff, case[c("mix1", "mix2", "p")])
    ref_probs <- pmix(mixref, case$q)
    ref_quants <- qmix(mixref, case$p)
    res_probs <- abs(ref_probs - test_probs)
    res_quants <- abs(ref_quants - test_quants)
    expect_true(all(res_probs < eps))
    expect_true(all(res_quants < eps))
}

test_that("Difference in beta variates evaluates correctly", mixdiff_cmp(beta))
test_that("Difference in beta mixture variates evaluates correctly", mixdiff_cmp(betaMix))

test_that("Difference in gamma variates evaluates correctly", mixdiff_cmp(gamma))
test_that("Difference in gamma mixture variates evaluates correctly", mixdiff_cmp(gammaMix))

test_that("Difference in normal variates evaluates correctly", mixdiff_cmp(norm))
test_that("Difference in normal mixture variates evaluates correctly", mixdiff_cmp(normMix))

## for the normal we can use exact analytical results
test_that("Difference in normal variates evaluates (analytically) correctly", mixdiff_cmp_norm(norm,norm_ref))
test_that("Difference in normal mixture variates evaluates (analytically) correctly", mixdiff_cmp_norm(normMix,normMix_ref))

## now test difference distributions on the link-transformed scales
## (the cannonical cases)
apply_link <- function(dists, link) {
    RBesT:::dlink(dists$mix1) <- link
    RBesT:::dlink(dists$mix2) <- link
    dists
}
 
## log-odds
test_that("Difference in beta variates with logit link evaluates correctly", mixdiff_cmp(apply_link(beta, RBesT:::logit_dlink)))
test_that("Difference in beta mixture variates with logit link evaluates correctly", mixdiff_cmp(apply_link(betaMix, RBesT:::logit_dlink)))

## relative risk
test_that("Difference in beta variates with log link evaluates correctly", mixdiff_cmp(apply_link(beta, RBesT:::log_dlink)))
test_that("Difference in beta mixture variates with log link evaluates correctly", mixdiff_cmp(apply_link(betaMix, RBesT:::log_dlink)))

## relative counts
test_that("Difference in gamma variates with log link evaluates correctly", mixdiff_cmp(apply_link(gamma, RBesT:::log_dlink)))
test_that("Difference in gamma mixture variates with log link evaluates correctly", mixdiff_cmp(apply_link(gammaMix, RBesT:::log_dlink)))

