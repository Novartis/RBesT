
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

test_that("Cumulative beta distribution function evaluates lower.tail correctly", { pmix_lower_tail_test(beta) })
test_that("Cumulative beta mixture distribution function evaluates lower.tail correctly", { pmix_lower_tail_test(betaMix) })

test_that("Cumulative normal distribution function evaluates lower.tail correctly", { pmix_lower_tail_test(norm) })
test_that("Cumulative normal mixture distribution function evaluates lower.tail correctly", { pmix_lower_tail_test(normMix) })

test_that("Cumulative gamma distribution function evaluates lower.tail correctly", { pmix_lower_tail_test(gamma) })
test_that("Cumulative gamma mixture distribution function evaluates lower.tail correctly", { pmix_lower_tail_test(gammaMix) })

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

test_that("Beta quantile function is correct", { mix_simul_test(beta, eps, c(0.1, 0.9)) })
test_that("Beta mixture quantile function is correct", { mix_simul_test(betaMix, eps, c(0.1, 0.9)) })

test_that("Normal quantile function is correct", { mix_simul_test(norm, eps, c(-1, 0)) })
test_that("Normal mixture quantile function is correct", { mix_simul_test(normMix, eps, c(4, 1)) })
test_that("Normal mixture with very weak component quantile function is correct", { mix_simul_test(normMixWeak, eps, c(4, 1)) })

test_that("Gamma quantile function is correct", { mix_simul_test(gamma, eps, c(2, 7)) })
test_that("Gamma mixture quantile function is correct", { mix_simul_test(gammaMix, eps, c(2, 7), ptest = seq(0.2, 0.8, by=0.1)) })

## problematic gamma (triggers internally a fallback to root finding)
gammaMix2 <- mixgamma(c(8.949227e-01, 7.051570e-01, 6.125121e-02),
                      c(1.049106e-01, 3.009986e-01, 5.169626e-04),
                      c(1.666667e-04, 1.836051e+04, 1.044005e-02))

test_that("Singular gamma mixture quantile function is correct", { mix_simul_test(gammaMix2, 10*eps, c(1, 1E3), ptest = seq(0.2, 0.8, by=0.1)) })


consistent_cdf <- function(mix, values) {
    dens <- dmix(mix, values)
    cdf <- pmix(mix, values)
    lcdf <- pmix(mix, values, log.p=TRUE)
    expect_true(all(diff(cdf) >= 0))
    expect_numeric(dens, lower=0, finite=TRUE, any.missing=FALSE)
    expect_numeric(cdf, lower=0, upper=1, finite=TRUE, any.missing=FALSE)
    expect_numeric(lcdf, upper=0, finite=FALSE, any.missing=FALSE)
}

consistent_ccdf <- function(mix, values) {
    dens <- dmix(mix, values)
    ccdf <- pmix(mix, values, FALSE)
    lccdf <- pmix(mix, values, FALSE, TRUE)
    expect_true(all(diff(ccdf) <= 0))
    expect_numeric(dens, lower=0, finite=TRUE, any.missing=FALSE)
    expect_numeric(ccdf, lower=0, upper=1, finite=TRUE, any.missing=FALSE)
    expect_numeric(lccdf, upper=0, finite=FALSE, any.missing=FALSE)
}

test_that("Beta CDF function is consistent", { consistent_cdf(beta, seq(0.1, 0.9, by=0.1)) })
test_that("Beta mixture CDF function is consistent", { consistent_cdf(betaMix, seq(0.1, 0.9, by=0.1)) })

test_that("Normal CDF is consistent", { consistent_cdf(norm, seq(-2, 2, by=0.1)) })
test_that("Normal mixture CDF is consistent", { consistent_cdf(norm, seq(-2, 2, by=0.1)) })

test_that("Gamma CDF function is consistent", { consistent_cdf(gamma, seq(2, 7, by=0.1)) })
test_that("Gamma mixture CDF function is consistent", { consistent_cdf(gammaMix, seq(2, 7, by=0.1)) })


## problematic beta which triggers that the cumulative of the
## predictive is not monotone (probably fixed with Stan 2.18, check
## again once 2.18 is out)

## problematic Beta density
bm1 <- mixbeta(c(1.0, 298.30333970, 146.75306521))
test_that("Problematic (1) BetaBinomial CDF function is consistent", { consistent_cdf(preddist(bm1, n=50), 0:50) })
test_that("Problematic (1) BetaBinomial CCDF function is consistent", { consistent_ccdf(preddist(bm1, n=50), 0:50) })
bm2 <- mixbeta(c(1.0, 3 + 1/3, 47 + 1/3))
test_that("Problematic (2) BetaBinomial CDF function is consistent", { consistent_cdf(preddist(bm2, n=50), 0:50) })
test_that("Problematic (2) BetaBinomial CCDF function is consistent", { consistent_ccdf(preddist(bm2, n=50), 0:50) })

## tests for the multivariate normal mixture density
p <- 4
Rho <- diag(p)
Rho[lower.tri(Rho)] <- c(0.3, 0.8, -0.2, 0.1, 0.5, -0.4)
Rho[upper.tri(Rho)] <- t(Rho)[upper.tri(Rho)]
s <- c(1, 2, 3, 4)
S <- diag(s, p) %*% Rho %*% diag(s, p)
rownames(S) <- colnames(S) <- 1:p

mvn_consistent_dimension <- function(mix, p) {
    s <- summary(mix)
    expect_numeric(s$mean, any.missing=FALSE, len=p)
    expect_matrix(s$cov, any.missing=FALSE, nrows=p, ncols=p)
}

test_that("Multivariate normal mixture has consistent dimensionality",
{
    for(i in 1:(nrow(S)-1)) {
        p_sub <- 4-i
        S_sub <- S[-c(1:i), -c(1:i), drop=FALSE]
        mvn_consistent_dimension(mixmvnorm(c(1, rep(0, p_sub), S_sub), sigma=S_sub), p_sub)
    }
})

test_that("Multivariate normal mixture has consistent dimension naming",
{
    for(i in 1:(nrow(S)-1)) {
        p_sub <- 4-i
        S_sub <- S[-c(1:i), -c(1:i), drop=FALSE]
        m_sub <- rep(0, p_sub)
        dim_labels <- letters[1:p_sub]
        names(m_sub) <- dim_labels
        test_mix <- mixmvnorm(c(1, m_sub, S_sub), sigma=S_sub)
        ## now test that names are used consistently
        expect_equal(rownames(sigma(test_mix)), dim_labels)
        expect_equal(colnames(sigma(test_mix)), dim_labels)
        expect_equal(names(summary(test_mix)$mean), dim_labels)
        expect_equal(rownames(summary(test_mix)$cov), dim_labels)
        expect_equal(colnames(summary(test_mix)$cov), dim_labels)
        expect_equal(colnames(rmix(test_mix, 1)), dim_labels)
    }
})

test_that("Multivariate normal mixture has consistent initialization",
{
    p <- nrow(S)
    mv1 <- mixmvnorm(c(1, rep(0, p), S), sigma=S, param="ms")
    mv2 <- mixmvnorm(c(1, rep(0, p), 1), sigma=S, param="mn")
    mv3 <- mixmvnorm(c(1, rep(0, p), 2), sigma=S, param="mn")

    expect_equal(summary(mv1)$cov, S, tolerance=eps_lower)
    expect_equal(summary(mv2)$cov, S, tolerance=eps_lower)
    expect_equal(summary(mv3)$cov, S/2, tolerance=eps_lower)
})

mvn_consistent_summaries <- function(mix, S=Nsamp_equant) {
    samp <- rmix(mix, S)
    m <- colMeans(samp)
    expect_equal(colMeans(samp), summary(mix)$mean, tolerance=eps)
    expect_equal(cov(samp), summary(mix)$cov, tolerance=eps)
}

test_that("Multivariate normal mixture has consistent summaries",
{
    p <- nrow(S)
    mv1 <- mixmvnorm(c(1, rep(0, p), S), sigma=S, param="ms")
    mv2 <- mixmvnorm(c(1, rep(0, p), 1), sigma=S, param="mn")
    mv3 <- mixmvnorm(c(0.2, rep(0, p), 2), c(0.8, rep(1, p), 6), sigma=S, param="mn")

    mvn_consistent_summaries(mv1)
    mvn_consistent_summaries(mv2)
    mvn_consistent_summaries(mv3)
})
