context("mixstanvar: Mixture Distribution brms Adapter")

stopifnot(requireNamespace("brms", quietly=TRUE))

## various tests around mixture distributions
set.seed(234534)

if(getOption("brms.backend", "not_set") == "not_set") {
    .brms_backend <- Sys.getenv("BRMS_BACKEND", "not_set")
    if(.brms_backend != "not_set") {
        options(brms.backend=.brms_backend)
    }
}
if(getOption("cmdstanr_write_stan_file_dir", "not_set") == "not_set") {
    .brms_cache_dir <- Sys.getenv("BRMS_CACHE_DIR", "not_set")
    if(.brms_cache_dir != "not_set") {
        options(cmdstanr_write_stan_file_dir=.brms_cache_dir)
    }
}

## sample based match requirements
eps <- 1E-1

## define the different test cases (univariate)
beta <- mixbeta(c(1, 11, 4))
betaMix <- mixbeta(c(0.8, 11, 4), c(0.2, 1, 1))

gamma <- mixgamma(c(1, 5, 10), param="mn")
gammaMix <- mixgamma(rob=c(0.25, 8, 0.5), inf=c(0.75, 8, 10), param="mn")

norm <- mixnorm(c(1, 0, sqrt(2)), sigma=1)

normMix <- mixnorm(c(0.2, 0, 2), c(0.8, 1, 2), sigma=1)
normMixWeak <- mixnorm(c(0.2, 0, 2), c(0.8, 1, 2), c(0, 0, 1), sigma=1)

## tests the quantile and distribution function against simulated
## samples when using brms prior sampling as reference
mixstanvar_simul_test <- function(mix,
                                  brms_args,
                                  eps, qtest, ptest = seq(0.2, 0.8, by=0.2)) {
    skip_on_cran()
    capture.output(brms_prior <- do.call(brms::brm, c(brms_args, list(seed=1423545, refresh=0, sample_prior="only", stanvars=mixstanvar(prior=mix)))))
    samp <- as.numeric(brms::as_draws_matrix(brms_prior, variable="b_Intercept")[,1])
    qtest_samp <- quantile(samp, ptest)
    qref_qmix <- qmix(mix, ptest)
    res_quants <- abs(qref_qmix - qtest_samp)
    expect_true(all(res_quants < eps))
    ptest_samp <- vapply(qtest, function(q) mean(samp < q), c(0.1))
    pref_pmix <- pmix(mix, qtest)
    res_probs <- abs(pref_pmix - ptest_samp)
    expect_true(all(res_probs < eps))
}

brms_beta_args <- list(formula=brms::bf(r | trials(n) ~ 1, family=brms::brmsfamily("binomial", link="identity"), center=FALSE),
                       data=data.frame(r=0, n=0),
                       prior=brms::prior(mixbeta(prior_w, prior_a, prior_b), coef=Intercept))

test_that("Beta quantiles are correct for brms sampled prior", mixstanvar_simul_test(beta, brms_beta_args, eps, c(0.1, 0.9)))
test_that("Beta mixture quantiles are correct for brms sampled prior", mixstanvar_simul_test(betaMix, brms_beta_args, eps, c(0.1, 0.9)))

brms_normal_args <- list(formula=brms::bf(y ~ 1, family=brms::brmsfamily("gaussian", link="identity"), center=FALSE),
                         data=data.frame(y=0),
                         prior=brms::prior(mixnorm(prior_w, prior_m, prior_s), coef=Intercept) + brms::prior(constant(1), class=sigma))
test_that("Normal quantiles are correct for brms sampled prior", mixstanvar_simul_test(norm, brms_normal_args, eps, c(-1, 0)))
test_that("Normal mixture quantiles are correct for brms sampled prior", mixstanvar_simul_test(normMix, brms_normal_args, eps, c(2, 1), ptest=c(0.3, 0.5, 0.7)))

brms_gamma_args <- list(formula=brms::bf(y ~ 1, family=brms::brmsfamily("gaussian", link="identity"), center=FALSE),
                        data=data.frame(y=1),
                        prior=brms::prior(mixgamma(prior_w, prior_a, prior_b), coef=Intercept) + brms::prior(constant(1), class=sigma))

test_that("Gamma quantiles are correct for brms sampled prior", mixstanvar_simul_test(gamma, brms_gamma_args, eps, c(2, 7)))
test_that("Gamma mixture quantile function is correct for brms sampled prior", mixstanvar_simul_test(gammaMix, brms_gamma_args, eps, c(2, 7), ptest = seq(0.2, 0.8, by=0.2)))

# Here we approximate the samples using a multi-variante normal via a
# moment based approxmation and compare this to the respective
# approximation of the multi variate normal mixture. While this is not
# correct per se, it is sufficient for testing as these
# approximationas are unique and they do test for marginal means and
# correlations. See for details
# https://statproofbook.github.io/P/mvn-kl.html
KLdiv_mvnorm <- function(m_1, sigma_1, m_2, sigma_2) {
    m_delta <- (m_2 - m_1)
    inv_sigma_2 <- solve(sigma_2)
    p <- length(m_1)
    0.5 * ( t(m_delta) %*% inv_sigma_2 %*% m_delta + sum(diag(inv_sigma_2 %*% sigma_1)) - log(det(sigma_1)) + log(det(sigma_2)) - p )
}

mixstanvar_simul_mv_test <- function(mvmix, brms_args, eps) {
    skip_on_cran()
    capture.output(brms_prior <- do.call(brms::brm, c(brms_args, list(seed=1423545, refresh=0, sample_prior="only", stanvars=mixstanvar(prior=mvmix)))))
    samp <- brms::as_draws_matrix(brms_prior, variable="^b_", regex=TRUE)
    samp_m <- colMeans(samp)
    samp_sigma <- cov(samp)    
    mix_m <- summary(mvmix)$mean
    mix_sigma <- summary(mvmix)$cov
    kl <- KLdiv_mvnorm(samp_m, samp_sigma, mix_m, mix_sigma)
    expect_true(abs(kl) < eps)
}

p <- 4
Rho <- diag(p)
Rho[lower.tri(Rho)] <- c(0.3, 0.8, -0.2, 0.1, 0.5, -0.4)
Rho[upper.tri(Rho)] <- t(Rho)[upper.tri(Rho)]
s <- c(1, 2, 3, 4)
S <- diag(s, p) %*% Rho %*% diag(s, p)
zero <- rep(0, p)
m1 <- 0:3
m2 <- 1:4

mvnorm_single_4 <- mixmvnorm(c(1, m1, 5),
                           param="mn", sigma=S)

mvnorm_heavy_4 <- mixmvnorm(c(0.5, m1, 0.25),
                          c(0.5, m2, 5),
                          param="mn", sigma=S)

brms_mvn_4_args <- list(formula=brms::bf(y ~ 1 + l1 + l2 + l3, family=brms::brmsfamily("gaussian", link="identity"), center=FALSE),
                        data=data.frame(y=1, l1=0, l2=0, l3=0),
                        prior=brms::prior(mixmvnorm(prior_w, prior_m, prior_sigma_L), class=b) + brms::prior(constant(1), class=sigma))

test_that("Multivariate normal (4D) is correct for brms sampled prior", mixstanvar_simul_mv_test(mvnorm_single_4, brms_mvn_4_args, eps))
test_that("Multivariate normal with heavy (4D) tails is correct for brms sampled prior", mixstanvar_simul_mv_test(mvnorm_heavy_4, brms_mvn_4_args, eps))

mvnorm_single_2 <- mixmvnorm(c(1, m1[1:2], 5),
                             param="mn", sigma=S[1:2,1:2])


mvnorm_heavy_2 <- mixmvnorm(c(0.5, m1[1:2], 0.25),
                            c(0.5, m2[1:2], 5),
                            param="mn", sigma=S[1:2,1:2])

brms_mvn_2_args <- list(formula=brms::bf(y ~ 1 + l1, family=brms::brmsfamily("gaussian", link="identity"), center=FALSE),
                        data=data.frame(y=1, l1=0),
                        prior=brms::prior(mixmvnorm(prior_w, prior_m, prior_sigma_L), class=b) + brms::prior(constant(1), class=sigma))

test_that("Multivariate normal (2D) is correct for brms sampled prior", mixstanvar_simul_mv_test(mvnorm_single_2, brms_mvn_2_args, eps))
test_that("Multivariate normal with heavy (2D) tails is correct for brms sampled prior", mixstanvar_simul_mv_test(mvnorm_heavy_2, brms_mvn_2_args, eps))
