skip_if_not_installed("brms")

## various tests around mixture distributions
set.seed(234534)

if (getOption("brms.backend", "not_set") == "not_set") {
  .brms_backend <- Sys.getenv("BRMS_BACKEND", "not_set")
  if (.brms_backend != "not_set") {
    options(brms.backend = .brms_backend)
  }
}
if (getOption("cmdstanr_write_stan_file_dir", "not_set") == "not_set") {
  .brms_cache_dir <- Sys.getenv("BRMS_CACHE_DIR", "not_set")
  if (.brms_cache_dir != "not_set") {
    options(cmdstanr_write_stan_file_dir = .brms_cache_dir)
  }
}

## sample based match requirements
eps <- 1E-1

## define the different test cases (univariate)
beta <- mixbeta(c(1, 11, 4))
betaMix <- mixbeta(c(0.8, 11, 4), c(0.2, 1, 1))

gamma <- mixgamma(c(1, 5, 10), param = "mn")
gammaMix <- mixgamma(rob = c(0.25, 8, 0.5), inf = c(0.75, 8, 10), param = "mn")

norm <- mixnorm(c(1, 0, sqrt(2)), sigma = 1)

normMix <- mixnorm(c(0.2, 0, 2), c(0.8, 1, 2), sigma = 1)
normMixWeak <- mixnorm(c(0.2, 0, 2), c(0.8, 1, 2), c(0, 0, 1), sigma = 1)

## tests the quantile and distribution function against simulated
## samples when using brms prior sampling as reference
mixstanvar_simul_test <- function(
  mix,
  brms_args,
  eps,
  qtest,
  ptest = seq(0.2, 0.8, by = 0.2)
) {
  skip_on_cran()
  skip_on_ci()
  capture.output(
    brms_prior <- do.call(
      brms::brm,
      c(
        brms_args,
        list(
          seed = 1423545,
          refresh = 0,
          sample_prior = "only",
          stanvars = mixstanvar(prior = mix)
        )
      )
    )
  )
  samp <- as.numeric(brms::as_draws_matrix(
    brms_prior,
    variable = "b_Intercept"
  )[, 1])
  qtest_samp <- quantile(samp, ptest)
  qref_qmix <- qmix(mix, ptest)
  res_quants <- abs(qref_qmix - qtest_samp)
  expect_true(all(res_quants < eps))
  ptest_samp <- vapply(qtest, function(q) mean(samp < q), c(0.1))
  pref_pmix <- pmix(mix, qtest)
  res_probs <- abs(pref_pmix - ptest_samp)
  expect_true(all(res_probs < eps))
}


mixstanvar_test <- function(mix, brms_args) {
  brms_prior_empty <- do.call(
    brms::brm,
    c(
      brms_args,
      list(
        seed = 1423545,
        refresh = 0,
        sample_prior = "only",
        stanvars = mixstanvar(prior = mix),
        empty = TRUE
      )
    )
  )
  mix_class <- gsub("Mix$", "", class(mix)[1])
  stan_dist_lpdf <- paste0("mix", mix_class, "_lpdf")
  stan_dist_lcdf <- paste0("mix", mix_class, "_lcdf")
  stan_dist_lccdf <- paste0("mix", mix_class, "_lccdf")
  stan_dist_cdf <- paste0("mix", mix_class, "_cdf")
  ## look for the declared density in Stan
  stan_code <- brms::stancode(brms_prior_empty)
  expect_true(
    grep(stan_dist_lpdf, stan_code) == 1,
    info = "Looking for declared Stan mixture density pdf in generated brms Stan code."
  )
  expect_true(
    grep(stan_dist_lcdf, stan_code) == 1,
    info = "Looking for declared Stan mixture density cdf in generated brms Stan code."
  )
  expect_true(
    grep(stan_dist_lccdf, stan_code) == 1,
    info = "Looking for declared Stan mixture density ccdf in generated brms Stan code."
  )
  expect_true(
    grep(stan_dist_cdf, stan_code) == 1,
    info = "Looking for declared Stan mixture density natural scale cdf in generated brms Stan code."
  )
  ## now check for the mixture being passed to Stan as data
  stan_data <- brms::standata(brms_prior_empty)
  for (i in 1:3) {
    param <- paste0("prior_", rownames(mix)[i])
    expect_true(all(unname(mix[i, ]) == unname(stan_data[[param]])))
  }
}

brms_beta_args <- list(
  formula = brms::bf(
    r | trials(n) ~ 1,
    family = brms::brmsfamily("binomial", link = "identity"),
    center = FALSE
  ),
  data = data.frame(r = 0, n = 0),
  prior = brms::prior(mixbeta(prior_w, prior_a, prior_b), coef = Intercept)
)

test_that("Beta quantiles are correct for brms sampled prior", {
  mixstanvar_simul_test(beta, brms_beta_args, eps, c(0.1, 0.9))
})
test_that("Beta prior is declared correctly in brms generated model and data", {
  mixstanvar_test(beta, brms_beta_args)
})
test_that("Beta mixture quantiles are correct for brms sampled prior", {
  mixstanvar_simul_test(betaMix, brms_beta_args, eps, c(0.1, 0.9))
})
test_that("Beta mixture prior is declared correctly in brms generated model and data", {
  mixstanvar_test(betaMix, brms_beta_args)
})

brms_beta_trunc_args <- list(
  formula = brms::bf(
    r | trials(n) ~ 1,
    family = brms::brmsfamily("binomial", link = "identity"),
    center = FALSE
  ),
  data = data.frame(r = 0, n = 0),
  prior = brms::prior(
    mixbeta(prior_w, prior_a, prior_b),
    class = b,
    lb = 0.1,
    ub = 0.9
  )
)
test_that("Beta (truncated) prior is declared correctly in brms generated model and data", {
  mixstanvar_test(beta, brms_beta_trunc_args)
})
test_that("Beta mixture (truncated) prior is declared correctly in brms generated model and data", {
  mixstanvar_test(betaMix, brms_beta_trunc_args)
})

brms_normal_args <- list(
  formula = brms::bf(
    y ~ 1,
    family = brms::brmsfamily("gaussian", link = "identity"),
    center = FALSE
  ),
  data = data.frame(y = 0),
  prior = brms::prior(mixnorm(prior_w, prior_m, prior_s), coef = Intercept) +
    brms::prior(constant(1), class = sigma)
)
test_that("Normal quantiles are correct for brms sampled prior", {
  mixstanvar_simul_test(norm, brms_normal_args, eps, c(-1, 0))
})
test_that("Normal prior is declared correctly in brms generated model and data", {
  mixstanvar_test(norm, brms_normal_args)
})
test_that("Normal mixture quantiles are correct for brms sampled prior", {
  mixstanvar_simul_test(
    normMix,
    brms_normal_args,
    eps,
    c(2, 1),
    ptest = c(0.3, 0.5, 0.7)
  )
})
test_that("Normal mixture prior is declared correctly in brms generated model and data", {
  mixstanvar_test(normMix, brms_normal_args)
})

brms_normal_trunc_args <- list(
  formula = brms::bf(
    y ~ 1,
    family = brms::brmsfamily("gaussian", link = "identity"),
    center = FALSE
  ),
  data = data.frame(y = 0),
  prior = brms::prior(
    mixnorm(prior_w, prior_m, prior_s),
    class = b,
    lb = -5,
    ub = 5
  ) +
    brms::prior(constant(1), class = sigma)
)
test_that("Normal (truncated) prior is declared correctly in brms generated model and data", {
  mixstanvar_test(norm, brms_normal_trunc_args)
})
test_that("Normal mixture (truncated) prior is declared correctly in brms generated model and data", {
  mixstanvar_test(normMix, brms_normal_trunc_args)
})

brms_gamma_args <- list(
  formula = brms::bf(
    y ~ 1,
    family = brms::brmsfamily("gaussian", link = "identity"),
    center = FALSE
  ),
  data = data.frame(y = 1),
  prior = brms::prior(mixgamma(prior_w, prior_a, prior_b), coef = Intercept) +
    brms::prior(constant(1), class = sigma)
)

test_that("Gamma quantiles are correct for brms sampled prior", {
  mixstanvar_simul_test(gamma, brms_gamma_args, eps, c(2, 7))
})
test_that("Gamma prior is declared correctly in brms generated model and data", {
  mixstanvar_test(gamma, brms_gamma_args)
})
test_that("Gamma mixture quantile function is correct for brms sampled prior", {
  mixstanvar_simul_test(
    gammaMix,
    brms_gamma_args,
    eps,
    c(2, 7),
    ptest = seq(0.2, 0.8, by = 0.2)
  )
})
test_that("Gamma mixture prior is declared correctly in brms generated model and data", {
  mixstanvar_test(gammaMix, brms_gamma_args)
})

brms_gamma_trunc_args <- list(
  formula = brms::bf(
    y ~ 1,
    family = brms::brmsfamily("gaussian", link = "identity"),
    center = FALSE
  ),
  data = data.frame(y = 1),
  prior = brms::prior(
    mixgamma(prior_w, prior_a, prior_b),
    class = b,
    lb = 0.1,
    ub = 10
  ) +
    brms::prior(constant(1), class = sigma)
)

test_that("Gamma (truncated) prior is declared correctly in brms generated model and data", {
  mixstanvar_test(gamma, brms_gamma_trunc_args)
})
test_that("Gamma mixture (truncated) prior is declared correctly in brms generated model and data", {
  mixstanvar_test(gammaMix, brms_gamma_trunc_args)
})

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
  0.5 *
    (t(m_delta) %*%
      inv_sigma_2 %*%
      m_delta +
      sum(diag(inv_sigma_2 %*% sigma_1)) -
      log(det(sigma_1)) +
      log(det(sigma_2)) -
      p)
}

mixstanvar_simul_mv_test <- function(mvmix, brms_args, eps) {
  skip_on_cran()
  skip_on_ci()
  capture.output(
    brms_prior <- do.call(
      brms::brm,
      c(
        brms_args,
        list(
          seed = 1423545,
          refresh = 0,
          sample_prior = "only",
          stanvars = mixstanvar(prior = mvmix)
        )
      )
    )
  )
  samp <- brms::as_draws_matrix(brms_prior, variable = "^b_", regex = TRUE)
  samp_m <- colMeans(samp)
  samp_sigma <- cov(samp)
  mix_m <- summary(mvmix)$mean
  mix_sigma <- summary(mvmix)$cov
  kl <- KLdiv_mvnorm(samp_m, samp_sigma, mix_m, mix_sigma)
  expect_true(abs(kl) < eps)
}

mixstanvar_test_mvnormMix <- function(mix, brms_args) {
  brms_prior_empty <- do.call(
    brms::brm,
    c(
      brms_args,
      list(
        seed = 1423545,
        refresh = 0,
        sample_prior = "only",
        stanvars = mixstanvar(prior = mix),
        empty = TRUE
      )
    )
  )
  stan_dist <- paste0("mix", gsub("Mix$", "", class(mix)[1]), "_lpdf")
  ## look for the declared density in Stan
  expect_true(
    grep(stan_dist, brms::stancode(brms_prior_empty)) == 1,
    info = "Looking for declared Stan mixture density in generated brms Stan code."
  )
  ## now check for the mixture being passed to Stan as data
  stan_data <- brms::standata(brms_prior_empty)
  ## number of mixture components
  Nc <- ncol(mix)
  expect_equal(stan_data$prior_Nc, Nc)
  ## dimensionality
  p <- length(summary(mix)$mean)
  expect_equal(stan_data$prior_p, p)
  ## weights per component
  expect_equal(stan_data$prior_w, array(unname(mix[1, ])))
  ## means per component
  expect_equal(unname(stan_data$prior_m), unname(t(mix[2:(p + 1), ])))
  ## covariance information
  for (i in seq_len(Nc)) {
    S_c <- matrix(stan_data$prior_sigma[i, , , drop = FALSE], p, p)
    expect_equal(sqrt(diag(S_c)), unname(mix[(p + 2):(1 + 2 * p), i]))
    Rho_c <- cov2cor(S_c)
    expect_equal(
      Rho_c[lower.tri(Rho_c)],
      unname(mix[(1 + 2 * p + 1):nrow(mix), i])
    )
  }
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

mvnorm_single_4 <- mixmvnorm(c(1, m1, 5), param = "mn", sigma = S)

mvnorm_heavy_4 <- mixmvnorm(
  c(0.5, m1, 0.25),
  c(0.5, m2, 5),
  param = "mn",
  sigma = S
)


brms_mvn_4_args <- list(
  formula = brms::bf(
    y ~ 1 + l1 + l2 + l3,
    family = brms::brmsfamily("gaussian", link = "identity"),
    center = FALSE
  ),
  data = data.frame(y = 1, l1 = 0, l2 = 0, l3 = 0),
  prior = brms::prior(mixmvnorm(prior_w, prior_m, prior_sigma_L), class = b) +
    brms::prior(constant(1), class = sigma)
)

test_that("Multivariate normal (4D) is correct for brms sampled prior", {
  mixstanvar_simul_mv_test(mvnorm_single_4, brms_mvn_4_args, eps)
})
test_that("Multivariate normal (4D) prior is declared correctly in brms generated model and data", {
  mixstanvar_test_mvnormMix(mvnorm_single_4, brms_mvn_4_args)
})
test_that("Multivariate normal with heavy (4D) tails is correct for brms sampled prior", {
  mixstanvar_simul_mv_test(mvnorm_heavy_4, brms_mvn_4_args, eps)
})
test_that("Multivariate normal with heavy (4D) tails is declared correctly in brms generated model and data", {
  mixstanvar_test_mvnormMix(mvnorm_heavy_4, brms_mvn_4_args)
})

mvnorm_single_2 <- mixmvnorm(
  c(1, m1[1:2], 5),
  param = "mn",
  sigma = S[1:2, 1:2]
)


mvnorm_heavy_2 <- mixmvnorm(
  c(0.5, m1[1:2], 0.25),
  c(0.5, m2[1:2], 5),
  param = "mn",
  sigma = S[1:2, 1:2]
)

brms_mvn_2_args <- list(
  formula = brms::bf(
    y ~ 1 + l1,
    family = brms::brmsfamily("gaussian", link = "identity"),
    center = FALSE
  ),
  data = data.frame(y = 1, l1 = 0),
  prior = brms::prior(mixmvnorm(prior_w, prior_m, prior_sigma_L), class = b) +
    brms::prior(constant(1), class = sigma)
)

test_that("Multivariate normal (2D) is correct for brms sampled prior", {
  mixstanvar_simul_mv_test(mvnorm_single_2, brms_mvn_2_args, eps)
})
test_that("Multivariate normal (2D) is declared correctly in brms generated model and data", {
  mixstanvar_test_mvnormMix(mvnorm_single_2, brms_mvn_2_args)
})
test_that("Multivariate normal with heavy (2D) tails is correct for brms sampled prior", {
  mixstanvar_simul_mv_test(mvnorm_heavy_2, brms_mvn_2_args, eps)
})
test_that("Multivariate normal with heavy (2D) is declared correctly in brms generated model and data", {
  mixstanvar_test_mvnormMix(mvnorm_heavy_2, brms_mvn_2_args)
})
