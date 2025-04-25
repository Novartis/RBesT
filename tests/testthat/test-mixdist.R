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

gamma <- mixgamma(c(1, 5, 10), param = "mn")
gammaMix <- mixgamma(rob = c(0.25, 8, 0.5), inf = c(0.75, 8, 10), param = "mn")

norm <- mixnorm(c(1, 0, sqrt(2)), sigma = 1)

normMix <- mixnorm(c(0.2, 0, 2), c(0.8, 2, 2), sigma = 1)
normMixWeak <- mixnorm(c(0.2, 0, 2), c(0.8, 2, 2), c(0, 0, 1), sigma = 1)

pmix_lower_tail_test <- function(mix, N = Nsamp_quant) {
  ## sample some random quantiles
  do_test <- function(mix) {
    q <- rmix(mix, N)
    pl <- pmix(mix, q, lower.tail = TRUE)
    pu <- pmix(mix, q, lower.tail = FALSE)
    res <- abs(pl - (1 - pu))
    expect_true(all(res < eps_lower))
  }
  ## now also test the respective predictive
  do_test(mix)
  do_test(preddist(mix, n = 100))
}

test_that("Cumulative beta distribution function evaluates lower.tail correctly", {
  pmix_lower_tail_test(beta)
})
test_that("Cumulative beta mixture distribution function evaluates lower.tail correctly", {
  pmix_lower_tail_test(betaMix)
})

test_that("Cumulative normal distribution function evaluates lower.tail correctly", {
  pmix_lower_tail_test(norm)
})
test_that("Cumulative normal mixture distribution function evaluates lower.tail correctly", {
  pmix_lower_tail_test(normMix)
})

test_that("Cumulative gamma distribution function evaluates lower.tail correctly", {
  pmix_lower_tail_test(gamma)
})
test_that("Cumulative gamma mixture distribution function evaluates lower.tail correctly", {
  pmix_lower_tail_test(gammaMix)
})

## tests the quantile and distribution function against simulated samples
mix_simul_test <- function(
  mix,
  eps,
  qtest,
  ptest = seq(0.1, 0.9, by = 0.1),
  S = Nsamp_equant
) {
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

test_that("Beta quantile function is correct", {
  mix_simul_test(beta, eps, c(0.1, 0.9))
})
test_that("Beta mixture quantile function is correct", {
  mix_simul_test(betaMix, eps, c(0.1, 0.9))
})

test_that("Normal quantile function is correct", {
  mix_simul_test(norm, eps, c(-1, 0))
})
test_that("Normal mixture quantile function is correct", {
  mix_simul_test(normMix, eps, c(4, 1))
})
test_that("Normal mixture with very weak component quantile function is correct", {
  mix_simul_test(normMixWeak, eps, c(4, 1))
})

test_that("Gamma quantile function is correct", {
  mix_simul_test(gamma, eps, c(2, 7))
})
test_that("Gamma mixture quantile function is correct", {
  mix_simul_test(gammaMix, eps, c(2, 7), ptest = seq(0.2, 0.8, by = 0.1))
})

## problematic gamma (triggers internally a fallback to root finding)
gammaMix2 <- mixgamma(
  c(8.949227e-01, 7.051570e-01, 6.125121e-02),
  c(1.049106e-01, 3.009986e-01, 5.169626e-04),
  c(1.666667e-04, 1.836051e+04, 1.044005e-02)
)

test_that("Singular gamma mixture quantile function is correct", {
  mix_simul_test(
    gammaMix2,
    10 * eps,
    c(1, 1E3),
    ptest = seq(0.2, 0.8, by = 0.1)
  )
})


consistent_cdf <- function(mix, values) {
  dens <- dmix(mix, values)
  cdf <- pmix(mix, values)
  lcdf <- pmix(mix, values, log.p = TRUE)
  expect_true(all(diff(cdf) >= 0))
  expect_numeric(dens, lower = 0, finite = TRUE, any.missing = FALSE)
  expect_numeric(cdf, lower = 0, upper = 1, finite = TRUE, any.missing = FALSE)
  expect_numeric(lcdf, upper = 0, finite = FALSE, any.missing = FALSE)
}

consistent_ccdf <- function(mix, values) {
  dens <- dmix(mix, values)
  ccdf <- pmix(mix, values, FALSE)
  lccdf <- pmix(mix, values, FALSE, TRUE)
  expect_true(all(diff(ccdf) <= 0))
  expect_numeric(dens, lower = 0, finite = TRUE, any.missing = FALSE)
  expect_numeric(ccdf, lower = 0, upper = 1, finite = TRUE, any.missing = FALSE)
  expect_numeric(lccdf, upper = 0, finite = FALSE, any.missing = FALSE)
}

test_that("Beta CDF function is consistent", {
  consistent_cdf(beta, seq(0.1, 0.9, by = 0.1))
})
test_that("Beta mixture CDF function is consistent", {
  consistent_cdf(betaMix, seq(0.1, 0.9, by = 0.1))
})

test_that("Normal CDF is consistent", {
  consistent_cdf(norm, seq(-2, 2, by = 0.1))
})
test_that("Normal mixture CDF is consistent", {
  consistent_cdf(norm, seq(-2, 2, by = 0.1))
})

test_that("Gamma CDF function is consistent", {
  consistent_cdf(gamma, seq(2, 7, by = 0.1))
})
test_that("Gamma mixture CDF function is consistent", {
  consistent_cdf(gammaMix, seq(2, 7, by = 0.1))
})


## problematic beta which triggers that the cumulative of the
## predictive is not monotone (probably fixed with Stan 2.18, check
## again once 2.18 is out)

## problematic Beta density
bm1 <- mixbeta(c(1.0, 298.30333970, 146.75306521))
test_that("Problematic (1) BetaBinomial CDF function is consistent", {
  consistent_cdf(preddist(bm1, n = 50), 0:50)
})
test_that("Problematic (1) BetaBinomial CCDF function is consistent", {
  consistent_ccdf(preddist(bm1, n = 50), 0:50)
})
bm2 <- mixbeta(c(1.0, 3 + 1 / 3, 47 + 1 / 3))
test_that("Problematic (2) BetaBinomial CDF function is consistent", {
  consistent_cdf(preddist(bm2, n = 50), 0:50)
})
test_that("Problematic (2) BetaBinomial CCDF function is consistent", {
  consistent_ccdf(preddist(bm2, n = 50), 0:50)
})

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
  expect_numeric(s$mean, any.missing = FALSE, len = p)
  expect_matrix(s$cov, any.missing = FALSE, nrows = p, ncols = p)
}

test_that("Multivariate normal mixture has consistent dimensionality", {
  for (i in 1:(nrow(S) - 1)) {
    p_sub <- 4 - i
    S_sub <- S[-c(1:i), -c(1:i), drop = FALSE]
    mvn_consistent_dimension(
      mixmvnorm(c(1, rep(0, p_sub), S_sub), sigma = S_sub),
      p_sub
    )
  }
})

test_that("Multivariate normal mixture has consistent dimension naming", {
  for (i in 1:(nrow(S) - 1)) {
    p_sub <- 4 - i
    S_sub <- S[-c(1:i), -c(1:i), drop = FALSE]
    m_sub <- rep(0, p_sub)
    dim_labels <- letters[1:p_sub]
    names(m_sub) <- dim_labels
    test_mix <- mixmvnorm(c(1, m_sub, S_sub), sigma = S_sub)
    ## now test that names are used consistently
    expect_equal(rownames(sigma(test_mix)), dim_labels)
    expect_equal(colnames(sigma(test_mix)), dim_labels)
    expect_equal(names(summary(test_mix)$mean), dim_labels)
    expect_equal(rownames(summary(test_mix)$cov), dim_labels)
    expect_equal(colnames(summary(test_mix)$cov), dim_labels)
    expect_equal(colnames(rmix(test_mix, 1)), dim_labels)
  }
})

test_that("Multivariate normal mixture has consistent initialization", {
  p <- nrow(S)
  mv1 <- mixmvnorm(c(1, rep(0, p), S), sigma = S, param = "ms")
  mv2 <- mixmvnorm(c(1, rep(0, p), 1), sigma = S, param = "mn")
  mv3 <- mixmvnorm(c(1, rep(0, p), 2), sigma = S, param = "mn")

  expect_equal(summary(mv1)$cov, S, tolerance = eps_lower)
  expect_equal(summary(mv2)$cov, S, tolerance = eps_lower)
  expect_equal(summary(mv3)$cov, S / 2, tolerance = eps_lower)
})

mvn_consistent_summaries <- function(mix, S = Nsamp_equant) {
  samp <- rmix(mix, S)
  m <- colMeans(samp)
  expect_equal(colMeans(samp), summary(mix)$mean, tolerance = eps)
  expect_equal(cov(samp), summary(mix)$cov, tolerance = eps)
}

test_that("Multivariate normal mixture has consistent summaries", {
  p <- nrow(S)
  mv1 <- mixmvnorm(c(1, rep(0, p), S), sigma = S, param = "ms")
  mv2 <- mixmvnorm(c(1, rep(0, p), 1), sigma = S, param = "mn")
  mv3 <- mixmvnorm(
    c(0.2, rep(0, p), 2),
    c(0.8, rep(1, p), 6),
    sigma = S,
    param = "mn"
  )

  mvn_consistent_summaries(mv1)
  mvn_consistent_summaries(mv2)
  mvn_consistent_summaries(mv3)
})


test_that("msr2mvnorm works with valid inputs", {
  m <- c(0, 1)
  s <- c(1, 2)
  r <- c(0.5)
  result <- msr2mvnorm(m, s, r)
  expect_true(is.numeric(result))
  expect_equal(length(result), length(m) + length(m)^2)
})

test_that("msr2mvnorm handles single variable case", {
  m <- c(0)
  s <- c(2)
  r <- numeric() # No correlations for a single variable
  result <- msr2mvnorm(m, s, r)
  expect_true(is.numeric(result))
  expect_equal(length(result), 2)
})

test_that("msr2mvnorm throws error for correlation in one-dimensional input", {
  m <- c(0)
  s <- c(1)
  r <- c(0.5) # Correlation provided for a single variable
  expect_error(msr2mvnorm(m, s, r), "Assertion on 'r' failed")
})

test_that("msr2mvnorm detects empty inputs", {
  m <- numeric(0)
  s <- numeric(0)
  r <- numeric(0)
  expect_error(msr2mvnorm(m, s, r), "Assertion on 'm' failed")
})

test_that("msr2mvnorm detects mismatched dimensions", {
  m <- c(0, 1)
  s <- c(1) # Length mismatch with ~m~
  r <- c(0.5)
  expect_error(msr2mvnorm(m, s, r), "Assertion on 's' failed")
})

test_that("msr2mvnorm detects invalid correlation length", {
  m <- c(0, 1, 2)
  s <- c(1, 2, 3)
  r <- c(0.5, 0.3) # Should be 3 correlations for 3 variables
  expect_error(msr2mvnorm(m, s, r), "Assertion on 'r' failed")
})

test_that("msr2mvnorm works with large inputs", {
  d <- 10
  m <- rnorm(d)
  s <- runif(d, 1, 2)
  r <- runif(d * (d - 1) / 2, -1, 1)
  result <- msr2mvnorm(m, s, r)
  expect_true(is.numeric(result))
  expect_equal(length(result), d + d^2)
})


test_that("mixmvnorm works with a single 3D component in msr format", {
  params <- c(1, 1:3, 4:6, 0.1, 0.2, 0.3)
  result <- mixmvnorm(params, param = "msr")

  # Check class
  expect_s3_class(result, c("mix", "mvnormMix"))

  # Check structure
  expect_true(is.matrix(result))
  expect_equal(ncol(result), 1) # Single component
  expect_equal(
    rownames(result),
    c(
      "w",
      "m[1]",
      "m[2]",
      "m[3]",
      "s[1]",
      "s[2]",
      "s[3]",
      "rho[2,1]",
      "rho[3,1]",
      "rho[3,2]"
    )
  )

  # Check values
  expect_equal(result["w", 1], 1)
  expect_equal(result["m[1]", 1], 1)
  expect_equal(result["m[2]", 1], 2)
  expect_equal(result["m[3]", 1], 3)
  expect_equal(result["s[1]", 1], 4)
  expect_equal(result["s[2]", 1], 5)
  expect_equal(result["s[3]", 1], 6)
  expect_equal(result["rho[2,1]", 1], 0.1)
  expect_equal(result["rho[3,1]", 1], 0.2)
  expect_equal(result["rho[3,2]", 1], 0.3)
})

test_that("mixmvnorm works with a two-component mixture in msr format", {
  params1 <- c(0.2, 1:3, 4:6, 0.1, 0.2, 0.3)
  params2 <- c(0.8, 1:3, 4:6, 0.1, 0.2, 0.3)
  result <- mixmvnorm(params1, params2, param = "msr")

  # Check class
  expect_s3_class(result, c("mix", "mvnormMix"))

  # Check structure
  expect_true(is.matrix(result))
  expect_equal(ncol(result), 2) # Two components
  expect_equal(
    rownames(result),
    c(
      "w",
      "m[1]",
      "m[2]",
      "m[3]",
      "s[1]",
      "s[2]",
      "s[3]",
      "rho[2,1]",
      "rho[3,1]",
      "rho[3,2]"
    )
  )

  # Check values for component 1
  expect_equal(result["w", 1], 0.2)
  expect_equal(result["m[1]", 1], 1)
  expect_equal(result["m[2]", 1], 2)
  expect_equal(result["m[3]", 1], 3)
  expect_equal(result["s[1]", 1], 4)
  expect_equal(result["s[2]", 1], 5)
  expect_equal(result["s[3]", 1], 6)
  expect_equal(result["rho[2,1]", 1], 0.1)
  expect_equal(result["rho[3,1]", 1], 0.2)
  expect_equal(result["rho[3,2]", 1], 0.3)

  # Check values for component 2
  expect_equal(result["w", 2], 0.8)
  expect_equal(result["m[1]", 2], 1)
  expect_equal(result["m[2]", 2], 2)
  expect_equal(result["m[3]", 2], 3)
  expect_equal(result["s[1]", 2], 4)
  expect_equal(result["s[2]", 2], 5)
  expect_equal(result["s[3]", 2], 6)
  expect_equal(result["rho[2,1]", 2], 0.1)
  expect_equal(result["rho[3,1]", 2], 0.2)
  expect_equal(result["rho[3,2]", 2], 0.3)
})

test_that("mixmvnorm throws error for invalid weights in msr format", {
  params1 <- c(-0.2, 1:3, 4:6, 0.1, 0.2, 0.3) # Negative weight
  params2 <- c(1.2, 1:3, 4:6, 0.1, 0.2, 0.3) # Weight > 1
  expect_error(
    mixmvnorm(params1, params2, param = "msr"),
    "Assertion on 'weights' failed"
  )
})

test_that("mixmvnorm throws error for mismatched dimensions in msr format", {
  params1 <- c(0.5, 1:3, 4:6, 0.1, 0.2, 0.3)
  params2 <- c(0.5, 1:4, 4:7, 0.1, 0.2, 0.3, 0.4) # Mismatched dimensions
  expect_error(
    mixmvnorm(params1, params2, param = "msr"),
    "All components must have equal number of parameters."
  )
})

test_that("mixmvnorm works with a large mixture in msr format", {
  d <- 5 # Dimensionality
  components <- 3 # Number of components
  params <- lapply(1:components, function(i) {
    weight <- 1 / components
    mean <- 1:d
    sd <- runif(d, 1, 2)
    corr <- runif(d * (d - 1) / 2, -1, 1)
    c(weight, mean, sd, corr)
  })
  result <- do.call(mixmvnorm, c(params, list(param = "msr")))

  # Check class
  expect_s3_class(result, c("mix", "mvnormMix"))

  # Check structure
  expect_true(is.matrix(result))
  expect_equal(ncol(result), components)
  expect_equal(sum(as.numeric(result["w", ])), 1, tolerance = 1e-6)
})


test_that("msr2mvnorm returns a flattened vector when unlist = TRUE", {
  m <- c(1, 2, 3)
  s <- c(4, 5, 6)
  r <- c(0.1, 0.2, 0.3)
  result <- msr2mvnorm(m, s, r, unlist = TRUE)

  # Check structure
  expect_true(is.numeric(result))
  expect_equal(length(result), length(m) + length(m)^2)
})

test_that("msr2mvnorm returns a list when unlist = FALSE", {
  m <- c(1, 2, 3)
  s <- c(4, 5, 6)
  r <- c(0.1, 0.2, 0.3)
  result <- msr2mvnorm(m, s, r, unlist = FALSE)

  # Check structure
  expect_true(is.list(result))
  expect_named(result, c("m", "s"))

  # Check values
  expect_equal(result$m, m)
  expect_true(is.matrix(result$s))
  expect_equal(dim(result$s), c(length(m), length(m)))
})

test_that("msr2mvnorm handles missing m when unlist = FALSE", {
  s <- c(4, 5, 6)
  r <- c(0.1, 0.2, 0.3)
  result <- msr2mvnorm(s = s, r = r, unlist = FALSE)

  # Check structure
  expect_true(is.list(result))
  expect_named(result, "s")

  # Check values
  expect_true(is.matrix(result$s))
  expect_equal(dim(result$s), c(length(s), length(s)))
})

test_that("msr2mvnorm throws error for invalid correlation length", {
  m <- c(1, 2, 3)
  s <- c(4, 5, 6)
  r <- c(0.1, 0.2) # Invalid length for 3D input
  expect_error(
    msr2mvnorm(m, s, r, unlist = TRUE),
    "Assertion on 'r' failed"
  )
})

test_that("msr2mvnorm works with single variable input and unlist = TRUE", {
  m <- c(1)
  s <- c(4)
  r <- NULL # No correlations for 1D input
  result <- msr2mvnorm(m, s, r, unlist = TRUE)

  # Check structure
  expect_true(is.numeric(result))
  expect_equal(length(result), 1 + 1^2)
})

test_that("msr2mvnorm works with single variable input and unlist = FALSE", {
  m <- c(1)
  s <- c(4)
  r <- NULL # No correlations for 1D input
  result <- msr2mvnorm(m, s, r, unlist = FALSE)

  # Check structure
  expect_true(is.list(result))
  expect_named(result, c("m", "s"))

  # Check values
  expect_equal(result$m, m)
  expect_true(is.matrix(result$s))
  expect_equal(dim(result$s), c(1, 1))
})

test_that("msr2mvnorm works with large input and unlist = FALSE", {
  withr::local_seed(765865)
  d <- 10
  m <- rnorm(d)
  s <- runif(d, 1, 2)
  r <- runif(d * (d - 1) / 2, -1, 1)
  result <- msr2mvnorm(m, s, r, unlist = FALSE)

  # Check structure
  expect_true(is.list(result))
  expect_named(result, c("m", "s"))

  # Check values
  expect_equal(result$m, m)
  expect_true(is.matrix(result$s))
  expect_equal(dim(result$s), c(d, d))
})


test_that("mixmvnorm handles zero standard deviations correctly", {
  # Define a covariance matrix with one dimension having zero variance
  S <- diag(c(0, 2)) %*% matrix(c(1, 0.5, 0.5, 1), 2, 2) %*% diag(c(0, 2))

  # Construct a multivariate normal mixture with zero standard deviation in one dimension
  expect_no_warning(
    mvnm <- mixmvnorm(
      c(1, c(0, 0), diag(c(0, 2)^2))
    )
  )

  # Check that the function does not error
  expect_s3_class(mvnm, c("mix", "mvnormMix"))

  # Check the summary output
  summary_mvnm <- summary(mvnm)

  # Verify that the standard deviations are correctly set to zero
  expect_equal(summary_mvnm$cov[1, 1], 0)
  expect_equal(summary_mvnm$cov[1, 2], 0)
  expect_equal(summary_mvnm$cov[2, 1], 0)
  expect_equal(summary_mvnm$cov[2, 2], 4)
})


test_that("mixmvnorm works for 1D mixtures with zero standard deviation", {
  mvnm <- mixmvnorm(
    rob = c(1, 0, 0) # Weight = 1, mean = 0, variance = 0
  )

  # Check structure
  expect_s3_class(mvnm, c("mix", "mvnormMix"))
  expect_equal(ncol(mvnm), 1) # Single component
  expect_equal(rownames(mvnm), c("w", "m[1]", "s[1]"))

  # Check values
  expect_equal(mvnm["w", 1], 1)
  expect_equal(mvnm["m[1]", 1], 0)
  expect_equal(mvnm["s[1]", 1], 0)
})

test_that("mixmvnorm works for 2D mixtures with one zero standard deviation", {
  S <- diag(c(0, 2)) %*% matrix(c(1, 0.5, 0.5, 1), 2, 2) %*% diag(c(0, 2))
  mvnm <- mixmvnorm(
    rob = c(0.5, c(0, 1), diag(c(0, 2)^2)),
    inf = c(0.5, c(1, 2), S / 4)
  )

  # Check structure
  expect_s3_class(mvnm, c("mix", "mvnormMix"))
  expect_equal(ncol(mvnm), 2) # Two components
  expect_equal(
    rownames(mvnm),
    c("w", "m[1]", "m[2]", "s[1]", "s[2]", "rho[2,1]")
  )

  # Check values for component 1
  expect_equal(mvnm["w", 1], 0.5)
  expect_equal(mvnm["m[1]", 1], 0)
  expect_equal(mvnm["m[2]", 1], 1)
  expect_equal(mvnm["s[1]", 1], 0)
  expect_equal(mvnm["s[2]", 1], 2)
  expect_equal(mvnm["rho[2,1]", 1], 0)

  # Check values for component 2
  expect_equal(mvnm["w", 2], 0.5)
  expect_equal(mvnm["m[1]", 2], 1)
  expect_equal(mvnm["m[2]", 2], 2)
  expect_equal(mvnm["s[1]", 2], 0)
  expect_equal(mvnm["s[2]", 2], 1)
  expect_equal(mvnm["rho[2,1]", 2], 0)
})


test_that("mixmvnorm works for 2D mixtures with one zero standard deviation (msr)", {
  mvnm <- mixmvnorm(
    c(0.5, 0, 1, 0, 2, 0.7), # Weight = 0.5, mean = (0, 1), SD = (0, 2), correlation = 0.7
    c(0.5, 1, 2, 0, 0.5, -0.3), # Weight = 0.5, mean = (1, 2), SD = (0, 0.5), correlation = -0.3
    param = "msr"
  )

  # Check structure
  expect_s3_class(mvnm, c("mix", "mvnormMix"))
  expect_equal(ncol(mvnm), 2) # Two components
  expect_equal(
    rownames(mvnm),
    c("w", "m[1]", "m[2]", "s[1]", "s[2]", "rho[2,1]")
  )

  # Check values for component 1
  expect_equal(mvnm["w", 1], 0.5)
  expect_equal(mvnm["m[1]", 1], 0)
  expect_equal(mvnm["m[2]", 1], 1)
  expect_equal(mvnm["s[1]", 1], 0)
  expect_equal(mvnm["s[2]", 1], 2)
  expect_equal(mvnm["rho[2,1]", 1], 0.7)

  # Check values for component 2
  expect_equal(mvnm["w", 2], 0.5)
  expect_equal(mvnm["m[1]", 2], 1)
  expect_equal(mvnm["m[2]", 2], 2)
  expect_equal(mvnm["s[1]", 2], 0)
  expect_equal(mvnm["s[2]", 2], 0.5)
  expect_equal(mvnm["rho[2,1]", 2], -0.3)
})
