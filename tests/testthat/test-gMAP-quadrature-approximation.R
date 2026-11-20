## Generalized quadrature-based sampler-scaling initialization.
##
## These tests are deterministic (no MCMC
## sampling); they cover `.gmap_quadrature_approximation()` directly and the
## `chains = 0` skeleton `gMAP()` produces from it. Sampled/geometry
## behaviour (endpoint identities, Jacobian exactness) is covered by
## `test-gMAP-s2z-center.R` and `test-gMAP-s2z.R`, which already exercise the
## `group_scale_guess`-aware Stan target via `helper-s2z.R`.

skip_on_cran()

test_that(".gmap_quadrature_approximation() returns finite, positive scales", {
  d <- RBesT::AS
  fit0 <- suppressMessages(gMAP(
    cbind(r, n - r) ~ 1 | study,
    data = d, family = binomial,
    tau.dist = "HalfNormal", tau.prior = 1, beta.prior = 2, chains = 0
  ))
  ins <- fit0$fit.data
  group.stratum <- rep(ins$tau_strata_pred, ins$n_groups)
  group.stratum[ins$group_index] <- ins$tau_strata_index

  res <- RBesT:::.gmap_quadrature_approximation(
    family = "binomial", y = ins$r, y.se = NULL, r = ins$r, n = ins$r_n,
    count = ins$count, log.offset = ins$log_offset, X = ins$X,
    beta.prior = ins$beta_prior, group.index = ins$group_index,
    group.stratum = group.stratum, tau.dist = "HalfNormal",
    tau.prior = ins$tau_prior, prior_PD = FALSE
  )

  expect_true(res$ok)
  expect_true(is.finite(res$log.tau.location))
  expect_true(is.finite(res$log.tau.scale) && res$log.tau.scale > 0)
  expect_length(res$beta.scale, ncol(ins$X))
  expect_true(all(is.finite(res$beta.scale)) && all(res$beta.scale > 0))
  expect_length(res$group.scale, ins$n_groups)
  expect_true(all(is.finite(res$group.scale)) && all(res$group.scale > 0))
  expect_identical(length(res$group.information), ins$n_groups)
  expect_true(all(res$group.information > 0)) # every AS study has data
})

test_that(".gmap_quadrature_approximation() is inert but ok for Fixed tau", {
  d <- RBesT::AS
  fit0 <- suppressMessages(gMAP(
    cbind(r, n - r) ~ 1 | study,
    data = d, family = binomial,
    tau.dist = "Fixed", tau.prior = 0.3, beta.prior = 2, chains = 0
  ))
  ins <- fit0$fit.data
  group.stratum <- rep(ins$tau_strata_pred, ins$n_groups)
  group.stratum[ins$group_index] <- ins$tau_strata_index

  res <- RBesT:::.gmap_quadrature_approximation(
    family = "binomial", y = ins$r, y.se = NULL, r = ins$r, n = ins$r_n,
    count = ins$count, log.offset = ins$log_offset, X = ins$X,
    beta.prior = ins$beta_prior, group.index = ins$group_index,
    group.stratum = group.stratum, tau.dist = "Fixed",
    tau.prior = ins$tau_prior, prior_PD = FALSE
  )

  ## tau_raw_guess is inert (ignored by the Stan model) when tau.dist ==
  ## "Fixed"; the quadrature approximation still succeeds (beta/group scaling
  ## remain well defined at the fixed heterogeneity value) but reports no
  ## log(tau) location/scale.
  expect_true(res$ok)
  expect_true(is.na(res$log.tau.location))
  expect_true(is.na(res$log.tau.scale))
  expect_true(all(is.finite(res$beta.scale)) && all(res$beta.scale > 0))
  expect_true(all(is.finite(res$group.scale)) && all(res$group.scale > 0))
})

test_that(".gmap_quadrature_approximation() is a no-op for prior_PD", {
  res <- RBesT:::.gmap_quadrature_approximation(
    family = "binomial", y = RBesT::AS$r, y.se = NULL, r = RBesT::AS$r,
    n = RBesT::AS$r + RBesT::AS$n, count = integer(8), log.offset = numeric(8),
    X = matrix(1, 8, 1), beta.prior = cbind(0, 2),
    group.index = seq_len(8), group.stratum = rep(1L, 8),
    tau.dist = "HalfNormal", tau.prior = cbind(0, 1), prior_PD = TRUE
  )
  expect_true(res$ok)
  expect_true(all(res$center == 0))
  expect_true(all(res$group.information == 0))
})

test_that(".gmap_quadrature_approximation() reports the documented guard reasons", {
  d4 <- transform(
    RBesT::AS,
    stratum = factor(rep(c("A", "B", "C", "D"), length.out = 8))
  )
  expect_warning(
    fit4b <- gMAP(
      cbind(r, n - r) ~ 1 | study,
      data = d4, tau.strata = stratum, family = binomial,
      tau.dist = "HalfNormal",
      tau.prior = matrix(c(0, 1, 0, 1, 0, 1, 0, 1), 4, 2, byrow = TRUE),
      beta.prior = 2, chains = 0
    ),
    "more than three observed non-fixed tau strata"
  )
  expect_identical(fit4b$fit.data$re_param, 1L)

  d.under <- transform(RBesT::AS, stratum = factor(c(rep("A", 7), "B")))
  expect_warning(
    fit.under <- gMAP(
      cbind(r, n - r) ~ 1 | study,
      data = d.under, tau.strata = stratum, family = binomial,
      tau.dist = "HalfNormal",
      tau.prior = matrix(c(0, 1, 0, 1), 2, 2, byrow = TRUE),
      beta.prior = 2, chains = 0
    ),
    "fewer than two informative groups"
  )
  expect_identical(fit.under$fit.data$re_param, 1L)
})

test_that("gMAP() populates the new per-stratum/per-group Stan data fields", {
  fit0 <- suppressMessages(gMAP(
    cbind(r, n - r) ~ 1 | study,
    data = RBesT::AS, family = binomial,
    tau.dist = "HalfNormal", tau.prior = 1, beta.prior = 2, chains = 0
  ))
  d <- fit0$fit.data

  ## tau_raw_guess is now a 2 x n_tau_strata matrix (location/scale per
  ## stratum), mirroring beta_raw_guess's existing per-coefficient shape.
  expect_identical(dim(d$tau_raw_guess), c(2L, 1L))
  expect_true(all(is.finite(d$tau_raw_guess)))
  expect_true(d$tau_raw_guess["scale", ] > 0)

  ## hybrid beta scaling: the pooled-GLM coefficient is kept as the location
  ## exactly; only the scale is quadrature-derived.
  pooled <- glm.fit(
    d$X, RBesT::AS$r / RBesT::AS$n,
    weights = as.vector(RBesT::AS$n),
    family = binomial()
  )$coefficients
  expect_equal(unname(d$beta_raw_guess["mean", ]), unname(pooled))
  expect_true(all(d$beta_raw_guess["sd", ] > 0))

  ## per-group physical-scale numeraire for partial centering.
  expect_identical(length(d$group_scale_guess), d$n_groups)
  expect_true(all(d$group_scale_guess > 0))

  ## per-group posterior location guess for the same affine map: consumed by
  ## both the conventional and s2z Stan partial-centering paths (see
  ## partial_center_loc() / the s2z branch in inst/stan/gMAP.stan).
  ## It must match the
  ## quadrature approximation's own `group.location` output exactly, i.e.
  ## gMAP() must actually consume it, not merely compute and discard it.
  expect_identical(length(d$group_location_guess), d$n_groups)
  expect_true(all(is.finite(d$group_location_guess)))
  group.stratum <- rep(d$tau_strata_pred, d$n_groups)
  group.stratum[d$group_index] <- d$tau_strata_index
  quad <- RBesT:::.gmap_quadrature_approximation(
    family = "binomial", y = d$r, y.se = NULL, r = d$r, n = d$r_n,
    count = d$count, log.offset = d$log_offset, X = d$X,
    beta.prior = d$beta_prior, group.index = d$group_index,
    group.stratum = group.stratum, tau.dist = "HalfNormal",
    tau.prior = d$tau_prior, prior_PD = FALSE
  )
  expect_true(quad$ok)
  expect_equal(as.numeric(d$group_location_guess), as.numeric(quad$group.location))
  ## a genuine posterior location guess is not all-zero for informative data
  expect_true(any(abs(d$group_location_guess) > 1e-6))

  ## At the fully centered S2Z endpoint, a zero raw coordinate must map to
  ## the estimated group location projected onto the zero-sum subspace.
  centered <- d
  centered$re_param <- 2L
  centered$re_center[] <- 1
  re.at.zero <- s2z_re(
    centered, numeric(centered$n_groups - 1L), tau1 = 0.5
  )
  expected.location <- centered$group_location_guess -
    mean(centered$group_location_guess)
  expect_equal(re.at.zero, expected.location)
})

test_that("RBesT.MC.rescale = FALSE unit-scales tau/beta/group guesses", {
  fit <- withr::with_options(
    list(RBesT.MC.rescale = FALSE),
    suppressMessages(gMAP(
      cbind(r, n - r) ~ 1 | study,
      data = RBesT::AS, family = binomial,
      tau.dist = "HalfNormal", tau.prior = 1, beta.prior = 2, chains = 0
    ))
  )
  d <- fit$fit.data
  expect_true(all(d$tau_raw_guess["scale", ] == 1))
  expect_true(all(d$beta_raw_guess["sd", ] == 1))
  expect_true(all(d$group_scale_guess == 1))
  ## rescale only touches scales, matching how beta_raw_guess's own location
  ## (the pooled-GLM mean) is left untouched; group_location_guess must
  ## still be the quadrature location, not reset to 0.
  expect_true(any(abs(d$group_location_guess) > 1e-6))
})

test_that("ncp = 2 selects the CP endpoint from quadrature fractions at the legacy-equivalent threshold", {
  ## Extreme per-group sample sizes push every fraction well above
  ## 400/401 (the fraction implied by the legacy `tau_guess / se > 20` rule),
  ## so ncp = 2 must select the centered endpoint.
  d.rich <- data.frame(
    r = c(5000L, 15000L, 25000L, 3000L),
    n = c(10000L, 20000L, 30000L, 4000L),
    study = factor(paste0("S", 1:4))
  )
  fit.rich <- withr::with_options(
    list(RBesT.MC.ncp = 2),
    suppressMessages(gMAP(
      cbind(r, n - r) ~ 1 | study,
      data = d.rich, family = binomial,
      tau.dist = "HalfNormal", tau.prior = 1, beta.prior = 2, chains = 0
    ))
  )
  expect_identical(fit.rich$fit.data$re_param, 0L)

  ## The legacy centered path requires an intercept. Automatic CP must use
  ## the exact all-centered partial endpoint for an intercept-free model.
  d.rich$x <- 1
  fit.no.intercept <- withr::with_options(
    list(RBesT.MC.ncp = 2),
    suppressMessages(gMAP(
      cbind(r, n - r) ~ 0 + x | study,
      data = d.rich, family = binomial,
      tau.dist = "HalfNormal", tau.prior = 1, beta.prior = 2, chains = 0
    ))
  )
  expect_identical(fit.no.intercept$fit.data$re_param, 2L)
  expect_true(all(fit.no.intercept$fit.data$re_center == 1))

  ## Sparse data must remain at the non-centered endpoint.
  d.sparse <- data.frame(
    r = c(1L, 0L, 2L, 1L), n = c(20L, 18L, 22L, 19L),
    study = factor(paste0("S", 1:4))
  )
  fit.sparse <- withr::with_options(
    list(RBesT.MC.ncp = 2),
    suppressMessages(gMAP(
      cbind(r, n - r) ~ 1 | study,
      data = d.sparse, family = binomial,
      tau.dist = "HalfNormal", tau.prior = 1, beta.prior = 2, chains = 0
    ))
  )
  expect_identical(fit.sparse$fit.data$re_param, 1L)
})

test_that("ncp = 2 falls back to the non-centered endpoint when quadrature fails", {
  d4 <- transform(
    RBesT::AS,
    stratum = factor(rep(c("A", "B", "C", "D"), length.out = 8))
  )
  fit <- withr::with_options(
    list(RBesT.MC.ncp = 2),
    suppressMessages(suppressWarnings(gMAP(
      cbind(r, n - r) ~ 1 | study,
      data = d4, tau.strata = stratum, family = binomial,
      tau.dist = "HalfNormal",
      tau.prior = matrix(c(0, 1, 0, 1, 0, 1, 0, 1), 4, 2, byrow = TRUE),
      beta.prior = 2, chains = 0
    )))
  )
  expect_identical(fit$fit.data$re_param, 1L)
})

test_that("ncp = 2 resolving to the centered endpoint on a no-intercept model rewrites to the all-one partial endpoint", {
  ## re_param = 0 (the literal CP fast path) requires an intercept (it zeroes
  ## X_param's first column, which must be the all-ones intercept column);
  ## ncp = 2 must therefore rewrite a CP resolution to the general partial
  ## machinery with re_center all equal to one, exactly like a literal
  ## `RBesT.MC.ncp = 0` request already does for a no-intercept model
  ## (see "centered no-intercept models use the exact partial endpoint" in
  ## test-gMAP-s2z-center.R). Extreme per-group sample sizes push every
  ## fraction above the 400/401 threshold, so ncp = 2 resolves to CP here.
  d.rich <- data.frame(
    r = c(50000L, 150000L, 250000L, 30000L),
    n = c(100000L, 200000L, 300000L, 40000L),
    study = factor(paste0("S", 1:4))
  )
  fit <- withr::with_options(
    list(RBesT.MC.ncp = 2),
    suppressMessages(suppressWarnings(gMAP(
      cbind(r, n - r) ~ -1 + study | study,
      data = d.rich, family = binomial,
      tau.dist = "HalfNormal", tau.prior = 1, beta.prior = 2, chains = 0
    )))
  )
  expect_identical(fit$fit.data$re_param, 2L)
  expect_true(all(fit$fit.data$re_center == 1))

  ## and this must actually be samplable (no "requires treatment contrast
  ## parametrization" Stan-level reject at the literal CP fast path)
  fit.sampled <- suppressMessages(suppressWarnings(withr::with_options(
    list(RBesT.MC.ncp = 2),
    gMAP(
      cbind(r, n - r) ~ -1 + study | study,
      data = d.rich, family = binomial,
      tau.dist = "HalfNormal", tau.prior = 1, beta.prior = 2,
      chains = 1, iter = 200, warmup = 100, thin = 1, cores = 1
    )
  )))
  expect_s3_class(fit.sampled, "gMAP")
  expect_identical(posterior::ndraws(fit.sampled$draws), 100L)
})
