## Tests for the sum-to-zero (S2Z) reparametrization of the gMAP Stan model.
##
## See design/design-sum-to-zero-gmap.md, section 7 (validation ladder).
##
## Three tiers, cheapest first:
##
##   1. Pure R arithmetic and the analytic derivation. CRAN-safe.
##   2. Deterministic checks against the compiled Stan model through
##      `chains = 0` skeletons and `rstan::constrain_pars()`. The structural
##      invariants hold for *any* parameter vector, so they are asserted on a
##      chosen grid of unconstrained vectors rather than on posterior draws:
##      exact, free of Monte Carlo error, and able to reach extreme `tau`.
##   3. Cached MCMC fixtures, used only where the assertion is genuinely about
##      a distribution. Built by `make -j4 test-fixtures` from the recipes in
##      `tests/testthat/fixtures-mcmc-src/`.

## ---------------------------------------------------------------------------
## 1. Arithmetic and derivation (pure R, CRAN-safe)
## ---------------------------------------------------------------------------

test_that("s2z recovery coefficients match the variance form", {
  ## The implementation uses sd_alpha = hypot(s1, tau/sqrt(J)) and
  ## r = (tau/sqrt(J)) / sd_alpha, so that
  ##   v / (s1^2 + v) = r^2  and  sqrt(s1^2 v / (s1^2 + v)) = s1 r
  ## with v = tau^2 / J. Verify against the direct variance form over a wide
  ## grid; this is exactly the algebra the Stan code performs. The grid spans
  ## the degenerate limits (tau -> 0, s1 -> 0, s1 -> Inf), which are spelled
  ## out separately below because they carry the interpretation.
  grid <- expand.grid(
    s1 = 10^seq(-3, 4, length.out = 15),
    tau = 10^seq(-6, log10(50), length.out = 15),
    J = c(1, 2, 3, 6, 11, 40)
  )

  v <- grid$tau^2 / grid$J
  mean_coef_var <- v / (grid$s1^2 + v)
  sd_coef_var <- sqrt(grid$s1^2 * v / (grid$s1^2 + v))

  r <- s2z_recovery(grid$s1, grid$tau, grid$J)$r

  expect_true(all(r > 0 & r < 1))
  expect_equal(r^2, mean_coef_var, tolerance = 1e-12)
  expect_equal(grid$s1 * r, sd_coef_var, tolerance = 1e-12)

  coefs <- function(s1, tau, J) {
    rr <- s2z_recovery(s1, tau, J)$r
    c(mean = rr^2, sd = s1 * rr)
  }

  ## tau -> 0: nothing to recover, the shift collapses to zero
  expect_equal(unname(coefs(2, 1e-12, 6)), c(0, 0), tolerance = 1e-12)

  ## s1 -> Inf: a vague intercept prior means alpha carries no information
  ## about the shift, so the mean coefficient vanishes while the shift keeps
  ## its full prior sd tau/sqrt(J)
  vague <- coefs(1e8, 1, 6)
  expect_lt(unname(vague["mean"]), 1e-12)
  expect_equal(unname(vague["sd"]), 1 / sqrt(6), tolerance = 1e-8)

  ## s1 -> 0: a point mass prior pins beta[1] = m1, i.e. abar = alpha - m1
  lim <- coefs(1e-12, 1, 6)
  expect_equal(unname(lim["mean"]), 1, tolerance = 1e-12)
  expect_lt(unname(lim["sd"]), 1e-11)
})

test_that("the Helmert basis spans the zero-sum subspace", {
  for (J in c(1, 2, 3, 6, 11, 40)) {
    Q <- s2z_helmert_basis(J)
    expect_equal(dim(Q), c(J, J - 1L))
    if (J > 1) {
      expect_equal(crossprod(Q), diag(J - 1), tolerance = 1e-12)
      expect_equal(colSums(Q), rep(0, J - 1), tolerance = 1e-12)
      ## QQ' is the centering projector
      expect_equal(
        Q %*% t(Q),
        diag(J) - matrix(1 / J, J, J),
        tolerance = 1e-12
      )
    }
  }
})

test_that("the s2z marginal likelihood equals the legacy one for every tau", {
  ## Closed-form, zero Monte Carlo error. For the normal endpoint the group
  ## effects and the intercept can be integrated out analytically under both
  ## parameterisations:
  ##
  ##   legacy: y ~ N(m1 1, s1^2 11' + tau^2 ZZ' + diag(y_se^2))
  ##   s2z:    y ~ N(m1 1, (s1^2 + tau^2/J) 11' + tau^2 ZQQ'Z' + diag(y_se^2))
  ##
  ## These are the same matrix because QQ' = I - 11'/J and Z1_J = 1_H, which
  ## is the algebraic identity the whole reparametrization rests on: the
  ## widened prior variance s1^2 + tau^2/J and the unscaled subspace effects
  ## are exactly what make the two agree, and either the (1 - 1/J)^-1/2
  ## rescale or a different widening would break the identity below.
  ##
  ## Scope, so this is not mistaken for more than it is: both matrices are
  ## built here in R, so this checks the *derivation*, not the Stan code. The
  ## implementation is held to it by the log_prob test below, which does read
  ## the compiled target.
  crohn <- RBesT::crohn
  y <- crohn$y
  y_se <- 88 / sqrt(crohn$n)
  J <- nrow(crohn)
  H <- J
  Z <- diag(J)[seq_len(H), , drop = FALSE]
  Q <- s2z_helmert_basis(J)
  one <- matrix(1, H, 1)
  m1 <- 5
  s1 <- 88

  dmvn <- function(x, mu, S) {
    ch <- chol(S)
    z <- backsolve(ch, x - mu, transpose = TRUE)
    -0.5 * length(x) * log(2 * pi) - sum(log(diag(ch))) - 0.5 * sum(z^2)
  }

  for (tau in c(1e-6, 0.01, 1, 5, 44, 200)) {
    S_legacy <- s1^2 * tcrossprod(one) +
      tau^2 * tcrossprod(Z) +
      diag(y_se^2)
    S_s2z <- (s1^2 + tau^2 / J) * tcrossprod(one) +
      tau^2 * tcrossprod(Z %*% Q) +
      diag(y_se^2)

    expect_lt(max(abs(S_legacy - S_s2z)), 1e-8)
    expect_equal(
      dmvn(y, rep(m1, H), S_legacy),
      dmvn(y, rep(m1, H), S_s2z),
      tolerance = 1e-12
    )
  }
})

test_that("the sampler control defaults follow the parametrization", {
  ## The default was lowered from 0.99 to 0.95 once the sum-to-zero
  ## reparametrization removed the beta[1]/mean(eps) ridge: benchmarked
  ## divergence free over 900 chains and 15 non-centered scenarios while
  ## sampling ~1.6x faster (design/s2z/benchmark-ncp-adapt-delta.R). Opting
  ## out restores the legacy 0.99, since the legacy geometry needs the more
  ## conservative target.
  ##
  ## This pins the values; that they actually reach the sampler is checked
  ## by the step-size test at the bottom of this file.
  expect_equal(.gmap_sampler_control(TRUE, list())$adapt_delta, 0.95)
  expect_equal(.gmap_sampler_control(FALSE, list())$adapt_delta, 0.99)

  expect_equal(.gmap_sampler_control(TRUE, list())$stepsize, 0.01)
  expect_equal(.gmap_sampler_control(TRUE, list())$max_treedepth, 20)

  ## an explicit user setting must win over either default, in both
  ## directions -- this is the documented escape hatch
  expect_equal(
    .gmap_sampler_control(TRUE, list(adapt_delta = 0.99))$adapt_delta,
    0.99
  )
  expect_equal(
    .gmap_sampler_control(FALSE, list(adapt_delta = 0.90))$adapt_delta,
    0.90
  )
  ## ... while leaving the other entries alone
  expect_equal(
    .gmap_sampler_control(TRUE, list(adapt_delta = 0.99))$max_treedepth,
    20
  )

  ## and the user settings are read from the documented option
  withr::with_options(
    list(RBesT.MC.control = list(adapt_delta = 0.8, stepsize = 0.5)),
    {
      expect_equal(.gmap_sampler_control(TRUE)$adapt_delta, 0.8)
      expect_equal(.gmap_sampler_control(TRUE)$stepsize, 0.5)
    }
  )
})

## ---------------------------------------------------------------------------
## 2. Deterministic checks against the compiled model (no sampling)
## ---------------------------------------------------------------------------

test_that("the S2Z path activates exactly as specified", {
  skip_on_cran()

  ## Both sides of the predicate in one sweep. The fallbacks are ordinary
  ## documented arguments, not opt-ins, so the negative side has to be
  ## asserted, not assumed; `optout` is the user escape hatch.
  ##
  ## The observable signature is the parameter layout itself, which the
  ## `chains = 0` model instance already carries: J - 1 free subspace
  ## coordinates plus the recovery direction under s2z, against J group
  ## effects and no recovery direction otherwise.
  scenarios <- s2z_test_scenarios()

  for (nm in names(scenarios)) {
    d <- s2z_test_skeleton(nm)$fit.data
    on <- scenarios[[nm]]$s2z

    expect_identical(s2z_active(s2z_test_skeleton(nm)), on, info = nm)

    n_up <- rstan::get_num_upars(s2z_stanfit(nm))
    cp <- s2z_constrained(nm, rep(0, n_up))

    n_re <- if (on) d$n_groups - 1L else d$n_groups
    n_rec <- if (on) 1L else 0L

    expect_identical(length(cp$xi_eta), as.integer(n_re), info = nm)
    expect_identical(length(cp$xi_abar), n_rec, info = nm)
    expect_identical(
      n_up,
      as.integer(d$mX + d$n_tau_strata + n_re + n_rec),
      info = nm
    )

    ## the reparametrization must not disturb the data-row bookkeeping, which
    ## is what `fitted()` reports -- in particular for J = 1, where the free
    ## coordinate vector is empty. The skeleton carries no draws, hence the
    ## expected "no posterior samples" warning.
    expect_identical(
      nrow(suppressWarnings(fitted(s2z_test_skeleton(nm)))),
      as.integer(d$H),
      info = nm
    )
  }
})

test_that("the s2z identities hold for every parameter vector", {
  skip_on_cran()

  ## The three structural invariants of design 3.2.1/3.2.2 in one sweep:
  ##
  ##   1'(Q xi) = 0                                  (zero-sum subspace)
  ##   theta = X beta[alpha] + (Q xi)[group]         (theta reconstruction)
  ##   beta[1] = alpha - r (r (alpha - m1) + s1 zeta)  (intercept recovery)
  ##
  ## None of these is a statement about the posterior -- they hold for any
  ## parameter vector -- so they are checked deterministically through
  ## `constrain_pars()` over a grid that sweeps tau across several orders of
  ## magnitude. That is both exact and sharper than sampled draws, which never
  ## visit the extremes where a wrong widening or a rescaled basis shows up.
  ##
  ## Note that this re-implements the recovery formula and so can only detect
  ## a change to it, never a wrong derivation: `abar` deliberately does not
  ## enter the target, which also makes it invisible to the log_prob check.
  ## The prior-marginal test below is the assertion with teeth there.
  scenarios <- s2z_test_scenarios()
  on_names <- names(scenarios)[vapply(scenarios, function(s) s$s2z, logical(1))]

  for (nm in on_names) {
    d <- s2z_test_skeleton(nm)$fit.data
    J <- d$n_groups
    m1 <- d$beta_prior[1, 1]
    s1 <- d$beta_prior[1, 2]

    for (upars in s2z_upars_grid(nm)) {
      cp <- s2z_constrained(nm, upars)

      beta <- as.numeric(cp$beta)
      theta <- as.numeric(cp$theta)
      tau1 <- as.numeric(cp$tau)[1]
      alpha <- s2z_alpha(d, as.numeric(cp$beta_raw)[1])
      re <- s2z_re(d, cp$xi_eta, tau1)

      ## a relative bound: tau ranges over orders of magnitude here, so an
      ## absolute one would be a statement about the scenario, not the algebra
      scale <- max(1, abs(alpha), tau1, max(abs(theta)), max(abs(beta)))
      tol <- 1e-10 * scale
      lbl <- paste0(nm, " (tau = ", signif(tau1, 3), ")")

      ## 1'(Q xi) = 0 to machine precision
      expect_lt(abs(sum(re)), tol, label = lbl)

      ## beta[1] is the recovered super-population intercept; theta is built
      ## from the sampled intercept alpha instead
      beta_alpha <- beta
      beta_alpha[1] <- alpha
      theta_hat <- as.numeric(d$X %*% beta_alpha) + re[d$group_index]
      expect_lt(max(abs(theta - theta_hat)), tol, label = lbl)

      ## the closed form of the recovered intercept
      r <- s2z_recovery(s1, tau1, J)$r
      abar <- r * (r * (alpha - m1) + s1 * as.numeric(cp$xi_abar))
      expect_lt(abs(beta[1] - (alpha - abar)), tol, label = lbl)

      ## J = 1: the free coordinate vector has length 0 and the random effect
      ## is entirely absorbed, so Stan must handle vector[0] and matrix[1,0]
      ## and theta must equal alpha exactly
      if (J == 1L) {
        expect_identical(length(cp$xi_eta), 0L, info = nm)
        expect_lt(max(abs(theta - alpha)), tol, label = lbl)
      }
    }
  }
})

test_that("the Stan s2z target matches an independent implementation", {
  skip_on_cran()

  ## Differences of the target are compared rather than absolute values,
  ## because Stan drops normalising constants that do not depend on any
  ## parameter. A *tau-dependent* constant -- the one real risk here, since the
  ## widened intercept prior has a tau-dependent sd -- does not cancel in such
  ## differences, so this is a sharp and completely noise-free detector.
  cases <- c("binomial_ncp", "binomial_cp", "fixed_tau", "nonzero_prior_mean")

  for (nm in cases) {
    d <- s2z_test_skeleton(nm)$fit.data
    sf <- s2z_stanfit(nm)
    ## vary tau over several orders of magnitude; that is what a dropped
    ## tau-dependent normalising constant would show up in
    grid <- s2z_upars_grid(nm, n_tau = 13L)

    stan_lp <- vapply(grid, function(u) rstan::log_prob(sf, u), numeric(1))
    ref_lp <- vapply(
      grid,
      function(u) s2z_reference_log_prob(d, u),
      numeric(1)
    )

    expect_true(all(is.finite(stan_lp)), info = nm)
    ## differences must agree to machine precision
    delta <- (stan_lp - stan_lp[1]) - (ref_lp - ref_lp[1])
    expect_lt(max(abs(delta)), 1e-8, label = nm)
  }
})

test_that("the draws skeleton keeps the legacy variable layout", {
  skip_on_cran()

  ## This is the legacy variable set, order and dimension. S2Z keeps beta in
  ## transformed parameters precisely so that this stays true; only lp__
  ## changes value (not position). The sampled counterpart is asserted in the
  ## summaries test below.
  legacy <- c(
    paste0("theta[", 1:8, "]"),
    "beta[1]",
    "tau[1]",
    "theta_pred",
    "theta_resp_pred",
    "lp__"
  )

  expect_identical(
    posterior::variables(s2z_test_skeleton("binomial_ncp")$draws),
    legacy
  )
  ## and opting out must not move anything
  expect_identical(
    posterior::variables(s2z_test_skeleton("optout")$draws),
    legacy
  )
})

test_that("a degenerate intercept prior is rejected rather than silently wrong", {
  skip_on_cran()

  ## s1 = 0 makes sd_alpha = 0 whenever tau is fixed at 0 as well. The guard
  ## lives in the Stan `transformed data` block, which already runs when the
  ## model is instantiated -- so the *identity* of the guard can be pinned
  ## without sampling at all. Instantiation reports the rejection on the
  ## message stream rather than raising, hence the capture.
  d <- s2z_test_skeleton("binomial_ncp")$fit.data
  d$beta_prior[1, 2] <- 0
  rejection <- capture.output(
    suppressWarnings(s2z_log_prob_fit(d)),
    type = "message"
  )
  expect_match(
    paste(rejection, collapse = "\n"),
    "s2z requires a strictly positive intercept prior sd"
  )

  ## and the user must see a hard failure rather than a silently wrong fit.
  ## gMAP() wraps the Stan rejection, so only the generic message is pinned
  ## here; the specific guard is the assertion above. Two iterations suffice,
  ## the rejection happens before the first gradient evaluation.
  expect_error(
    suppressMessages(suppressWarnings(gMAP(
      cbind(r, n - r) ~ 1 | study,
      data = AS,
      family = binomial,
      tau.dist = "Fixed",
      tau.prior = 0.25,
      beta.prior = rbind(c(0, 0)),
      warmup = 1,
      iter = 2,
      chains = 1
    ))),
    regexp = "Stan sampler did not run successfully"
  )
})

test_that("the S2Z option is validated", {
  skip_on_cran()

  ## `chains = 0` never reaches the sampler, so this stays cheap; the flag is
  ## validated during model setup.
  bad_option <- function(value) {
    withr::with_options(list(RBesT.MC.s2z = value), {
      suppressMessages(suppressWarnings(gMAP(
        cbind(r, n - r) ~ 1 | study,
        data = AS,
        family = binomial,
        tau.dist = "HalfNormal",
        tau.prior = 0.5,
        beta.prior = 2,
        chains = 0
      )))
    })
  }

  ## the message is pinned to the flag being validated, so this cannot start
  ## passing because of some unrelated setup failure
  expect_error(bad_option("yes"), regexp = "use_s2z")
  expect_error(bad_option(NA), regexp = "use_s2z")
})

## ---------------------------------------------------------------------------
## 3. Distributional assertions (cached MCMC fixtures)
## ---------------------------------------------------------------------------

test_that("the recovery direction zeta stays exactly standard normal", {
  ## zeta enters no density statement other than its own std_normal() prior,
  ## and beta[1] enters neither the likelihood nor any prior. The posterior
  ## therefore factorises and zeta must be an exact N(0,1) draw. This is the
  ## cheapest detector for beta[1] leaking back into the target, e.g. by
  ## re-adding a vectorised `beta ~ normal(...)` under s2z.
  ##
  ## The claim is scenario-independent by construction, so it is checked on
  ## both parametrizations and on a second endpoint rather than on every
  ## scenario in the table.
  for (nm in c("zeta_binomial_ncp", "zeta_binomial_cp", "zeta_normal_ncp")) {
    fit <- s2z_fixture_fit(nm)
    expect_true(s2z_active(fit), info = nm)

    zeta_draws <- posterior::subset_draws(fit$draws, variable = "xi_abar[1]")
    zeta <- as.numeric(zeta_draws)

    ## distribution: the actual correctness detector, sharp but not flaky
    expect_gt(
      suppressWarnings(stats::ks.test(zeta, "pnorm")$p.value),
      1e-3,
      label = nm
    )
    expect_equal(mean(zeta), 0, tolerance = 0.1, info = nm)
    expect_equal(stats::sd(zeta), 1, tolerance = 0.1, info = nm)

    ## efficiency: zeta is independent in the target but shares a step size
    ## and mass matrix with the rest, so it can sit somewhat below N. The band
    ## is deliberately loose: a genuine leak of beta[1] into the target would
    ## collapse this far below 0.5, while the sharp detectors above are the KS
    ## and moment checks.
    ess_rel <- posterior::ess_bulk(zeta_draws) / length(zeta)
    expect_gt(ess_rel, 0.5, label = nm)
  }
})

test_that("the recovered intercept has exactly its prior marginal", {
  ## The deterministic identity test above re-implements the recovery formula
  ## and so can only detect a change to it, never a wrong derivation. This is
  ## the assertion with teeth. Under prior_PD the intercept marginal is known
  ## in closed form whatever tau does,
  ##
  ##   beta[1] ~ Normal(m1, s1),
  ##
  ## because the widening of the alpha prior and the recovery of abar must
  ## compose back exactly. Splitting the recovery as r instead of r^2, or
  ## mismatching the widened variance, breaks the sd while leaving every
  ## other invariant in this file intact.
  fit <- s2z_fixture_fit("prior_pd")
  expect_true(s2z_active(fit))

  ## read the prior back off the fit, so the fixture recipe and the test
  ## cannot drift apart
  m1 <- fit$fit.data$beta_prior[1, 1]
  s1 <- fit$fit.data$beta_prior[1, 2]
  tau_sd <- fit$fit.data$tau_prior[1, 2]

  beta1 <- posterior::extract_variable_array(fit$draws, "beta[1]")

  ## compare on the Monte Carlo scale, so the tolerance follows the actual
  ## precision of the run rather than a hand-picked constant
  expect_lt(abs(mean(beta1) - m1) / posterior::mcse_mean(beta1), 4)
  expect_lt(abs(sd(beta1) - s1) / posterior::mcse_sd(beta1), 4)

  ## tau is untouched by the reparametrization and must keep its own marginal
  tau <- posterior::extract_variable_array(fit$draws, "tau[1]")
  expect_lt(
    abs(mean(tau) - tau_sd * sqrt(2 / pi)) / posterior::mcse_mean(tau),
    4
  )
  expect_lt(
    abs(sd(tau) - tau_sd * sqrt(1 - 2 / pi)) / posterior::mcse_sd(tau),
    4
  )
})

test_that("neither xi_abar nor the raw parameters pollute the summaries", {
  ## verbose fits keep xi_abar and the raw parameters in the draws; the
  ## downstream selectors match base names and must not pick them up
  fit <- s2z_fixture_fit("zeta_binomial_ncp")
  expect_true("xi_abar[1]" %in% posterior::variables(fit$draws))

  ## the sampled layout, once the raw parameters are dropped, is the legacy
  ## one -- i.e. exactly what a non-verbose fit and the draw-free skeleton
  ## report
  expect_identical(
    s2z_reported_variables(fit),
    posterior::variables(s2z_test_skeleton("binomial_ncp")$draws)
  )

  s <- summary(fit)
  ## the legacy contract: beta rows are the design matrix columns, `theta` is
  ## the one-row super-population mean and `tau` carries no tau_raw
  expect_identical(rownames(s$beta), colnames(fit$X))
  expect_identical(rownames(s$tau), "tau[1]")
  expect_identical(rownames(s$theta), "theta_resp")
  expect_identical(rownames(summary(fit, type = "link")$theta), "theta")
  expect_identical(nrow(fitted(fit)), nrow(AS))

  expect_no_error(plot(fit))
  expect_no_error(forest_plot(fit))
})

test_that("switching S2Z off leaves the reported quantities unchanged", {
  ## the switch changes the sampling geometry, not the model, so the two
  ## posteriors must agree up to Monte Carlo error
  on_fit <- s2z_fixture_fit("switch_on")
  off_fit <- s2z_fixture_fit("switch_off")

  expect_true(s2z_active(on_fit))
  expect_false(s2z_active(off_fit))

  for (v in c("beta[1]", "tau[1]", "theta_pred")) {
    on_draws <- posterior::extract_variable(on_fit$draws, v)
    off_draws <- posterior::extract_variable(off_fit$draws, v)
    se <- sqrt(
      posterior::mcse_mean(on_draws)^2 + posterior::mcse_mean(off_draws)^2
    )
    expect_lt(abs(mean(on_draws) - mean(off_draws)) / se, 4, label = v)
  }

  ## and the reported draws layout is identical either way
  expect_identical(
    s2z_reported_variables(on_fit),
    s2z_reported_variables(off_fit)
  )
})

test_that("the adapt_delta default reaches the sampler", {
  ## `.gmap_sampler_control()` is unit tested above; this checks the wiring,
  ## so it fails if the value stops reaching the sampler: a higher target
  ## acceptance rate must produce a smaller adapted step size. Step sizes are
  ## only comparable within one parametrization, so both fixtures are s2z
  ## fits differing solely in adapt_delta.
  step_size <- function(nm) {
    mean(s2z_fixture_fit(nm)$draws_diag[, , "stepsize__"])
  }

  expect_lt(step_size("stepsize_ad99"), step_size("stepsize_ad95"))
})
