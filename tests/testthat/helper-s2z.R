## Helpers for the sum-to-zero (S2Z) reparametrization tests.
##
## These helpers are intentionally self-contained (they do not source
## anything from outside of the package build).
##
## The structural invariants of the reparametrization hold for *any* parameter
## vector, not only for posterior draws. They are therefore checked on
## `chains = 0` skeletons via `rstan::constrain_pars()`: deterministic, exact,
## and able to reach extreme `tau` values that no posterior ever visits.
## Sampling is reserved for the few assertions that are genuinely about a
## distribution; those are backed by cached MCMC fixtures.

#' Orthonormal Helmert basis of the zero-sum subspace.
#'
#' Independent R implementation of the Stan `zero_sum_basis()` function, so the
#' tests verify the Stan code rather than restate it.
#'
#' @param J number of groups
#' @return `J x (J-1)` matrix with `Q'Q = I`, `1'Q = 0`
s2z_helmert_basis <- function(J) {
  Q <- matrix(0, nrow = J, ncol = J - 1)
  for (k in seq_len(J - 1)) {
    s <- 1 / sqrt(k * (k + 1))
    Q[seq_len(k), k] <- s
    Q[k + 1, k] <- -k * s
  }
  Q
}

s2z_strata_data <- function() {
  transform(
    RBesT::AS,
    stratum = factor(rep(c("A", "B"), length.out = nrow(RBesT::AS)))
  )
}

s2z_country_data <- function() {
  transform(
    RBesT::transplant,
    country = cut(1:11, c(0, 5, 8, Inf), c("CH", "US", "DE"))
  )
}

#' Scenarios exercised by the S2Z invariant tests.
#'
#' Every entry costs one `chains = 0` model instantiation, not an MCMC run, so
#' the table can cover the full contract: all three endpoints, both
#' parametrizations, fixed tau, covariates, a non-zero prior mean, the `J = 1`
#' and `J = 3` edges, and the documented fallbacks. The wide *sampled* scenario
#' sweep (extreme response, no-rescale, legacy equivalence) is not part of
#' the test suite.
#'
#' `s2z` records whether the scenario is expected to activate the S2Z path.
#' The negative entries are ordinary documented arguments, not opt-ins, so the
#' fallback has to be asserted rather than assumed. Dropping `re_dist == 0`
#' would be a statistical error rather than a slow sampler -- a vector of iid
#' Student-t values is not spherically symmetric, so the common shift is not
#' data-free and marginalising it is simply wrong.
s2z_test_scenarios <- function() {
  binomial_args <- function(...) {
    ## NB: plain `[<-` rather than utils::modifyList(), which recurses into
    ## the data.frame and would splice columns instead of replacing it
    base <- list(
      formula = cbind(r, n - r) ~ 1 | study,
      data = RBesT::AS,
      family = binomial,
      tau.dist = "HalfNormal",
      tau.prior = 0.5,
      beta.prior = 2
    )
    override <- list(...)
    base[names(override)] <- override
    base
  }

  list(
    binomial_ncp = list(
      s2z = TRUE,
      opts = list(RBesT.MC.ncp = 1),
      args = function() binomial_args()
    ),
    binomial_cp = list(
      s2z = TRUE,
      opts = list(RBesT.MC.ncp = 0),
      args = function() binomial_args()
    ),
    normal_ncp = list(
      s2z = TRUE,
      opts = list(RBesT.MC.ncp = 1),
      args = function() {
        crohn <- RBesT::crohn
        crohn$y.se <- 88 / sqrt(crohn$n)
        list(
          formula = cbind(y, y.se) ~ 1 | study,
          data = crohn,
          family = gaussian,
          tau.dist = "HalfNormal",
          tau.prior = 44,
          beta.prior = 88
        )
      }
    ),
    poisson_ncp = list(
      s2z = TRUE,
      opts = list(RBesT.MC.ncp = 1),
      args = function() {
        list(
          formula = y ~ 1 + offset(log(n)) | study,
          data = data.frame(
            y = c(11L, 4L, 8L, 21L, 6L, 15L),
            n = c(120, 55, 92, 210, 71, 160),
            study = factor(paste0("Study ", 1:6))
          ),
          family = poisson,
          tau.dist = "HalfNormal",
          tau.prior = 0.5,
          beta.prior = 2
        )
      }
    ),
    fixed_tau = list(
      s2z = TRUE,
      opts = list(),
      args = function() binomial_args(tau.dist = "Fixed", tau.prior = 0.25)
    ),
    covariate = list(
      s2z = TRUE,
      opts = list(),
      args = function() {
        list(
          formula = cbind(r, n - r) ~ 1 + country | study,
          data = s2z_country_data(),
          family = binomial,
          tau.dist = "HalfNormal",
          tau.prior = 1,
          beta.prior = rbind(c(0, 2), c(0, 1), c(0, 1))
        )
      }
    ),
    nonzero_prior_mean = list(
      s2z = TRUE,
      opts = list(),
      args = function() binomial_args(beta.prior = rbind(c(-1.5, 2)))
    ),
    single_group = list(
      s2z = TRUE,
      opts = list(),
      args = function() binomial_args(data = RBesT::AS[1, ])
    ),
    three_groups = list(
      s2z = TRUE,
      opts = list(),
      args = function() binomial_args(data = RBesT::AS[1:3, ])
    ),
    ## documented fallbacks: the negative side of the S2Z predicate
    fallback_student_t = list(
      s2z = FALSE,
      opts = list(),
      args = function() binomial_args(REdist = "t", t.df = 5)
    ),
    fallback_tau_strata = list(
      s2z = FALSE,
      opts = list(),
      args = function() {
        binomial_args(
          data = s2z_strata_data(),
          tau.strata = quote(stratum),
          tau.prior = c(0.5, 0.5)
        )
      }
    ),
    fallback_no_intercept = list(
      s2z = FALSE,
      opts = list(),
      args = function() {
        list(
          formula = cbind(r, n - r) ~ 0 + country | study,
          data = s2z_country_data(),
          family = binomial,
          tau.dist = "HalfNormal",
          tau.prior = 1,
          beta.prior = rbind(c(0, 2), c(0, 2), c(0, 2))
        )
      }
    ),
    ## the user escape hatch; not a correctness condition, but the layout must
    ## match the legacy one exactly
    optout = list(
      s2z = FALSE,
      opts = list(RBesT.MC.s2z = FALSE),
      args = function() binomial_args()
    ),
    optout_normal = list(
      s2z = FALSE,
      opts = list(RBesT.MC.s2z = FALSE),
      args = function() {
        crohn <- RBesT::crohn
        crohn$y.se <- 88 / sqrt(crohn$n)
        list(
          formula = cbind(y, y.se) ~ 1 | study,
          data = crohn,
          family = gaussian,
          tau.dist = "HalfNormal",
          tau.prior = 44,
          beta.prior = 88
        )
      }
    )
  )
}

.s2z_cache <- new.env(parent = emptyenv())

.s2z_memoise <- function(key, value) {
  if (is.null(.s2z_cache[[key]])) {
    .s2z_cache[[key]] <- value
  }
  .s2z_cache[[key]]
}

#' Draw-free `gMAP()` skeleton for a scenario, memoised.
#'
#' `chains = 0` runs the whole model setup -- the design matrix, the prior
#' encoding and the S2Z predicate -- without any sampling, so `fit.data` is the
#' very data list a real fit would use.
s2z_test_skeleton <- function(name) {
  .s2z_memoise(paste0("skeleton_", name), {
    scenario <- s2z_test_scenarios()[[name]]
    stopifnot(!is.null(scenario))
    withr::with_options(scenario$opts, {
      suppressMessages(suppressWarnings(do.call(
        RBesT::gMAP,
        c(scenario$args(), list(chains = 0))
      )))
    })
  })
}

#' Build a sampler-free Stan fit object usable with `rstan::log_prob()`.
#'
#' `gMAP()` deliberately does not return the `stanfit`, so the target is
#' accessed by re-instantiating the compiled model with the very same data list
#' the fit used.
s2z_log_prob_fit <- function(fit_data) {
  suppressMessages(rstan::sampling(
    RBesT:::stanmodels$gMAP,
    data = fit_data,
    chains = 0
  ))
}

#' Compiled Stan model instance for a scenario, memoised.
s2z_stanfit <- function(name) {
  .s2z_memoise(
    paste0("stanfit_", name),
    s2z_log_prob_fit(s2z_test_skeleton(name)$fit.data)
  )
}

#' Constrain an unconstrained parameter vector through the Stan model.
#'
#' Returns the raw parameters (`beta_raw`, `tau_raw`, `xi_eta`) and the
#' transformed ones (`theta`, `beta_param`, `tau`) regardless of
#' `RBesT.verbose`, which only controls what a *sampled* fit retains.
#'
#' Note that the *super-population* `beta` is a generated quantity -- the
#' recovery of the common shift is an RNG draw -- and `constrain_pars()` does
#' not run `generated quantities`. Neither `beta` nor the RNG driven
#' `theta_pred` / `theta_resp_pred` can be asserted on here; under s2z
#' `beta_param[1]` is the sampled intercept alpha.
s2z_constrained <- function(name, upars) {
  rstan::constrain_pars(s2z_stanfit(name), upars)
}

#' A deterministic sweep of unconstrained parameter vectors for a data list.
#'
#' The unconstrained layout is `c(beta_raw, tau_raw, xi_eta)`, where `xi_eta`
#' has `J - 1` entries under s2z and `J` otherwise. `tau` is swept over several
#' orders of magnitude, which is exactly where a wrong widening of the
#' intercept prior, a rescaled subspace basis or a missing tau-dependent
#' normaliser would show up, and is unreachable from posterior draws.
gmap_upars_grid <- function(d, fit, n_tau = 7L, seed = 9871L) {
  n_up <- rstan::get_num_upars(fit)
  tau_slots <- d$mX + seq_len(d$n_tau_strata)

  withr::with_seed(seed, {
    grid <- lapply(seq(-3, 3, length.out = n_tau), function(tr) {
      u <- stats::rnorm(n_up, 0, 0.5)
      u[tau_slots] <- tr
      u
    })
    c(list(rep(0, n_up)), grid)
  })
}

#' The same sweep for one of the named scenarios.
s2z_upars_grid <- function(name, n_tau = 7L, seed = 9871L) {
  gmap_upars_grid(
    s2z_test_skeleton(name)$fit.data,
    s2z_stanfit(name),
    n_tau = n_tau,
    seed = seed
  )
}

#' Load a cached MCMC fixture for one of the sampled S2Z tests.
#'
#' Built from the committed recipes in `tests/testthat/fixtures-mcmc-src/` by
#' `make -j4 test-fixtures`. The switch fixtures are sampled with
#' `RBesT.verbose = TRUE` so that the raw parameters stay in the draws.
s2z_fixture_fit <- function(name) {
  load_gmap_fixture(paste0("gmap_s2z_", name), type = "mcmc")
}

#' R implementation of the S2Z target, used to check the Stan target.
#'
#' Mirrors the `transformed parameters` and `model` blocks of the s2z branch.
#' It follows the same branch structure and the same Helmert formula, so it is
#' not an independent derivation and cannot vouch for the derivation itself --
#' see the prior-marginal test for that. What it does do, and what it exists
#' for, is pin down the normalising constants: it includes all of them, Stan
#' drops the parameter-independent ones, so comparisons must be made on
#' *differences* of the target across parameter vectors, which is exactly what
#' exposes a dropped tau-dependent constant.
#'
#' @param d the Stan data list (`fit$fit.data`)
#' @param upars unconstrained parameter vector `c(beta_raw, tau_raw, xi_eta)`
s2z_reference_log_prob <- function(d, upars) {
  mX <- d$mX
  J <- d$n_groups
  n_strata <- d$n_tau_strata
  n_re <- J - 1L

  beta_raw <- upars[seq_len(mX)]
  tau_raw <- upars[mX + seq_len(n_strata)]
  xi_eta <- upars[mX + n_strata + seq_len(n_re)]

  g1 <- d$beta_raw_guess[1, ]
  g2 <- d$beta_raw_guess[2, ]
  beta <- g1 + g2 * beta_raw

  tau <- if (d$tau_prior_dist == -1) {
    d$tau_prior[, 1]
  } else {
    exp(d$tau_raw_guess[1, ] + d$tau_raw_guess[2, ] * tau_raw)
  }

  m1 <- d$beta_prior[1, 1]
  s1 <- d$beta_prior[1, 2]
  alpha <- beta[1]
  sd_alpha <- s2z_recovery(s1, tau[1], J)$sd_alpha

  ## partial-centering model parametrization; an all-zero `re_center` is the
  ## plain non-centered/centered pair the shipped program used before the
  ## option existed. Mirrors the `transformed parameters` block of gMAP.stan.
  center <- if (is.null(d$re_center)) rep(0, J) else as.numeric(d$re_center)
  plain <- max(center) == 0

  qxi <- as.numeric(s2z_helmert_basis(J) %*% xi_eta)
  if (plain) {
    re <- if (d$re_param == 1) tau[1] * qxi else g2[1] * qxi
  } else {
    group_g <- d$group_scale_guess
    group_loc <- if (is.null(d$group_location_guess)) {
      rep(0, J)
    } else {
      d$group_location_guess
    }
    sc <- 1 - center + center * (tau[1] / group_g)
    center_location <- -center * (group_loc / group_g)
    w <- tau[1] * (qxi - center_location) / sc
    re <- w - mean(w)
  }

  beta_alpha <- beta
  beta_alpha[1] <- alpha
  theta <- as.numeric(d$X %*% beta_alpha) + re[d$group_index]

  if (plain) {
    lp <- sum(stats::dnorm(
      xi_eta,
      0,
      if (d$re_param == 1) 1 else tau[1] / g2[1],
      log = TRUE
    ))
  } else {
    ## under partial centering the density is stated on the physical group
    ## effects and the model parametrization enters through the log-Jacobian,
    ## whose
    ## determinant is det(Q' diag(d) Q) = (prod d_j)(sum 1/d_j)/J, computed
    ## directly from mean(sc) (exact for heterogeneous sc, i.e. a per-group
    ## group_scale_guess)
    lp <- -0.5 * sum(re^2) / tau[1]^2 -
      sum(log(sc)) +
      log(mean(sc)) -
      0.5 * n_re * log(2 * pi)
  }
  ## the widened intercept prior; its sd depends on tau, so its normalising
  ## constant must not be dropped
  lp <- lp + stats::dnorm(alpha, m1, sd_alpha, log = TRUE)
  if (mX > 1) {
    lp <- lp + sum(stats::dnorm(
      beta[-1],
      d$beta_prior[-1, 1],
      d$beta_prior[-1, 2],
      log = TRUE
    ))
  }

  ## tau prior and the Jacobian of the log transform (unchanged by s2z; only
  ## the families actually exercised here are implemented)
  if (d$tau_prior_dist == -1) {
    lp <- lp + sum(stats::dnorm(tau_raw, 0, 1, log = TRUE))
  } else if (d$tau_prior_dist == 0) {
    lp <- lp + sum(stats::dnorm(tau, 0, d$tau_prior[, 2], log = TRUE))
    lp <- lp + sum(d$tau_raw_guess[2, ] * tau_raw)
  } else {
    stop("tau prior family not implemented in the reference target")
  }

  if (d$prior_PD == 0) {
    lp <- lp +
      switch(
        d$link,
        sum(stats::dnorm(d$y, theta, d$y_se, log = TRUE)),
        sum(stats::dbinom(d$r, d$r_n, stats::plogis(theta), log = TRUE)),
        sum(stats::dpois(d$count, exp(d$log_offset + theta), log = TRUE))
      )
  }

  lp
}

#' Group index to tau stratum mapping, as Stan builds it in `transformed data`.
#'
#' Groups without data keep the prediction stratum.
gmap_tau_strata_gindex <- function(d) {
  idx <- rep(d$tau_strata_pred, d$n_groups)
  for (i in seq_len(d$H)) {
    idx[d$group_index[i]] <- d$tau_strata_index[i]
  }
  idx
}

#' Location-scale log density of the random effect distribution.
gmap_re_ldens <- function(x, location, scale, re_dist, t_df) {
  if (re_dist == 0) {
    stats::dnorm(x, location, scale, log = TRUE)
  } else {
    stats::dt((x - location) / scale, df = t_df, log = TRUE) - log(scale)
  }
}

#' R reference target for the *conventional* (non-S2Z) parametrization.
#'
#' Deliberately written in the physical coordinates, which is a different
#' formulation from the Stan program: here the group effects
#' \eqn{eps_j \sim D(0, \tau_j)} carry the density and the sampled coordinates
#' enter only through the model parametrization
#'
#'   eps_j = tau_j (xi_j - m_j) / sc_j,
#'   sc_j  = (1 - c_j) + c_j tau_j / g,   m_j = c_j (beta[1] - anchor) / g,
#'
#' plus its log-Jacobian `sum log(tau_j / sc_j)`. gMAP.stan instead states the
#' location-scale density `xi_j ~ D(m_j, sc_j)` directly on the sampled
#' coordinate, so the two agree only if the Stan location, scale *and* the
#' implied Jacobian are all right. `c = 0` reduces to the shipped non-centered
#' model parametrization and `c = 1` to the shipped centered model
#' parametrization, so the same reference pins down the two fast paths as well.
#'
#' As with `s2z_reference_log_prob()` all normalising constants are kept and
#' Stan drops the parameter-independent ones, so only *differences* of the
#' target across parameter vectors are comparable.
#'
#' @param d the Stan data list (`fit$fit.data`)
#' @param upars unconstrained parameter vector `c(beta_raw, tau_raw, xi_eta)`
#' @param drop deliberately break one ingredient of the parametrization.
#'   Exists so the tests can demonstrate that they are sensitive to the
#'   location (`"loc"`), to the scale (`"scale"`) and to the log-Jacobian
#'   (`"jacobian"`); a passing comparison against `"none"` alone would not
#'   prove that.
gmap_conv_reference_log_prob <- function(
  d,
  upars,
  drop = c("none", "loc", "scale", "jacobian")
) {
  drop <- match.arg(drop)
  mX <- d$mX
  J <- d$n_groups
  n_strata <- d$n_tau_strata

  beta_raw <- upars[seq_len(mX)]
  tau_raw <- upars[mX + seq_len(n_strata)]
  xi_eta <- upars[mX + n_strata + seq_len(J)]

  anchor <- d$beta_raw_guess[1, ]
  g <- d$beta_raw_guess[2, ]
  beta <- anchor + g * beta_raw

  tau <- if (d$tau_prior_dist == -1) {
    d$tau_prior[, 1]
  } else {
    exp(d$tau_raw_guess[1, ] + d$tau_raw_guess[2, ] * tau_raw)
  }
  tau_group <- tau[gmap_tau_strata_gindex(d)]

  center <- if (is.null(d$re_center)) rep(0, J) else as.numeric(d$re_center)
  group_g <- d$group_scale_guess
  group_loc <- if (is.null(d$group_location_guess)) rep(0, J) else d$group_location_guess
  sc <- 1 - center + center * (tau_group / group_g)
  location <- if (d$has_intercept == 1) beta[1] else 0
  center.anchor <- if (d$has_intercept == 1) anchor[1] else 0
  m <- center * ((location - center.anchor - group_loc) / group_g)
  if (drop == "scale") {
    sc <- rep(1, J)
  }
  if (drop == "loc") {
    m <- rep(0, J)
  }

  eps <- tau_group * (xi_eta - m) / sc
  theta <- as.numeric(d$X %*% beta) + eps[d$group_index]

  ## the random effects carry the density; the parametrization contributes its
  ## Jacobian
  lp <- sum(gmap_re_ldens(eps, 0, tau_group, d$re_dist, d$re_dist_t_df))
  if (drop != "jacobian") {
    lp <- lp + sum(log(tau_group / sc))
  }

  lp <- lp +
    sum(stats::dnorm(beta, d$beta_prior[, 1], d$beta_prior[, 2], log = TRUE))

  if (d$tau_prior_dist == -1) {
    lp <- lp + sum(stats::dnorm(tau_raw, 0, 1, log = TRUE))
  } else if (d$tau_prior_dist == 0) {
    lp <- lp + sum(stats::dnorm(tau, 0, d$tau_prior[, 2], log = TRUE))
    lp <- lp + sum(d$tau_raw_guess[2, ] * tau_raw)
  } else {
    stop("tau prior family not implemented in the reference target")
  }

  if (d$prior_PD == 0) {
    lp <- lp +
      switch(
        d$link,
        sum(stats::dnorm(d$y, theta, d$y_se, log = TRUE)),
        sum(stats::dbinom(d$r, d$r_n, stats::plogis(theta), log = TRUE)),
        sum(stats::dpois(d$count, exp(d$log_offset + theta), log = TRUE))
      )
  }

  lp
}

#' Is a fit on the S2Z code path?
s2z_active <- function(fit) {
  d <- fit$fit.data
  ## use_s2z is absent from data lists built before the opt-out switch
  ## existed, where the path was unconditional
  use_s2z <- if (is.null(d[["use_s2z"]])) 1L else d[["use_s2z"]]
  use_s2z == 1 && d$re_dist == 0 && d$n_tau_strata == 1 && d$has_intercept == 1
}

#' Recovery coefficients of the marginalised common shift.
#'
#' The implementation uses `sd_alpha = hypot(s1, tau/sqrt(J))` and
#' `r = (tau/sqrt(J)) / sd_alpha`, so that `r^2 = v / (s1^2 + v)` and
#' `s1 * r = sqrt(s1^2 v / (s1^2 + v))` with `v = tau^2 / J`.
s2z_recovery <- function(s1, tau, J) {
  sd_a <- tau / sqrt(J)
  sd_alpha <- sqrt(s1^2 + sd_a^2)
  list(r = sd_a / sd_alpha, sd_alpha = sd_alpha)
}

#' Reconstruct the sampled intercept `alpha = beta[1] + mean(eps)`.
s2z_alpha <- function(d, beta_raw1) {
  g <- d$beta_raw_guess
  as.numeric(g[1, 1] + g[2, 1] * beta_raw1)
}

#' Reconstruct the group random effects `re = Q %*% xi`.
s2z_re <- function(d, xi_eta, tau1) {
  J <- d$n_groups
  g <- d$beta_raw_guess[2, 1]
  center <- if (is.null(d$re_center)) rep(0, J) else as.numeric(d$re_center)
  qxi <- as.numeric(
    s2z_helmert_basis(J) %*% matrix(as.numeric(xi_eta), ncol = 1)
  )
  if (max(center) == 0) {
    return(qxi * if (d$re_param == 1) tau1 else g)
  }
  ## partial-centering parametrization, mirroring gMAP.stan: per-group scale
  ## numeraire and posterior location guess, then re-project onto the
  ## zero-sum subspace
  group_g <- d$group_scale_guess
  group_loc <- if (is.null(d$group_location_guess)) rep(0, J) else d$group_location_guess
  sc <- 1 - center + center * (tau1 / group_g)
  center_location <- -center * (group_loc / group_g)
  w <- tau1 * (qxi - center_location) / sc
  w - mean(w)
}

#' The variables a non-verbose fit would report.
#'
#' Verbose fits keep the sampled parameters (including `beta_param`, whose
#' first entry is the s2z intercept alpha) in the draws; dropping them must
#' leave exactly the reported variable set, in the reported order.
s2z_reported_variables <- function(fit) {
  ## s2z_quad/s2z_log_det are declared in transformed parameters (the s2z
  ## partial-centering log-Jacobian pieces) and are excluded from a
  ## non-verbose fit's draws via the same `exclude_pars` list as
  ## xi_eta/beta_raw/tau_raw/beta_param (see R/gMAP.R); a verbose fixture
  ## retains them, so they must be filtered here too. Pre-existing gap,
  ## unrelated to group_location_guess: verbose fixtures were never
  ## available to exercise this helper until now.
  grep(
    "^(xi_eta|beta_raw|tau_raw|beta_param|s2z_quad|s2z_log_det)",
    posterior::variables(fit$draws),
    value = TRUE,
    invert = TRUE
  )
}
