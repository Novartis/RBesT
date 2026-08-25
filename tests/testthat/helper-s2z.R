## Helpers for the sum-to-zero (S2Z) reparametrization tests.
##
## See design/design-sum-to-zero-gmap.md. These helpers are intentionally
## self-contained (they do not source anything from design/, which is not part
## of the package build).
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
#' sweep (extreme response, no-rescale, legacy equivalence) lives in the
#' equivalence harness under `design/s2z/`.
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
#' Returns the raw parameters (`beta_raw`, `tau_raw`, `xi_eta`, `xi_abar`) and
#' the transformed ones (`theta`, `beta`, `tau`) regardless of
#' `RBesT.verbose`, which only controls what a *sampled* fit retains. The
#' generated quantities `theta_pred` / `theta_resp_pred` are RNG driven and
#' must not be asserted on.
s2z_constrained <- function(name, upars) {
  rstan::constrain_pars(s2z_stanfit(name), upars)
}

#' A deterministic sweep of unconstrained parameter vectors for a scenario.
#'
#' The unconstrained layout is `c(beta_raw, tau_raw, xi_eta, xi_abar)`. `tau`
#' is swept over several orders of magnitude, which is exactly where a wrong
#' widening of the intercept prior or a rescaled subspace basis would show up,
#' and is unreachable from posterior draws.
s2z_upars_grid <- function(name, n_tau = 7L, seed = 9871L) {
  d <- s2z_test_skeleton(name)$fit.data
  n_up <- rstan::get_num_upars(s2z_stanfit(name))
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

#' Load a cached MCMC fixture for one of the sampled S2Z tests.
#'
#' Built from the committed recipes in `tests/testthat/fixtures-mcmc-src/` by
#' `make -j4 test-fixtures`. The zeta and switch fixtures are sampled with
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
#' @param upars unconstrained parameter vector
#'   `c(beta_raw, tau_raw, xi_eta, xi_abar)`
s2z_reference_log_prob <- function(d, upars) {
  mX <- d$mX
  J <- d$n_groups
  n_strata <- d$n_tau_strata
  n_re <- J - 1L

  beta_raw <- upars[seq_len(mX)]
  tau_raw <- upars[mX + seq_len(n_strata)]
  xi_eta <- upars[mX + n_strata + seq_len(n_re)]
  xi_abar <- upars[mX + n_strata + n_re + 1L]

  g1 <- d$beta_raw_guess[1, ]
  g2 <- d$beta_raw_guess[2, ]
  beta <- g1 + g2 * beta_raw

  tau <- if (d$tau_prior_dist == -1) {
    d$tau_prior[, 1]
  } else {
    exp(d$tau_raw_guess[1] + d$tau_raw_guess[2] * tau_raw)
  }

  m1 <- d$beta_prior[1, 1]
  s1 <- d$beta_prior[1, 2]
  alpha <- beta[1]
  sd_alpha <- s2z_recovery(s1, tau[1], J)$sd_alpha

  xi <- if (d$ncp == 1) tau[1] * xi_eta else g2[1] * xi_eta
  re <- as.numeric(s2z_helmert_basis(J) %*% xi)

  beta_alpha <- beta
  beta_alpha[1] <- alpha
  theta <- as.numeric(d$X %*% beta_alpha) + re[d$group_index]

  lp <- stats::dnorm(xi_abar, 0, 1, log = TRUE)
  lp <- lp + sum(stats::dnorm(
    xi_eta,
    0,
    if (d$ncp == 1) 1 else tau[1] / g2[1],
    log = TRUE
  ))
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
    lp <- lp + sum(d$tau_raw_guess[2] * tau_raw)
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
  scale <- if (d$ncp == 1) tau1 else d$beta_raw_guess[2, 1]
  xi <- matrix(as.numeric(xi_eta) * scale, ncol = 1)
  as.numeric(s2z_helmert_basis(d$n_groups) %*% xi)
}

#' The variables a non-verbose fit would report.
#'
#' Verbose fits keep the raw parameters in the draws; dropping them must leave
#' exactly the legacy variable set, in the legacy order.
s2z_reported_variables <- function(fit) {
  grep(
    "^(xi_eta|xi_abar|beta_raw|tau_raw)",
    posterior::variables(fit$draws),
    value = TRUE,
    invert = TRUE
  )
}
