# MAP Prior Derivation with `gMAP`

Use when the objective is to turn historical data into a prior: binary, normal,
Poisson/count, negative-binomial recurrent events, or variances.

## Follow-ups to settle first

1. The model formula, including the study/group identifier after `|`.
2. The endpoint columns and how the sample size is encoded.
3. `family` and, for normal endpoints, the known sampling σ.
4. The heterogeneity prior (`tau.dist`, `tau.prior`) — offer the defaults below.
5. Whether all studies are equally exchangeable, or need `tau.strata`.
6. For counts: whether the historical inputs are rates or counts, i.e. whether an
   exposure `offset` is needed.

## Data encoding by family

| Family | Response | Notes |
|---|---|---|
| binary | `cbind(r, n - r) ~ 1 \| study` | `family = binomial`, logit link; `n` is encoded in the response |
| normal | `cbind(y, y.se) ~ 1 \| study` | `family = gaussian`, identity link; pass `weights = n` |
| Poisson | `r ~ 1 + offset(log(exposure)) \| study` | `family = poisson`, log link; offset must be on the log scale |

For normal endpoints `weights` is *recommended*, not optional in practice: it is
what lets `gMAP` estimate the reference σ that the derived mixture carries. Omit
it and `sigma()` of the mixture is unset, so every downstream `postmix`, `oc*`
and `pos*` call has to be given `sigma` explicitly. If the summaries come without
per-study sample sizes, set it by hand afterwards with `sigma(mix) <- ...`.

## Prior defaults to offer

- **Binary, response rate roughly 20–80%**: `tau.dist = "HalfNormal"` with
  `tau.prior = 1` is the very conservative choice, `tau.prior = 1/2` a
  conservative one.
- **Normal**: `HalfNormal(0, σ/2)` is the very conservative choice, `σ/4` a
  conservative one.
- **Log variance**: sampling SD of a log variance is √2, so `HalfNormal(0, √2/2)`
  is very conservative and `√2/4` conservative.
- With few historical studies the data say little about τ, so this prior
  largely decides how much is borrowed: show the sensitivity. With a single
  study τ is not estimable at all and the prior needs an external rationale.
- **`beta.prior`**: a single number is the SD with mean 0; `cbind(mean, SD)` sets
  both. For the intercept a unit-information prior — precision of a single
  observation — is the usual recommendation.

`tau.prior` shapes: for one stratum and a distribution with a default location
(e.g. `HalfNormal`, whose μ is fixed at 0) a **scalar** is the dispersion. A
two-column form `cbind(mean, SD)` is equally valid. With several strata, pass a
vector of length = number of strata. Watch the collision: with a **single**
stratum a length-two vector is read as the pair (a, b), not as two strata.

## Binary

```r
library(RBesT)
set.seed(34563)
map_bin <- gMAP(cbind(r, n - r) ~ 1 | study,
  family = binomial, data = AS,
  tau.dist = "HalfNormal", tau.prior = 1,
  beta.prior = 2
)
print(map_bin)
```

## Normal

`crohn` reports no standard error, so build it from the known σ:

```r
dat <- crohn
crohn_sigma <- 88
dat$y.se <- crohn_sigma / sqrt(dat$n)

set.seed(689654)
map_norm <- gMAP(cbind(y, y.se) ~ 1 | study,
  weights = n, data = dat, family = gaussian,
  beta.prior = cbind(0, crohn_sigma),
  tau.dist = "HalfNormal", tau.prior = cbind(0, crohn_sigma / 2)
)
print(map_norm)
```

## Count / negative binomial recurrent events

Whether an offset is needed depends on how the history is reported. `asthma`
already carries yearly **log-rates**, so it is modelled as a normal endpoint on
the log-rate scale with no offset. Had the inputs been log mean *counts*
`log(λ·d)`, `offset = log(d)` would recover the log-rate.

```r
asthma_ph3 <- subset(asthma, phase == "phase III")

## per-patient sampling SD on the log-rate scale, sqrt((1 + kappa*lambda*t)/(lambda*t))
lambda_ctrl <- 1.8 # placebo rate, events per patient-year
kappa_true <- 1.9 # overdispersion
followup <- 1 * 0.85 + 0.5 * 0.15 # mean exposure in years, allowing for dropout
sigma_ctrl <- round(
  sqrt((1 + kappa_true * lambda_ctrl * followup) / (lambda_ctrl * followup)), 2
)

set.seed(46456)
map_rate <- gMAP(cbind(log_mu_hat, se_log_mu_hat) ~ 1 | study,
  data = asthma_ph3, family = gaussian,
  tau.dist = "HalfNormal", tau.prior = sigma_ctrl / 4,
  beta.prior = 2
)
map_rate_mix <- automixfit(map_rate)
sigma(map_rate_mix) <- sigma_ctrl
print(map_rate_mix)
```

A genuine Poisson fit instead uses the log-offset form:
`gMAP(events ~ 1 + offset(log(exposure_y)) | study, family = poisson, ...)`.

## Variance MAP

RBesT has no Γ likelihood, so the documented route approximates the log-variance
as normal. With `s²` on `df` degrees of freedom, `α = df/2`, `β = α/s²`:

```r
hdata <- data.frame(
  study = 1:6,
  sd = c(12.11, 10.97, 10.94, 9.41, 10.97, 10.95),
  df = c(597, 60, 548, 307, 906, 903)
)
hdata$alpha <- hdata$df / 2
hdata$beta <- hdata$alpha / hdata$sd^2
hdata$logvar_mean <- log(hdata$sd^2 * hdata$alpha) - digamma(hdata$alpha)
hdata$logvar_var <- psigamma(hdata$alpha, 1)

set.seed(45356)
map_var <- gMAP(cbind(logvar_mean, sqrt(logvar_var)) ~ 1 | study,
  data = hdata,
  tau.dist = "HalfNormal", tau.prior = sqrt(2) / 2,
  beta.prior = cbind(4.8, 100)
)
print(map_var)
```

The normal approximation is useful from about 10 degrees of freedom and very
good beyond 20. Below that, escalate.

## Differential discounting (`tau.strata`)

When some studies are less relevant (observational vs randomized, say), give
them their own heterogeneity parameter:

- `tau.strata` assigns each study to a stratum.
- `tau.prior` then needs one entry per stratum, e.g. `c(0.5, 1)` for two strata.
- `tau.strata.pred` (default `1`) picks the stratum the MAP prediction targets.

## Covariates (meta-regression)

`gMAP` accepts study-level covariates on the right-hand side of the formula,
e.g. `cbind(r, n - r) ~ 1 + country | study`; a runnable example with the
prediction for a new study is in [analysis](./analysis.md). Use it when a
characteristic that drives the control outcome is known for **every**
historical study and for the new trial. Limits to state:

- The regression coefficients are common to all studies; the only random
  effect is the study intercept.
- Each coefficient needs its own row in `beta.prior`.
- With summary data a covariate describes a whole study or arm (a region, a
  mean age), not individual patients.

Anything beyond this — adjustment for patient-level baseline characteristics,
subgroup summaries reported by only some studies, a trend or random-walk model
for extrapolation — is an advanced topic. Recommend the worked case studies at
<https://opensource.nibr.com/bamdd/>, where such models are demonstrated, and
[escalate](./escalation.md) rather than constructing one yourself.

## Diagnostics before anything else

```r
pl <- plot(map_bin)
names(pl)
## guards the element names quoted below
stopifnot("forest_model" %in% names(pl), !("tau_forest" %in% names(pl)))
```

`plot.gMAP` returns `traceBeta`, `traceTau`, `traceLogTau`, `densityBeta`,
`densityTau`, `densityLogTau`, `densityThetaStar`, `densityThetaStarLink` and
`forest_model`. There is **no** `tau_forest` element — use `densityTau` or
`traceLogTau` for the heterogeneity, and `forest_model` for the shrinkage vs
stratified comparison. `forest_plot()` is the standalone, customizable version.

Check, in order: no divergences and `Rhat` near 1 in `print(map_bin)`; the τ
posterior against its prior; then the forest plot for shrinkage.

## Approximate with a mixture

```r
map_mix <- automixfit(map_bin)
print(map_mix)
plot(map_mix)$mix
```

`automixfit` searches `Nc = seq(1, 4)` by AIC with penalty `k = 6`. Always look
at `plot(...)$mix` — the density overlay on the sample histogram. If it fits
poorly, override with `mixfit(map_bin, type = "beta", Nc = 3)`.

## Sensitivity and refits

`update()` works on a `gMAP` object through `stats::update.default` and the
stored `$call` — there is no `update.gMAP` method. It **re-runs the full MCMC**,
so it is a refit, not a cheap tweak:

- `update(map_bin, tau.prior = 1/2)` — τ-prior sensitivity.
- `update(map_bin, data = all_data)` — the MAC step when co-data arrives.

## Sampling behaviour in 1.12-0

Defaults changed; the model did not. Relevant when diagnosing a fit:

- Group effects use a **sum-to-zero** parametrization whenever the model has an
  intercept, normal random effects and a single τ stratum. Student-t effects,
  several strata and intercept-free models keep the old parametrization.
  `options(RBesT.MC.s2z = FALSE)` restores the conventional one.
- **Automatic partial centering** is the default (`options(RBesT.MC.ncp = 3)`);
  quadrature estimates per-group centering fractions and `gMAP` warns and falls
  back to the exact non-centered model if its numerical checks disagree.
- Default `adapt_delta` is now **0.95** (was 0.99). On divergent transitions the
  package tells you to raise it back to 0.99 — that is the first thing to try.
- Posterior summaries use modern `posterior` functions, so `Rhat` and friends can
  differ numerically from pre-1.10 runs.

## Next steps

[prior-tuning](./prior-tuning.md) for ESS and robustification,
[analysis](./analysis.md) for using the prior on data,
[design](./design.md) for operating characteristics.
