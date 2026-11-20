# Analysing Trial Data

Use when data exist — final or interim — and the task is to update a prior, form
a predictive distribution, compare two arms, or extract shrinkage estimates.

## Conjugate updating with `postmix`

`postmix(priormix, data, ...)` returns the posterior mixture. Pass either the raw
data vector or summary statistics:

| Class | Summary form | Notes |
|---|---|---|
| `betaMix` | `postmix(prior, r = , n = )` | responders and total |
| `normMix` | `postmix(prior, m = , se = )` | mean and its standard error |
| `normMix` | `postmix(prior, m = , n = )` | `se` derived from `sigma(prior)` |
| `gammaMix` | `postmix(prior, m = , n = )` | mean count and sample size |

```r
library(RBesT)
map_bin <- mixbeta(c(0.7, 18, 133), c(0.3, 3.5, 15))
rob <- robustify(map_bin, weight = 0.2, mean = 0.5)
post_ctrl <- postmix(rob, r = 14, n = 60)
print(post_ctrl)
summary(post_ctrl)
```

For normal endpoints the reference scale matters. Supplying both `n` and `se`
updates `sigma(prior)`; supplying only `n` uses the prior's existing `sigma`.

```r
prior_norm <- mixnorm(c(1, -50, 1), param = "mn", sigma = 88)
postmix(prior_norm, m = -35, se = 88 / sqrt(40))
```

## Predictive distribution

`preddist(mix, n = 1, ...)` turns a parameter distribution into a distribution
for a future observation of size `n` — beta-binomial, normal, gamma-Poisson.
Normal mixtures take `sigma`.

```r
pred_ctrl <- preddist(post_ctrl, n = 40)
summary(pred_ctrl)
```

Use it for predictive power and for "what could the next cohort look like",
never as a parameter distribution.

## Two-arm comparison

The difference of two mixtures has its own d/p/q/r family, all with the same
argument order — `mix1`, `mix2`, then the quantile/probability:

```r
post_trt <- postmix(mixbeta(c(1, 1, 1)), r = 27, n = 60)
pmixdiff(post_trt, post_ctrl, 0, lower.tail = FALSE) # Pr(theta_trt - theta_ctrl > 0)
qmixdiff(post_trt, post_ctrl, c(0.025, 0.975)) # credible interval for the difference
```

`lower.tail = FALSE` gives `Pr(θ₁ − θ₂ > q)`; `TRUE` gives `Pr(θ₁ − θ₂ ≤ q)`. Fix
the direction with the user before writing the call — this is the single most
common source of silently wrong answers.

For the full decision-theoretic machinery — thresholds, boundaries, operating
characteristics — go to [design](./design.md).

## Prediction and shrinkage from a `gMAP` fit

`predict()` on a `gMAP` object returns the posterior for each input group
(shrinkage estimates). Passing new data returns the *predictive* distribution as
if those were new studies.

```r
trans_cov <- transform(
  transplant,
  country = cut(1:11, c(0, 5, 8, Inf), c("CH", "US", "DE"))
)
set.seed(34246)
map_cov <- gMAP(cbind(r, n - r) ~ 1 + country | study,
  data = trans_cov, family = binomial,
  tau.dist = "HalfNormal", tau.prior = 1,
  ## regression coefficients are contrasts, so they need tighter priors
  beta.prior = rbind(c(0, 2), c(0, 1), c(0, 1))
)

shrunk <- predict(map_cov) # one row per input study
summary(shrunk)

pred_new <- predict(map_cov, data.frame(country = "CH", study = 12))
print(pred_new)
```

Notes:

- Predictions are on the **response scale** by default; `summary()` returns a
  matrix, `as.matrix()` the underlying draws.
- New data must carry every covariate in the formula, including the grouping
  column, with levels matching the fit.
- `mixfit()` accepts a prediction object directly, so a per-country MAP prior is
  `automixfit(pred_new)`.

## Sanity checks

1. Did the posterior move where the data said it should, and by a plausible amount?
2. For a robustified prior, how much weight does the vague component retain after
   updating? Large residual weight is the signal of prior-data conflict.
3. Is the reported scale the one the user asked for — response or link?

## Next steps

[design](./design.md) for decisions and operating characteristics,
[prior-tuning](./prior-tuning.md) for the ESS of the updated prior.
