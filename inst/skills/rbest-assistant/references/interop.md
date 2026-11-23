# Interoperability: brms, Posterior Draws, JSON, Multivariate

Use when the prior has to leave RBesT, when raw draws are needed, or when the
endpoint is multivariate.

## Mixture priors in brms via `mixstanvar`

`mixstanvar(..., verbose = FALSE)` emits a `stanvars` object carrying the mixture
parameters plus the matching Stan density. Four densities are generated, and the
**argument names matter** because they form the variable names:

| Density | Arguments |
|---|---|
| `mixbeta(w, a, b)` | weights, shapes, shapes |
| `mixgamma(w, a, b)` | weights, shapes, inverse scales |
| `mixnorm(w, m, s)` | weights, means, standard deviations |
| `mixmvnorm(w, m, sigma_L)` | weights, means, Cholesky factors |

Each mixture passed to `mixstanvar` gets a label from its argument name, and the
prior is declared as `label_argument`. So `mixstanvar(map = map_AS)` is referred
to in a prior as `mixnorm(map_w, map_m, map_s)`.

Two rules that are easy to miss:

1. brms regressions work on the **link scale**, so a beta mixture must be moved
   to the logit scale first — sample, transform, refit.
2. brms centers covariates by default, which shifts the intercept and silently
   invalidates an intercept prior. Pass `center = FALSE` to `brms::bf()`.

```r
library(RBesT)
map_AS_beta <- mixbeta(c(0.62, 19.2, 57.8), c(0.38, 3.5, 9.4))
set.seed(3453)
map_AS <- mixfit(logit(rmix(map_AS_beta, 1e4)), type = "norm", Nc = 2)
print(map_AS)
```

The brms half needs `brms` and `glue` and compiles a Stan model, so it is not
executed by the snippet verifier:

```r
# verify: skip
library(brms)
trial <- data.frame(
  n = c(6, 24), r = c(1, 15),
  arm = factor(c("placebo", "secukinumab"))
)
model <- bf(r | trials(n) ~ 1 + arm, family = binomial, center = FALSE)
model_prior <- prior(mixnorm(map_w, map_m, map_s), coef = Intercept) +
  prior(normal(0, 2), class = b)
analysis <- brm(model,
  data = trial, prior = model_prior,
  stanvars = mixstanvar(map = map_AS),
  seed = 365634, refresh = 0
)
hypothesis(analysis, "armsecukinumab > 0")
```

Use `get_prior(model, data = trial)` to find the exact parameter names before
declaring priors. Truncated mixture priors work: the adapter also generates the
distribution functions brms needs for truncation.

## Posterior draws from a `gMAP` fit

Since 1.10-0 a `gMAP` object stores `posterior` draws objects with aligned
sampler diagnostics. The raw `rstan` `stanfit` `fit` slot **was removed** — code
doing `map$fit` or `rstan::extract(map$fit, ...)` no longer works. Use the
`as_draws*()` accessors:

```r
set.seed(34563)
map_bin <- gMAP(cbind(r, n - r) ~ 1 | study,
  family = binomial, data = AS,
  tau.dist = "HalfNormal", tau.prior = 1, beta.prior = 2
)

posterior::variables(as_draws(map_bin))
tau_draws <- posterior::as_draws_matrix(map_bin, variable = "tau")
str(as_draws_rvars(map_bin))
```

`as_draws`, `as_draws_matrix`, `as_draws_array`, `as_draws_df`, `as_draws_list`
and `as_draws_rvars` are all available. Warmup draws are dropped by default
(`options(RBesT.MC.save_warmup = TRUE)` keeps them, and then
`inc_warmup = TRUE` retrieves them). Since 1.12-0 the variable order is `theta`,
`tau`, `beta`, `theta_pred`, `theta_resp_pred`, `lp__` — names and meanings are
unchanged, so only code indexing by column position is affected.

`as.matrix()` still works on `gMAP` and prediction objects for a plain matrix of
draws.

## Saving and reloading mixtures

```r
f <- tempfile(fileext = ".json")
write_mix_json(map_AS_beta, f, digits = 8)
identical(class(read_mix_json(f)), class(map_AS_beta))
unlink(f)
```

`write_mix_json` defaults to 4 digits and says so; it warns when that rounds a
component weight to zero, and errors if every weight rounds to zero. Set
`digits` explicitly. Reading keeps zero-weight components.

The file is the record of the pre-specified prior. From then on work with the
mixture read back from it — for the design evaluation as much as for the final
analysis — and not with the object it was written from.

## Multivariate normal mixtures

`mixmvnorm` supports a mean/SD/correlation parametrization through
`msr2mvnorm`, named dimensions, and components with zero standard deviations:

```r
bi <- mixmvnorm(
  c(1, c(0, 0), c(1, 2), 0.5),
  sigma = msr2mvnorm(s = c(1, 2), r = 0.5, unlist = FALSE)$s,
  param = "msr"
)
print(bi)
```

Support for multivariate workflows is thinner than for the univariate ones — if
the request goes beyond the documented examples, see
[escalation](./escalation.md).

## Next steps

[version-notes](./version-notes.md) when old code breaks,
[api-gotchas](./api-gotchas.md) for signatures.
