# Tuning a Prior: Robustification, ESS, Combination, Elicitation

Use when a prior already exists — MAP-derived or hand-specified — and the task is
to characterize it, protect it, combine it, or build it from elicited summaries.

## Robustification

`robustify(priormix, weight, mean, n = 1, ...)` mixes in a vague component to
protect against prior-data conflict. The `weight` is the prior probability given
to that vague component and is a substantive modelling choice, not a default to
pick silently. Ask for it, or show a sensitivity grid.

Starting points to offer, for the user to confirm — they are examples, not
defaults:

- **10–20%** when the historical data are strongly relevant: same disease,
  similar population and study design.
- **30–50%** when relevance is moderate: a different population (adult to
  pediatric, say) or a related disease.

The weight is the *prior* probability of the vague component; the trial data
re-weight the components. `1 - weight` is therefore **not** the amount
borrowed, which depends on the trial data; the operating characteristics show
it. Plan a
sensitivity analysis that moves the weight from the chosen value towards 1 (no
borrowing); the weight at which the conclusion changes is the tipping point.

Set `mean` explicitly — it should sit at the null hypothesis. The defaults are
deliberately unexciting and `robustify` tells you when it falls back on them:
`1/2` for beta mixtures, the prior mean for gamma and normal ones. For normal
mixtures also pass `sigma`.

The vague component must be weakly informative, not flat: about one
observation's worth of information (the default `n = 1`, a unit-information
prior), with its mean chosen explicitly. A normal component with a very large
variance gives practically **no** robustification — the data favour it over the informative
part only under extreme conflict, far too late (Lindley's paradox). Do not widen it "to be safe".

```r
library(RBesT)
map_bin <- mixbeta(c(0.7, 18, 133), c(0.3, 3.5, 15))
rob <- robustify(map_bin, weight = 0.2, mean = 0.5)
print(rob)
```

The beta method uses `n + 1` rather than `n` as the robust component's sample
size, so the default robust part is the uniform `Beta(1, 1)` instead of the
strictly unit-information `Beta(1/2, 1/2)`.

## Effective sample size

The ESS measures how informative a prior is. It is **not** a number of patients
to subtract from, or add to, a sample size or power calculation: what the prior
contributes depends on the control data the trial observes — see
[borrowing-design](./borrowing-design.md).

`ess(mix, method = c("elir", "moment", "morita"))` — `elir` is the default and
the recommended choice. `moment` is conservative; `morita` can be optimistic for
mixtures with many components. Quote the method whenever you quote a number.

```r
sapply(c("elir", "moment", "morita"), function(m) ess(map_bin, method = m))
ess(rob) # robustification always costs ESS
```

A robustification weight grid makes the trade-off concrete:

```r
weights <- c(0, 0.1, 0.2, 0.5)
data.frame(
  weight = weights,
  ess = sapply(weights, function(w) ess(robustify(map_bin, weight = w, mean = 0.5)))
)
```

### ESS on a link scale

Normal mixtures get a `family` argument, so the ESS of a mixture living on a
transformed scale is evaluated in the exponential family it actually belongs to.
This is what you want for a logit-scale prior destined for a brms model:

```r
map_logit <- mixnorm(c(1, logit(0.25), 0.4), sigma = 2)
ess(map_logit, family = binomial)
```

For normal mixtures on their own scale, `sigma` is taken from the mixture; pass
`sigma = ` explicitly when the mixture carries none:

```r
map_rate <- mixnorm(c(1, log(1.8), 0.15), sigma = 1.58)
ess(map_rate, sigma = 1.58)
```

## Combining mixtures

`mixcombine(..., weight, rescale = TRUE)` takes mixtures of the same class:

```r
inf <- mixbeta(c(1, 18, 133))
vague <- mixbeta(c(1, 1, 1))
mixcombine(inf, vague, weight = c(0.8, 0.2))
```

`robustify` is exactly this call with a canonical vague component.

## Building a prior from elicited summaries

Constructors accept alternative parametrizations rather than requiring natural
parameters:

```r
mixbeta(c(1, 0.25, 0.1), param = "ms") # mean and standard deviation
mixbeta(c(1, 0.25, 30), param = "mn") # mean and number of observations
mixnorm(c(1, -50, 20), param = "mn", sigma = 88) # n via s = sigma/sqrt(n)
```

The standalone converters return the natural parameters:

```r
ms2beta(m = 0.25, s = 0.1)
mn2beta(m = 0.25, n = 30)
ms2gamma(m = 1.8, s = 0.4)
mn2gamma(m = 1.8, n = 20, likelihood = "poisson")
```

`msr2mvnorm` does the same for multivariate normals from means, SDs and
correlations.

## Reference scale

Normal mixtures carry a reference σ, and most downstream functions need it:

```r
sigma(map_rate)
sigma(map_rate) <- 1.6
```

A mixture derived from a `gMAP` fit inherits σ only when the fit could estimate
it — for gaussian models that means `weights` was supplied. Otherwise set it by
hand, as above, before any `postmix`, `oc*` or `pos*` call.

## Sanity checks

1. Does `plot(mix)$mix` overlay the sample well? (For `mixfit`/`automixfit` output.)
2. Is the ESS plausible against the number of historical patients and the
   heterogeneity between the studies?
3. Does the robust component sit at the null, and is its weight defensible?

## Next steps

[analysis](./analysis.md) to use the prior on data, [design](./design.md) for
operating characteristics, [interop](./interop.md) to hand it to brms.
