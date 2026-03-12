# Effective Sample Size for a Conjugate Prior

Calculates the Effective Sample Size (ESS) for a mixture prior. The ESS
indicates how many experimental units the prior is roughly equivalent
to.

## Usage

``` r
ess(mix, method = c("elir", "moment", "morita"), ...)

# S3 method for class 'betaMix'
ess(mix, method = c("elir", "moment", "morita"), ..., s = 100)

# S3 method for class 'gammaMix'
ess(mix, method = c("elir", "moment", "morita"), ..., s = 100, eps = 1e-04)

# S3 method for class 'normMix'
ess(
  mix,
  method = c("elir", "moment", "morita"),
  ...,
  family = gaussian,
  sigma,
  s = 100
)
```

## Arguments

- mix:

  Prior (mixture of conjugate distributions).

- method:

  Selects the used method. Can be either `elir` (default), `moment` or
  `morita`.

- ...:

  Optional arguments applicable to specific methods.

- s:

  For `morita` method large constant to ensure that the prior scaled by
  this value is vague (default 100); see Morita et al. (2008) for
  details.

- eps:

  Probability mass left out from the numerical integration of the
  expected information for the Poisson-Gamma case of Morita method
  (defaults to 1E-4).

- family:

  defines data likelihood and link function (`binomial`, `gaussian`, or
  `poisson`).

- sigma:

  reference scale.

## Value

Returns the ESS of the prior as floating point number.

## Details

The ESS is calculated using either the expected local information ratio
(elir) *Neuenschwander et al. (2020)*, the moments approach or the
method by *Morita et al. (2008)*.

The elir approach measures effective sample size in terms of the average
curvature of the prior in relation to the Fisher information. Informally
this corresponds to the average peakiness of the prior in relation to
the information content of a single observation. The elir approach is
the only ESS which fulfills predictive consistency. The predictive
consistency of the ESS requires that the ESS of a prior is consistent
when considering an averaged posterior ESS of additional data
distributed according to the predictive distribution of the prior. The
expectation of the posterior ESS is taken wrt to the prior predictive
distribution and the averaged posterior ESS corresponds to the sum of
the prior ESS and the number of forward simulated data items. The elir
approach results in ESS estimates which are neither conservative nor
liberal whereas the moments method yields conservative and the morita
method liberal results. See the example section for a demonstration of
predictive consistency.

For the moments method the mean and standard deviation of the mixture
are calculated and then approximated by the conjugate distribution with
the same mean and standard deviation. For conjugate distributions, the
ESS is well defined. See the examples for a step-wise calculation in the
beta mixture case.

The Morita method used here evaluates the mixture prior at the mode
instead of the mean as proposed originally by Morita. The method may
lead to very optimistic ESS values, especially if the mixture contains
many components. The calculation of the Morita approach here follows the
approach presented in Neuenschwander B. et all (2019) which avoids the
need for a minimization and does not restrict the ESS to be an integer.

The arguments `sigma` and `family` are specific for normal mixture
densities. These specify the sampling standard deviation for a
`gaussian` family (the default) while also allowing to consider the ESS
of standard one-parameter exponential families, i.e. `binomial` or
`poisson`. The function supports non-gaussian families with unit
dispersion only.

## Methods (by class)

- `ess(betaMix)`: ESS for beta mixtures.

- `ess(gammaMix)`: ESS for gamma mixtures.

- `ess(normMix)`: ESS for normal mixtures.

## Supported Conjugate Prior-Likelihood Pairs

|                     |                             |                             |                |
|---------------------|-----------------------------|-----------------------------|----------------|
| **Prior/Posterior** | **Likelihood**              | **Predictive**              | **Summaries**  |
| Beta                | Binomial                    | Beta-Binomial               | `n`, `r`       |
| Normal              | Normal (*fixed \\\sigma\\*) | Normal                      | `n`, `m`, `se` |
| Gamma               | Poisson                     | Gamma-Poisson               | `n`, `m`       |
| Gamma               | Exponential                 | Gamma-Exp (*not supported*) | `n`, `m`       |

## References

Morita S, Thall PF, Mueller P. Determining the effective sample size of
a parametric prior. *Biometrics* 2008;64(2):595-602.

Neuenschwander B., Weber S., Schmidli H., O’Hagan A. (2020).
Predictively consistent prior effective sample sizes. *Biometrics*,
76(2), 578–587. https://doi.org/10.1111/biom.13252

## Examples

``` r
# Conjugate Beta example
a <- 5
b <- 15
prior <- mixbeta(c(1, a, b))

ess(prior)
#> [1] 20
(a + b)
#> [1] 20

# Beta mixture example
bmix <- mixbeta(rob = c(0.2, 1, 1), inf = c(0.8, 10, 2))

ess(bmix, "elir")
#> [1] 7.65152

ess(bmix, "moment")
#> [1] 3.161034
# moments method is equivalent to
# first calculate moments
bmix_sum <- summary(bmix)
# then calculate a and b of a matching beta
ab_matched <- ms2beta(bmix_sum["mean"], bmix_sum["sd"])
# finally take the sum of a and b which are equivalent
# to number of responders/non-responders respectivley
round(sum(ab_matched))
#> [1] 3

ess(bmix, method = "morita")
#> [1] 8.487603

# One may also calculate the ESS on the logit scale, which
# gives slightly different results due to the parameter
# transformation, e.g.:
prior_logit <- mixnorm(c(1, log(5 / 15), sqrt(1 / 5 + 1 / 15)))
ess(prior_logit, family = binomial)
#> [1] 21.78289

bmix_logit <- mixnorm(
  rob = c(0.2, 0, 2),
  inf = c(0.8, log(10 / 2), sqrt(1 / 10 + 1 / 2))
)
ess(bmix_logit, family = binomial)
#> [1] 10.12276

# Predictive consistency of elir
n_forward <- 1E1
bmixPred <- preddist(bmix, n = n_forward)
pred_samp <- rmix(bmixPred, 1E2)
# use more samples here for greater accuracy, e.g.
# pred_samp <- rmix(bmixPred, 1E3)
pred_ess <- sapply(pred_samp, function(r) {
  ess(postmix(bmix, r = r, n = n_forward), "elir")
})
ess(bmix, "elir")
#> [1] 7.65152
mean(pred_ess) - n_forward
#> [1] 7.54457

# Normal mixture example
nmix <- mixnorm(rob = c(0.5, 0, 2), inf = c(0.5, 3, 4), sigma = 10)

ess(nmix, "elir")
#> Using default prior reference scale 10
#> [1] 10.82796

ess(nmix, "moment")
#> Using default prior reference scale 10
#> [1] 8.163265

# the reference scale determines the ESS
sigma(nmix) <- 20
ess(nmix)
#> Using default prior reference scale 20
#> [1] 43.31185

# we may also interpret normal mixtures as densities assigned to
# parameters of a logit transformed response rate of a binomial
nmix_logit <- mixnorm(c(1, logit(1 / 4), 2 / sqrt(10)))
ess(nmix_logit, family = binomial)
#> [1] 15.17836

# Gamma mixture example
gmix <- mixgamma(rob = c(0.3, 20, 4), inf = c(0.7, 50, 10))

ess(gmix) ## interpreted as appropriate for a Poisson likelihood (default)
#> [1] 7.159378

likelihood(gmix) <- "exp"
ess(gmix) ## interpreted as appropriate for an exponential likelihood
#> [1] 34.93388
```
