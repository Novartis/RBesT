# Conjugate Posterior Analysis

Calculates the posterior distribution for data `data` given a prior
`priormix`, where the prior is a mixture of conjugate distributions. The
posterior is then also a mixture of conjugate distributions.

## Usage

``` r
postmix(priormix, data, ...)

# S3 method for class 'betaMix'
postmix(priormix, data, n, r, ...)

# S3 method for class 'normMix'
postmix(priormix, data, n, m, se, ...)

# S3 method for class 'gammaMix'
postmix(priormix, data, n, m, ...)
```

## Arguments

- priormix:

  prior (mixture of conjugate distributions).

- data:

  individual data. If the individual data is not given, then summary
  data has to be provided (see below).

- ...:

  includes arguments which depend on the specific case, see description
  below.

- n:

  sample size.

- r:

  Number of successes.

- m:

  Sample mean.

- se:

  Sample standard error.

## Details

A conjugate prior-likelihood pair has the convenient property that the
posterior is in the same distributional class as the prior. This
property also applies to mixtures of conjugate priors. Let

\$\$p(\theta;\mathbf{w},\mathbf{a},\mathbf{b})\$\$

denote a conjugate mixture prior density for data

\$\$y\|\theta \sim f(y\|\theta),\$\$

where \\f(y\|\theta)\\ is the likelihood. Then the posterior is again a
mixture with each component \\k\\ equal to the respective posterior of
the \\k\\th prior component and updated weights \\w'\_k\\,

\$\$p(\theta;\mathbf{w'},\mathbf{a'},\mathbf{b'}\|y) = \sum\_{k=1}^K
w'\_k \\ p_k(\theta;a'\_k,b'\_k\|y).\$\$

The weight \\w'\_k\\ for \\k\\th component is determined by the marginal
likelihood of the new data \\y\\ under the \\k\\th prior distribution
which is given by the predictive distribution of the \\k\\th component,

\$\$w'\_k \propto w_k \\ \int p_k(\theta;a_k,b_k) \\ f(y\|\theta) \\
d\theta \equiv w^\ast_k .\$\$

The final weight \\w'\_k\\ is then given by appropriate normalization,
\\w'\_k = w^\ast_k / \sum\_{k=1}^K w^\ast_k\\. In other words, the
weight of component \\k\\ is proportional to the likelihood that data
\\y\\ is generated from the respective component, i.e. the marginal
probability; for details, see for example *Schmidli et al., 2015*.

*Note:* The prior weights \\w_k\\ are fixed, but the posterior weights
\\w'\_k \neq w_k\\ still change due to the changing normalization.

The data \\y\\ can either be given as individual data or as summary data
(sufficient statistics). See below for details for the implemented
conjugate mixture prior densities.

## Methods (by class)

- `postmix(betaMix)`: Calculates the posterior beta mixture
  distribution. The individual data vector is expected to be a vector of
  0 and 1, i.e. a series of Bernoulli experiments. Alternatively, the
  sufficient statistics `n` and `r` can be given, i.e. number of trials
  and successes, respectively.

- `postmix(normMix)`: Calculates the posterior normal mixture
  distribution with the sampling likelihood being a normal with fixed
  standard deviation. Either an individual data vector `data` can be
  given or the sufficient statistics which are the standard error `se`
  and sample mean `m`. If the sample size `n` is used instead of the
  sample standard error, then the reference scale of the prior is used
  to calculate the standard error. Should standard error `se` and sample
  size `n` be given, then the reference scale of the prior is updated;
  however it is recommended to use the command
  [`sigma()`](https://opensource.nibr.com/RBesT/reference/mixnorm.md) to
  set the reference standard deviation.

- `postmix(gammaMix)`: Calculates the posterior gamma mixture
  distribution for Poisson and exponential likelihoods. Only the Poisson
  case is supported in this version.

## Supported Conjugate Prior-Likelihood Pairs

|                     |                             |                             |                |
|---------------------|-----------------------------|-----------------------------|----------------|
| **Prior/Posterior** | **Likelihood**              | **Predictive**              | **Summaries**  |
| Beta                | Binomial                    | Beta-Binomial               | `n`, `r`       |
| Normal              | Normal (*fixed \\\sigma\\*) | Normal                      | `n`, `m`, `se` |
| Gamma               | Poisson                     | Gamma-Poisson               | `n`, `m`       |
| Gamma               | Exponential                 | Gamma-Exp (*not supported*) | `n`, `m`       |

## References

Schmidli H, Gsteiger S, Roychoudhury S, O'Hagan A, Spiegelhalter D,
Neuenschwander B. Robust meta-analytic-predictive priors in clinical
trials with historical control information. *Biometrics*
2014;70(4):1023-1032.

## Examples

``` r
# binary example with individual data (1=event,0=no event), uniform prior
prior.unif <- mixbeta(c(1, 1, 1))
data.indiv <- c(1, 0, 1, 1, 0, 1)
posterior.indiv <- postmix(prior.unif, data.indiv)
print(posterior.indiv)
#> Univariate beta mixture
#> Mixture Components:
#>   comp1
#> w 1    
#> a 5    
#> b 3    
# or with summary data (number of events and number of patients)
r <- sum(data.indiv)
n <- length(data.indiv)
posterior.sum <- postmix(prior.unif, n = n, r = r)
print(posterior.sum)
#> Univariate beta mixture
#> Mixture Components:
#>   comp1
#> w 1    
#> a 5    
#> b 3    

# binary example with robust informative prior and conflicting data
prior.rob <- mixbeta(c(0.5, 4, 10), c(0.5, 1, 1))
posterior.rob <- postmix(prior.rob, n = 20, r = 18)
print(posterior.rob)
#> Univariate beta mixture
#> Mixture Components:
#>   comp1        comp2       
#> w  0.002672948  0.997327052
#> a 22.000000000 19.000000000
#> b 12.000000000  3.000000000

# normal example with individual data
sigma <- 88
prior.mean <- -49
prior.se <- sigma / sqrt(20)
prior <- mixnorm(c(1, prior.mean, prior.se), sigma = sigma)
data.indiv <- c(-46, -227, 41, -65, -103, -22, 7, -169, -69, 90)
posterior.indiv <- postmix(prior, data.indiv)
# or with summary data (mean and number of patients)
mn <- mean(data.indiv)
n <- length(data.indiv)
posterior.sum <- postmix(prior, m = mn, n = n)
#> Using default prior reference scale 88
print(posterior.sum)
#> Univariate normal mixture
#> Reference scale: 88
#> Mixture Components:
#>   comp1    
#> w   1.00000
#> m -51.43333
#> s  16.06653
```
