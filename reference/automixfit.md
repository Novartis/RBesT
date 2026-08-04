# Automatic Fitting of Mixtures of Conjugate Distributions to a Sample

Fitting a series of mixtures of conjugate distributions to a `sample`,
using Expectation-Maximization (EM). The number of mixture components is
specified by the vector `Nc`. First a `Nc[1]` component mixture is
fitted, then a `Nc[2]` component mixture, and so on. The mixture
providing the best AIC value is then selected.

## Usage

``` r
automixfit(sample, Nc = seq(1, 4), k = 6, thresh = -Inf, verbose = FALSE, ...)
```

## Arguments

- sample:

  Sample to be fitted by a mixture distribution.

- Nc:

  Vector of mixture components to try out (default `seq(1,4)`).

- k:

  Penalty parameter for AIC calculation (default 6)

- thresh:

  The procedure stops if the difference of subsequent AIC values is
  smaller than this threshold (default -Inf). Setting the threshold to 0
  stops `automixfit` once the AIC becomes worse.

- verbose:

  Enable verbose logging.

- ...:

  Further arguments passed to
  [`mixfit()`](https://opensource.nibr.com/RBesT/reference/mixfit.md),
  including `type`.

## Value

As result the best fitting mixture model is returned, i.e. the model
with lowest AIC. All other models are saved in the attribute `models`.

## Details

The `type` argument specifies the distribution of the mixture
components, and can be a normal, beta or gamma distribution.

The penalty parameter `k` is 2 for the standard AIC definition. *Collet
(2003)* suggested to use values in the range from 2 to 6, where larger
values of `k` penalize more complex models. To favor mixtures with fewer
components a value of 6 is used as default.

## References

Collett D (2003). *Modelling Survival Data in Medical Research*, 2nd
edition. Chapman and Hall/CRC.

## Examples

``` r
# random sample of size 1000 from a mixture of 2 beta components
bm <- mixbeta(beta1 = c(0.4, 20, 90), beta2 = c(0.6, 35, 65))
bmSamp <- rmix(bm, 1000)

# fit with EM mixture models with up to 10 components and stop if
# AIC increases
bmFit <- automixfit(bmSamp, Nc = 1:10, thresh = 0, type = "beta")
bmFit
#> EM for Beta Mixture Model
#> Log-Likelihood = 1130.983
#> 
#> Univariate beta mixture
#> Mixture Components:
#>   comp1      comp2     
#> w  0.6098402  0.3901598
#> a 34.8965697 20.8913411
#> b 64.8756213 93.3172690

# advanced usage: find out about all discarded models
bmFitAll <- attr(bmFit, "models")

sapply(bmFitAll, AIC, k = 6)
#>         2         3         1 
#> -2231.965 -2212.242 -1921.889 
```
