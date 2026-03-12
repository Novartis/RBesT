# Ankylosing Spondylitis.

Data set containing historical information for placebo for a phase II
trial of ankylosing spondylitis patients. The primary efficacy endpoint
was the percentage of patients with a 20% response according to the
Assessment of SpondyloArthritis international Society criteria for
improvement (ASAS20) at week 6.

## Usage

``` r
AS
```

## Format

A data frame with 8 rows and 3 variables:

- study:

  study

- n:

  study size

- r:

  number of events

## References

Baeten D. et. al, *The Lancet*, 2013, (382), 9906, p 1705

## Examples

``` r
## Setting up dummy sampling for fast execution of example
## Please use 4 chains and 20x more warmup & iter in practice
.user_mc_options <- options(RBesT.MC.warmup=50, RBesT.MC.iter=100,
                            RBesT.MC.chains=2, RBesT.MC.thin=1)

set.seed(34563)
map_AS <- gMAP(cbind(r, n - r) ~ 1 | study,
  family = binomial,
  data = AS,
  tau.dist = "HalfNormal", tau.prior = 1,
  beta.prior = 2
)
#> Assuming default prior location   for beta: 0
#> Warning: The largest R-hat is 1.31, indicating chains have not mixed.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#r-hat
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
#> Warning: Maximal Rhat > 1.1. Consider increasing RBesT.MC.warmup MCMC parameter.
#> Final MCMC sample equivalent to less than 1000 independent draws.
#> Please consider increasing the MCMC simulation size.
## Recover user set sampling defaults
options(.user_mc_options)
```
