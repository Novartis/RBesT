# Scrambles the order of a mcmc array object for usage as a mcmc sample. It is advisable to set order once per mcmc run, otherwise correlations in the mcmc sample will be lost.

Scrambles the order of a mcmc array object for usage as a mcmc sample.
It is advisable to set order once per mcmc run, otherwise correlations
in the mcmc sample will be lost.

## Usage

``` r
chains2sample(chains, order, drop = TRUE)
```
