# Asthma exacerbation recurrent event data.

Data set containing historical information for placebo arms of relevant
trials for the treatment of asthma. The primary outcome is the rate of
asthma exacerbations, a recurrent event modelled with a negative
binomial distribution. The full data set as published in Holzhauer, Wang
& Schmidli (2018) summarizes ten historical placebo arms by their
back-calculated log mean event rate and dispersion together with the
associated standard errors on the log scale.

## Usage

``` r
asthma
```

## Format

A data frame with 10 rows and 11 variables:

- study:

  study label

- NCT:

  ClinicalTrials.gov / ISRCTN registry identifier(s)

- d:

  follow-up (exposure) duration in years

- n:

  study size

- mu_hat:

  estimated mean event rate

- log_mu_hat:

  log mean event rate

- se_log_mu_hat:

  standard error of the log mean event rate

- kappa_hat:

  estimated dispersion parameter

- log_kappa_hat:

  log dispersion parameter

- se_log_kappa_hat:

  standard error of the log dispersion parameter

- phase:

  development phase of the trial

## References

Holzhauer B., Wang C. and Schmidli H. *Statistics in Medicine*, 2018,
37(10):1640-1657

## Examples

``` r
## Setting up dummy sampling for fast execution of example
## Please use 4 chains and 20x more warmup & iter in practice
.user_mc_options <- options(RBesT.MC.warmup=50, RBesT.MC.iter=100,
                            RBesT.MC.chains=2, RBesT.MC.thin=1)

set.seed(34563)
asthma_ph3 <- subset(asthma, phase == "phase III")
map_asthma <- gMAP(cbind(log_mu_hat, se_log_mu_hat) ~ 1 + offset(log(d)) | study,
  family = gaussian,
  data = asthma_ph3,
  tau.dist = "HalfNormal", tau.prior = 0.5,
  beta.prior = cbind(0, 2)
)
#> Warning: The largest R-hat is 1.12, indicating chains have not mixed.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#r-hat
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
#> Final MCMC sample equivalent to less than 1000 independent draws.
#> Please consider increasing the MCMC simulation size.
## Recover user set sampling defaults
options(.user_mc_options)
```
