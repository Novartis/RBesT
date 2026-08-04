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

Holzhauer B, Wang C, Schmidli H (2018). “Historical control information
for clinical trials with a recurrent event endpoint.” *Statistics in
Medicine*, **37**(10), 1640–1657.

## Examples

``` r
.user_mc_options <- options()

set.seed(34563)
asthma_ph3 <- subset(asthma, phase == "phase III")
map_asthma <- gMAP(cbind(log_mu_hat, se_log_mu_hat) ~ 1 + offset(log(d)) | study,
  family = gaussian,
  data = asthma_ph3,
  tau.dist = "HalfNormal", tau.prior = 0.5,
  beta.prior = cbind(0, 2)
)
## Recover user set sampling defaults
options(.user_mc_options)
```
