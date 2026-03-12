# internal function used for integration of densities which appears to be much more stable from -Inf to +Inf in the logit space while the density to be integrated recieves inputs from 0 to 1 such that the inverse distribution function must be used. The integral solved is int_x dmix(mix,x) integrand(x) where integrand must be given as log and we integrate over the support of mix.

integrate density in logit space and split by component such that the
quantile function of each component is used. This ensures that the R
implementation of the quantile function is always used.

## Usage

``` r
integrate_density_log(
  log_integrand,
  mix,
  Lplower = -Inf,
  Lpupper = Inf,
  eps = getOption("RBesT.integrate_prob_eps", 1e-06)
)
```

## Arguments

- log_integrand:

  function to integrate over which must return the log(f)

- mix:

  density over which to integrate

- Lplower:

  logit of lower cumulative density

- Lpupper:

  logit of upper cumulative density
