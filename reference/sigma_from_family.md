# Compute sigma as a function of eta for a GLM family

For a GLM family with variance function V(mu), link function g, and
dispersion parameter phi, the sampling standard deviation of the MLE on
the link scale is: sigma(eta) = sqrt(phi \* V(mu)) / \|dmu/deta\| where
mu = linkinv(eta + offset).

## Usage

``` r
sigma_from_family(eta, family, offset = 0, phi = 1)
```

## Arguments

- eta:

  Numeric vector of values on the link scale (parameter of interest).

- family:

  An R family object (e.g., binomial(), MASS::negative.binomial(0.5)).

- offset:

  Numeric scalar. Shift applied before evaluating the variance function.
  For NB/Poisson with log link, this is log(exposure). Default 0.

- phi:

  Numeric scalar. Dispersion parameter. For gaussian family this is
  sigma^2 (user-supplied). For binomial/NB/Poisson this is 1 (the
  default).

## Value

Numeric vector of sigma values (same length as eta).
