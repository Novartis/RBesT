# Integrate a log-space function against a mixture density

Computes \\\int \exp(\text{log\\integrand}(x)) f\_{\text{mix}}(x) dx\\.
Uses logit-space transformation for numerical stability with adaptive
integration. S3 generic dispatching on the mixture type.

## Usage

``` r
# S3 method for class 'normMix'
integrate_density_log(mix, log_integrand, ...)

# S3 method for class 'betaMix'
integrate_density_log(mix, log_integrand, ...)

# S3 method for class 'gammaMix'
integrate_density_log(mix, log_integrand, ...)

integrate_density_log(mix, log_integrand, ...)

# Default S3 method
integrate_density_log(
  mix,
  log_integrand,
  Lplower = -Inf,
  Lpupper = Inf,
  eps = getOption("RBesT.integrate_prob_eps", 1e-06),
  ...
)
```

## Arguments

- mix:

  mixture density to integrate over (dispatch argument)

- log_integrand:

  function returning `log(g(x))`

- ...:

  additional arguments passed to methods

## Methods (by class)

- `integrate_density_log(normMix)`: Gauss-Hermite method for normMix

- `integrate_density_log(betaMix)`: Gauss-Jacobi method for betaMix

- `integrate_density_log(gammaMix)`: Gauss-Laguerre method for gammaMix

- `integrate_density_log(default)`: Default method using adaptive
  logit-space integration
