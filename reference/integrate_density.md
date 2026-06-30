# Integrate a function against a mixture density

Computes \\\int g(x) f\_{\text{mix}}(x) dx\\ where \\f\_{\text{mix}}\\
is a mixture density. S3 generic dispatching on the mixture type.

## Usage

``` r
# S3 method for class 'normMix'
integrate_density(mix, integrand, ...)

integrate_density(mix, integrand, ...)

# Default S3 method
integrate_density(
  mix,
  integrand,
  Lplower = -Inf,
  Lpupper = Inf,
  eps = getOption("RBesT.integrate_prob_eps", 1e-06),
  ...
)
```

## Arguments

- mix:

  mixture density to integrate over (dispatch argument)

- integrand:

  function to integrate

- ...:

  additional arguments passed to methods

## Methods (by class)

- `integrate_density(normMix)`: Gauss-Hermite method for normMix

- `integrate_density(default)`: Default method using adaptive
  integration
