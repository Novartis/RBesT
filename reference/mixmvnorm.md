# Multivariate Normal Mixture Density

The multivariate normal mixture density and auxiliary functions.

## Usage

``` r
mixmvnorm(..., sigma, param = c("ms", "mn", "msr"))

msr2mvnorm(m, s, r, unlist = TRUE)

# S3 method for class 'mvnormMix'
print(x, ...)

# S3 method for class 'mvnormMix'
summary(object, ...)

# S3 method for class 'mvnormMix'
sigma(object, ...)
```

## Arguments

- ...:

  List of mixture components.

- sigma:

  Reference covariance.

- param:

  Determines how the parameters in the list are interpreted. See
  details.

- m:

  Mean vector.

- s:

  Standard deviation vector.

- r:

  Vector of correlations in column-major format of the lower triangle of
  the correlation matrix.

- unlist:

  Logical. Controls whether the result is a flattened vector (`TRUE`) or
  a list with mean `m` and covariance `s` (`FALSE`). Defaults to `TRUE`.

- object, x:

  Multivariate normal mixture object.

## Value

Returns a multivariate normal mixture with the specified mixture
components.

## Details

Each entry in the `...` argument list is a numeric vector defining one
component of the mixture multivariate normal distribution. The first
entry of the component defining vector is the weight of the mixture
component followed by the vector of means in each dimension and finally
a specification of the covariance matrix, which depends on the chosen
parametrization. The covariance matrix is expected to be given as
numeric vector in a column-major format, which is standard conversion
applied to matrices by the vector concatenation function
[`base::c()`](https://rdrr.io/r/base/c.html). Please refer to the
examples section below.

Each component defining vector can be specified in different ways as
determined by the `param` option:

- ms:

  Mean vector and covariance matrix `s`. Default.

- mn:

  Mean vector and number of observations. `n` determines the covariance
  for each component via the relation \\\Sigma/n\\ with \\\Sigma\\ being
  the known reference covariance.

- msr:

  Mean vector, standard deviations and correlations in column-major
  format (corresponds to order when printing multi-variate normal
  mixtures).

The reference covariance \\\Sigma\\ is the known covariance in the
normal-normal model (observation covariance). The function `sigma` can
be used to query the reference covariance and may also be used to assign
a new reference covariance, see examples below. In case `sigma` is not
specified, the user has to supply `sigma` as argument to functions which
require a reference covariance.

## See also

Other mixdist:
[`mix`](https://opensource.nibr.com/RBesT/reference/mix.md),
[`mixbeta()`](https://opensource.nibr.com/RBesT/reference/mixbeta.md),
[`mixcombine()`](https://opensource.nibr.com/RBesT/reference/mixcombine.md),
[`mixgamma()`](https://opensource.nibr.com/RBesT/reference/mixgamma.md),
[`mixjson`](https://opensource.nibr.com/RBesT/reference/mixjson.md),
[`mixnorm()`](https://opensource.nibr.com/RBesT/reference/mixnorm.md),
[`mixplot`](https://opensource.nibr.com/RBesT/reference/mixplot.md)

## Examples

``` r
# default mean & covariance parametrization
S <- diag(c(1, 2)) %*% matrix(c(1, 0.5, 0.5, 1), 2, 2) %*% diag(c(1, 2))
mvnm1 <- mixmvnorm(
  rob = c(0.2, c(0, 0), diag(c(2, 2)^2)),
  inf = c(0.8, c(0.5, 1), S / 4), sigma = S
)

print(mvnm1)
#> Multivariate normal mixture
#> Outcome dimension: 2
#> Reference covariance:
#>   1 2
#> 1 1 1
#> 2 1 4
#> Mixture Components:
#>          rob inf
#> w        0.2 0.8
#> m[1]     0.0 0.5
#> m[2]     0.0 1.0
#> s[1]     2.0 0.5
#> s[2]     2.0 1.0
#> rho[2,1] 0.0 0.5
summary(mvnm1)
#> $mean
#>   1   2 
#> 0.4 0.8 
#> 
#> $cov
#>      1    2
#> 1 1.04 0.28
#> 2 0.28 1.76
#> 

set.seed(657846)
mixSamp1 <- rmix(mvnm1, 500)
colMeans(mixSamp1)
#>         1         2 
#> 0.3709829 0.8039847 

# alternative mean, sd and correlation parametrization
mvnm1_alt <- mixmvnorm(
  rob = c(0.2, c(0, 0), c(2, 2), 0.0),
  inf = c(0.8, c(0.5, 1), c(1, 2) / 2, 0.5),
  sigma = msr2mvnorm(s = c(1, 2), r = 0.5, unlist = FALSE)$s,
  param = "msr"
)

print(mvnm1_alt)
#> Multivariate normal mixture
#> Outcome dimension: 2
#> Reference covariance:
#>   1 2
#> 1 1 1
#> 2 1 4
#> Mixture Components:
#>          rob inf
#> w        0.2 0.8
#> m[1]     0.0 0.5
#> m[2]     0.0 1.0
#> s[1]     2.0 0.5
#> s[2]     2.0 1.0
#> rho[2,1] 0.0 0.5
```
