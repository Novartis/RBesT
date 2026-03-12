# The Gamma Mixture Distribution

The gamma mixture density and auxiliary functions.

## Usage

``` r
mixgamma(..., param = c("ab", "ms", "mn"), likelihood = c("poisson", "exp"))

ms2gamma(m, s, drop = TRUE)

mn2gamma(m, n, likelihood = c("poisson", "exp"), drop = TRUE)

# S3 method for class 'gammaMix'
print(x, ...)

# S3 method for class 'gammaPoissonMix'
print(x, ...)

# S3 method for class 'gammaExpMix'
print(x, ...)

# S3 method for class 'gammaMix'
summary(object, probs = c(0.025, 0.5, 0.975), ...)

# S3 method for class 'gammaPoissonMix'
summary(object, probs = c(0.025, 0.5, 0.975), ...)
```

## Arguments

- ...:

  List of mixture components.

- param:

  Determines how the parameters in the list are interpreted. See
  details.

- likelihood:

  Defines with what likelihood the Gamma density is used (Poisson or
  Exp). Defaults to `poisson`.

- m:

  Vector of means of the Gamma mixture components

- s:

  Vector of standard deviations of the gamma mixture components,

- drop:

  Delete the dimensions of an array which have only one level.

- n:

  Vector of sample sizes of the Gamma mixture components.

- x:

  The mixture to print

- object:

  Gamma mixture object.

- probs:

  Quantiles reported by the `summary` function.

## Value

`mixgamma` returns a gamma mixture with the specified mixture
components. `ms2gamma` and `mn2gamma` return the equivalent natural `a`
and `b` parametrization given parameters `m`, `s`, or `n`.

## Details

Each entry in the `...` argument list is expected to be a triplet of
numbers which defines the weight \\w_k\\, first and second parameter of
the mixture component \\k\\. A triplet can optionally be named which
will be used appropriately.

The first and second parameter can be given in different
parametrizations which is set by the `param` option:

- ab:

  Natural parametrization of Gamma density (`a`=shape and `b`=rate).
  Default.

- ms:

  Mean and standard deviation, \\m=a/b\\ and \\s=\sqrt{a}/b\\.

- mn:

  Mean and number of observations. Translation to natural parameter
  depends on the `likelihood` argument. For a Poisson likelihood \\n=b\\
  (and \\a=m \cdot n\\), for an Exp likelihood \\n=a\\ (and \\b=n/m\\).

## See also

Other mixdist:
[`mix`](https://opensource.nibr.com/RBesT/reference/mix.md),
[`mixbeta()`](https://opensource.nibr.com/RBesT/reference/mixbeta.md),
[`mixcombine()`](https://opensource.nibr.com/RBesT/reference/mixcombine.md),
[`mixjson`](https://opensource.nibr.com/RBesT/reference/mixjson.md),
[`mixmvnorm()`](https://opensource.nibr.com/RBesT/reference/mixmvnorm.md),
[`mixnorm()`](https://opensource.nibr.com/RBesT/reference/mixnorm.md),
[`mixplot`](https://opensource.nibr.com/RBesT/reference/mixplot.md)

## Examples

``` r
# Gamma mixture with robust and informative component
gmix <- mixgamma(rob = c(0.3, 20, 4), inf = c(0.7, 50, 10))

# objects can be printed
gmix
#> Univariate Gamma mixture
#> Mixture Components:
#>   rob  inf 
#> w  0.3  0.7
#> a 20.0 50.0
#> b  4.0 10.0
# or explicitly
print(gmix)
#> Univariate Gamma mixture
#> Mixture Components:
#>   rob  inf 
#> w  0.3  0.7
#> a 20.0 50.0
#> b  4.0 10.0

# summaries are defined
summary(gmix)
#>      mean        sd      2.5%     50.0%     97.5% 
#> 5.0000000 0.8514693 3.4362134 4.9560695 6.8210139 

# sub-components may be extracted
# by component number
gmix[[2]]
#> Univariate Gamma mixture
#> Mixture Components:
#>   inf 
#> w  0.7
#> a 50.0
#> b 10.0
# or component name
gmix[["inf"]]
#> Univariate Gamma mixture
#> Mixture Components:
#>   inf 
#> w  0.7
#> a 50.0
#> b 10.0

# alternative mean and standard deviation parametrization
gmsMix <- mixgamma(rob = c(0.5, 8, 0.5), inf = c(0.5, 9, 2), param = "ms")

# or mean and number of observations parametrization
gmnMix <- mixgamma(rob = c(0.2, 2, 1), inf = c(0.8, 2, 5), param = "mn")

# and mixed parametrizations are also possible
gfmix <- mixgamma(rob1 = c(0.15, mn2gamma(2, 1)),
                  rob2 = c(0.15, ms2gamma(2, 5)),
                  inf = c(0.7, 50, 10))
```
